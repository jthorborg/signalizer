/*************************************************************************************

	Signalizer - cross-platform audio visualization plugin - v. 0.x.y

	Copyright (C) 2016 Janus Lynggaard Thorborg (www.jthorborg.com)

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.

	See \licenses\ for additional details on licenses associated with this program.

**************************************************************************************

	file:COscilloscopeDSP.cpp

		Implementation of signal processing code for the oscilloscope.

*************************************************************************************/

#ifndef SIGNALIZER_OSCILLOSCOPEDSP_INL
#define SIGNALIZER_OSCILLOSCOPEDSP_INL

#include "Oscilloscope.h"
#include "StreamPreprocessing.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include <cpl/simd.h>
#include <cpl/lib/variable_array.h>

namespace Signalizer
{

	template<typename ISA, typename Eval>
	void Oscilloscope::analyseAndSetupState(const EvaluatorParams& params)
	{
		calculateFundamentalPeriod<ISA, Eval>(params);
		calculateTriggeringOffset<ISA, Eval>(params);

		resizeAudioStorage();

		if (state.envelopeMode == EnvelopeModes::PeakDecay)
		{
			runPeakFilter<ISA>();
		}
		else if (state.envelopeMode == EnvelopeModes::None)
		{
			shared.autoGainEnvelope.store(1, std::memory_order_release);
		}
	}

	template<typename ISA, typename Eval>
	void Oscilloscope::calculateFundamentalPeriod(const EvaluatorParams& params)
	{
#ifdef PHASE_VOCODER
		auto const TransformSize = OscilloscopeContent::LookaheadSize >> 1;
#else
		auto const TransformSize = OscilloscopeContent::LookaheadSize;
#endif

		if (state.customTrigger)
		{
			auto const normalizedFrequency = state.customTriggerFrequency / audioStream.getAudioHistorySamplerate();
			triggerState.record = BinRecord{ 0, 1, normalizedFrequency * TransformSize};

			auto fundamental = state.customTriggerFrequency;

			triggerState.fundamental = state.customTriggerFrequency;
			triggerState.cycleSamples = audioStream.getAudioHistorySamplerate() / fundamental;
		}
		else if(state.triggerMode == OscilloscopeContent::TriggeringMode::Spectral)
		{
			Eval eval(params);

			if (!eval.isWellDefined())
				return;

			// we will try to analyse the points closest to the sync point (latest in time)
			auto offset = std::max<std::size_t>(std::ceil(state.effectiveWindowSize), OscilloscopeContent::LookaheadSize);

			eval.startFrom(-static_cast<cpl::ssize_t>(offset));

			for (std::size_t i = 0; i < OscilloscopeContent::LookaheadSize; ++i)
			{
				transformBuffer[i] = eval.evaluateSampleInc();
			}


			signaldust::DustFFT_fwdDa(reinterpret_cast<double*>(transformBuffer.data()), TransformSize);
#ifdef PHASE_VOCODER
			signaldust::DustFFT_fwdDa(reinterpret_cast<double*>(transformBuffer.data() + TransformSize), TransformSize);
#endif
			// estimates the true frequency by calculating a bin offset to the current bin w
			auto quadDelta = [&](auto w) {
#if PHASE_VOCODER
				auto deltaBinOffset = std::arg(transformBuffer[w]) - std::arg(transformBuffer[TransformSize + w]);

				deltaBinOffset -= w * 2 * M_PI;

				while (deltaBinOffset < 2 * M_PI)
					deltaBinOffset += 2 * M_PI;

				while (deltaBinOffset > 2 * M_PI)
					deltaBinOffset -= 2 * M_PI;

				return deltaBinOffset / (2 * M_PI);
#else
				const auto
					x0 = transformBuffer[w],
					x1 = transformBuffer[w + 1],
					xm1 = transformBuffer[w == 0 ? 1 : w - 1];

				const auto denom = x0 * 2.0 - xm1 - x1;

				return (denom.real() + denom.imag()) != 0 ? std::real((xm1 - x1) / denom) : 0;
#endif
			};

			const double quarterSemitone = std::pow(2, 0.25 / 12.0) - 1;
			const double threshold = content->triggerThreshold.getTransformedValue();
			const double hysteresis = content->triggerHysteresis.getTransformedValue();
			const auto invHysteresis = 1 - hysteresis;

			// Reduce to 1/4 (to slip through "vastly better case" + half transform size
			BinRecord max{ 1, std::max(threshold * TransformSize / 6.0, std::abs(transformBuffer[1])), quadDelta(1) };

			for (std::size_t i = 2; i < (TransformSize >> 1); ++i)
			{
				BinRecord current{ i, std::abs(transformBuffer[i]) };

				// candidate must be vastly better
				if (invHysteresis * current.value > max.value * 2)
				{
					// weird parabolas
					if (max.omega() > 0)
					{
						// check if it is somewhat harmonically related, in which case we discard the candidate
						current.offset = quadDelta(i);

						// harmonic relationship
						auto factor = current.omega() / max.omega();

						auto sensivity = current.value / max.value;

						// shortcut if the value is 20 times bigger
						if (invHysteresis * sensivity > 20)
						{
							max = current;
							continue;
						}

						// the same value, just a better estimate, from another bin
						// TODO: fix this case by polynomially interpolate the value as well
						if (std::abs(1 - factor) < invHysteresis * quarterSemitone)
						{
							max = current;
							continue;
						}

						auto multipleDeviation = std::abs(factor - std::floor(factor + 0.5));

						// check if the harmonic series is more than half a semi-tone away, in which case we take the candidate
						if (invHysteresis * std::abs(multipleDeviation) > quarterSemitone)
						{
							max = current;
						}
					}
					else
					{
						max = current;
						max.offset = quadDelta(max.index);
					}

				}
			}

			// copy old filter
			auto localMedian = medianTriggerFilter;

			// store new data
			medianTriggerFilter[medianPos].record = max;

			medianPos++;
			medianPos &= (MedianData::FilterSize - 1);

			const auto middle = (MedianData::FilterSize >> 1);
			std::nth_element(
				localMedian.begin(),
				localMedian.begin() + middle,
				localMedian.end(),
				[](const auto & a, const auto & b)
				{
					return a.record.index < b.record.index;
				}
			);

			auto & oldMedianBin = localMedian[middle];

			// check to discard (temporarily) much higher frequencies through a median filter
			if (oldMedianBin.record.index != -1 && std::abs(max.omega() - (oldMedianBin.record.omega())) > 0.5)
			{
				max = oldMedianBin.record;
			}

			triggerState.record = max;

			// interpolate peak position quadratically

			auto fundamental = audioStream.getAudioHistorySamplerate() * (max.omega()) / TransformSize;

			triggerState.fundamental = fundamental = std::max(5.0, fundamental);
			triggerState.cycleSamples = audioStream.getAudioHistorySamplerate() / fundamental;
		}
	}

	inline double Oscilloscope::getGain()
	{
		return (std::isfinite(state.autoGain) ? state.autoGain : 1) * state.manualGain;
	}

	template<typename ISA, typename Eval>
	void Oscilloscope::calculateTriggeringOffset(const EvaluatorParams& params)
	{
		if (state.triggerMode == OscilloscopeContent::TriggeringMode::EnvelopeHold || state.triggerMode == OscilloscopeContent::TriggeringMode::ZeroCrossing)
		{
			triggerState.cycleSamples = 0;
			triggerState.sampleOffset = 0;
			// -1 - (state.effectiveWindowSize - (int)state.effectiveWindowSize)
			triggerState.sampleOffset = (state.effectiveWindowSize * 0.5 - (int)(state.effectiveWindowSize * 0.5)) - 1.5;
			return;
		}
		else if (state.triggerMode != OscilloscopeContent::TriggeringMode::Spectral)
		{
			triggerState.sampleOffset = 0;
			triggerState.cycleSamples = 0;

			return;
		}

#ifdef PHASE_VOCODER
		auto const TransformSize = OscilloscopeContent::LookaheadSize >> 1;
#else
		auto const TransformSize = OscilloscopeContent::LookaheadSize;
#endif

		Eval eval(params);

		if (!eval.isWellDefined())
			return;

		const auto tau = cpl::simd::consts<double>::tau;
		const auto radians = tau * (triggerState.record.omega()) / TransformSize;

		// we will try to analyse the points closest to the sync point (latest in time)
		auto offsetReal = std::max<double>(OscilloscopeContent::LookaheadSize, state.effectiveWindowSize + triggerState.cycleSamples);
		auto const offset = (std::size_t)std::ceil(offsetReal);

		auto sampleDifference = offset - (state.effectiveWindowSize + triggerState.cycleSamples);

		eval.startFrom(-static_cast<cpl::ssize_t>(offset));

		for (std::size_t i = 0; i < OscilloscopeContent::LookaheadSize; ++i)
		{
			temporaryBuffer[i] = eval.evaluateSampleInc();
		}

		// get the complex sinusoid phase
		auto z = cpl::dsp::goertzel(temporaryBuffer, OscilloscopeContent::LookaheadSize, radians);

		// rotate the sinusoid by the fractional difference in samples
		// a basic identity of the DFT, if the time domain is moved by k,
		// it is a complex rotation in the frequency domain by
		// -i*2*pi*k*n/N
		auto rotation = -sampleDifference * radians;

		z *= std::complex<double>(std::cos(rotation), -std::sin(rotation));

		// invert unit circle as the display is "backwards" in time
		auto phase = tau - std::arg(z);
		// correct phase by delta
		phase += triggerState.record.offset * tau;
		// phase correct to sines
		phase -= cpl::simd::consts<double>::pi_half;
		// add user-defined phase offset
		phase += tau * content->triggerPhaseOffset.getParameterView().getValueTransformed() / 360;
		// since we can't go back in time, travel around the unit circle until the phase is positive.
		phase = std::fmod(phase, tau);
		while (phase < 0)
			phase += tau;

		triggerState.phase = phase;

		// normalize phase to cycles
		auto cycles = phase / tau;

		// convert to samples
		triggerState.sampleOffset = cycles * audioStream.getAudioHistorySamplerate() / (triggerState.fundamental) - 1;

	}

	inline void Oscilloscope::resizeAudioStorage()
	{
		std::size_t requiredSampleBufferSize = 0;
		// TODO: Add
		// std::size_t additionalSamples = state.sampleInterpolation == SubSampleInterpolation::Lanczos ? OscilloscopeContent::InterpolationKernelSize : 0;

		if (state.triggerMode == OscilloscopeContent::TriggeringMode::Spectral)
		{
			// buffer size = length of detected freq in samples + display window size + lookahead
			requiredSampleBufferSize = std::max(static_cast<std::size_t>(0.5 + triggerState.cycleSamples + std::ceil(state.effectiveWindowSize)), OscilloscopeContent::LookaheadSize);
		}
		else
		{
			//requiredSampleBufferSize = static_cast<std::size_t>(0.5 + triggerState.cycleSamples + std::ceil(state.effectiveWindowSize) * 2) + OscilloscopeContent::LookaheadSize;
			requiredSampleBufferSize = static_cast<std::size_t>(std::ceil(state.effectiveWindowSize + 1));
		}
		channelData.front.resizeStorage(requiredSampleBufferSize, std::max(requiredSampleBufferSize, audioStream.getAudioHistoryCapacity()));
		channelData.back.resizeStorage(requiredSampleBufferSize, std::max(requiredSampleBufferSize, audioStream.getAudioHistoryCapacity()));
	}


	template<typename ISA, class Analyzer>
	void Oscilloscope::executeSamplingWindows(AFloat ** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		if (numChannels < 1)
			return;

		Analyzer ana(numChannels, numSamples, audioStream.getASyncPlayhead().getSteadyClock(), *triggerState.triggeringProcessor);

		auto mode = content->channelConfiguration.param.getAsTEnum<OscChannels>();
		auto triggeringChannel = std::min(numChannels, static_cast<std::size_t>(content->triggeringChannel.getTransformedValue())) - 1;

		if (numChannels == 1)
			mode = OscChannels::Left;

		if (numChannels == 1)
		{
			for (std::size_t n = 0; n < numSamples; ++n)
			{
				ana.process(buffer[0][n]);
			}
		}
		else
		{
			std::size_t offset = 0;
			switch (mode)
			{
			case OscChannels::Right: offset = 1;
			case OscChannels::Left:
			{
				const auto channel = buffer[offset];
				for (std::size_t n = 0; n < numSamples; ++n)
				{
					ana.process(channel[n]);
				}
				break;
			}
			case OscChannels::Side:
				for (std::size_t n = 0; n < numSamples; ++n)
				{
					ana.process(0.5 * (buffer[0][n] - buffer[1][n]));
				}
				break;
			case OscChannels::Separate:
				for (std::size_t n = 0; n < numSamples; ++n)
				{
					ana.process(buffer[triggeringChannel][n]);
				}

				break;
			case OscChannels::Mid: case OscChannels::MidSide:
				for (std::size_t n = 0; n < numSamples; ++n)
				{
					ana.process(0.5 * (buffer[0][n] + buffer[1][n]));
				}
				break;
			}
		}


	}

	template<typename ISA>
	void Oscilloscope::preAnalyseAudio(AFloat ** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		const bool isPeakHolder = state.triggerMode == OscilloscopeContent::TriggeringMode::EnvelopeHold;
		const bool isZeroCrossingPeakSampler = state.triggerMode == OscilloscopeContent::TriggeringMode::ZeroCrossing;

		if (isZeroCrossingPeakSampler)
		{
			executeSamplingWindows<ISA, ZeroCrossingProcessor<ISA>>(buffer, numChannels, numSamples);
		}
		else if (isPeakHolder)
		{
			executeSamplingWindows<ISA, PeakHoldProcessor<ISA>>(buffer, numChannels, numSamples);
		}
	}

	template<typename ISA>
	void Oscilloscope::audioEntryPoint(AFloat ** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		if (numSamples == 0 || numChannels == 0)
			return;

		cpl::CMutex scopedLock(bufferLock);

		channelData.resizeChannels(numChannels);

		// get channel names for legend display
		channelNames = audioStream.getAsyncChannelNames();

		cpl::variable_array<float*> localBuffers(buffer, buffer + numChannels);

		triggerState.triggeringProcessor->update(audioStream.getASyncPlayhead().getSteadyClock());
		preAnalyseAudio<ISA>(localBuffers.data(), numChannels, numSamples);

		if (state.triggerMode != OscilloscopeContent::TriggeringMode::EnvelopeHold && state.triggerMode != OscilloscopeContent::TriggeringMode::ZeroCrossing)
		{
			audioProcessing<ISA>(localBuffers.data(), numChannels, numSamples, channelData.front);
		}
		else
		{
			triggerState.triggeringProcessor->processMutating<ISA>(*this, localBuffers.data(), numChannels, numSamples);
		}
	}

	template<typename ISA>
		void Oscilloscope::audioProcessing(AFloat ** buffer, const std::size_t numChannels, const std::size_t numSamples, ChannelData::Buffer & target)
		{
			using namespace cpl::simd;
			typedef AFloat T;

			auto const sampleRate = audioStream.getAudioHistorySamplerate();

			channelData.tuneColourSmoothing(content->colourSmoothing.getTransformedValue(), sampleRate);
			channelData.tuneCrossOver(300, 3000, sampleRate);

			if (target.defaultChannel().audioData.getSize() < 1)
				return;

			auto colourSmoothPole = channelData.smoothFilterPole;

			// (division by zero is well-defined)
			const auto envelopeCoeff = std::exp(-1.0 / (content->envelopeWindow.getNormalizedValue() * audioStream.getAudioHistorySamplerate()));

			cpl::variable_array<AFloat> filterEnv(numChannels);

			for (std::size_t i = 0; i < numChannels; ++i)
				filterEnv[i] = channelData.filterStates.channels[i].envelope;

			auto colourArray = [](const auto & colourParam)
			{
				return std::array<float, 3> {static_cast<AFloat>(colourParam.r.getValue()), static_cast<AFloat>(colourParam.g.getValue()), static_cast<AFloat>(colourParam.b.getValue())};
			};

			auto filterStates = [=](const auto & bands, auto & states)
			{
				for (std::size_t i = 0; i < ChannelData::Bands; ++i)
				{
					auto const input = bands[i] * bands[i];
					states[i] = input + colourSmoothPole * (states[i] - input); // can also rectify
				}

			};

			decltype(colourArray(content->lowColour)) colours[] = { colourArray(content->lowColour), colourArray(content->midColour), colourArray(content->highColour) };

			auto accumulateColour = [&colours](const auto & state, ChannelData::PixelType key, float blend)
			{
				typedef ChannelData::PixelType::ComponentType C;
				ChannelData::PixelType ret;
				constexpr auto PixelMax = static_cast<AFloat>(std::numeric_limits<C>::max());
				AFloat red(0), green(0), blue(0);

				for (std::size_t i = 0; i < ChannelData::Bands; ++i)
				{
					red += state[i] * colours[i][0];
					green += state[i] * colours[i][1];
					blue += state[i] * colours[i][2];
				}

				auto invMax = PixelMax / std::max(red, std::max(blue, green));

				ret.pixel.a = std::numeric_limits<C>::max();
				ret.pixel.r = static_cast<C>(red * invMax);
				ret.pixel.g = static_cast<C>(green * invMax);
				ret.pixel.b = static_cast<C>(blue * invMax);

				return ret.lerp(key, blend);
			};

			using fs = ChannelData;

			auto mode = content->channelConfiguration.param.getAsTEnum<OscChannels>();

			if (numChannels == 1)
				mode = OscChannels::Left;

			const auto blend = 1 - content->frequencyColouringBlend.parameter.getValue();

			if (numChannels >= 2)
			{

				std::size_t offset = 0;

				// process envelopes so we're not thrashing the icache
				if (state.envelopeMode != EnvelopeModes::None)
				{
					switch (mode)
					{
					case OscChannels::Right: offset = 1;
					case OscChannels::Left:
					{
						const auto channel = buffer[offset];

						for (std::size_t n = 0; n < numSamples; ++n)
						{
							const auto sample = cpl::Math::square(channel[n]);
							filterEnv[fs::Left] = sample + envelopeCoeff * (filterEnv[fs::Left] - sample);
						}

						for (std::size_t c = 1; c < numChannels; ++c)
							filterEnv[c] = filterEnv[fs::Left];

						break;
					}

					case OscChannels::Mid:
						for (std::size_t n = 0; n < numSamples; ++n)
						{
							const auto mid = 0.5 * (buffer[fs::Left][n] + buffer[fs::Right][n]);
							const auto sample = cpl::Math::square(mid);
							filterEnv[fs::Left] = sample + envelopeCoeff * (filterEnv[fs::Left] - sample);
						}

						for (std::size_t c = 1; c < numChannels; ++c)
							filterEnv[c] = filterEnv[fs::Left];

						break;
					case OscChannels::Side:
						for (std::size_t n = 0; n < numSamples; ++n)
						{
							const auto side = 0.5 * (buffer[fs::Left][n] - buffer[fs::Right][n]);
							const auto sample = cpl::Math::square(side);
							filterEnv[fs::Left] = sample + envelopeCoeff * (filterEnv[fs::Left] - sample);
						}

						for (std::size_t c = 1; c < numChannels; ++c)
							filterEnv[c] = filterEnv[fs::Left];

						break;

					case OscChannels::Separate:
						for (std::size_t n = 0; n < numSamples; ++n)
						{
							for (std::size_t c = 0; c < numChannels; ++c)
							{
								const auto sample = cpl::Math::square(buffer[c][n]);
								filterEnv[c] = sample + envelopeCoeff * (filterEnv[c] - sample);
							} 
						}

						break;

					case OscChannels::MidSide:
						for (std::size_t n = 0; n < numSamples; ++n)
						{
							const auto left = buffer[fs::Left][n], right = buffer[fs::Right][n];
							const auto mid = 0.5 * cpl::Math::square(left + right), side = 0.5 * cpl::Math::square(left - right);

							filterEnv[fs::Left] = mid + envelopeCoeff * (filterEnv[fs::Left] - mid);
							filterEnv[fs::Right] = side + envelopeCoeff * (filterEnv[fs::Right] - side);
						}

						for (std::size_t c = 2; c < numChannels; ++c)
							filterEnv[c] = filterEnv[fs::Right];

						break;
					}

				}

				auto sideSignal = [](const auto & left, const auto & right)
				{
					unq_typeof(left) ret;
					for (std::size_t i = 0; i < ChannelData::Bands; ++i)
						ret[i] = left[i] - right[i];
					return ret;
				};

				auto midSignal = [](const auto & left, const auto & right)
				{
					unq_typeof(left) ret;
					for (std::size_t i = 0; i < ChannelData::Bands; ++i)
						ret[i] = left[i] + right[i];
					return ret;
				};

				typedef ChannelData::Crossover::BandArray SmoothState;
				typedef unq_typeof(target.channels[fs::Left].colourData.createWriter()) ColourWriter;

				cpl::variable_array<SmoothState> smoothStates(numChannels, [&](std::size_t i) { return channelData.filterStates.channels[i].smoothFilters; });
				cpl::variable_array<ColourWriter> colourWriters(numChannels, [&](std::size_t i) { return target.channels[i].colourData.createWriter(); });
				cpl::variable_array<ChannelData::Crossover*> networks(numChannels, [&](std::size_t i) { return &channelData.filterStates.channels[i].network; });
				cpl::variable_array<ChannelData::PixelType> colours(numChannels, [&](auto i) { return content->getColour(i).getAsJuceColour(); });

				auto
					midSmoothState = channelData.filterStates.midSideSmoothsFilters[0],
					sideSmoothState = channelData.filterStates.midSideSmoothsFilters[1];

				auto &&
					mw = target.midSideColour[0].createWriter(),
					sw = target.midSideColour[1].createWriter();


				for (std::size_t n = 0; n < numSamples; ++n)
				{

					const auto 
						left = buffer[fs::Left][n],
						right = buffer[fs::Right][n];

					// split signal into bands:
					auto leftBands = networks[0]->process(left, channelData.networkCoeffs);
					auto rightBands = networks[1]->process(right, channelData.networkCoeffs);

					filterStates(leftBands, smoothStates[0]);
					filterStates(rightBands, smoothStates[1]);

					// magnitude doesn't matter for these, as we normalize the data anyway -
					filterStates(midSignal(leftBands, rightBands), midSmoothState);
					filterStates(sideSignal(leftBands, rightBands), sideSmoothState);

					for (std::size_t c = 2; c < numChannels; ++c)
					{
						auto bands = networks[c]->process(buffer[c][n], channelData.networkCoeffs);
						filterStates(rightBands, smoothStates[c]);
					}

					mw.setHeadAndAdvance(accumulateColour(midSmoothState, colours[0], blend));
					sw.setHeadAndAdvance(accumulateColour(sideSmoothState, colours[1], blend));

					for (std::size_t c = 0; c < numChannels; ++c)
					{
						colourWriters[c].setHeadAndAdvance(accumulateColour(smoothStates[c], colours[c], blend));
					}

				}

				for (std::size_t c = 0; c < numChannels; ++c)
				{
					channelData.filterStates.channels[c].smoothFilters = smoothStates[c];
				}


				channelData.filterStates.midSideSmoothsFilters[0] = midSmoothState;
				channelData.filterStates.midSideSmoothsFilters[1] = sideSmoothState;

			}
			else if (numChannels == 1)
			{
				ChannelData::PixelType colour(content->getColour(0).getAsJuceColour());

				filterEnv[1] = 0;
				auto && lw = target.channels[fs::Left].colourData.createWriter();
				auto leftSmoothState = channelData.filterStates.channels[fs::Left].smoothFilters;

				for (std::size_t n = 0; n < numSamples; n++)
				{
					const auto left = buffer[fs::Left][n];
					const auto lSquared = left * left;

					// average envelope
					filterEnv[fs::Left] = lSquared + envelopeCoeff * (filterEnv[fs::Left] - lSquared);
					// get peak sample

					// split signal into bands:
					auto leftBands = channelData.filterStates.channels[fs::Left].network.process(left, channelData.networkCoeffs);

					filterStates(leftBands, leftSmoothState);
					lw.setHeadAndAdvance(accumulateColour(leftSmoothState, colour, blend));
				}

				channelData.filterStates.channels[fs::Left].smoothFilters = leftSmoothState;

			}

			// store calculated envelope
			if (state.envelopeMode == EnvelopeModes::RMS)
			{
				// we end up calculating envelopes even though its not possibly needed, but the overhead
				// is negligible

				auto start = std::sqrt(filterEnv[0]);

				for (std::size_t c = 0; c < numChannels; ++c)
					start = std::max(start, std::sqrt(filterEnv[c]));

				shared.autoGainEnvelope.store(1.0 / start, std::memory_order_release);

				// only update filters if this mode is on.
				smoothEnvelopeState(ChannelData::Left) = filterEnv[0];
				smoothEnvelopeState(ChannelData::Right) = filterEnv[1];

			}

			// save audio data
			for(std::size_t c = 0; c < target.channels.size(); ++c)
				target.channels[c].audioData.createWriter().copyIntoHead(buffer[c], numSamples);

			state.transportPosition = audioStream.getASyncPlayhead().getPositionInSamples() + numSamples;

		}

	template<typename ISA>
		void Oscilloscope::runPeakFilter()
		{
			// there is a number of optimisations we can do here, mostly that we actually don't care about
			// timing, we are only interested in the current largest value in the set.
			using namespace cpl;
			using namespace cpl::simd;
			using cpl::simd::load;

			typedef typename ISA::V V;

			V
				vLMax = zero<V>(),
				vRMax = zero<V>(),
				vSign = consts<V>::sign_mask;

			auto const loopIncrement = elements_of<V>::value;

			auto const vHalf = consts<V>::half;

			auto const numChannels = channelData.front.channels.size();

			auto const channelMode = numChannels == 1 ? OscChannels::Left : state.channelMode;

			auto const numSamples = channelData.front.channels.begin()->audioData.getSize();

			auto stop = numSamples - (numSamples & (loopIncrement - 1));
			if (stop <= 0)
				stop = 0;

			// since this runs in every frame, we need to scale the coefficient by how often this function runs
			// (and the amount of samples)
			double power = numSamples * (avgFps.getAverage() / juce::Time::getHighResolutionTicksPerSecond());
			double coeff = std::pow(std::exp(-static_cast<double>(loopIncrement) / (content->envelopeWindow.getNormalizedValue() * audioStream.getAudioHistorySamplerate())), power);

			if (channelMode <= OscChannels::OffsetForMono)
			{
				auto && view = channelData.front.channels[0].audioData.createProxyView();

				const auto * leftBuffer = view.begin();

				if (channelMode == OscChannels::Left)
				{
					for (std::size_t i = 0; i < stop; i += loopIncrement)
					{
						auto const vLInput = load<V>(leftBuffer + i);
						vLMax = max(vand(vLInput, vSign), vLMax);
					}
				}
				else
				{
					auto && rightView = channelData.front.channels[1].audioData.createProxyView();
					const auto * rightBuffer = rightView.begin();

					switch (channelMode)
					{
					case OscChannels::Right:

						for (std::size_t i = 0; i < stop; i += loopIncrement)
						{
							auto const vLInput = load<V>(rightBuffer + i);
							vLMax = max(vand(vLInput, vSign), vLMax);
						}
						break;
					case OscChannels::Mid:
						for (std::size_t i = 0; i < stop; i += loopIncrement)
						{
							auto const vInput = load<V>(leftBuffer + i) + load<V>(rightBuffer + i);
							vLMax = max(vand(vInput * vHalf, vSign), vLMax);
						}
						break;
					case OscChannels::Side:

						for (std::size_t i = 0; i < stop; i += loopIncrement)
						{
							auto const vInput = load<V>(leftBuffer + i) - load<V>(rightBuffer + i);
							vLMax = max(vand(vInput * vHalf, vSign), vLMax);
						}
						break;
					}
				}

				vRMax = vLMax;

				// TODO: remainder?

				suitable_container<V> lmax = vLMax, rmax = vRMax;

				double highestLeft = *std::max_element(lmax.begin(), lmax.end());
				double highestRight = *std::max_element(rmax.begin(), rmax.end());

				smoothEnvelopeState(0) = std::max(smoothEnvelopeState(0) * coeff, highestLeft  * highestLeft);
				smoothEnvelopeState(1) = std::max(smoothEnvelopeState(1) * coeff, highestRight * highestRight);

				for (std::size_t c = 2; c < numChannels; ++c)
				{
					smoothEnvelopeState(c) = smoothEnvelopeState(1);
				}

			}
			else
			{
				switch (channelMode)
				{
					case OscChannels::Separate:
					{
	
						for (std::size_t c = 0; c < numChannels; ++c)
						{
							auto&& view = channelData.front.channels[c].audioData.createProxyView();
							auto buffer = view.begin();

							for (std::size_t i = 0; i < stop; i += loopIncrement)
							{
								auto const vInput = load<V>(buffer + i);
								vLMax = max(vand(vInput, vSign), vLMax);
							}

							suitable_container<V> lmax = vLMax;
							auto highestValue = *std::max_element(lmax.begin(), lmax.end());

							smoothEnvelopeState(c) = std::max<float>(smoothEnvelopeState(c) * coeff, highestValue * highestValue);
						}

						break;
					}

					case OscChannels::MidSide:
					{
						ChannelData::AudioBuffer::ProxyView views[2] = { channelData.front.channels[0].audioData.createProxyView(), channelData.front.channels[1].audioData.createProxyView() };

						const auto * leftBuffer = views[0].begin();
						const auto * rightBuffer = views[1].begin();

						for (std::size_t i = 0; i < stop; i += loopIncrement)
						{
							auto const vLInput = load<V>(leftBuffer + i);
							auto const vRInput = load<V>(rightBuffer + i);

							auto const a = vLInput + vRInput;
							auto const b = vLInput - vRInput;

							vLMax = max(vand(a * vHalf, vSign), vLMax);
							vRMax = max(vand(b * vHalf, vSign), vRMax);
						}

						suitable_container<V> lmax = vLMax, rmax = vRMax;

						double highestLeft = *std::max_element(lmax.begin(), lmax.end());
						double highestRight = *std::max_element(rmax.begin(), rmax.end());

						smoothEnvelopeState(ChannelData::Left) = std::max(smoothEnvelopeState(0) * coeff, highestLeft  * highestLeft);
						smoothEnvelopeState(ChannelData::Right) = std::max(smoothEnvelopeState(1) * coeff, highestRight * highestRight);

						for (std::size_t c = 2; c < numChannels; ++c)
						{
							smoothEnvelopeState(c) = smoothEnvelopeState(1);
						}

						break;
					}

				}

			}


			auto start = std::sqrt(smoothEnvelopeState(0));

			for (std::size_t c = 0; c < numChannels; ++c)
			{
				start = std::max(start, std::sqrt(smoothEnvelopeState(c)));
			}

			shared.autoGainEnvelope.store(1.0 / start, std::memory_order_release);
		}
};
#endif