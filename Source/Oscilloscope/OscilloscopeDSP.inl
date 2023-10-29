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

	file: OscilloscopeDSP.inl

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
	void Oscilloscope::analyseAndSetupState(const EvaluatorParams& params, Oscilloscope::StreamState& cs)
	{
		calculateFundamentalPeriod<ISA, Eval>(params);
		calculateTriggeringOffset<ISA, Eval>(params);

		cs.channelData.resizeAudioStorage(cs.triggerMode, state.effectiveWindowSize, triggerState.cycleSamples, cs.historyCapacity);

		if (cs.envelopeMode == EnvelopeModes::PeakDecay)
		{
			runPeakFilter<ISA>(cs.channelData);
		}
		else if (cs.envelopeMode == EnvelopeModes::None)
		{
			state.autoGain = 1;
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
			auto const normalizedFrequency = state.customTriggerFrequency / state.sampleRate;
			triggerState.record = BinRecord{ 0, 1, normalizedFrequency * TransformSize};

			auto fundamental = state.customTriggerFrequency;

			triggerState.fundamental = state.customTriggerFrequency;
			triggerState.cycleSamples = state.sampleRate / fundamental;
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

			triggerFFT.forward(transformBuffer, transformBuffer, triggerWork);

#ifdef PHASE_VOCODER
			triggerFFT.forward(cpl::as_uarray(transformBuffer).slice(TransformSize, TransformSize), cpl::as_uarray(transformBuffer).slice(0, TransformSize), triggerWork);
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

				const auto denom = x0 * ChannelFloat(2) - xm1 - x1;

				return (denom.real() + denom.imag()) != 0 ? std::real((xm1 - x1) / denom) : 0;
#endif
			};

			const double quarterSemitone = std::pow(2, 0.25 / 12.0) - 1;
			const double threshold = content->triggerThreshold.getTransformedValue();
			const double hysteresis = content->triggerHysteresis.getTransformedValue();
			const auto invHysteresis = 1 - hysteresis;

			// Reduce to 1/4 (to slip through "vastly better case" + half transform size
			BinRecord max{ 1, std::max<double>(threshold * TransformSize / 6.0, std::abs(transformBuffer[1])), quadDelta(1) };

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

			auto fundamental = state.sampleRate * (max.omega()) / TransformSize;

			triggerState.fundamental = fundamental = std::max(5.0, fundamental);
			triggerState.cycleSamples = state.sampleRate / fundamental;
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
			transformBuffer[i] = eval.evaluateSampleInc();
		}

		// get the complex sinusoid phase
		auto z = cpl::dsp::goertzel(transformBuffer, OscilloscopeContent::LookaheadSize, radians);

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
		triggerState.sampleOffset = cycles * state.sampleRate / (triggerState.fundamental) - 1;

	}

	template<typename ISA, class Analyzer>
	void Oscilloscope::StreamState::executeSamplingWindows(AudioStream::ListenerContext& ctx, AFloat ** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		if (numChannels < 2)
			return;

		Analyzer ana(numChannels, numSamples, ctx.getPlayhead().getSteadyClock(), *triggeringProcessor);

		CPL_RUNTIME_ASSERTION(numChannels % 2 == 0);

		auto localMode = channelMode;

		std::size_t triggerSeparate, triggerPair;
		content->calculateTriggerIndices(numChannels, triggerSeparate, triggerPair);

		if (numChannels == 1)
			localMode = OscChannels::Left;

		if (numChannels == 1)
		{
			// currently unreachable
			for (std::size_t n = 0; n < numSamples; ++n)
			{
				ana.process(buffer[0][n]);
			}
		}
		else
		{
			if (localMode == OscChannels::MidSide)
			{
				// Mid/side mode user can choose either from a pair, depending on the odd bit.
				if (triggerSeparate & 0x1)
				{
					localMode = OscChannels::Mid;
				}
				else
				{
					localMode = OscChannels::Mid;
				}

				triggerPair = triggerSeparate & ~0x1;
			}

			switch (localMode)
			{
			case OscChannels::Right: triggerPair++;
			case OscChannels::Left:
			{
				const auto channel = buffer[triggerPair];
				for (std::size_t n = 0; n < numSamples; ++n)
				{
					ana.process(channel[n]);
				}
				break;
			}
			case OscChannels::Separate:
				for (std::size_t n = 0; n < numSamples; ++n)
				{
					ana.process(buffer[triggerSeparate][n]);
				}
				break;
			case OscChannels::Mid:
				for (std::size_t n = 0; n < numSamples; ++n)
				{
					ana.process(0.5f * (buffer[triggerPair + 0][n] + buffer[triggerPair + 1][n]));
				}
				break;
			case OscChannels::Side:
				for (std::size_t n = 0; n < numSamples; ++n)
				{
					ana.process(0.5f * (buffer[triggerPair + 0][n] - buffer[triggerPair + 1][n]));
				}
				break;
			}
		}
	}

	template<typename ISA>
	void Oscilloscope::StreamState::preAnalyseAudio(AudioStream::ListenerContext& ctx, AFloat ** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		const bool isPeakHolder = triggerMode == OscilloscopeContent::TriggeringMode::EnvelopeHold;
		const bool isZeroCrossingPeakSampler = triggerMode == OscilloscopeContent::TriggeringMode::ZeroCrossing;

		if (isZeroCrossingPeakSampler)
		{
			executeSamplingWindows<ISA, ZeroCrossingProcessor<ISA>>(ctx, buffer, numChannels, numSamples);
		}
		else if (isPeakHolder)
		{
			executeSamplingWindows<ISA, PeakHoldProcessor<ISA>>(ctx, buffer, numChannels, numSamples);
		}
	}

	template<typename ISA>
	void Oscilloscope::StreamState::audioEntryPoint(AudioStream::ListenerContext& ctx, AudioStream::DataType** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		if (numSamples == 0 || numChannels == 0)
			return;

		channelData.resizeChannels(numChannels);

		cpl::variable_array<float*> localBuffers(buffer, buffer + numChannels);

		triggeringProcessor->update(ctx.getPlayhead().getSteadyClock());
		preAnalyseAudio<ISA>(ctx, localBuffers.data(), numChannels, numSamples);

		if (triggerMode != OscilloscopeContent::TriggeringMode::EnvelopeHold && triggerMode != OscilloscopeContent::TriggeringMode::ZeroCrossing)
		{
			audioProcessing<ISA>(ctx.getInfo(), ctx.getPlayhead(), localBuffers.data(), numChannels, numSamples, channelData.front);
		}
		else
		{
			triggeringProcessor->processMutating<ISA>(*this, ctx, localBuffers.data(), numChannels, numSamples);
		}
	}

	template<typename ISA>
		void Oscilloscope::StreamState::audioProcessing(
			const AudioStream::Info& info, 
			const AudioStream::Playhead& playhead, 
			const AudioStream::DataType* const* buffer,
			const std::size_t numChannels, 
			const std::size_t numSamples, 
			ChannelData::Buffer& target
		)
		{
			using namespace cpl::simd;
			typedef AFloat T;

			channelData.tuneColourSmoothing(content->colourSmoothing.getTransformedValue(), info.sampleRate);
			channelData.tuneCrossOver(300, 3000, info.sampleRate);

			if (target.defaultChannel().audioData.getSize() < 1)
				return;

			auto colourSmoothPole = channelData.smoothFilterPole;

			// (division by zero is well-defined)
			const auto envelopeCoeff = static_cast<Signalizer::AFloat>(std::exp(-1.0 / (content->envelopeWindow.getNormalizedValue() * info.sampleRate)));

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

			auto mode = channelMode;

			if (numChannels == 1)
				mode = OscChannels::Left;

			const auto blend = 1 - static_cast<Signalizer::AFloat>(content->frequencyColouringBlend.parameter.getValue());

			if (numChannels >= 2)
			{

				std::size_t offset = 0;

				// process envelopes so we're not thrashing the icache
				if (envelopeMode != EnvelopeModes::None)
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
							const auto mid = 0.5f * (buffer[fs::Left][n] + buffer[fs::Right][n]);
							const auto sample = cpl::Math::square(mid);
							filterEnv[fs::Left] = sample + envelopeCoeff * (filterEnv[fs::Left] - sample);
						}

						for (std::size_t c = 1; c < numChannels; ++c)
							filterEnv[c] = filterEnv[fs::Left];

						break;
					case OscChannels::Side:
						for (std::size_t n = 0; n < numSamples; ++n)
						{
							const auto side = 0.5f * (buffer[fs::Left][n] - buffer[fs::Right][n]);
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
							const auto 
								mid = 0.5f * cpl::Math::square(left + right), 
								side = 0.5f * cpl::Math::square(left - right);

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

				CPL_RUNTIME_ASSERTION((numChannels % 2) == 0);

				// Process colours for each channel pair
				for (std::size_t channelPair = 0; channelPair < numChannels; channelPair += 2)
				{
					auto& leftMid = channelData.filterStates.channels[channelPair + fs::Left];
					auto& rightSide = channelData.filterStates.channels[channelPair + fs::Right];

					auto 
						cwLeft = target.channels[channelPair + fs::Left].colourData.createWriter(),
						cwRight = target.channels[channelPair + fs::Right].colourData.createWriter(),
						cwMid = target.channels[channelPair + fs::Mid].auxColourData.createWriter(),
						cwSide = target.channels[channelPair + fs::Side].auxColourData.createWriter();

					auto& netLeft = leftMid.network;
					auto& netRight = rightSide.network;

					auto
						&smLeft = leftMid.smoothFilters,
						&smRight = rightSide.smoothFilters,
						&smMid = leftMid.auxSmoothFilter,
						&smSide = rightSide.auxSmoothFilter;

					const ChannelData::PixelType leftColour = leftMid.defaultKey;
					const ChannelData::PixelType rightColour = rightSide.defaultKey;

					for (std::size_t n = 0; n < numSamples; ++n)
					{
						auto leftBands = netLeft.process(buffer[channelPair + fs::Left][n], channelData.networkCoeffs);
						auto rightBands = netRight.process(buffer[channelPair + fs::Right][n], channelData.networkCoeffs);

						filterStates(leftBands, smLeft);
						filterStates(rightBands, smRight);

						// magnitude doesn't matter for these, as we normalize the data anyway -
						filterStates(midSignal(leftBands, rightBands), smMid);
						filterStates(sideSignal(leftBands, rightBands), smSide);

						cwLeft.setHeadAndAdvance(accumulateColour(smLeft, leftColour, blend));
						cwRight.setHeadAndAdvance(accumulateColour(smRight, rightColour, blend));

						cwMid.setHeadAndAdvance(accumulateColour(smMid, leftColour, blend));
						cwSide.setHeadAndAdvance(accumulateColour(smSide, rightColour, blend));

					}
				}

			}
			else if (numChannels == 1)
			{
				ChannelData::PixelType colour(channelData.filterStates.channels[0].defaultKey);

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
			if (envelopeMode == EnvelopeModes::RMS)
			{
				// we end up calculating envelopes even though its not possibly needed, but the overhead
				// is negligible

				auto start = std::sqrt(filterEnv[0]);

				for (std::size_t c = 0; c < numChannels; ++c)
					start = std::max(start, std::sqrt(filterEnv[c]));

				envelopeGain = 1.0 / start;

				// only update filters if this mode is on.
				channelData.filterStates.channels[ChannelData::Left].envelope = filterEnv[0];
				channelData.filterStates.channels[ChannelData::Right].envelope = filterEnv[1];

			}

			// save audio data
			for(std::size_t c = 0; c < target.channels.size(); ++c)
				target.channels[c].audioData.createWriter().copyIntoHead(buffer[c], numSamples);

			transportPosition = playhead.getPositionInSamples() + numSamples;
			bpm = playhead.getBPM();
			sampleRate = info.sampleRate;

		}

	template<typename ISA>
		void Oscilloscope::runPeakFilter(ChannelData& data)
		{
			// there is a number of optimisations we can do here, mostly that we actually don't care about
			// timing, we are only interested in the current largest value in the set.
			using namespace cpl;
			using namespace cpl::simd;
			using cpl::simd::load;

			typedef typename ISA::V V;

			auto smoothEnvelopeState = [&](std::size_t i) -> AFloat& { return data.filterStates.channels[i].envelope; };

			V
				vLMax = zero<V>(),
				vRMax = zero<V>(),
				vSign = consts<V>::sign_mask;

			auto const loopIncrement = elements_of<V>::value;

			auto const vHalf = consts<V>::half;

			auto const numChannels = data.front.channels.size();

			auto const channelMode = numChannels == 1 ? OscChannels::Left : shared.channelMode.load();

			auto const numSamples = data.front.channels.begin()->audioData.getSize();

			auto stop = numSamples - (numSamples & (loopIncrement - 1));
			if (stop <= 0)
				stop = 0;

			// since this runs in every frame, we need to scale the coefficient by how often this function runs
			// (and the amount of samples)
			double power = numSamples * openGLDeltaTime();
			double coeff = std::pow(std::exp(-static_cast<double>(loopIncrement) / (content->envelopeWindow.getNormalizedValue() * state.sampleRate)), power);

			if (channelMode <= OscChannels::OffsetForMono)
			{
				auto && view = data.front.channels[0].audioData.createProxyView();

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
					auto && rightView = data.front.channels[1].audioData.createProxyView();
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
							auto&& view = data.front.channels[c].audioData.createProxyView();
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
						ChannelData::AudioBuffer::ProxyView views[2] = { data.front.channels[0].audioData.createProxyView(), data.front.channels[1].audioData.createProxyView() };

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

			state.autoGain = 1.0 / start;
		}
};
#endif
