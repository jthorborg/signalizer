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


#include "COscilloscope.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include <cpl/simd.h>
#include <cpl/CDBMeterGraph.h>

namespace Signalizer
{

	void COscilloscope::calculateFundamentalPeriod()
	{
		switch (state.triggerMode)
		{
			case OscilloscopeContent::TriggeringMode::Spectral:
			{

#ifdef PHASE_VOCODER
				auto const TransformSize = OscilloscopeContent::LookaheadSize >> 1;
#else
				auto const TransformSize = OscilloscopeContent::LookaheadSize;
#endif

				const auto && view = lifoStream.createProxyView();

				cpl::ssize_t
					cursor = view.cursorPosition(),
					totalBufferSamples = view.size();
				
				if (totalBufferSamples < OscilloscopeContent::LookaheadSize)
					return;

				// we will try to analyse the points closest to the sync point (latest in time)
				auto offset = std::max<std::size_t>(std::ceil(state.effectiveWindowSize), OscilloscopeContent::LookaheadSize);

				auto pointer = view.begin() + (cursor - offset);
				while (pointer < view.begin())
					pointer += totalBufferSamples;

				auto end = view.end();

				for (cpl::ssize_t i = 0; i < OscilloscopeContent::LookaheadSize; ++i)
				{
					transformBuffer[i] = *pointer++;

					if (pointer == end)
						pointer -= totalBufferSamples;
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

				BinRecord max{ 1, std::abs(transformBuffer[1]), quadDelta(1) };

				for (std::size_t i = 2; i < (TransformSize >> 1); ++i)
				{
					BinRecord current{ i, std::abs(transformBuffer[i]) };

					// candidate must be vastly better
					if (current.value > max.value * 2)
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
							if (sensivity > 20)
							{
								max = current;
								continue;
							}

							// the same value, just a better estimate, from another bin
							// TODO: fix this case by polynomially interpolate the value as well
							if (std::abs(1 - factor) < quarterSemitone)
							{
								max = current;
								continue;
							}

							auto multipleDeviation = std::abs(factor - std::floor(factor + 0.5));

							// check if the harmonic series is more than half a semi-tone away, in which case we take the candidate
							if (std::abs(multipleDeviation) > quarterSemitone)
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
	}

	void COscilloscope::calculateTriggeringOffset()
	{
		if (state.triggerMode == OscilloscopeContent::TriggeringMode::None)
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

		auto && view = lifoStream.createProxyView();
		cpl::ssize_t
			cursor = view.cursorPosition(),
			totalBufferSamples = view.size();

		if (totalBufferSamples < OscilloscopeContent::LookaheadSize)
			return;

		const auto tau = cpl::simd::consts<double>::tau;
		const auto radians = tau * (triggerState.record.omega()) / TransformSize;

		// we will try to analyse the points closest to the sync point (latest in time)
		auto offsetReal = std::max<double>(OscilloscopeContent::LookaheadSize, state.effectiveWindowSize + triggerState.cycleSamples);
		auto const offset = (std::size_t)std::ceil(offsetReal);

		auto sampleDifference = offset - (state.effectiveWindowSize + triggerState.cycleSamples);

		auto pointer = view.begin() + (cursor - offset);
		while (pointer < view.begin())
			pointer += totalBufferSamples;

		auto const end = view.end();

		for (cpl::ssize_t i = 0; i < OscilloscopeContent::LookaheadSize; ++i)
		{
			auto sample = *pointer;

			temporaryBuffer[i] = sample;
			pointer++;
			if (pointer == end)
				pointer -= totalBufferSamples;
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
		triggerState.sampleOffset = cycles * audioStream.getAudioHistorySamplerate() / (triggerState.fundamental);

	}

	void COscilloscope::resizeAudioStorage()
	{
		std::size_t requiredSampleBufferSize = 0;
		if (state.timeMode != OscilloscopeContent::TimeMode::Beats)
		{
			// buffer size = length of detected freq in samples + display window size + lookahead
			requiredSampleBufferSize = static_cast<std::size_t>(0.5 + triggerState.cycleSamples + std::ceil(state.effectiveWindowSize)) + OscilloscopeContent::LookaheadSize;
		}
		else
		{
			requiredSampleBufferSize = static_cast<std::size_t>(0.5 + triggerState.cycleSamples + std::ceil(state.effectiveWindowSize) * 2) + OscilloscopeContent::LookaheadSize;
		}
		lifoStream.setStorageRequirements(requiredSampleBufferSize, std::max(requiredSampleBufferSize, audioStream.getAudioHistoryCapacity() + OscilloscopeContent::LookaheadSize));
	}

	bool COscilloscope::onAsyncAudio(const AudioStream & source, AudioStream::DataType ** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		switch (cpl::simd::max_vector_capacity<float>())
		{
		case 32:
		case 16:
		case 8:
#ifdef CPL_COMPILER_SUPPORTS_AVX
			audioProcessing<cpl::Types::v8sf>(buffer, numChannels, numSamples);
			break;
#endif
		case 4:
			audioProcessing<cpl::Types::v4sf>(buffer, numChannels, numSamples);
			break;
		default:
			audioProcessing<float>(buffer, numChannels, numSamples);
			break;
		}
		return false;
	}
	
	template<typename V>
		void COscilloscope::audioProcessing(typename cpl::simd::scalar_of<V>::type ** buffer, std::size_t numChannels, std::size_t numSamples)
		{
			using namespace cpl::simd;
			typedef typename scalar_of<V>::type T;
			if (numChannels != 2)
				return;

			cpl::CMutex scopedLock(bufferLock);

			T filterEnv[2] = { filters.envelope[0], filters.envelope[1] };

			for (std::size_t n = 0; n < numSamples; n++)
			{
				using fs = FilterStates;
				// collect squared inputs (really just a cheap abs)
				const auto lSquared = buffer[fs::Left][n] * buffer[fs::Left][n];
				const auto rSquared = buffer[fs::Right][n] * buffer[fs::Right][n];

				// average envelope
				filterEnv[fs::Left] = lSquared + state.envelopeCoeff * (filterEnv[fs::Left] - lSquared);
				filterEnv[fs::Right] = rSquared + state.envelopeCoeff * (filterEnv[fs::Right] - rSquared);



			}
			// store calculated envelope
			if (state.envelopeMode == EnvelopeModes::RMS && state.normalizeGain)
			{
				// we end up calculating envelopes even though its not possibly needed, but the overhead
				// is negligible
				double currentEnvelope = 1.0 / (2 * std::max(std::sqrt(filterEnv[0]), std::sqrt(filterEnv[1])));

				// only update filters if this mode is on.
				filters.envelope[0] = filterEnv[0];
				filters.envelope[1] = filterEnv[1];

				if (std::isnormal(currentEnvelope))
				{
					content->inputGain.getParameterView().updateFromProcessorTransformed(
						currentEnvelope,
						cpl::Parameters::UpdateFlags::All & ~cpl::Parameters::UpdateFlags::RealTimeSubSystem
					);
				}
			}
			// save data
			lifoStream.createWriter().copyIntoHead(buffer[0], numSamples);
			state.transportPosition = audioStream.getASyncPlayhead().getPositionInSamples() + numSamples;
		}

};
