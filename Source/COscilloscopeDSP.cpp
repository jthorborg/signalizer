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

	double COscilloscope::getTriggeringOffset()
	{
		switch (cpl::enum_cast<OscilloscopeContent::TriggeringMode>(content->triggerMode.param.getTransformedValue()))
		{
			case OscilloscopeContent::TriggeringMode::Spectral:
			{
#ifdef PHASE_VOCODER
				auto const TransformSize = OscilloscopeContent::LookaheadSize >> 1;
#else
				auto const TransformSize = OscilloscopeContent::LookaheadSize;
#endif
				transformBuffer.resize(OscilloscopeContent::LookaheadSize);
				temporaryBuffer.resize(OscilloscopeContent::LookaheadSize);

				auto && view = lifoStream.createProxyView();

				cpl::ssize_t
					cursor = view.cursorPosition(),
					size = view.size();

				if (size < OscilloscopeContent::LookaheadSize)
					return 0;
				// we will try to analyse the points closest to the sync point (latest in time)
				auto offset = std::max<std::size_t>(std::ceil(state.effectiveWindowSize), OscilloscopeContent::LookaheadSize);
				auto sampleDifference = std::max<double>(state.effectiveWindowSize, OscilloscopeContent::LookaheadSize) - offset;


				//offset = OscilloscopeContent::LookaheadSize;
				auto pointer = view.begin() + (cursor - offset);
				while (pointer < view.begin())
					pointer += size;

				auto end = view.end();

				for (cpl::ssize_t i = 0; i < OscilloscopeContent::LookaheadSize; ++i)
				{
					auto sample = *pointer;

					temporaryBuffer[i] = sample;
					transformBuffer[i] = sample;
					pointer++;
					if (pointer == end)
						pointer -= size;
				}



				signaldust::DustFFT_fwdDa(reinterpret_cast<double*>(transformBuffer.data()), TransformSize);
#ifdef PHASE_VOCODER
				signaldust::DustFFT_fwdDa(reinterpret_cast<double*>(transformBuffer.data() + TransformSize), TransformSize);
#endif

				auto quadDelta = [&](auto w) {
#if PHASE_VOCODER
					auto delta = std::arg(transformBuffer[w]) - std::arg(transformBuffer[TransformSize + w]);

					delta -= w * 2 * M_PI;

					while (delta < 2 * M_PI)
						delta += 2 * M_PI;

					while (delta > 2 * M_PI)
						delta -= 2 * M_PI;


					return delta / (2 * M_PI);


#else
					auto x0 = transformBuffer[w], x1 = transformBuffer[w + 1], xm1 = transformBuffer[w == 0 ? 1 : w - 1];
					return std::real((xm1 - x1) / ((x0 * 2.0) - xm1 - x1));
#endif
				};

				std::size_t maxIndex = 1;
				double maxValue = std::abs(transformBuffer[1]);
				double delta = quadDelta(1);

				const double halfSemitone = std::pow(2, 0.25 / 12.0) - 1;

				for (std::size_t i = 2; i < (TransformSize >> 1); ++i)
				{
					auto const newValue = std::abs(transformBuffer[i]);
					// candidate must be vastly better
					if (newValue > maxValue * 2)
					{
						auto oldOmega = (maxIndex + delta);

						// weird parabolas
						if (oldOmega > 0)
						{
							// check if it is somewhat harmonically related, in which case we discard the candidate
							auto const newDelta = quadDelta(i);
							auto newOmega = (i + newDelta);

							auto factor = newOmega / oldOmega;

							auto sensivity = newValue / maxValue;

							// shortcut if the value is 20 times bigger
							if (sensivity > 20)
							{
								maxValue = newValue;
								maxIndex = i;
								delta = newDelta;
								continue;
							}

							// the same value, just a better estimate, from another bin
							// TODO: fix this case by polynomially interpolate the value as well
							if (std::abs(1 - factor) < halfSemitone)
							{
								maxValue = newValue;
								maxIndex = i;
								delta = newDelta;
								continue;
							}

							auto multipleDeviation = std::abs(factor - std::floor(factor + 0.5));

							// check if the harmonic series is more than half a semi-tone away, in which case we take the candidate
							if (std::abs(multipleDeviation) > halfSemitone)
							{
								maxValue = newValue;
								maxIndex = i;
								delta = newDelta;
							}
						}
						else
						{
							maxValue = newValue;
							maxIndex = i;
							delta = quadDelta(i);
						}

					}
				}

				// copy old filter
				std::array<MedianData, MedianFilterSize> localMedian = medianTriggerFilter;

				// store new data
				medianTriggerFilter[medianPos].bin = maxIndex;
				medianTriggerFilter[medianPos].delta = delta;

				medianPos++;
				medianPos %= medianTriggerFilter.size();

				auto middle = (MedianFilterSize >> 1);
				std::nth_element(localMedian.begin(), localMedian.begin() + middle, localMedian.end(), [](const auto & a, const auto & b) { return a.bin < b.bin; });
				auto oldMedianBin = localMedian[middle];

				// check to discard (temporarily) much higher frequencies through a median filter
				if (std::abs(maxIndex + delta - (oldMedianBin.bin + oldMedianBin.delta)) > 0.5)
				{
					maxIndex = oldMedianBin.bin;
					delta = oldMedianBin.delta;
				}

				quantizedFreq = audioStream.getAudioHistorySamplerate() * double(maxIndex) / TransformSize;

				// interpolate peak position quadratically
				delta = quadDelta(maxIndex);
				detectedFreq = audioStream.getAudioHistorySamplerate() * (maxIndex + delta) / TransformSize;

				detectedFreq = std::max(1.0, detectedFreq);
				cycleLength = audioStream.getAudioHistorySamplerate() / detectedFreq;

				// we will try to analyse the points closest to the sync point (latest in time)
				auto offsetReal = std::max<double>(OscilloscopeContent::LookaheadSize, state.effectiveWindowSize + cycleLength);
				offset = (std::size_t)std::ceil(offsetReal);

				sampleDifference = offset - (state.effectiveWindowSize + cycleLength);

				auto z2 = cpl::dsp::goertzel(temporaryBuffer, temporaryBuffer.size(), cpl::simd::consts<double>::tau * (maxIndex + delta) / TransformSize);

				pointer = view.begin() + (cursor - offset);
				while (pointer < view.begin())
					pointer += size;

				end = view.end();

				auto count2 = OscilloscopeContent::LookaheadSize;
				//count2 = OscilloscopeContent::LookaheadSize;
				for (cpl::ssize_t i = 0; i < count2; ++i)
				{
					auto sample = *pointer;

					temporaryBuffer[i] = sample;
					pointer++;
					if (pointer == end)
						pointer -= size;
				}

				// get the complex sinusoid phase
				auto z = cpl::dsp::goertzel(temporaryBuffer, count2, cpl::simd::consts<double>::tau * (maxIndex + delta) / TransformSize);

				auto argz = -sampleDifference * cpl::simd::consts<double>::tau * (maxIndex + delta) / TransformSize;

				auto matrix = std::complex<double>(std::cos(argz), -std::sin(argz));

				z *= matrix;

				// invert unit circle as the display is "backwards" in time
				auto phase = cpl::simd::consts<double>::tau - std::arg(z);
				// correct phase by delta
				phase += delta * cpl::simd::consts<double>::tau;
				// phase correct to sines
				phase -= cpl::simd::consts<double>::pi_half;
				// add user-defined phase offset
				phase += cpl::simd::consts<double>::tau * content->triggerPhaseOffset.getParameterView().getValueTransformed() / 360;
				// since we can't go back in time, travel around the unit circle until the phase is positive.
				while (phase < 0)
					phase += cpl::simd::consts<double>::tau;

				while (phase > cpl::simd::consts<double>::tau)
					phase -= cpl::simd::consts<double>::tau;

				detectedPhase = phase;

				// normalize phase to cycles
				auto cycles = (phase / (cpl::simd::consts<double>::tau));


				// invert unit circle as the display is "backwards" in time
				auto phase2 = cpl::simd::consts<double>::tau - std::arg(z2);
				// correct phase by delta
				phase2 += delta * cpl::simd::consts<double>::tau;
				// phase correct to sines
				phase2 -= cpl::simd::consts<double>::pi_half;
				// add user-defined phase offset
				phase2 += cpl::simd::consts<double>::tau * content->triggerPhaseOffset.getParameterView().getValueTransformed() / 360;
				// since we can't go back in time, travel around the unit circle until the phase is positive.
				while (phase2 < 0)
					phase2 += cpl::simd::consts<double>::tau;

				// normalize phase to cycles
				auto cycles2 = (phase2 / (cpl::simd::consts<double>::tau));

				temp8kOffset = cycles2 * audioStream.getAudioHistorySamplerate() / (detectedFreq);

				// convert to samples
				return cycles * audioStream.getAudioHistorySamplerate() / (detectedFreq);

			}


		}

		return 0;
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
		}

};
