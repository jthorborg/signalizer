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

	double COscilloscope::getTriggeringOffset(const AudioStream::AudioBufferAccess & buf)
	{
		if (buf.getNumChannels() < 1)
			return 0;
		switch (cpl::enum_cast<OscilloscopeContent::TriggeringMode>(content->triggerMode.param.getTransformedValue()))
		{
			case OscilloscopeContent::TriggeringMode::Spectral:
			{
				auto && view = buf.getView(0);
					
				if (view.size() < OscilloscopeContent::LookaheadSize)
					return 0;

				transformBuffer.resize(OscilloscopeContent::LookaheadSize);
				buffer.resize(OscilloscopeContent::LookaheadSize);

				for (std::size_t i = 0; i < transformBuffer.size(); ++i)
				{
					transformBuffer[i] = view[i];
				}

				signaldust::DustFFT_fwdDa(reinterpret_cast<double*>(transformBuffer.data()), transformBuffer.size());

				auto quadDelta = [&](auto w) {
					auto x0 = transformBuffer[w], x1 = transformBuffer[w + 1], xm1 = transformBuffer[w == 0 ? 1 : w - 1];
					return std::real((xm1 - x1) / ((x0 * 2.0) - xm1 - x1));
				};

				std::size_t maxIndex = 1;
				double maxValue = std::abs(transformBuffer[1]);
				double delta = quadDelta(1);

				const double halfSemitone = std::pow(2, 0.5 / 12.0) - 1;

				for (std::size_t i = 2; i < transformBuffer.size() >> 1; ++i)
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

							//auto sensivity = newValue / maxValue;

							// the same value, just a better estimate.
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

				quantizedFreq = audioStream.getAudioHistorySamplerate() * double(maxIndex) / transformBuffer.size();

				// interpolate peak position quadratically
				delta = quadDelta(maxIndex);
				detectedFreq = audioStream.getAudioHistorySamplerate() * (maxIndex + delta) / transformBuffer.size();

				// set up goertzel buffer
				for (std::size_t i = 0; i < transformBuffer.size(); ++i)
				{
					buffer[i] = view[i];
				}

				// get the complex sinusoid phase
				auto z = cpl::dsp::goertzel(buffer, buffer.size(), cpl::simd::consts<double>::tau * (maxIndex + delta) / buffer.size());

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

				// normalize phase to cycles
				auto cycles = (phase / (cpl::simd::consts<double>::tau));

				// convert to samples
				return cycles * audioStream.getAudioHistorySamplerate() / (detectedFreq);

			}


		}

		return 0;
	}

	
};
