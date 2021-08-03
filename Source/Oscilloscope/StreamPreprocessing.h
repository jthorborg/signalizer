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

	file:OscilloscopeStreamProcessing.h

		Functions and classes for preprocessing audio streams before they
		reach the oscilloscope.

*************************************************************************************/

#ifndef OSCILLOSCOPE_STREAM_PROCESSING_H
	#define OSCILLOSCOPE_STREAM_PROCESSING_H

	#include "Signalizer.h"
	#include "Oscilloscope.h"
	#include <queue>

	namespace Signalizer
	{
		class TriggeringProcessor
		{
		public:

			template<typename ISA> friend class SignalStreamBaseProcessor;

			void setSettings(OscilloscopeContent::TriggeringMode triggerMode, double newWindowSize, double valueThreshold, double newHysteresis)
			{
				triggerType = triggerMode;
				windowChanged = std::ceil(newWindowSize) != std::ceil(windowSize);
				windowSize = newWindowSize;
				hysteresis = newHysteresis;
				threshold = valueThreshold;
			}

			void update(std::uint64_t currentSteadyClock)
			{
				steadyClock = currentSteadyClock;
				if (windowChanged)
				{
					windowChanged = false;

					if (isWorkingOnPeak)
					{
						peaks.pop();
					}

					while (peaks.size() && peaks.front() < currentSteadyClock)
						peaks.pop();

					bufferedSamples = currentPeak = oldPeak = 0;
					frontOrigin = currentSteadyClock;

					isWorkingOnPeak = false;

				}
			}

			template<typename ISA>
			void processMutating(Oscilloscope & o, AFloat ** localPointers, std::size_t numChannels, std::size_t numSamples)
			{
				if (frontOrigin + bufferedSamples < steadyClock)
				{
					frontOrigin = steadyClock;
					bufferedSamples = 0;
				}

				auto ceilingSize = std::ceil(windowSize);
				auto halfSize = ceilingSize / 2;

				auto processIntoBackBuffer = [&](auto samples)
				{
					o.audioProcessing<ISA>(localPointers, numChannels, samples, o.channelData.back);
					numSamples -= samples;
					for (std::size_t c = 0; c < numChannels; ++c)
						localPointers[c] += samples;

					auto oldSamples = bufferedSamples;
					steadyClock += samples;
					bufferedSamples += samples;
					bufferedSamples = std::min<std::uint64_t>(bufferedSamples, ceilingSize + 1);

					frontOrigin += (oldSamples + samples) - bufferedSamples;
				};

				if (ceilingSize == 0 && peaks.size())
				{
					std::queue<std::uint64_t>().swap(peaks);
				}

				while (numSamples != 0)
				{
					if (!peaks.size())
					{
						processIntoBackBuffer(numSamples);
						break;
					}
					else if (!isWorkingOnPeak)
					{

						isWorkingOnPeak = true;

						auto nextPeak = peaks.front();
						if (nextPeak >= steadyClock)
						{
							// select closest peak that has a full buffer
							auto deltaToPeak = nextPeak - steadyClock;
							auto numSamplesToProcess = std::min<std::uint64_t>(numSamples, deltaToPeak + halfSize);
							processIntoBackBuffer(numSamplesToProcess);
							currentPeak = nextPeak;
						}
						else
						{
							// select closest peak that has a full buffer
							currentPeak = nextPeak;

						}
					}

					std::uint64_t windowEnd;
					bool isPeakOutsideOfWindow = false;
					bool readyForBufferSwap = false;

					if (currentPeak - oldPeak < halfSize)
					{
						windowEnd = oldPeak + halfSize;
					}
					else
					{
						isPeakOutsideOfWindow = true;
						windowEnd = frontOrigin + bufferedSamples;

					}

					std::uint64_t peakWindowEnd = (isPeakOutsideOfWindow ? 1 : 0) + currentPeak + halfSize;
					std::uint64_t numSamplesToProcess = 0;

					auto missingBufferSamples = peakWindowEnd - std::min(peakWindowEnd, windowEnd);

					auto neededPreSamples = std::min<std::uint64_t>(halfSize, std::max<double>((currentPeak - oldPeak), halfSize) - halfSize);

					if (isPeakOutsideOfWindow)
					{
						numSamplesToProcess = std::min<std::uint64_t>(numSamples, missingBufferSamples);

						if (numSamplesToProcess > 0)
							processIntoBackBuffer(numSamplesToProcess);

						readyForBufferSwap = missingBufferSamples == numSamplesToProcess;
					}
					else
					{
						if (bufferedSamples >= missingBufferSamples)
						{
							readyForBufferSwap = true;
						}
						else
						{
							auto numRemaining = missingBufferSamples - bufferedSamples;

							numSamplesToProcess = std::min<std::uint64_t>(numSamples, numRemaining);

							if (numSamplesToProcess > 0)
								processIntoBackBuffer(numSamplesToProcess);

							readyForBufferSwap = numRemaining == numSamplesToProcess;
						}

					}

					if (readyForBufferSwap)
					{
						auto amount = (isPeakOutsideOfWindow ? halfSize : missingBufferSamples) + neededPreSamples;

						auto cappedSize = std::min<std::size_t>(bufferedSamples, std::ceil(amount + 1));
						// 1
						o.channelData.swapBuffers(cappedSize, -(cpl::ssize_t)bufferedSamples);
						// 2
						bufferedSamples -= std::min<std::uint64_t>(bufferedSamples, cappedSize);
						// 3
						frontOrigin += cappedSize;
						oldPeak = currentPeak;
						isWorkingOnPeak = false;
						peaks.pop();
					}
				}
			}

		private:

			double hysteresis;
			double threshold;
			double windowSize;
			bool windowChanged;
			double state;
			OscilloscopeContent::TriggeringMode triggerType;
			std::uint64_t crossOrigin;
			bool isPeakHold;
			std::uint64_t oldPeak;
			std::uint64_t currentPeak, bufferedSamples;
			std::uint64_t frontOrigin;
			std::uint64_t steadyClock;
			std::queue<std::uint64_t> peaks;
			bool isWorkingOnPeak;


		};

		template<typename ISA>
		class SignalStreamBaseProcessor
		{
		public:

			SignalStreamBaseProcessor(std::size_t numChannels, std::size_t numSamples, std::uint64_t steadyClock, TriggeringProcessor & outsideState)
				: numSamples(numSamples)
				, outsideState(outsideState)
				, state(outsideState.state)
				, numChannels(numChannels)
				, steadyClock(steadyClock)
				, isPeakHolding(outsideState.isPeakHold)
				, threshold(outsideState.threshold)
				, hysteresis(outsideState.hysteresis)
				, peaks(outsideState.peaks)
				, crossOrigin(outsideState.crossOrigin)
			{

			}

			~SignalStreamBaseProcessor()
			{
				outsideState.state = state;
				outsideState.isPeakHold = isPeakHolding;
				outsideState.crossOrigin = crossOrigin;
			}

		protected:

			std::size_t numSamples;
			const std::uint64_t steadyClock;
			std::uint64_t crossOrigin;
			TriggeringProcessor & outsideState;
			std::queue<std::uint64_t> & peaks;
			bool isPeakHolding;
			const double threshold, hysteresis;
			double state;
			std::size_t processedSamples;
			const std::size_t numChannels;
		};

		template<typename ISA>
		class PeakHoldProcessor : SignalStreamBaseProcessor<ISA>
		{
		    using SignalStreamBaseProcessor<ISA>::state;
		    using SignalStreamBaseProcessor<ISA>::threshold;
		    using SignalStreamBaseProcessor<ISA>::peaks;
		    using SignalStreamBaseProcessor<ISA>::isPeakHolding;
		    using SignalStreamBaseProcessor<ISA>::hysteresis;
		    using SignalStreamBaseProcessor<ISA>::steadyClock;

		public:
			std::size_t count = 0;

			using SignalStreamBaseProcessor<ISA>::SignalStreamBaseProcessor;

			inline void process(double sample) noexcept
			{
				sample *= sample;
				auto delta = sample - state;

				if (delta < 0)
				{
					state *= 0.9999;

					state = std::max(threshold * threshold, state);

					if (isPeakHolding)
					{
						// minus one since this is the first sample that doesn't qualify as a (rising) peak
						peaks.push(steadyClock + count - 1);
						isPeakHolding = false;
					}
				}
				else
				{
					if (delta > hysteresis * state)
						isPeakHolding = true;

					state = sample;
				}

				count++;
			}

		};

		template<typename ISA>
		class ZeroCrossingProcessor : SignalStreamBaseProcessor<ISA>
		{
		    using SignalStreamBaseProcessor<ISA>::state;
		    using SignalStreamBaseProcessor<ISA>::threshold;
		    using SignalStreamBaseProcessor<ISA>::peaks;
		    using SignalStreamBaseProcessor<ISA>::isPeakHolding;
		    using SignalStreamBaseProcessor<ISA>::hysteresis;
		    using SignalStreamBaseProcessor<ISA>::steadyClock;
		    using SignalStreamBaseProcessor<ISA>::crossOrigin;

		public:
			std::size_t count = 0;

			using SignalStreamBaseProcessor<ISA>::SignalStreamBaseProcessor;

			inline void process(double sample) noexcept
			{
				if (sample > 0 && state < 0)
				{
					isPeakHolding = true;
					crossOrigin = steadyClock + count;
				}

				if (isPeakHolding && sample > threshold)
				{
					isPeakHolding = false;
					peaks.push(crossOrigin);
				}

				state = sample;
				count++;
			}

		};

	};



#endif
