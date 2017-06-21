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

		Interface for the Oscilloscope view

*************************************************************************************/

#ifndef OSCILLOSCOPE_STREAM_PROCESSING_H
	#define OSCILLOSCOPE_STREAM_PROCESSING_H

	#include "Signalizer.h"
	#include <queue>

	namespace Signalizer
	{
		struct PreprocessingTriggerState
		{
			double hysterisis;
			double timeConstant;
			double windowSize;
			bool windowChanged;
			double state;
			double lastOffset, currentOffset;
			bool isPeakHold;
			std::uint64_t oldPeak;
			std::uint64_t currentPeak, bufferedSamples;
			std::uint64_t frontOrigin, backOrigin;
			std::queue<std::uint64_t> peaks;
			bool isWorkingOnPeak;
			volatile bool block;

			void handleStateChanges(std::uint64_t currentTime)
			{
				if (windowChanged)
				{
					windowChanged = false;

					if (isWorkingOnPeak)
					{
						peaks.pop();
					}

					while (peaks.size() && peaks.front() < currentTime)
						peaks.pop();

					currentPeak = oldPeak = 0;
					backOrigin = bufferedSamples = 0;
					frontOrigin = currentTime;

					isWorkingOnPeak = false;

				}
			}

		};

		template<typename ISA>
		class SignalStreamBaseProcessor
		{
		public:

			SignalStreamBaseProcessor(AFloat ** buffer, std::size_t numChannels, std::size_t & numSamples, std::uint64_t steadyClock, PreprocessingTriggerState & outsideState)
				: buffer(buffer)
				, numSamples(numSamples)
				, outsideState(outsideState)
				, state(outsideState.state)
				, lastOffset(outsideState.lastOffset)
				, currentOffset(outsideState.currentOffset)
				, windowSize(outsideState.windowSize)
				, processedSamples(0)
				, numChannels(numChannels)
				, steadyClock(steadyClock)
				, isPeakHolding(outsideState.isPeakHold)
			{

			}

			~SignalStreamBaseProcessor()
			{
				outsideState.state = state;
				outsideState.lastOffset = lastOffset;
				outsideState.currentOffset = currentOffset;
				outsideState.isPeakHold = isPeakHolding;
			}

		protected:
			AFloat ** buffer;
			std::size_t & numSamples;
			std::uint64_t steadyClock;
			PreprocessingTriggerState & outsideState;
			bool isPeakHolding;
			double state, lastOffset, currentOffset, windowSize;
			std::size_t processedSamples, numChannels;
		};

		template<typename ISA>
		class PeakHoldProcessor : SignalStreamBaseProcessor<ISA>
		{
		public:
			std::size_t count = 0;

			using SignalStreamBaseProcessor<ISA>::SignalStreamBaseProcessor;

			inline void process(double sample) noexcept 
			{
				sample *= sample;

				if (sample < state)
				{
					state *= 0.9999;
					
					if (isPeakHolding)
					{
						// minus one since this is the first sample that doesn't qualify as a (rising) peak
						outsideState.peaks.push(steadyClock + count - 1);
						isPeakHolding = false;
					}
				}
				else
				{
					state = sample;
					isPeakHolding = true;
				}

				count++;
			}

			~PeakHoldProcessor()
			{
				/*auto skippedSamples = numSamples - processedSamples;

				for (std::size_t c = 0; c < numChannels; ++c)
				{
					buffer[c] += skippedSamples;
				}

				numSamples = processedSamples; */
			}

		};

		template<typename ISA>
		class ZeroCrossingProcessor : SignalStreamBaseProcessor<ISA>
		{
		public:

			using SignalStreamBaseProcessor<ISA>::SignalStreamBaseProcessor;

			inline void process(double sample) noexcept 
			{
				if (sample > 0 && state < 0)
				{
					state = sample;
					lastOffset = currentOffset;
					currentOffset = 0;
				}
				else
				{
					lastOffset += 1;
					currentOffset += 1;
				}

				state = sample;
			}

		};

	};



#endif
