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

	namespace Signalizer
	{
		struct PreprocessingTriggerState
		{
			double hysterisis;
			double timeConstant;
			double windowSize;
			double state;
			double lastOffset, currentOffset;
		};

		template<typename ISA>
		class SignalStreamBaseProcessor
		{
		public:

			SignalStreamBaseProcessor(AFloat ** buffer, std::size_t numChannels, std::size_t & numSamples, PreprocessingTriggerState & outsideState)
				: buffer(buffer)
				, numSamples(numSamples)
				, outsideState(outsideState)
				, state(outsideState.state)
				, lastOffset(outsideState.lastOffset)
				, currentOffset(outsideState.currentOffset)
				, windowSize(outsideState.windowSize)
				, processedSamples(0)
				, numChannels(numChannels)
			{

			}

			~SignalStreamBaseProcessor()
			{
				outsideState.state = state;
				outsideState.lastOffset = lastOffset;
				outsideState.currentOffset = currentOffset;
			}

		protected:
			AFloat ** buffer;
			std::size_t & numSamples;
			
			PreprocessingTriggerState & outsideState;

			double state, lastOffset, currentOffset, windowSize;
			std::size_t processedSamples, numChannels;
		};

		template<typename ISA>
		class PeakHoldProcessor : SignalStreamBaseProcessor<ISA>
		{
		public:

			using SignalStreamBaseProcessor<ISA>::SignalStreamBaseProcessor;

			inline void process(double sample) noexcept 
			{
				sample *= sample;
				if (sample < state)
				{
					state *= 0.9999;

					if (2 * windowSize > currentOffset)
					{
						lastOffset += 1;
						currentOffset += 1;
						processedSamples += 1;
					}

				}
				else
				{
					state = sample;
					lastOffset = currentOffset;
					currentOffset = 0;
				}
			}

			~PeakHoldProcessor()
			{
				auto skippedSamples = numSamples - processedSamples;

				for (std::size_t c = 0; c < numChannels; ++c)
				{
					buffer[c] += skippedSamples;
				}

				numSamples = processedSamples;
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
