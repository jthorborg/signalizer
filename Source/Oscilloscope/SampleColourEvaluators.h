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
 
	file:SampleColuorEvaluators.h

		Implementation of evaluators
 
*************************************************************************************/

#ifndef SIGNALIZER_SAMPLECOLOUREVALUATORS_H
	#define SIGNALIZER_SAMPLECOLOUREVALUATORS_H

	#include "Oscilloscope.h"
	#include <functional>

	namespace Signalizer
	{

		class Oscilloscope::DefaultKey
		{
		public:
			DefaultKey(ChannelData & data, std::size_t index)
				: defaultKey(data.filterStates.channels.at(index).defaultKey)
			{

			}

			ChannelData::PixelType getDefaultKey() const noexcept
			{
				return defaultKey;
			}

		private:
			ChannelData::PixelType defaultKey;
		};

		class Oscilloscope::DynamicChannelEvaluator : public Oscilloscope::SampleColourEvaluatorBase, public Oscilloscope::DefaultKey
		{
		public:

			DynamicChannelEvaluator(const EvaluatorParams& params)
				: DefaultKey(params.data, params.colourIndex)
				, audioView(params.data.front.channels.at(params.channelIndex).audioData.createProxyView())
				, colourView(params.data.front.channels.at(params.channelIndex).colourData.createProxyView())
			{

			}

			inline bool isWellDefined() const noexcept
			{
				return audioView.size() > 0 && colourView.size() > 0;
			}

			void startFrom(cpl::ssize_t offset)
			{
				startFrom(offset, offset);
			}

			void startFrom(cpl::ssize_t audioOffset, cpl::ssize_t colourOffset)
			{
				audioPointer = audioView.begin() + audioView.cursorPosition() + audioOffset;

				while (audioPointer < audioView.begin())
					audioPointer += audioView.size();

				while (audioPointer >= audioView.end())
					audioPointer -= audioView.size();

				colourPointer = colourView.begin() + colourView.cursorPosition() + colourOffset;

				while (colourPointer < colourView.begin())
					colourPointer += colourView.size();

				while (colourPointer >= colourView.end())
					colourPointer -= colourView.size();

			}

			inline void inc() noexcept
			{
				audioPointer++, colourPointer++;

				if (audioPointer == audioView.end())
					audioPointer -= audioView.size();

				if (colourPointer == colourView.end())
					colourPointer -= colourView.size();
			}

			inline std::pair<std::ptrdiff_t, std::ptrdiff_t> distance() const noexcept
			{
				return std::make_pair(
					std::distance<AudioIt>(audioPointer, audioView.begin()),
					std::distance<ColourIt>(colourPointer, colourView.begin())
				);
			}

			inline std::pair<AudioT, ColourT> evaluate() const noexcept
			{
				return { *audioPointer, *colourPointer };
			}

			AudioT evaluateSample() const noexcept
			{
				return *audioPointer;
			}

			ColourT evaluateColour() const noexcept
			{
				return *colourPointer;
			}

			AudioT evaluateSampleInc() noexcept
			{
				auto ret = *audioPointer++;

				if (audioPointer == audioView.end())
					audioPointer -= audioView.size();

				return ret;
			}

			ColourT evaluateColourInc() noexcept
			{
				auto ret = *colourPointer++;

				if (colourPointer == colourView.end())
					colourPointer -= colourView.size();

				return ret;
			}

		private:

			juce::Colour defaultKey;

			ChannelData::AudioBuffer::ProxyView audioView;
			ChannelData::ColourBuffer::ProxyView colourView;

			AudioIt audioPointer {};
			ColourIt colourPointer {};
		};

		template<std::size_t ChannelIndex, std::size_t ColourIndex>
		class Oscilloscope::SimpleChannelEvaluator : public Oscilloscope::DynamicChannelEvaluator
		{
		public:
			SimpleChannelEvaluator(const EvaluatorParams& params)
				: Oscilloscope::DynamicChannelEvaluator({ params.data, ChannelIndex, ColourIndex })
			{

			}
		};

		template<std::size_t ChannelIndex, std::size_t ColourIndex, typename BinaryFunction>
			class Oscilloscope::MidSideEvaluatorBase : public Oscilloscope::SampleColourEvaluatorBase, public Oscilloscope::DefaultKey
			{
			public:

				MidSideEvaluatorBase(const EvaluatorParams& params)
					: DefaultKey(params.data, ColourIndex)
					, audioViewLeft(params.data.front.channels.at(0).audioData.createProxyView())
					, audioViewRight(params.data.front.channels.at(1).audioData.createProxyView())
					, colourView(params.data.front.midSideColour[ChannelIndex].createProxyView())
				{

				}

				inline bool isWellDefined() const noexcept
				{
					return audioViewLeft.size() > 0 && audioViewRight.size() > 0 && audioViewLeft.size() == audioViewRight.size() && colourView.size() > 0;
				}

				void startFrom(cpl::ssize_t offset)
				{
					startFrom(offset, offset);
				}

				void startFrom(cpl::ssize_t audioOffset, cpl::ssize_t colourOffset)
				{
					audioOffset += audioViewLeft.cursorPosition();
					audioPointerLeft = audioViewLeft.begin() + audioOffset;
					audioPointerRight = audioViewRight.begin() + audioOffset;

					while (audioPointerLeft < audioViewLeft.begin())
					{
						audioPointerLeft += audioViewLeft.size();
						audioPointerRight += audioViewLeft.size();
					}


					while (audioPointerLeft >= audioViewLeft.end())
					{
						audioPointerLeft -= audioViewLeft.size();
						audioPointerRight -= audioViewLeft.size();
					}

					colourPointer = colourView.begin() + colourView.cursorPosition() + colourOffset;

					while (colourPointer < colourView.begin())
						colourPointer += colourView.size();

					while (colourPointer >= colourView.end())
						colourPointer -= colourView.size();

				}

				inline void inc() noexcept
				{
					audioPointerLeft++, audioPointerRight++, colourPointer++;

					if (audioPointerLeft == audioViewLeft.end())
					{
						audioPointerLeft -= audioViewLeft.size();
						audioPointerRight -= audioViewLeft.size();
					}

					if (colourPointer == colourView.end())
						colourPointer -= colourView.size();
				}

				inline std::pair<std::ptrdiff_t, std::ptrdiff_t> distance() const noexcept
				{
					return std::make_pair(
						std::distance<AudioIt>(audioPointerLeft, audioViewLeft.begin()),
						std::distance<ColourIt>(colourPointer, colourView.begin())
					);
				}

				inline std::pair<AudioT, ColourT> evaluate() const noexcept
				{
					return{ evaluateSample(), evaluateColour() };
				}

				AudioT evaluateSample() const noexcept
				{
					return static_cast<AudioT>(0.5) * BinaryFunction()(*audioPointerLeft, *audioPointerRight);
				}

				ColourT evaluateColour() const noexcept
				{
					return *colourPointer;
				}

				AudioT evaluateSampleInc() noexcept
				{
					auto ret = evaluateSample();

					audioPointerLeft++, audioPointerRight++;

					if (audioPointerLeft == audioViewLeft.end())
					{
						audioPointerLeft -= audioViewLeft.size();
						audioPointerRight -= audioViewLeft.size();
					}

					return ret;
				}

				ColourT evaluateColourInc() noexcept
				{
					auto ret = *colourPointer++;

					if (colourPointer == colourView.end())
						colourPointer -= colourView.size();

					return ret;
				}

			private:

				ChannelData::AudioBuffer::ProxyView audioViewLeft, audioViewRight;
				ChannelData::ColourBuffer::ProxyView colourView;

				AudioIt audioPointerLeft {}, audioPointerRight {};
				ColourIt colourPointer{};
			};

		template<>
			class Oscilloscope::SampleColourEvaluator<OscChannels::Left, 0> : public SimpleChannelEvaluator<0, 0>
			{
				using SimpleChannelEvaluator<0, 0>::SimpleChannelEvaluator;
			};

		template<>
			class Oscilloscope::SampleColourEvaluator<OscChannels::Left, 1> : public SimpleChannelEvaluator<0, 1>
			{
				using SimpleChannelEvaluator<0, 1>::SimpleChannelEvaluator;
			};

		template<>
			class Oscilloscope::SampleColourEvaluator<OscChannels::Right, 0> : public SimpleChannelEvaluator<1, 0>
			{
				using SimpleChannelEvaluator<1, 0>::SimpleChannelEvaluator;
			};

		template<>
			class Oscilloscope::SampleColourEvaluator<OscChannels::Right, 1> : public SimpleChannelEvaluator<1, 1>
			{
				using SimpleChannelEvaluator<1, 1>::SimpleChannelEvaluator;
			};


		template<>
			class Oscilloscope::SampleColourEvaluator<OscChannels::Mid, 0> : public MidSideEvaluatorBase<0, 0, std::plus<>>
			{
				using MidSideEvaluatorBase<0, 0, std::plus<>>::MidSideEvaluatorBase;
			};

		template<>
			class Oscilloscope::SampleColourEvaluator<OscChannels::Mid, 1> : public MidSideEvaluatorBase<0, 1, std::plus<>>
			{
				using MidSideEvaluatorBase<0, 1, std::plus<>>::MidSideEvaluatorBase;
			};

		template<>
			class Oscilloscope::SampleColourEvaluator<OscChannels::Side, 0> : public MidSideEvaluatorBase<1, 0, std::minus<>>
			{
				using MidSideEvaluatorBase<1, 0, std::minus<>>::MidSideEvaluatorBase;
			};

		template<>
			class Oscilloscope::SampleColourEvaluator<OscChannels::Side, 1> : public MidSideEvaluatorBase<1, 1, std::minus<>>
			{
				using MidSideEvaluatorBase<1, 1, std::minus<>>::MidSideEvaluatorBase;
			};
	};

#endif