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

	file:ChannelData.h

		Management of the I/O channel data for the oscillscope

*************************************************************************************/

#ifndef SIGNALIZER_CHANNEL_DATA_H
	#define SIGNALIZER_CHANNEL_DATA_H

	#include "Signalizer.h"
	#include <cpl/simd.h>
	#include <cpl/dsp/LinkwitzRileyNetwork.h>
	#include <cpl/dsp/SmoothedParameterState.h>

	namespace Signalizer
	{
		struct ChannelData
		{
			static const std::size_t Bands = 3;
			typedef cpl::dsp::LinkwitzRileyNetwork<AFloat, Bands> Crossover;
			typedef cpl::GraphicsND::UPixel<cpl::GraphicsND::ComponentOrder::OpenGL> PixelType;
			typedef cpl::CLIFOStream<AFloat, 32> AudioBuffer;
			typedef cpl::CLIFOStream<PixelType, 32> ColourBuffer;

			enum Entry
			{
				Slow = 0,
				Left = 0,
				Fast = 1,
				Right = 1
			};

			struct Channel
			{
				AudioBuffer audioData;
				ColourBuffer colourData;
			};

			struct FilterStates
			{
				struct ChannelState
				{
					Crossover::BandArray smoothFilters{};
					Crossover network;
					juce::Colour defaultKey;
					AFloat envelope;
				};

				std::vector<ChannelState> channels;

				Crossover::BandArray midSideSmoothsFilters[2];
			};

			struct Buffer
			{
				std::vector<Channel> channels{ 0};
				ColourBuffer midSideColour[2];

				Channel & defaultChannel()
				{
					return channels[0];
				}

				void resizeStorage(std::size_t samples, std::size_t capacity = -1)
				{
					if (capacity == static_cast<std::size_t>(-1))
						capacity = cpl::Math::nextPow2Inc(samples);

					for (auto & c : channels)
					{
						c.audioData.setStorageRequirements(samples, capacity);
						c.colourData.setStorageRequirements(samples, capacity);
					}

					for (auto & c : midSideColour)
					{
						c.setStorageRequirements(samples, capacity);
					}
				}
			};

			void resizeChannels(std::size_t newChannels)
			{
				for (auto buffer : { &back, &front })
				{
					auto& channelsForBuffer = buffer->channels;

					bool alreadyHasData = channelsForBuffer.size() > 0;

					channelsForBuffer.resize(newChannels);

					if (alreadyHasData)
						buffer->resizeStorage(channelsForBuffer.front().audioData.getSize(), channelsForBuffer.front().audioData.getCapacity());
				}

				filterStates.channels.resize(std::max(filterStates.channels.size(), newChannels));

			}

			void swapBuffers(std::size_t historySize, cpl::ssize_t offset)
			{
				auto swapBuf = [historySize, offset](const auto & inBuf, auto & outBuf)
				{
					outBuf.createWriter().copyIntoHead(inBuf.createProxyView(), historySize, offset);
				};

				for (std::size_t i = 0; i < std::extent<decltype(Buffer::midSideColour)>::value; ++i)
				{
					swapBuf(back.midSideColour[i], front.midSideColour[i]);
				}

				for (std::size_t i = 0; i < back.channels.size(); ++i)
				{
					swapBuf(back.channels[i].audioData, front.channels[i].audioData);
					swapBuf(back.channels[i].colourData, front.channels[i].colourData);
				}
			}

			void tuneCrossOver(double lowCrossover, double highCrossover, double sampleRate)
			{
				networkCoeffs = Crossover::Coefficients::design({ static_cast<AFloat>(lowCrossover / sampleRate), static_cast<AFloat>(highCrossover / sampleRate) });
			}

			void tuneColourSmoothing(double milliseconds, double sampleRate)
			{
				smoothFilterPole = cpl::dsp::SmoothedParameterState<AFloat, 1>::design(milliseconds, sampleRate);
			}

			std::size_t numChannels() const noexcept { return front.channels.size(); }
			bool empty() const noexcept { return numChannels() > 0; }

			Crossover::Coefficients networkCoeffs;
			cpl::dsp::SmoothedParameterState<AFloat, 1>::PoleState smoothFilterPole;

			FilterStates filterStates;
			Buffer back, front;

		};
	};

#endif
