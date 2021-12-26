/*************************************************************************************

	Signalizer - cross-platform audio visualization plugin - v. 0.x.y

	Copyright (C) 2021 Janus Lynggaard Thorborg (www.jthorborg.com)

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

	file:MixGraphListener.h

		Runtime mixer instance servicing the views

*************************************************************************************/


#ifndef SIGNALIZER_MIXGRAPHLISTENER_H
	#define SIGNALIZER_MIXGRAPHLISTENER_H

	#include <cpl/Common.h>
	#include <cpl/state/Serialization.h>
	#include <cpl/CAudioStream.h>
	#include <mutex>
	#include <shared_mutex>
	#include <map>
	#include <set>
	#include <string>
	#include <vector>
	#include <optional>
	#include "SignalizerConfiguration.h"

	namespace Signalizer
	{
		typedef std::make_signed<std::size_t>::type ssize_t;
		
		template<class T>
		using vector_set = std::set<T>;

		typedef std::int32_t PinInt;

		struct DirectedPortPair
		{
			PinInt Source;
			PinInt Destination;

			inline friend bool operator < (DirectedPortPair a, DirectedPortPair b)
			{
				return std::tie(a.Destination, a.Source) < std::tie(b.Destination, b.Source);
			}
		};

		// TODO: Refactor into CPL, share between ape and such.
		class AuxMatrix
		{
		public:

			void resizeChannels(std::size_t length)
			{
				auxBuffers.resize(length);
			}

			void softBufferResize(std::size_t length)
			{
				auto newSize = length * auxBuffers.size();
				auxData.resize(std::max(newSize, auxData.size()));

				for (std::size_t i = 0; i < auxBuffers.size(); ++i)
					auxBuffers[i] = auxData.data() + length * i;

				bufferLength = length;
			}

			void copy(const float* const* buffers, std::size_t channelMatrixOffset, std::size_t numBuffers)
			{
				for (std::size_t i = 0; i < numBuffers; ++i)
				{
					std::memcpy(auxBuffers[i + channelMatrixOffset], buffers[i], bufferLength * sizeof(float));
				}
			}

			void accumulate(const float* const* buffers, std::size_t index, std::size_t numBuffers, float start, float end)
			{
				const auto delta = end - start;

				for (std::size_t i = 0; i < numBuffers; ++i)
				{
					for (std::size_t n = 0; n < bufferLength; ++n)
					{
						const float progress = n / float(bufferLength - 1);
						auxBuffers[i + index][n] += buffers[i][n] * (start + progress * delta);
					}
				}
			}

			void clear(std::size_t index, std::size_t numBuffers)
			{
				for (std::size_t i = 0; i < numBuffers; ++i)
				{
					std::memset(auxBuffers[i + index], 0, bufferLength * sizeof(float));
				}
			}

			void copyResample(const float* buffer, std::size_t index, std::size_t numSamples)
			{
				if (numSamples == bufferLength)
					return copy(&buffer, index, 1);

				auto ratio = (double)numSamples / bufferLength;
				double x = 0;

				cpl::Types::fsint_t samples = static_cast<cpl::Types::fsint_t>(numSamples);

				for (std::size_t i = 0; i < bufferLength; ++i)
				{
					auxBuffers[index][i] = cpl::dsp::linearFilter<float>(buffer, samples, x);
					x += ratio;
				}
			}

			float* operator [] (std::size_t index) const { return auxBuffers[index]; }
			float** data() noexcept { return auxBuffers.data(); }
			std::size_t size() const noexcept { return auxBuffers.size(); }

		private:
			std::size_t bufferLength = 0;
			std::vector<float> auxData;
			std::vector<float*> auxBuffers;
		};

		
		class MixGraphListener
		{
		public:
			friend class HostGraph;

			MixGraphListener(std::shared_ptr<AudioStream> realtimeStream, AudioStream& presentation);

			void connect(std::shared_ptr<AudioStream>& stream, DirectedPortPair pair);
			void disconnect(std::shared_ptr<AudioStream>& stream, DirectedPortPair pair);
			void setExpectedBufferSize(std::size_t bufferSize);

		private:
			
			typedef cpl::CLIFOStream<AFloat> Buffer;

			struct ConnectionCommand
			{
				std::shared_ptr<AudioStream> stream;
				DirectedPortPair pair;
				bool isConnection;
			};

			struct State : AudioStream::Listener
			{
				struct Channel
				{
					cpl::CLIFOStream<AFloat> buffer;
				};

				std::map<DirectedPortPair, Channel> channelQueues;
				std::atomic_size_t current {};
				std::int64_t globalPosition{};
				std::int64_t endpoint{};
				std::int64_t offset{};
				std::int32_t refCount;
				MixGraphListener& parent;
				// retain the source.
				std::shared_ptr<AudioStream> source;

				State(const State& other) = delete;
				State(State&& other) = default;
				State& operator = (State&& other) = default;

				State(MixGraphListener& parent, std::shared_ptr<AudioStream>& source)
					: parent(parent), source(source), refCount(1)
				{
					listenToSource(*source.get(), false, 1000);
				}

				virtual void onAsyncChangedProperties(const Stream& changedSource, const typename Stream::AudioStreamInfo& before)
				{
					parent.asyncPropertiesChanged(*this, changedSource);
				}

				virtual bool onAsyncAudio(const AudioStream& source, AFloat** buffer, std::size_t numChannels, std::size_t numSamples) override
				{
					parent.onAsyncAudio(*this, buffer, numChannels, numSamples, source.getASyncPlayhead().getPositionInSamples());
					return false;
				}

				~State()
				{
					detachFromSource();
				}
			};

			auto emplace(std::shared_ptr<AudioStream>& stream);

			void handleStructuralChange(double sampleRate, std::size_t numSamples);
			void deliver(std::size_t numSamples);
			void asyncPropertiesChanged(State& s, const AudioStream& source);
			void onAsyncAudio(State& s, AFloat** buffer, std::size_t numChannels, std::size_t numSamples, std::int64_t globalPosition);

			void updateTopologyCommands();

			std::map<const AudioStream*, State> graph;
			std::size_t initialSampleCapacity {};
			std::shared_ptr<AudioStream> realtime;
			AudioStream& presentation;
			State* self;
			std::shared_mutex dataMutex;

			std::mutex connectDisconnectMutex;
			std::vector<ConnectionCommand> connectionCommands;

			AuxMatrix matrix;
			std::atomic_bool structuralChange;
			std::atomic_bool enabled;
			std::atomic_size_t maximumLatency;
		};
	}

#endif
