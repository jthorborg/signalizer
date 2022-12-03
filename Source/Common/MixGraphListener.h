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
	#include <memory>
	#include "ConcurrentConfig.h"
	#include "CommonSignalizer.h"

	namespace Signalizer
	{
		class AudioProcessor;

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

			void clear()
			{
				std::memset(auxData.data(), 0, auxData.size() * sizeof(float));
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

		
		class MixGraphListener : public AudioStream::Listener, public std::enable_shared_from_this<MixGraphListener>
		{
		public:
			friend class HostGraph;

			/// <summary>
			/// While in scope, keeps a mix graph alive pumping from inputs specified from a host graph into a presentation output.
			/// </summary>
			class Handle
			{
				friend class MixGraphListener;
				friend class HostGraph;
				std::shared_ptr<MixGraphListener> listener;

				Handle(std::shared_ptr<MixGraphListener> ret) : listener(std::move(ret)) {}

			public:

				Handle() = default;
				Handle(Handle&&) = default;
				Handle& operator = (Handle&&) = default;

				~Handle()
				{
					if (listener)
						listener->close();
				}
			};


			static std::pair<Handle, std::shared_ptr<AudioStream::Output>> create(AudioProcessor& processor);

			void connect(std::shared_ptr<AudioStream::Output>& stream, DirectedPortPair pair, const std::string& name);
			void disconnect(std::shared_ptr<AudioStream::Output>& stream, DirectedPortPair pair);


			~MixGraphListener();

		private:

			void close();

			MixGraphListener::MixGraphListener(AudioProcessor& p, AudioStream::IO&& presentation);

			typedef cpl::CLIFOStream<AFloat> Buffer;

			struct ConnectionCommand
			{
				std::shared_ptr<AudioStream::Output> stream;
				std::string name;
				DirectedPortPair pair;
				bool isConnection;
			};

			struct State
			{
				struct Channel
				{
					cpl::CLIFOStream<AFloat> buffer;
					std::string originName;
				};

				std::map<DirectedPortPair, Channel> channelQueues;
				std::atomic_size_t current {};
				std::int64_t globalPosition{};
				std::int64_t endpoint{};
				std::int64_t offset{};
				std::int32_t refCount {};

				std::weak_ptr<AudioStream::Output> source;

				State(const State& other) = delete;
				State(State&& other) noexcept
					: channelQueues(std::move(other.channelQueues))
					, current(other.current.load(std::memory_order_acquire))
					, globalPosition(other.globalPosition)
					, endpoint(other.endpoint)
					, offset(other.offset)
					, refCount(other.refCount)
					, source(std::move(other.source))
				{
				
				}
				State& operator = (State&& other) = delete;

				static State empty() { return {}; }

			private:
				State() {}
			};

			auto emplace(std::shared_ptr<AudioStream::Output>& stream);
			void remove(std::map<AudioStream::Handle, State>::iterator position);

			void onStreamPropertiesChanged(AudioStream::ListenerContext& changedSource, const AudioStream::AudioStreamInfo& before) override final;
			void onStreamAudio(AudioStream::ListenerContext& source, AFloat** buffer, std::size_t numChannels, std::size_t numSamples) override final;
			void onStreamDied(AudioStream::ListenerContext& dyingSource) override final;

			void handleStructuralChange(AudioStream::ListenerContext&, std::size_t numSamples);
			void deliver(AudioStream::ListenerContext& ctx, std::size_t numSamples);

			void updateTopologyCommands();
			void assignSelf();

			int id;
			std::map<AudioStream::Handle, State> graph;
			std::size_t initialSampleCapacity {};
			std::shared_ptr<AudioStream::Output> realtime;
			AudioStream::Input presentationInput;
			AudioStream::Handle presentationOutput;
			std::weak_ptr<AudioStream::Output> weakPresentationOutput;
			State* self;
			std::shared_mutex dataMutex;

			std::mutex connectDisconnectMutex;
			std::vector<ConnectionCommand> connectionCommands;

			AuxMatrix matrix;
			std::atomic_bool structuralChange;
			std::atomic_bool enabled;
			cpl::weak_atomic<std::size_t> maximumLatency;
			// TODO: should be shared?
			ConcurrentConfig& concurrentConfig;
		};
	}

#endif
