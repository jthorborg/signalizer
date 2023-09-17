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
	#include <cpl/Core.h>
	#include <cpl/state/Serialization.h>
	#include <cpl/dsp/ChannelMatrix.h>
	#include <mutex>
	#include <cpl/shared_mutex.h>
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

				void setSuspended(bool shouldBeSuspended)
				{
					listener->enabled = !shouldBeSuspended;
				}

				~Handle()
				{
					if (listener)
						listener->close();
				}
			};


			static std::pair<Handle, std::shared_ptr<AudioStream::Output>> create(AudioProcessor& processor);

			void connect(std::shared_ptr<AudioStream::Output>& stream, DirectedPortPair pair, const std::string& name);
			void disconnect(std::shared_ptr<AudioStream::Output>& stream, DirectedPortPair pair);
			std::size_t reportLatency() const noexcept;
			bool reportSynchronized() const noexcept;

			~MixGraphListener();

		private:

			void close();

			MixGraphListener(AudioProcessor& p, AudioStream::IO&& presentation);

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
				/// <summary>
				/// The amount of samples contained leading up to <see cref="endpoint"/>
				/// </summary>
				cpl::weak_atomic<std::int64_t> containedSamples {};
				/// <summary>
				/// Reported discontinuity, ie. if the last <see cref="globalPosition"/> didn't update the current
				/// <see cref="endpoint"/>
				/// </summary>
				std::int64_t discontinuity{};
				/// <summary>
				/// The last timing difference between this state and the host.
				/// </summary>
				std::int64_t typicalOffset{};
				/// <summary>
				/// The last reported time stamp of this current graph (not including the samples leading up to the endpoint).
				/// </summary>
				std::int64_t globalPosition{};
				/// <summary>
				/// The end of this current time line, ie. the global position plus the amount of samples this state received
				/// last time around.
				/// </summary>
				cpl::weak_atomic<std::int64_t> endpoint{};
				std::int32_t refCount {};

				std::weak_ptr<AudioStream::Output> source;

				State(const State& other) = delete;
				State(State&& other) noexcept
					: channelQueues(std::move(other.channelQueues))
					, containedSamples(other.containedSamples)
					, globalPosition(other.globalPosition)
					, endpoint(other.endpoint)
					, refCount(other.refCount)
					, source(std::move(other.source))
				{
				
				}
				State& operator = (State&& other) = delete;

				static State empty() { return {}; }

			private:
				State() {}
			};

			auto emplace(std::shared_ptr<AudioStream::Output> stream);
			void remove(std::map<AudioStream::Handle, State>::iterator position);

			void onStreamPropertiesChanged(AudioStream::ListenerContext& changedSource, const AudioStream::AudioStreamInfo& before) override final;
			void onStreamAudio(AudioStream::ListenerContext& source, AFloat** buffer, std::size_t numChannels, std::size_t numSamples) override final;
			void onStreamDied(AudioStream::ListenerContext& dyingSource) override final;

			void handleStructuralChange(AudioStream::ListenerContext&, std::size_t numSamples, cpl::unique_lock<cpl::shared_mutex>& lock);
			void deliver(AudioStream::ListenerContext& ctx, std::size_t numSamples);

			void updateTopologyCommands();
			void assignSelf();

			const int id;
			std::map<AudioStream::Handle, State> graph;
			std::size_t initialSampleCapacity {};
			std::shared_ptr<AudioStream::Output> realtime;
			AudioStream::Input presentationInput;
			AudioStream::Handle presentationOutput;
			std::weak_ptr<AudioStream::Output> weakPresentationOutput;
			State* self;
			cpl::shared_mutex dataMutex;

			std::mutex connectDisconnectMutex;
			std::vector<ConnectionCommand> connectionCommands;

			cpl::ChannelMatrix<AudioStream::DataType> matrix;
			std::atomic_bool structuralChange;
			std::atomic_bool enabled;
			cpl::weak_atomic<std::size_t> maximumLatency;
			cpl::relaxed_atomic<std::size_t> currentLatency;
			cpl::relaxed_atomic<bool> isSynchronized;

			ConcurrentConfig& concurrentConfig;

			JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MixGraphListener);

		};
	}

#endif
