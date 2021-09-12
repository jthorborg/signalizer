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

	file:HostGraph.h

		Live & serializable represenation of the running instances connected in the 
		host graph.

*************************************************************************************/


#ifndef SIGNALIZER_HOSTGRAPH_H
	#define SIGNALIZER_HOSTGRAPH_H

	#include <cpl/Common.h>
	#include <cpl/state/Serialization.h>
	#include <cpl/CAudioStream.h>
	#include <mutex>
	#include <shared_mutex>
	#include <map>
	#include <set>
	#include "SignalizerConfiguration.h"

	namespace Signalizer
	{
		typedef std::make_signed<std::size_t>::type ssize_t;

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

			MixGraphListener(AudioStream& realtime, AudioStream& presentation, bool isFirst) 
				: realtime(realtime)
				, presentation(presentation)
				, structuralChange(false)
				, isFirst(isFirst)
			{
				setExpectedBufferSize(realtime.getInfo().anticipatedSize);
				auto link = std::make_unique<State>(*this, realtime);
				self = link.get();
				queue.emplace_back(std::make_pair(Command::Connect, 0));
				connectionCommands.emplace_back(std::move(link));
			}

			void connect(AudioStream& stream)
			{
				// mix graph listener can't be destroyed while this is happening.
				std::lock_guard<std::mutex> lock(connectDisconnectMutex);
				// this will start calling back but only operate on itself, 
				// since it hasn't been added to the graph yet
				queue.emplace_back(std::make_pair(Command::Connect, connectionCommands.size()));
				connectionCommands.emplace_back(std::make_unique<State>(*this, stream));
			}

			void disconnect(AudioStream& stream)
			{
				// mix graph listener can't be destroyed while this is happening.
				std::lock_guard<std::mutex> lock(connectDisconnectMutex);
				// it's fine source itself is already dead, we only operate on local copies.
				queue.emplace_back(std::make_pair(Command::Disconnect, disconnectionCommands.size()));

				disconnectionCommands.emplace_back(&stream);
			}

			void setExpectedBufferSize(std::size_t bufferSize)
			{
				maximumLatency.store(std::max<std::size_t>(128, bufferSize * 2), std::memory_order_release);
			}

			AudioStream& realtimeStream() { return realtime; }

		private:

			typedef cpl::CLIFOStream<AFloat> Buffer;

			struct State : AudioStream::Listener
			{
				struct Channel
				{
					cpl::CLIFOStream<AFloat> buffer;
					std::int16_t channelIndex {};
				};

				std::vector<Channel> channelQueues;
				std::atomic_size_t current {};
				std::int64_t globalPosition{};
				std::int64_t endpoint{};
				std::int64_t offset{};
				MixGraphListener& parent;
				AudioStream& source;

				State(MixGraphListener& parent, AudioStream& source)
					: parent(parent), source(source)
				{
					listenToSource(source, false, 1000);
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

			void handleStructuralChange(double sampleRate, std::size_t numSamples)
			{
				if (structuralChange)
				{
					std::size_t channels = 0;

					for (auto& g : graph)
					{
						channels += g.second->channelQueues.size();
					}

					matrix.resizeChannels(channels);

					auto info = presentation.getInfo();
					info.anticipatedChannels = channels;
					info.sampleRate = sampleRate;

					presentation.initializeInfo(info);

					structuralChange = false;
				}

				matrix.softBufferResize(numSamples);

			}

			void deliver(std::size_t numSamples)
			{
				std::unique_lock<std::shared_mutex> graphLayoutAndDataLock(dataMutex);

				handleStructuralChange(realtime.getInfo().sampleRate, numSamples);
				std::size_t channels = 0;

				const auto hostOrigin = self->endpoint;
				const auto hostSamples = self->current.load(std::memory_order_acquire);

				for (auto& g : graph)
				{
					auto& state = *g.second;
					auto currentG = state.current.load(std::memory_order_acquire);

					if (currentG != 0 && hostOrigin != 0)
					{
						const auto sampleDifference = static_cast<std::int64_t>(currentG) - static_cast<std::int64_t>(hostSamples);
						const auto tlDifference = state.endpoint - hostOrigin;
						auto difference = sampleDifference - tlDifference;
						//difference -= state.offset;

						if (difference != 0)
						{
							if (difference > 0)
							{
								// this dependency is ahead of the origin.
								currentG -= std::min(currentG, static_cast<std::size_t>(difference));
							}
							else
							{
								currentG += static_cast<std::size_t>(-difference);
							}

							//state.globalPosition += difference;
						}


					}



					for(std::size_t c = 0; c < state.channelQueues.size(); ++c, ++channels)
					{
						if (currentG >= numSamples || currentG != 0)
						{
							auto reader = g.second->channelQueues[c].buffer.createProxyView();
							reader.offset(-static_cast<ssize_t>(currentG));

							reader.copyFromHead(matrix[channels], numSamples);
						}
						else
						{
							// nothing to deliver... zero it out.
							matrix.clear(channels, 1);
						}

					}

					currentG -= std::min(currentG, numSamples);
					
					state.current.store(currentG, std::memory_order_release);
				}

				// TODO: Run this in place without async thread.
				presentation.processIncomingRTAudio(matrix.data(), channels, numSamples, realtime.getASyncPlayhead());
			}

			void asyncPropertiesChanged(State& s, const AudioStream& source)
			{
				if (&source == &realtime)
				{
					setExpectedBufferSize(realtime.getInfo().anticipatedSize);
					structuralChange = true;
				}
			}

			void onAsyncAudio(State& s, AFloat** buffer, std::size_t numChannels, std::size_t numSamples, std::int64_t globalPosition)
			{
				if (!isFirst)
					return;

				const auto localMaxLatency = maximumLatency.load(std::memory_order_acquire);
				const auto maxBufferSize = localMaxLatency * 8;
				{
					std::shared_lock<std::shared_mutex> lock(dataMutex);

					s.globalPosition = globalPosition;
					s.endpoint = globalPosition + numSamples;
					if (s.channelQueues.size() != numChannels)
					{
						s.channelQueues.resize(numChannels);
						structuralChange = true;
					}

					for (auto& q : s.channelQueues)
					{
						// clamp to max of 16 * buffersize.
						const auto neededSize = std::min(maxBufferSize, s.current + numSamples);
						if (q.buffer.getSize() < neededSize)
							q.buffer.setStorageRequirements(neededSize, cpl::Math::nextPow2(neededSize), true);
					}

					for (std::size_t i = 0; i < numChannels; ++i)
					{
						s
							.channelQueues[i]
							.buffer
							.createWriter()
							.copyIntoHead(buffer[i], numSamples);
					}

					const auto currentContained = s.current.load(std::memory_order_acquire) + numSamples;

					s.current.store(std::min(currentContained, maxBufferSize), std::memory_order_release);
				}


				if (&s == self)
				{
					decltype(queue) localQueue;
					decltype(connectionCommands) localNewConnections;
					decltype(disconnectionCommands) localNewDisconnections;

					{
						std::lock_guard<std::mutex> cdlock(connectDisconnectMutex);
						std::swap(localQueue, queue);
						std::swap(localNewConnections, connectionCommands);
						std::swap(localNewDisconnections, disconnectionCommands);
					}

					// eager delivery. only "we" can alter the graph, so we don't need protection.

					if(!localQueue.empty())
						structuralChange = true;

					for (auto command : localQueue)
					{
						if (command.first == Command::Connect)
						{
							auto source = &localNewConnections[command.second]->source;
							graph[source] = std::move(localNewConnections[command.second]);
						}
						else
						{
							graph.erase(localNewDisconnections[command.second]);
						}
					}

					if (graph.empty())
						return;

					const auto hostSamples = s.current.load(std::memory_order_acquire);
					auto min = hostSamples;

					if (min == 0)
					{
						CPL_RUNTIME_ASSERTION(numSamples == 0);
						return;
					}

					for (auto& g : graph)
					{
						const auto currentG = g.second->current.load(std::memory_order_acquire);

						if (currentG != 0)
						{
							// we can keep progressing.
							min = std::min(min, currentG);
						}
						else if (hostSamples > localMaxLatency)
						{
							// this link in the chain is still zero, but we are spilling over our max latency which is unexpected 
							// even for a reverse dependency.
							// is something on the way, at least?
							// TODO: requires shared_ptr
							if (g.second->source.getApproximateInFlightPackets() > 0)
							{
								// OK: we will get notified at a later stage.
								return;
							}
							else
							{
								// ignore this dependency - it's not processing for some reason
								continue;
							}
						}
						else
						{
							// Not spilling over, but there's globally not > 0 samples so just early out.
							return;
						}
					}

					if (min > 0)
					{
						deliver(min);
					}
				}
			}



			std::size_t initialSampleCapacity {};
			AudioStream& realtime;
			AudioStream& presentation;
			State* self;
			std::shared_mutex dataMutex;
			std::mutex connectDisconnectMutex;

			enum class Command { Connect, Disconnect };

			std::vector<std::pair<Command, std::size_t>> queue;
			std::vector<std::unique_ptr<State>> connectionCommands;
			std::vector<AudioStream*> disconnectionCommands;

			AuxMatrix matrix;
			std::map<const AudioStream*, std::unique_ptr<State>> graph;
			std::atomic_bool structuralChange;
			std::atomic_size_t maximumLatency;
			bool isFirst;
		};

		class HostGraph : public cpl::CSerializer::Serializable
		{
		public:

			HostGraph(AudioStream& realtime, AudioStream& presentation);
			~HostGraph();

			std::shared_ptr<juce::Component> createEditor();


			void serialize(cpl::CSerializer::Archiver& ar, cpl::Version version) override;
			void deserialize(cpl::CSerializer::Builder& ar, cpl::Version version) override;

		private:

			struct SerializedHandle
			{
				juce::Uuid Handle; // todo: create on demand

				friend inline bool operator < (const HostGraph::SerializedHandle& a, const HostGraph::SerializedHandle& b)
				{
					return std::memcmp(a.Handle.getRawData(), b.Handle.getRawData(), sizeof(juce::Uuid)) < 0;
				}
			};

			typedef HostGraph* HHandle;
			typedef std::int32_t PinInt;
			typedef std::pair<PinInt, SerializedHandle> SerializedEdge;
			typedef std::map<SerializedHandle, PinInt> ExpectedEdgeMap;
			static constexpr PinInt InvalidPin = -1;

			class Editor : public juce::Component, public juce::AsyncUpdater
			{
				// Inherited via AsyncUpdater
				virtual void handleAsyncUpdate() override;
			};

			SerializedHandle serializeReference(HHandle h);
			HostGraph* resolve(HHandle h);

			void broadcastCreate();
			void broadcastDestruct();

			void onNodeCreated(HHandle n);
			void onNodeDestroyed(HHandle n);


			static std::mutex staticMutex;
			static std::set<HHandle> staticSet;

			SerializedHandle nodeID;
			std::string name;
			std::vector<HHandle> pinInputs;
			ExpectedEdgeMap expectedEdges;
			std::weak_ptr<Editor> editor;
			MixGraphListener mix;
			bool isFirst;
		};

	}

#endif
