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

			MixGraphListener(AudioStream& hh, std::size_t bufferSize) : host(hh), bufferSize(bufferSize), structuralChange(false)
			{
				connect(host);
				self = graph[&hh].get();
			}

			void setBufferSize(std::size_t samples)
			{
				std::lock_guard<std::mutex> lock(mutex);
				bufferSize = samples;
			}

			void connect(AudioStream& stream)
			{
				// can deadlock with audio thread.
				std::lock_guard<std::mutex> lock(mutex);
				graph[&stream] = std::make_unique<State>(*this, bufferSize);
				structuralChange = true;
			}

			AudioStream& localStream() { return host; }

		private:

			typedef cpl::CLIFOStream<AFloat> Buffer;

			struct State : AudioStream::Listener
			{
				struct Channel
				{
					cpl::CLIFOStream<AFloat> buffer;
					std::int16_t channelIndex;
				};

				std::vector<Channel> channelQueues;
				std::size_t current{}, initialSize;
				std::int64_t globalPosition{};
				MixGraphListener& parent;

				State(MixGraphListener& parent, std::size_t initialSize)
					: parent(parent), initialSize(initialSize)
				{
					listenToSource(parent.host, false, 1000);
				}


				virtual bool onAsyncAudio(const AudioStream& source, AFloat** buffer, std::size_t numChannels, std::size_t numSamples) override
				{
					globalPosition = source.getASyncPlayhead().getPositionInSamples();
					if (channelQueues.size() != numChannels)
					{
						channelQueues.resize(numChannels);
					}

					for(auto & q : channelQueues)
					{
						const auto neededSize = current + numSamples;
						if (q.buffer.getSize() < neededSize)
							q.buffer.setStorageRequirements(neededSize, cpl::Math::nextPow2(neededSize), true);
					}

					return parent.onAsyncAudio(*this, buffer, numChannels, numSamples);
				}

				~State()
				{
					detachFromSource();
				}
			};

			void handleStructuralChange(std::size_t numSamples)
			{
				if (structuralChange)
				{
					std::size_t channels = 0;

					for (auto& g : graph)
					{
						channels += g.second->channelQueues.size();
					}

					matrix.resizeChannels(channels);

					structuralChange = false;
				}

				matrix.softBufferResize(numSamples);

			}

			void deliver(std::size_t numSamples)
			{
				handleStructuralChange(numSamples);
				std::size_t channels = 0;

				for (auto& g : graph)
				{
					auto& state = *g.second;

					for(std::size_t c = 0; c < state.channelQueues.size(); ++c, ++channels)
					{
						auto reader = g.second->channelQueues[c].buffer.createProxyView();
						reader.offset(-static_cast<ssize_t>(state.current));

						reader.copyFromHead(matrix[channels], numSamples);
					}
					
					state.current -= numSamples;
				}
			}


			bool onAsyncAudio(State& s, AFloat** buffer, std::size_t numChannels, std::size_t numSamples)
			{
				std::lock_guard<std::mutex> lock(mutex);

				if (graph.size() == 1)
					return true;
				
				for (std::size_t i = 0; i < numChannels; ++i)
				{
					s
						.channelQueues[i]
						.buffer
						.createWriter()
						.copyIntoHead(buffer[i], numSamples);
				}

				s.current += numSamples;

				if (&s == self)
				{
					// eager delivery.

					auto min = s.current;

					for (auto& g : graph)
					{
						// issue
						min = std::min(min, g.second->current);
					}

					if (min > 0)
					{
						deliver(min);
					}
				}

				return false;
			}



			std::size_t bufferSize{};
			AudioStream& host;
			State* self;
			std::mutex mutex;
			AuxMatrix matrix;
			std::map<const AudioStream*, std::unique_ptr<State>> graph;
			bool structuralChange;
		};

		class HostGraph : public cpl::CSerializer::Serializable
		{
		public:

			HostGraph(AudioStream& stream);
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
