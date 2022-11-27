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

	file:MixGraphListener.cpp

		Implementation of MixGraphListener.h

*************************************************************************************/

#include "MixGraphListener.h"
#include <atomic>
#include <memory>

namespace Signalizer
{

	auto MixGraphListener::emplace(std::shared_ptr<AudioStream::Output>& stream)
	{
		auto its = graph.emplace(
			std::piecewise_construct,
			std::forward_as_tuple(stream->getHandle()),
			std::forward_as_tuple( State{})
		);
		if (its.second)
		{
			its.first->second.source = stream;
			stream->addListener(shared_from_this());
		}
		else
		{
			CPL_RUNTIME_ASSERTION(false && "Graph already contained stream?");
		}

		return its.first;
	}

	void MixGraphListener::close()
	{
		std::unique_lock<std::shared_mutex> graphLayoutAndDataLock(dataMutex);

		// there might be further processing, but from now on we're in a dead state.
		enabled = false;

		for (auto& g : graph)
		{
			g.second.source->removeListener(shared_from_this());
		}

		graph.clear();
	}

	void MixGraphListener::assignSelf()
	{
		self = &emplace(realtime)->second;
		presentationOutput->addListener(shared_from_this());
	}

	std::shared_ptr<AudioStream::Output>& MixGraphListener::getPresentationOutput()
	{
		return presentationOutput;
	}


	const ConcurrentConfig& MixGraphListener::getConcurrentConfig() const noexcept
	{
		return concurrentConfig;
	}


	void MixGraphListener::onStreamPropertiesChanged(AudioStream::ListenerContext& changedSource, const AudioStream::AudioStreamInfo&)
	{
		const auto& info = changedSource.getInfo();

		if (changedSource.getHandle() == realtime->getHandle())
		{
			maximumLatency = std::max<std::size_t>(128, info.anticipatedSize * 2);
			structuralChange = true;
			concurrentConfig.bpm = changedSource.getPlayhead().getBPM();
			concurrentConfig.numChannels = info.channels;
			concurrentConfig.sampleRate = info.sampleRate;
		}
		else if (changedSource.getHandle() == presentationOutput->getHandle())
		{
			concurrentConfig.historyCapacity = info.audioHistoryCapacity;
			concurrentConfig.historySize = info.audioHistorySize;
		}
	}

	void MixGraphListener::onStreamDied(AudioStream::ListenerContext& dyingSource)
	{
	}

	MixGraphListener::MixGraphListener(std::shared_ptr<AudioStream::Output> realtimeStream, AudioStream::IO&& presentation)
		: realtime(std::move(realtimeStream))
		, presentationInput(std::move(std::get<0>(presentation)))
		, presentationOutput(std::move(std::get<1>(presentation)))
		, structuralChange(false)
		, enabled(true)
		, self(nullptr)
	{
	}

	MixGraphListener::~MixGraphListener()
	{
		// solely here as a checkpoint.
		if (graph.size() != 0)
		{
			CPL_DEBUGOUT("Releasing mix graph listener that hasn't been cleaned up (but it isn't really a problem)");
			CPL_BREAKIFDEBUGGED();
		}
	}

	std::shared_ptr<MixGraphListener> MixGraphListener::create(std::shared_ptr<AudioStream::Output> realtimeStream)
	{
		auto io = AudioStream::create(false);

		auto mixGraph = std::shared_ptr<MixGraphListener>(new MixGraphListener(std::move(realtimeStream), std::move(io)));
		mixGraph->assignSelf();
		return mixGraph;
	}

	void MixGraphListener::connect(std::shared_ptr<AudioStream::Output>& other, DirectedPortPair pair)
	{
		// mix graph listener can't be destroyed while this is happening.
		std::lock_guard<std::mutex> lock(connectDisconnectMutex);
		// this will start calling back but only operate on itself, 
		// since it hasn't been added to the graph yet
		connectionCommands.emplace_back(ConnectionCommand{ other, pair, true });
	}

	void MixGraphListener::disconnect(std::shared_ptr<AudioStream::Output>& other, DirectedPortPair pair)
	{
		// mix graph listener can't be destroyed while this is happening.
		std::lock_guard<std::mutex> lock(connectDisconnectMutex);
		// it's fine source itself is already dead, we only operate on local copies.
		connectionCommands.emplace_back(ConnectionCommand{ other, pair, false });
	}

	void MixGraphListener::handleStructuralChange(AudioStream::ListenerContext& ctx, std::size_t numSamples)
	{
		auto& realInfo = ctx.getInfo();

		if (structuralChange)
		{
			PinInt maxDestinationPort = -1;

			for (auto& g : graph)
			{
				for (auto& entry : g.second.channelQueues)
				{
					maxDestinationPort = std::max(maxDestinationPort, entry.first.Destination);
				}
			}

			// from first-indexing to zero-based
			maxDestinationPort++;

			// round up to next multiple of two, so there's always pairs of channels
			maxDestinationPort += maxDestinationPort % 2;

			matrix.resizeChannels(maxDestinationPort);

			presentationInput.initializeInfo(
				[&](AudioStream::ProducerInfo & info)
				{
					info.channels = maxDestinationPort;
					info.sampleRate = realInfo.sampleRate;
					info.anticipatedSize = static_cast<std::uint32_t>(numSamples);
				}
			);

			presentationOutput->modifyConsumerInfo(
				[&](AudioStream::ConsumerInfo& info)
				{
					info.storeAudioHistory = true;
				}
			);


			structuralChange = false;
		}

		matrix.softBufferResize(numSamples);
	}

	void MixGraphListener::deliver(AudioStream::ListenerContext& ctx, std::size_t numSamples)
	{
		std::unique_lock<std::shared_mutex> graphLayoutAndDataLock(dataMutex);

		handleStructuralChange(ctx, numSamples);

		// no channels to show
		if (matrix.size() < 1)
			return;

		// clear the matrix for additive / empty slots - might not be needed?
		matrix.clear();

		const auto hostOrigin = self->endpoint;
		const auto hostSamples = self->current.load(std::memory_order_acquire);

		for (auto& g : graph)
		{
			auto& state = g.second;
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

			for (auto& q : state.channelQueues)
			{
				if (currentG >= numSamples || currentG != 0)
				{
					auto reader = q.second.buffer.createProxyView();
					reader.offset(-static_cast<ssize_t>(currentG));

					// TODO: should be additive
					reader.copyFromHead(matrix[q.first.Destination], numSamples);
				}
				else
				{
					// nothing to deliver... zero it out. TODO don't?
					matrix.clear(q.first.Destination, 1);
				}
			}

			currentG -= std::min(currentG, numSamples);
			state.current.store(currentG, std::memory_order_release);
		}

		// TODO: don't copy into a matrix, rig a provider from the read heads instead
		presentationInput.processIncomingRTAudio(matrix.data(), matrix.size(), numSamples, ctx.getPlayhead());
	}

	void MixGraphListener::onStreamAudio(AudioStream::ListenerContext& ctx, AFloat** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		const auto&& handle = ctx.getHandle();

		// We listen to the presentation output to get synchronized callbacks for changes in the stream properties.
		// We however also get callbacks for audio (that we emitted ourselves), which is a deadlock we shall ignore.
		if (handle == presentationOutput->getHandle())
			return;

		const auto localMaxLatency = maximumLatency.load();
		const auto maxBufferSize = localMaxLatency * 8;
		const auto globalPosition = ctx.getPlayhead().getPositionInSamples();

		if (enabled)
		{
			std::shared_lock<std::shared_mutex> lock(dataMutex);

			auto& s = graph[handle];

			s.globalPosition = globalPosition;
			s.endpoint = globalPosition + numSamples;

			const auto neededSize = std::min(maxBufferSize, s.current + numSamples);

			for (auto& q : s.channelQueues)
			{
				auto& channelQ = q.second.buffer;
				// clamp to max of 16 * buffersize.
				if (channelQ.getSize() < neededSize)
					channelQ.setStorageRequirements(neededSize, cpl::Math::nextPow2(neededSize), true);

				channelQ
					.createWriter()
					.copyIntoHead(buffer[q.first.Source], numSamples);
			}

			const auto currentContained = s.current.load(std::memory_order_acquire) + numSamples;
			s.current.store(std::min(currentContained, maxBufferSize), std::memory_order_release);
		}


		if (handle == realtime->getHandle())
		{
			updateTopologyCommands();

			if (!enabled || graph.empty())
				return;

			const auto hostSamples = self->current.load(std::memory_order_acquire);
			auto min = hostSamples;

			if (min == 0)
			{
				// TODO: Happens when changing self (resetting the enqueued count)
				/* CPL_RUNTIME_ASSERTION(numSamples == 0); */
				return;
			}

			for (auto& g : graph)
			{
				const auto currentG = g.second.current.load(std::memory_order_acquire);

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
					if (g.second.source->getApproximateInFlightPackets() > 0)
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
				deliver(ctx, min);
			}
		}
	}

	void MixGraphListener::updateTopologyCommands()
	{
		decltype(connectionCommands) localNewToplogy;

		{
			std::lock_guard<std::mutex> cdlock(connectDisconnectMutex);
			std::swap(localNewToplogy, connectionCommands);
		}

		// eager delivery. only "we" can alter the graph, so we don't need protection.
		if (!localNewToplogy.empty())
			structuralChange = true;
		else
			return;

		// so we can alter mapping tables. can be done without.
		std::unique_lock<std::shared_mutex> lock(dataMutex);

		for (auto& command : localNewToplogy)
		{
			const auto key = command.stream.get();
			auto it = graph.find(key->getHandle());

			if (it == graph.end())
			{
				CPL_RUNTIME_ASSERTION(command.isConnection);
				CPL_RUNTIME_ASSERTION(command.stream.get() != self->source.get());
				it = emplace(std::move(command.stream));
			}

			auto& source = it->second;

			if (command.isConnection)
			{
				source.refCount++;

				source.channelQueues[command.pair]; // just ensure it exists.

				// TODO: Pauses other channels from this source (and everything)
				source.current = 0;
			}
			else
			{
				source.channelQueues.erase(command.pair);

				if (--source.refCount == 0)
					graph.erase(it);
			}
		}
	}

}
