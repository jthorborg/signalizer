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
#include "../Processor/PluginProcessor.h"

namespace Signalizer
{
	static std::atomic_int mgCounter;

	auto MixGraphListener::emplace(std::shared_ptr<AudioStream::Output>& stream)
	{
		auto its = graph.emplace(
			std::piecewise_construct,
			std::forward_as_tuple(stream->getHandle()),
			std::forward_as_tuple( State::empty())
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

	void MixGraphListener::remove(std::map<AudioStream::Handle, State>::iterator position)
	{
		CPL_RUNTIME_ASSERTION(position->second.refCount == 0);
		if (auto sh = position->second.source.lock())
		{
			sh->removeListener(shared_from_this());
		}

		graph.erase(position);
	}


	void MixGraphListener::close()
	{
		std::unique_lock<std::shared_mutex> graphLayoutAndDataLock(dataMutex);

		// there might be further processing, but from now on we're in a dead state.
		enabled = false;

		for (auto& g : graph)
		{
			if (auto sh = g.second.source.lock())
				sh->removeListener(shared_from_this());
		}

		if (auto sh = weakPresentationOutput.lock())
			sh->removeListener(shared_from_this());

		graph.clear();
	}

	void MixGraphListener::assignSelf()
	{
		self = &emplace(realtime)->second;
		// always keep ourselves alive.
		self->refCount++;

		if(auto sh = weakPresentationOutput.lock())
			sh->addListener(shared_from_this());
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
		else if (changedSource.getHandle() == presentationOutput)
		{
			concurrentConfig.historyCapacity = info.audioHistoryCapacity;
			concurrentConfig.historySize = info.audioHistorySize;
		}
	}

	void MixGraphListener::onStreamDied(AudioStream::ListenerContext& dyingSource)
	{
	}

	MixGraphListener::MixGraphListener(AudioProcessor& p, AudioStream::IO&& presentation)
		: realtime(p.getRealtimeOutput())
		, presentationInput(std::move(std::get<0>(presentation)))
		, weakPresentationOutput(std::get<1>(presentation))
		, presentationOutput(std::get<1>(presentation)->getHandle())
		, structuralChange(false)
		, enabled(true)
		, self(nullptr)
		, concurrentConfig(*p.config)
		, id(mgCounter.fetch_add(1))
	{
		auto& output = std::get<1>(presentation);

		output->modifyConsumerInfo(
			[&](auto&& info)
			{
				// Here we restore initial parameters recovered from serialized storage in the engine,
				// previously (< 0.3.5) this was only recovered once the editor appears.
				info.storeAudioHistory = true;
				info.audioHistorySize = concurrentConfig.historySize;
				info.audioHistoryCapacity = concurrentConfig.historyCapacity;
			}
		);
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

	std::pair<MixGraphListener::Handle, std::shared_ptr<AudioStream::Output>> MixGraphListener::create(AudioProcessor& p)
	{
		auto io = AudioStream::create(false);
		auto presentationOutput = std::get<1>(io);
		auto mixGraph = std::shared_ptr<MixGraphListener>(new MixGraphListener(p, std::move(io)));
		mixGraph->assignSelf();
		return { Handle(mixGraph), presentationOutput };
	}

	void MixGraphListener::connect(std::shared_ptr<AudioStream::Output>& other, DirectedPortPair pair, const std::string& name)
	{
		// mix graph listener can't be destroyed while this is happening.
		std::lock_guard<std::mutex> lock(connectDisconnectMutex);
		auto portName = name + "[" + std::to_string(pair.Source) + "]";

		connectionCommands.emplace_back(ConnectionCommand{ other, std::move(portName), pair, true});
	}

	void MixGraphListener::disconnect(std::shared_ptr<AudioStream::Output>& other, DirectedPortPair pair)
	{
		// mix graph listener can't be destroyed while this is happening.
		std::lock_guard<std::mutex> lock(connectDisconnectMutex);
		// it's fine source itself is already dead, we only operate on local copies.
		connectionCommands.emplace_back(ConnectionCommand{ other, {}, pair, false });
	}

	std::size_t MixGraphListener::reportLatency() const noexcept
	{
		return currentLatency;
	}

	bool MixGraphListener::reportSynchronized() const noexcept
	{
		return isSynchronized;
	}

	void MixGraphListener::handleStructuralChange(AudioStream::ListenerContext& ctx, std::size_t numSamples, std::unique_lock<std::shared_mutex>& lock)
	{
		auto& realInfo = ctx.getInfo();

		if (structuralChange)
		{
			PinInt maxDestinationPort = -1;
			std::int64_t setPorts = 0;
			for (auto& g : graph)
			{
				for (auto& entry : g.second.channelQueues)
				{
					maxDestinationPort = std::max(maxDestinationPort, entry.first.Destination);
					setPorts |= 1ll << (std::int64_t)entry.first.Destination;

					auto nameCopy = entry.second.originName;
					presentationInput.enqueueChannelName(entry.first.Destination, std::move(nameCopy));
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


			for (std::int64_t i = 0; i < maxDestinationPort; ++i)
			{
				if ((setPorts & (1ll << i)) == 0)
				{
					presentationInput.enqueueChannelName(i, "nothing");
				}
			}


			structuralChange = false;
		}

		matrix.softBufferResize(numSamples);
	}

	void MixGraphListener::deliver(AudioStream::ListenerContext& ctx, std::size_t numSamples)
	{
		std::unique_lock<std::shared_mutex> graphLayoutAndDataLock(dataMutex);

		handleStructuralChange(ctx, numSamples, graphLayoutAndDataLock);

		// no channels to show
		if (matrix.size() < 1)
			return;

		// clear the matrix for additive / empty slots - might not be needed?
		matrix.clear();

		const auto hostEndpoint = self->endpoint.load();
		const auto hostSamples = static_cast<std::int64_t>(self->containedSamples.load());

		bool seeminglySynchronized = true;

		for (auto& g : graph)
		{
			auto& state = g.second;
			auto containedInState = state.containedSamples.load();

			if (containedInState != 0 && ctx.getPlayhead().isPlaying())
			{
				const auto sampleDifference = containedInState - hostSamples;
				const auto tlDifference = state.endpoint - hostEndpoint;
				const auto difference = sampleDifference - tlDifference;

				const auto hasDiscontinuity = self->discontinuity || state.discontinuity;

				if (difference != 0)
				{
					seeminglySynchronized = false;
					// twice we had the same delay? try to adjust
					if (!hasDiscontinuity && state.typicalOffset == difference)
					{
						if (difference > 0)
						{
							// erase history, we are too far behind...
							containedInState -= std::min(containedInState, difference);
							state.containedSamples = containedInState;
						}
						else
						{
							// add history, insert silence...
							for (auto& q : state.channelQueues)
							{
								q.second.buffer.createWriter().advance(static_cast<std::size_t>(-difference));
							}

							containedInState += -difference;
							state.containedSamples = containedInState;
						}

						state.typicalOffset = 0;
					}
					else
					{
						state.typicalOffset = difference;
					}
				}
			}

			// if we don't have enough to deliver, do nothing (matrix cleared at beginning)
			if (containedInState >= static_cast<std::int64_t>(numSamples) || containedInState != 0)
			{
				for (auto& q : state.channelQueues)
				{
					auto reader = q.second.buffer.createProxyView();
					reader.offset(-containedInState);

					reader.copyFromHead<true>(matrix[q.first.Destination], numSamples);
				}
			}

			// seems ultra wrong but makes them sync
			// something to do with this exceeding max buffer size back in onStreamAudio
			auto current = state.containedSamples.load();
			auto toSubtract = std::min(current, static_cast<std::int64_t>(numSamples));
			state.containedSamples.store(current - toSubtract);
		}

		this->isSynchronized = seeminglySynchronized;

		// TODO: don't copy into a matrix, rig a provider from the read heads instead?
		presentationInput.processIncomingRTAudio(matrix.data(), matrix.size(), numSamples, ctx.getPlayhead());
	}

	void MixGraphListener::onStreamAudio(AudioStream::ListenerContext& ctx, AFloat** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		const auto&& handle = ctx.getHandle();

		// We listen to the presentation output to get synchronized callbacks for changes in the stream properties.
		// We however also get callbacks for audio (that we emitted ourselves), which is a deadlock we shall ignore.
		if (handle == presentationOutput)
			return;

		const auto localMaxLatency = maximumLatency.load();
		const auto maxBufferSize = localMaxLatency * 8;
		const auto globalPosition = ctx.getPlayhead().getPositionInSamples();

		const bool isSelf = handle == realtime->getHandle();

		if (enabled)
		{
			std::shared_lock<std::shared_mutex> lock(dataMutex);

			// certain conditions can cause callbacks to temporarily appear, even though we deregistrered from this source and no longer know it.
			auto it = graph.find(handle);
			if (it == graph.end())
				return;

			auto& s = it->second;

			s.discontinuity = s.endpoint - globalPosition;

			s.globalPosition = globalPosition;
			s.endpoint = globalPosition + numSamples;

			const auto neededSize = std::min<std::size_t>(maxBufferSize, s.containedSamples + numSamples);

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

			const auto currentContained = s.containedSamples.load() + numSamples;
			s.containedSamples.store(std::min<std::int64_t>(currentContained, maxBufferSize));
		}


		if (isSelf)
		{
			updateTopologyCommands();

			if (!enabled || graph.empty())
				return;

			const auto hostSamples = self->containedSamples.load();
			const auto hostOrigin = self->endpoint - hostSamples;

			auto min = hostSamples;
			auto latency = min;
			
			if (min == 0)
			{
				// Happens when changing self (resetting the enqueued count)
				return;
			}

			for (auto& g : graph)
			{
				if (min <= 0)
					break;

				auto containedInState = g.second.containedSamples.load();
				latency = std::max(containedInState, latency);

				// when playing, we want to clamp the available amount of samples with respect to the alignment
				// of the timeline compared to the host. this ensures the inner chomping doesn't play "catch up"
				if (ctx.getPlayhead().isPlaying() && !self->discontinuity && !g.second.discontinuity)
				{
					const auto available = g.second.endpoint - hostOrigin;

					if(available < containedInState)
						containedInState = std::min(containedInState, std::max(0ll, available));
				}

				if (containedInState > 0)
				{
					// we can keep progressing.
					min = std::min(min, containedInState);
				}
				else if (hostSamples > static_cast<std::int64_t>(localMaxLatency))
				{
					// this link in the chain is still zero, but we are spilling over our max latency which is unexpected 
					// even for a reverse dependency.
					// is something on the way, at least?

					if(auto sh = g.second.source.lock(); sh && sh->getApproximateInFlightPackets() > 0)
					{
						// OK: we will get notified at a later stage.
						return;
					}
					else
					{
						// super special case for *current having been updated inbetween these statements (easy to provoke with a debugger)
						const auto newAvailable = g.second.containedSamples.load();
						if (newAvailable == containedInState)
						{
							// ignore this dependency - it's not processing yet for some reason or deleted itself
							continue;
						}
						else
						{
							min = std::min(min, newAvailable);
							continue;
						}
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
				// make sure when we're processing more than a stereo pair to at least present a sizable buffer,
				// since multithreading is going to be involved.
				constexpr int minMultiChannel = 64;

				if (matrix.size() > 2 && min < minMultiChannel)
					return;

				currentLatency = latency;
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

		isSynchronized = false;

		// so we can alter mapping tables. can be done without.
		std::unique_lock<std::shared_mutex> lock(dataMutex);

		for (auto& command : localNewToplogy)
		{
			const auto key = command.stream.get();
			auto it = graph.find(key->getHandle());

			if (it == graph.end())
			{
				CPL_RUNTIME_ASSERTION(command.isConnection);
				CPL_RUNTIME_ASSERTION(command.stream.get() != nullptr);
				CPL_RUNTIME_ASSERTION(command.stream.get() != self->source.lock().get());
				it = emplace(std::move(command.stream));
			}

			auto& source = it->second;

			if (command.isConnection)
			{
				source.refCount++;

				// also ensures it exists.
				source.channelQueues[command.pair].originName = std::move(command.name); 

				// TODO: Pauses other channels from this source (and everything)
				source.containedSamples = 0;
				
			}
			else
			{
				source.channelQueues.erase(command.pair);

				if (--source.refCount == 0)
					remove(it);
			}
		}
	}

}
