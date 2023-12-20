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

	file:HostGraph.cpp

		Implementation of HostGraph.h

*************************************************************************************/

#include "HostGraph.h"
#include <atomic>
#include "MixGraphListener.h"

namespace Signalizer
{
	std::mutex HostGraph::staticMutex;
	std::set<HostGraph::HHandle> HostGraph::staticSet;
	static std::atomic_int globalVersion;

	constexpr const char* kSerializationControlKey = "serialization-control";
	constexpr const char* kTopologyKey = "topology-data";

	HostGraph::HostGraph(std::shared_ptr<AudioStream::Output> realtimeOutput)
		: name("unnamed")
		, realtime(std::move(realtimeOutput))
	{
		broadcastCreate(GraphLock(staticMutex));
	}


	HostGraph::~HostGraph()
	{
		const GraphLock lock(staticMutex);
		broadcastDestruct(lock);
		resurrectNextAlias(lock);
	}

	void HostGraph::addModelListener(std::weak_ptr<juce::AsyncUpdater> shared)
	{
		modelChangedCallback = std::move(shared);
	}

	void HostGraph::serialize(cpl::CSerializer::Archiver& ar, cpl::Version version)
	{
		GraphLock lock(staticMutex);

		ar[kSerializationControlKey] << serializationControl;

		ar << name;

		if (!serializationControl.shouldSerializeGraph())
			return;

		auto serializeTopology = [this](cpl::CSerializer::Builder& serializer)
		{
			serializer << cpl::CSerializer::OptionalWrapper(nodeID);

			std::uint32_t counter = std::accumulate(
				topology.begin(),
				topology.end(),
				0u,
				[&](const std::uint32_t& c, const auto& pair) { return c + static_cast<std::uint32_t>(pair.second.inputs.size()); }
			);

			serializer << counter;

			for (const auto& pairRelation : topology)
			{
				for (const auto& pair : pairRelation.second.inputs)
				{
					serializer << std::make_pair(pairRelation.first, pair);
				}
			}
		};

		serializeTopology(ar[kTopologyKey]);
	}

	void HostGraph::deserialize(cpl::CSerializer::Builder& builder, cpl::Version version)
	{		
		TriggerModelUpdateOnExit exit{ this };
		std::lock_guard<std::mutex> lock(staticMutex);

		auto& serializedControl = builder[kSerializationControlKey];
		if(!serializedControl.isEmpty())
			serializedControl >> serializationControl;

		auto oldName = name;
		builder >> name;

		// Prior to version 0.4.2, topology was stored inline at this offset.
		if (version >= cpl::Version(0, 4, 2))
		{
			// And this is how we conditionally load it now.
			auto& serializedTopology = builder[kTopologyKey];
			if (!serializedTopology.isEmpty())
				deserializeTopology(serializedTopology, version, lock);
		}
		else
		{		
			// Before, we would always store/restore the data, so mimic this behaviour.
			deserializeTopology(builder, version, lock);
		}

		if (oldName != name)
			broadcastDetailChange(DetailChange::Rename, lock);
	}

	void HostGraph::deserializeTopology(cpl::CSerializer::Builder& ar, cpl::Version version, const GraphLock& lock)
	{
		hadTopologyDeserialized = true;
		decltype(nodeID) copyID;
		
		if (version >= cpl::Version(0, 4, 2))
		{
			ar >> cpl::CSerializer::OptionalWrapper(copyID);
		}
		else
		{
			ar >> cpl::CSerializer::OptionalWrapper(copyID, cpl::CSerializer::DeprecatedBinaryDeserialization {});
		}

		clearTopology(lock);

		uint32_t count;
		ar >> count;

		SerializedEdge copy{ serializeReference(this, lock), DirectedPortPair{} };

		for (uint32_t i = 0; i < count; ++i)
		{
			ar >> copy;
			topology[copy.first].inputs.insert(copy.second);
		}

		expectedNodesToResurrect = topology.size();

		changeIdentity(copyID, lock);

		if (count || expectedNodesToResurrect)
		{
			for (auto h : staticSet)
			{
				if (h->hasSerializedRepresentation())
					tryRebuildTopology(resolve(h), lock, false);
			}
		}
	}


	void HostGraph::changeIdentity(const std::optional<HostGraph::SerializedHandle>& newIdentity, const GraphLock& lock)
	{
		/*
			1. we don't have a name and we're not getting one (do nothing)
			2. we don't have a name but we're getting one (broadcast)
			3. we do have a name but we're erasing it ("destroy self then resurrect nameless")
			4. we do have a name and we're replacing it (broadcast)
			5. we are getting a name that already exists and become an alias (destroy self then orphan into existing, potentially to be resurrected)
			6. we already were an alias, in combination with any other
		*/

		// Ideally.. serialized references live in scene serialization and not in presets.
		// how can we detect that difference? ie. don't adopt names from presets,
		// only from scenes.
		// solution: either have archive::isUserRequest() or - make sure clients of the host graph
		// only ever deserialize this from direct host DAW request (not presets).
		// thus, this is serialized "outside" of the game.

		if (isAlias)
		{
			// 6
			CPL_RUNTIME_ASSERTION(nodeID.has_value());
			CPL_RUNTIME_ASSERTION(aliases.empty());

			if (auto collision = lookupForeign(*nodeID, lock))
			{
				CPL_RUNTIME_ASSERTION(collision != this);

				auto thisShared = shared_from_this();

				auto pos = std::find_if(collision->aliases.begin(), collision->aliases.end(), [&] (auto w) { return w.lock() == thisShared; });
				CPL_RUNTIME_ASSERTION(pos != collision->aliases.end());
				collision->aliases.erase(pos);
			}

			isAlias = false;
		}

		broadcastDestruct(lock);
		if (newIdentity.has_value())
		{
			if (auto collision = lookupForeign(*newIdentity, lock); collision && collision != this)
			{
				// 5
				isAlias = true;
				collision->aliases.push_back(weak_from_this());
			}
		}

		nodeID = newIdentity;
		broadcastCreate(lock);

		if (aliases.size())
		{
			resurrectNextAlias(lock);
		}
	}

	void HostGraph::resurrectNextAlias(const GraphLock& lock)
	{
		if (aliases.empty())
			return;

		auto deadAliases = std::remove_if(aliases.begin(), aliases.end(), [] (auto w) { return w.expired(); });
		aliases.erase(deadAliases, aliases.end());

		if (aliases.size() > 0)
		{
			auto front = aliases.front().lock();
			aliases.erase(aliases.begin());

			front->isAlias = false;
			std::swap(front->aliases, aliases);
			front->broadcastDetailChange(DetailChange::Reidentified, lock);
		}
	}

	void HostGraph::assumeNonAliasedIdentity()
	{
		TriggerModelUpdateOnExit exit{ this };

		const GraphLock lock(staticMutex);
		CPL_RUNTIME_ASSERTION(isAlias);

		changeIdentity(SerializedHandle::generateUnique(), lock);
	}

	bool HostGraph::computeIsDefaultLayout(const TriggerModelUpdateOnExit& lock) const noexcept
	{
		if (!hasSerializedRepresentation())
			return true;

		if (topology.size() != 1)
			return false;

		if (auto it = topology.find(*nodeID); it != topology.end())
		{
			if (it->first != *nodeID)
				return false;

			if (it->second.inputs.size() != 2)
				return false;

			if (it->second.inputs.count({ 0, 0 }) == 0)
				return false;

			if (it->second.inputs.count({ 1, 1 }) == 0)
				return false;
		}
		else
		{
			return false;
		}

		return true;
	}


	HostGraph::Model HostGraph::getModel()
	{
		Model m;
		m.isAlias = isAlias;

		GraphLock lock(staticMutex);

		for (auto n : staticSet)
		{
			if (n->isAlias)
			{
				continue;
			}

			const auto offset = m.connections.size();
			const auto& serialized = serializeReference(n, lock);

			if (auto it = topology.find(serialized); it != topology.end())
			{
				for (const auto& pair : it->second.inputs)
				{
					m.connections.emplace_back(pair);
				}
			}

			m.nodes.emplace_back(
				Model::NodeView
				{
					serialized,
					n->name,
					static_cast<int>(offset),
					static_cast<int>(m.connections.size() - offset),
					static_cast<int>(n->getNumChannels()),
					n->version,
					false
				}
			);

			if (n == this)
			{
				m.hostIndex = 0;
				std::swap(m.nodes[m.nodes.size() - 1], m.nodes[0]);
			}
		}

		// add missing nodes, rare but slow intersection
		if (m.nodes.size() != topology.size())
		{
			for (auto& entry : topology)
			{
				if (auto it = std::find_if(m.nodes.begin(), m.nodes.end(), [&](const auto& n) { return n.node == entry.first; }); it == m.nodes.end())
				{
					const auto offset = m.connections.size();

					for (const auto& pair : entry.second.inputs)
					{
						m.connections.emplace_back(pair);
					}

					auto peer = entry.second.liveReference;

					m.nodes.emplace_back(
						Model::NodeView
						{
							entry.first,
							peer ? peer->name : "<missing>",
							static_cast<int>(offset),
							static_cast<int>(m.connections.size() - offset),
							static_cast<int>(peer ? peer->getNumChannels() : 2),
							peer ? peer->version : 0,
							true
						}
					);
				}
			}
		}

		return m;
	}

	void HostGraph::updateModel(HostGraph::Model& m)
	{
		auto current = globalVersion.load();
		m = getModel();
		m.previousVersion = current; // bad.
	}

	void HostGraph::setName(std::string_view newName)
	{
		name = newName;
		broadcastDetailChange(DetailChange::Rename, GraphLock(staticMutex));
	}

	bool HostGraph::connect(const SerializedHandle& input, DirectedPortPair pair)
	{
		TriggerModelUpdateOnExit exit{ this };

		GraphLock lock(staticMutex);

		return internalConnect(input, pair, lock);
	}

	bool HostGraph::internalConnect(const SerializedHandle& input, DirectedPortPair pair, const GraphLock& lock)
	{
		bool known = topology.count(input);
		auto& relation = topology[input];

		if (relation.inputs.count(pair) != 0)
			return false;

		relation.inputs.insert(pair);

		if (auto h = resolve(input, lock))
		{
			submitConnect(h, pair, lock);
		}
		else if (!known)
		{
			// (only expect this once)
			expectedNodesToResurrect++;
		}

		return true;
	}

	bool HostGraph::disconnect(const SerializedHandle& input, DirectedPortPair pair)
	{		
		TriggerModelUpdateOnExit exit{ this };

		GraphLock lock(staticMutex);

		return internalDisconnect(input, pair, lock);
	}

	bool HostGraph::toggleSet(const std::vector<SerializedHandle>& handles)
	{
		TriggerModelUpdateOnExit exit{ this };

		GraphLock lock(staticMutex);

		if (handles.empty())
			return true;

		bool wasEverythingConnected = true;
		std::size_t disconnections = 0;

		// check whether everything is connected, and disconnect everything on the way
		for (const auto& h : handles)
		{
			bool didDisconnections = resetInstancedTopologyFor(h, lock, true);
			wasEverythingConnected = wasEverythingConnected && didDisconnections;
			
			if (didDisconnections)
				disconnections++;
		}

		// if everything was connected, we've now disconnected everything
		if (wasEverythingConnected || disconnections >= MaxInputChannels / 2)
			return true;

		// otherwise, let's reconstruct a new set of topology 
		const auto maxInputs = MaxInputChannels;
		std::bitset<32> connectedPorts;
		CPL_RUNTIME_ASSERTION(connectedPorts.size() > MaxInputChannels);

		// compile list of free ports
		for (const auto& rel : topology)
		{
			for (const auto& c : rel.second.inputs)
			{
				connectedPorts[c.Destination] = true;
			}
		}
			
		for (const auto& h : handles)
		{
			auto other = lookupForeign(h, lock);

			if (!other)
				return false;

			auto numChannels = other->getNumChannels();

			// search for an offset where we fit
			for(std::size_t offset = 0; offset < connectedPorts.size(); offset += 2) // (stereo offsets)
			{
				bool locallyFound = true;
				for (std::size_t z = 0; locallyFound && z < numChannels; ++z)
				{
					if (z + offset >= maxInputs || connectedPorts[z + offset])
						locallyFound = false;
				}

				if (locallyFound)
				{
					for (PinInt i = 0; i < numChannels; ++i)
					{
						auto hostPin = static_cast<PinInt>(offset) + i;
						connectedPorts[hostPin] = true;
						internalConnect(h, { i, hostPin }, lock);
					}
					break;
				}				
			}
		}

		return true;
	}

	bool HostGraph::internalDisconnect(const SerializedHandle& input, DirectedPortPair pair, const GraphLock& lock)
	{
		bool known = topology.count(input);
		auto& relation = topology[input];

		if (relation.inputs.count(pair) == 0)
			return false;

		relation.inputs.erase(pair);

		if (auto h = resolve(input, lock))
		{
			submitDisconnect(h, pair, lock);
		}
		else if (!known)
		{
			expectedNodesToResurrect++;
		}

		return true;
	}

	void HostGraph::setMixGraph(MixGraphListener::Handle& handle)
	{
		GraphLock lock(staticMutex);

		mix = handle.listener;

		if (!mix)
			return;

		expectedNodesToResurrect = topology.size();

		if (expectedNodesToResurrect)
		{
			for (auto h : staticSet)
			{
				if(h->hasSerializedRepresentation())
					tryRebuildTopology(resolve(h), lock, true);
			}
		}
	}

	void HostGraph::applyDefaultLayoutFromRuntime()
	{
		GraphLock g(staticMutex);

		if (!NONTERMINAL_ASSUMPTION(topology.size() == 0))
			return;

		PinInt numChannels = getNumChannels();

		if (numChannels > 0)
		{
			TriggerModelUpdateOnExit exit { this };

			const auto& ref = serializeReference(this, g);

			for (PinInt i = 0; i < numChannels; ++i)
			{
				internalConnect(ref, { i, i }, g);
			}
		}

		isDefaultRuntimeLayout = true;
	}


	int HostGraph::getNumChannels() const
	{
		return 2; // TODO: Fix
		//return mix.realtime->getInfo().anticipatedChannels;
	}

	bool HostGraph::hasSerializedRepresentation() const
	{
		return !isAlias && nodeID.has_value();
	}

	void HostGraph::submitConnect(HHandle h, DirectedPortPair pair, const GraphLock&)
	{
		if (!mix)
			return;

		auto& other = *resolve(h);

		mix->connect(other.realtime, pair, other.name);
	}

	void HostGraph::submitDisconnect(HHandle h, DirectedPortPair pair, const GraphLock&)
	{
		if (!mix)
			return;

		auto& other = *resolve(h);

		mix->disconnect(
			other.realtime, pair
		);
	}

	HostGraph::HHandle HostGraph::lookupForeign(const SerializedHandle& h, const GraphLock&)
	{
		for (auto g : staticSet)
		{
			if (g->nodeID.has_value() && *g->nodeID == h && !g->isAlias)
			{
				return g;
			}
		}

		return nullptr;
	}

	HostGraph::HHandle HostGraph::resolve(const SerializedHandle& h, const GraphLock& g)
	{
		if (auto it = topology.find(h); it != topology.end())
		{
			if (it->second.liveReference == nullptr)
				it->second.liveReference = lookupForeign(h, g);

			return it->second.liveReference;
		}

		return nullptr;
	}


	HostGraph::HHandle HostGraph::lookupPotentiallyForeign(const SerializedHandle& h, const GraphLock& g)
	{
		auto attempt = resolve(h, g);

		if (attempt)
			return attempt;

		return lookupForeign(h, g);
	}

	void HostGraph::clearTopology(const GraphLock& g)
	{
		for (auto& serialized : topology)
			resetInstancedTopologyFor(serialized.first, g, false);

		topology.clear();
	}

	void HostGraph::tryRebuildTopology(HostGraph* other, const GraphLock& g, bool rebuildPreviouslySeen)
	{
		CPL_RUNTIME_ASSERTION(other);
		CPL_RUNTIME_ASSERTION(other->hasSerializedRepresentation());

		auto serialized = serializeReference(other, g);

		if (auto it = topology.find(serialized); it != topology.end() && (rebuildPreviouslySeen || it->second.liveReference == nullptr))
		{
			CPL_RUNTIME_ASSERTION(expectedNodesToResurrect > 0);

			it->second.liveReference = other;
			expectedNodesToResurrect--;

			for (const auto& pair : it->second.inputs)
			{
				submitConnect(other, pair, g);
			}
		}
	}

	bool HostGraph::resetInstancedTopologyFor(const SerializedHandle& serialized, const GraphLock& g, bool eraseSerializedInfo)
	{
		if (auto it = topology.find(serialized); it != topology.end())
		{
			if (it->second.liveReference != nullptr)
			{
				// expect to see this node again
				// TODO ^: What if we're calling this twice? do we expect this node twice?
				if(!eraseSerializedInfo)
					expectedNodesToResurrect++;

				for (const auto& pair : it->second.inputs)
				{
					// TODO: What if this is a common quick create-destroy-create sequence?
					submitDisconnect(it->second.liveReference, pair, g);
				}

				it->second.liveReference = nullptr;
			}

			if(eraseSerializedInfo)
				topology.erase(serialized);

			return true;
		}

		return false;
	}

	HostGraph::SerializedHandle HostGraph::serializeReference(HHandle h, const GraphLock&)
	{
		// We should only serialize references to aliases in the edge case where it's ourselves.
		CPL_RUNTIME_ASSERTION(!h->isAlias || resolve(h) == this);

		if (h->nodeID.has_value())
			return *h->nodeID;
		else
			return h->nodeID.emplace(SerializedHandle::generateUnique());
	}

	HostGraph* HostGraph::resolve(HHandle h)
	{
		return h;
	}

	void HostGraph::broadcastDetailChange(HostGraph::DetailChange change, const GraphLock& g)
	{
		version = globalVersion.fetch_add(1);

		for (auto n : staticSet)
			n->onDetailChange(this, change, g);
	}

	void HostGraph::broadcastCreate(const GraphLock& g)
	{
		for (auto n : staticSet)
			n->onNodeCreated(this, g);

		staticSet.insert(this);
	}

	void HostGraph::broadcastDestruct(const GraphLock& g)
	{
		for (auto n : staticSet)
			n->onNodeDestroyed(this, g);

		staticSet.erase(this);
	}

	// ----- reactive

	void HostGraph::onNodeCreated(HostGraph::HHandle n, const GraphLock& g)
	{
		TriggerModelUpdateOnExit exit { this };

		if (expectedNodesToResurrect == 0)
			return;

		auto node = resolve(n);

		if (!node->hasSerializedRepresentation())
			return;

		tryRebuildTopology(node, g, false);
	}

	void HostGraph::onDetailChange(HostGraph::HHandle n, HostGraph::DetailChange change, const GraphLock& g)
	{
		TriggerModelUpdateOnExit exit { this };

		auto node = resolve(n);

		// TODO: Handle nodes that change reference...
		if (change != DetailChange::Reidentified || expectedNodesToResurrect == 0 || !node->hasSerializedRepresentation())
			return;

		tryRebuildTopology(node, g, false);
	}

	void HostGraph::onNodeDestroyed(HostGraph::HHandle n, const GraphLock& g)
	{
		TriggerModelUpdateOnExit exit { this };

		auto node = resolve(n);

		if (!node->hasSerializedRepresentation())
			return;

		const auto& serialized = serializeReference(node, g);

		// TODO: when to delete old connections? show them to user or accept transactions with dead plugins be forgotten?
		if (!resetInstancedTopologyFor(serialized, g /*, node->aliases.empty()*/))
			return;
	}
}
