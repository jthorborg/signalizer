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

namespace Signalizer
{
	std::mutex HostGraph::staticMutex;
	std::set<HostGraph::HHandle> HostGraph::staticSet;

	static int isTheFirst = 0;

	HostGraph::HostGraph(AudioStream& realtime, AudioStream& presentation)
		: mix(realtime, presentation, isTheFirst++ < 2)
	{
		isFirst = isTheFirst < 3;
		//isTheFirst = true;


		broadcastCreate();
	}

	HostGraph::~HostGraph()
	{
		broadcastDestruct();




	}

	std::shared_ptr<juce::Component> HostGraph::createEditor()
	{
		auto shared = std::make_shared<Editor>();
		editor = shared;

		return std::move(shared);
	}

	void HostGraph::serialize(cpl::CSerializer::Archiver& ar, cpl::Version version)
	{
		// TODO: No lock?
		std::lock_guard<std::mutex> lock(staticMutex);

		ar << name;
		ar << nodeID;

		PinInt counter = std::count_if(pinInputs.begin(), pinInputs.end(), [&](auto el) { return el != nullptr; });

		ar << counter;

		for(PinInt i = 0; i < pinInputs.size(); ++i)
		{
			if (pinInputs[i] != nullptr)
				ar << std::make_pair(i, serializeReference(pinInputs[i]));
		}
	}

	void HostGraph::deserialize(cpl::CSerializer::Builder& ar, cpl::Version version)
	{		
		// TODO: No lock?
		std::lock_guard<std::mutex> lock(staticMutex);

		ar >> name;
		ar >> nodeID;

		std::set<SerializedEdge> connections;

		PinInt count;
		ar >> count;

		for (PinInt i = 0; i < count; ++i)
		{
			SerializedEdge edge;
			ar >> edge;
			connections.insert(edge);
		}
	}

	HostGraph::SerializedHandle HostGraph::serializeReference(HHandle h)
	{
		return h->nodeID;
	}

	HostGraph* HostGraph::resolve(HHandle h)
	{
		return h;
	}

	void HostGraph::broadcastCreate()
	{
		std::lock_guard<std::mutex> lock(staticMutex);

		for (auto n : staticSet)
			n->onNodeCreated(this);

		staticSet.insert(this);
	}

	void HostGraph::broadcastDestruct()
	{
		std::lock_guard<std::mutex> lock(staticMutex);

		for (auto n : staticSet)
			n->onNodeDestroyed(this);

		staticSet.erase(this);
	}

	void HostGraph::onNodeCreated(HostGraph::HHandle n)
	{
		auto node = resolve(n);
		auto serialized = serializeReference(n);
		if (auto it = expectedEdges.find(serialized); it != expectedEdges.end())
		{
			pinInputs[it->second] = node;
			expectedEdges.erase(serialized);
		}

		if(auto ed = editor.lock())
		{
			ed->triggerAsyncUpdate();
		}

		if (isFirst)
			mix.connect(node->mix.realtimeStream());
		
	}

	void HostGraph::onNodeDestroyed(HostGraph::HHandle n)
	{
		auto node = resolve(n);

		// TODO: Scan and remove.
		if (auto ed = editor.lock())
		{
			ed->triggerAsyncUpdate();
		}

		if (isFirst)
			mix.disconnect(node->mix.realtimeStream());
	}

	void HostGraph::Editor::handleAsyncUpdate()
	{
	}
}
