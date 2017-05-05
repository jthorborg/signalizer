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
 
	file:SentientViewState.h

		A state object that automatically controls creation, serialization, caching
		and deserialization of views and their editors.
 
*************************************************************************************/

#ifndef SIGNALIZER_SENTIENTVIEWSTATE_H
#define SIGNALIZER_SENTIENTVIEWSTATE_H

#include "CommonSignalizer.h"
#include <memory>
#include <cpl/Utility.h>	

namespace Signalizer
{
	class SentientViewState
	{
	public:

		SentientViewState(SentientViewState && other)
			: state(other.state)
			, name(std::move(other.name))
			, editor(std::move(other.editor))
			, view(std::move(other.view))
		{
			editor.replaceGenerator([this] { return state.createEditor(); });
		}

		SentientViewState & operator = (SentientViewState && other) = delete;
		SentientViewState(const SentientViewState & other) = delete;

		SentientViewState(const std::string & name, ProcessorState & viewProcessorState, DecoupledStateObject<cpl::CSubView>::FGenerator generator)
			: state(viewProcessorState)
			, name(name)
			, editor(
				[this] { return state.createEditor(); }, 
				[](StateEditor & editor, cpl::CSerializer & sz, cpl::Version v) { editor.getEditorSO().serializeObject(sz, v); },
				[](StateEditor & editor, cpl::CSerializer & sz, cpl::Version v) { editor.getEditorSO().deserializeObject(sz, v); }
			)
			, view(
				generator,
				[](cpl::CSubView & view, cpl::CSerializer & sz, cpl::Version v) { view.serializeObject(sz, v); },
				[](cpl::CSubView & view, cpl::CSerializer & sz, cpl::Version v) { view.deserializeObject(sz, v); }
			)
		{

		}

		DecoupledStateObject<StateEditor> & getEditorDSO() noexcept { return editor; }
		DecoupledStateObject<cpl::CSubView> & getViewDSO() noexcept { return view; }
		const std::string & getName() const noexcept { return name; }
		ProcessorState & getProcessorState() noexcept { return state; }

	private:

		ProcessorState & state;
		std::string name;
		DecoupledStateObject<StateEditor> editor;
		DecoupledStateObject<cpl::CSubView> view;
	};

};


#endif