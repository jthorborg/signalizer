/*************************************************************************************

	Signalizer - cross-platform audio visualization plugin - v. 0.x.y

	Copyright (C) 2026 Janus Lynggaard Thorborg (www.jthorborg.com)

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

	file:ProfilerWindow.cpp

		Implementation of ProfilerWindow.cpp

*************************************************************************************/

#include <memory>

#if CPL_PROFILING
#include <cpl/profiling/JuceModelDisplay.h>
#endif

#include "ProfilerWindow.h"
#include "MainEditor.h"
#include "../Common/SharedBehaviour.h"

namespace Signalizer
{
	ProfilerWindow::ProfilerWindow(MainEditor* editor, const std::vector<std::shared_ptr<cpl::Profiling::Lane>>& lanes)
		: juce::DocumentWindow("Profiler", juce::Colours::black, juce::DocumentWindow::TitleBarButtons::allButtons)
		, editor(editor)
	{
		setUsingNativeTitleBar(true);
		setResizable(true, false);
#if CPL_PROFILING
		setContentOwned(new cpl::Profiling::EWMAProfilerComponent(lanes), false);
#endif
		centreWithSize(800, 500);
		setVisible(true);
	}

	ProfilerWindow::~ProfilerWindow()
	{
		if (editor)
			editor->profilerWindowDied();
	}

	void ProfilerWindow::mainEditorDied()
	{
		editor = nullptr;
		delete this;
	}

	void ProfilerWindow::closeButtonPressed()
	{
		delete this;
	}
}

