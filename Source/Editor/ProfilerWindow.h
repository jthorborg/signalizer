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

	file:ProfilerWindow.h

		A diagnostics window visualizing the profiling lanes of this Signalizer.

*************************************************************************************/

#ifndef SIGNALIZER_PROFILERWINDOW_H
	#define SIGNALIZER_PROFILERWINDOW_H

	#include <memory>
	#include <vector>

	namespace cpl
	{
		namespace Profiling
		{
			struct Lane;
		}
	}

	namespace Signalizer
	{
		class MainEditor;

		class ProfilerWindow : public juce::DocumentWindow
		{
		public:

			ProfilerWindow(MainEditor* editor, const std::vector<std::shared_ptr<cpl::Profiling::Lane>>& lanes);
			~ProfilerWindow();

			void mainEditorDied();
			void closeButtonPressed() override;

		private:

			MainEditor* editor;
		};
	}

#endif
