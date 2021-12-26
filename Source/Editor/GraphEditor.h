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

	file:GraphEditor.h

		An UI for editing the host graph connections feeding into this Signalizer.

*************************************************************************************/

#ifndef SIGNALIZER_GRAPHEDITOR_H
	#define SIGNALIZER_GRAPHEDITOR_H

	#include "../Common/HostGraph.h"
	#include <cpl/Common.h>
	#include <cpl/gui/GUI.h>
	#include <vector>
	#include <memory>

	namespace Signalizer
	{
		class AudioProcessor;
		class MainEditor;
		class Content;

		class GraphEditor : public juce::DocumentWindow
		{
		public:

			GraphEditor(MainEditor* editor, HostGraph& h);
			~GraphEditor();

			void mainEditorDied();
			void closeButtonPressed() override;

		private:
			std::shared_ptr<Content> content;

			MainEditor* editor;
			HostGraph& host;
		};
	}

#endif
