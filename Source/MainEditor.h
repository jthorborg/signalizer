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
 
	file:MainEditor.h

		Definitions for the main plugin editor, that handles all control logic
		and delegations to other UI views, as well as settings and communication
		between processor and UI.
 
*************************************************************************************/

#ifndef _MAINEDITOR_H
	#define _MAINEDITOR_H

	#include <cpl/Common.h>
	#include "PluginProcessor.h"
	#include <cpl/gui/gui.h>
	#include <map>
	#include <stack>
	#include "SignalizerDesign.h"
	#include <array>

	namespace Signalizer
	{

		class MainEditor
		: 
			public		AudioProcessorEditor, 
			private		juce::Timer, 
			private		juce::HighResolutionTimer,
			protected	cpl::CBaseControl::Listener,
			private		cpl::CBaseControl::ValueFormatter,
			public		cpl::CTopView, 
			private		cpl::COpenGLView::OpenGLEventListener,
			protected	cpl::CTextTabBar<>::CTabBarListener,
			private		juce::ComponentBoundsConstrainer,
			protected	juce::KeyListener,
			public		juce::ComponentListener
		{

		public:

			MainEditor(AudioProcessor * e, ParameterSet * params);
			~MainEditor();

			// CView overrides
			juce::Component * getWindow() override { return this; }
			std::unique_ptr<juce::Component> createEditor() override;
			void panelOpened(cpl::CTextTabBar<> * obj) override;
			void panelClosed(cpl::CTextTabBar<> * obj) override;
			void tabSelected(cpl::CTextTabBar<> * obj, int index) override;
			void activeTabClicked(cpl::CTextTabBar<> * obj, int index) override;
			void suspend() override;
			void resume() override;

			// Components and listeners.
			virtual bool keyPressed(const KeyPress &key, Component *originatingComponent) override;
			void paint(Graphics& g) override;
			void resized() override;
			void resizeEnd() override;
			void resizeStart() override;
			void focusGained(FocusChangeType cause) override;
			void focusLost(FocusChangeType cause) override;
			virtual void mouseDown(const MouseEvent& event) override;
			virtual void mouseUp(const MouseEvent& event) override;

			void componentMovedOrResized(Component& component, bool wasMoved, bool wasResized) override;
			void componentParentHierarchyChanged(Component& component) override;

			// timers
			void timerCallback() override;
			void hiResTimerCallback() override;

			// functionality
			void setRefreshRate(int rateInMs);

			// no parameter = fetch antialiasing from UI combo box
			void setAntialiasing(int multiSamplingLevel = -1);
			// doesn't actually change anything - only updates the selected value in the preset list.
			void setSelectedPreset(juce::File newPreset);

			void showAboutBox();
		
		protected:

			static const int elementSize = 25;
			// border around all elements, from which the background shines through'
			static const int elementBorder = 1;

			int getViewTopCoordinate() const noexcept;
			void onOGLRendering(cpl::COpenGLView * view) noexcept override;
			void onOGLContextCreation(cpl::COpenGLView * view) noexcept override;
			void onOGLContextDestruction(cpl::COpenGLView * view) noexcept override;

			// cpl::CBaseControl interface
			void valueChanged(const cpl::CBaseControl * cbc) override;
			bool stringToValue(const cpl::CBaseControl * ctrl, const std::string & valString, cpl::iCtrlPrec_t & val) override;
			bool valueToString(const cpl::CBaseControl * ctrl, std::string & valString, cpl::iCtrlPrec_t val) override;
			void onObjectDestruction(const cpl::Utility::DestructionServer<cpl::CBaseControl>::ObjectProxy & destroyedObject) override;
		private:

			typedef std::vector<std::unique_ptr<juce::Component>>::iterator EditorIterator;
		
			void deserialize(cpl::CSerializer & se, cpl::Version version) override;
			void serialize(cpl::CSerializer & se, cpl::Version version) override;

			struct Flags
			{
				cpl::ABoolFlag
					/// <summary>
					/// Set this to alter the swap interval
					/// </summary>
					swapIntervalChanged,
					/// <summary>
					/// Set this to alter whether the context is repainting continuously
					/// </summary>
					continuousRepaint;
			} mtFlags;

			struct NewChanges
			{
				std::atomic<int> swapInterval;
				std::atomic<bool> repaintContinuously;
			} newc;


			// the z-ordering system ensures this is basically a FIFO system
			void pushEditor(juce::Component * editor);
			void pushEditor(std::unique_ptr<juce::Component> editor);
			template<typename Pred>
				bool removeAnyEditor(Pred pred);
			juce::Component * getTopEditor() const;
			void popEditor();
			void deleteEditor(EditorIterator i);
			void clearEditors();
			cpl::CView * viewFromIndex(std::size_t index);
			void addTab(const std::string & name);
			void restoreTab();
			int getRenderEngine();
			void initUI();
			void suspendView(cpl::CView * view);
			void initiateView(cpl::CView * view, bool spawnNewEditor = false);
			void enterFullscreenIfNeeded(juce::Point<int> where);
			void enterFullscreenIfNeeded();
			void exitFullscreen();
			void setPreferredKioskCoords(juce::Point<int>) noexcept;

			// Relations
			AudioProcessor * engine;

			// Constant UI
			cpl::CTextTabBar<> tabs;
			cpl::CSVGButton ksettings, kfreeze, khelp, kkiosk;

			// Editor controls
			cpl::CButton kstableFps, kvsync, krefreshState, kidle;
			cpl::CInputControl kmaxHistorySize;
			cpl::CKnobSlider krefreshRate, kswapInterval;
			cpl::CComboBox krenderEngine, kantialias;
			cpl::CPresetWidget kpresets;
			std::array<cpl::CColourControl, cpl::CLookAndFeel_CPL::numColours>  colourControls;

			// state variables.
			int refreshRate;
			int viewTopCoord;
			std::int32_t selTab;
			cpl::iCtrlPrec_t oldRefreshRate;
			// TODO: refactor to not use. Do not trust.
			bool isEditorVisible;
			bool unFocused, idleInBack, firstKioskMode, hasAnyTabBeenSelected;
			juce::Point<int> kioskCoords;
			juce::Rectangle<int> preFullScreenSize;
			// View related data
			Signalizer::CDefaultView defaultView;
			juce::OpenGLContext oglc;

			struct ViewWithSerializedFlag
			{
				std::unique_ptr<cpl::CSubView> view;
				bool hasBeenRestored;
			};

			// .first = whether the view has been serialized.
			std::map<std::string, ViewWithSerializedFlag> views;
			std::vector<std::unique_ptr<juce::Component>> editorStack;
			cpl::CView * currentView;
			ResizableCornerComponent rcc;
			cpl::CSerializer viewSettings;
			ParameterSet * params;
			//cpl::CMessageSystem messageSystem;
		};
	};

#endif  // MainEditor
