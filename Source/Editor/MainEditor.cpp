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

	file:MainEditor.cpp

		Implementation of MainEditor.h

*************************************************************************************/

#include "MainEditor.h"
#include "../Processor/PluginProcessor.h"
#include "../Vectorscope/Vectorscope.h"
#include "../Oscilloscope/Oscilloscope.h"
#include "../Spectrum/Spectrum.h"
#include <cpl/CPresetManager.h>
#include <cpl/LexicalConversion.h>
#include "version.h"
#include <cpl/Mathext.h>

namespace cpl
{
	const ProgramInfo programInfo
	{
		"Signalizer",
		cpl::Version::fromParts(SIGNALIZER_MAJOR, SIGNALIZER_MINOR, SIGNALIZER_BUILD),
		"Janus Thorborg",
		"sgn",
		false,
		nullptr,
		SIGNALIZER_BUILD_INFO
	};

};

namespace Signalizer
{
	const static int kdefaultMaxSkippedFrames = 10;

	const static int kdefaultLength = 700, kdefaultHeight = 480;
	const static std::vector<std::string> RenderingEnginesList = { "Software", "OpenGL" };

	const static juce::String MainEditorName = "Main Editor Settings";

	std::string MainPresetName = "main";
	std::string DefaultPresetName = "default";

	enum class ViewTypes
	{
		Vectorscope,
		Oscilloscope,
		Spectrum,
		end
	};

	std::array<const char *, static_cast<std::size_t>(ViewTypes::end)> ViewIndexToMap =
	{
		"Vectorscope",
		"Oscilloscope",
		"Spectrum"
		//"Statistics"
	};

	template<typename T>
	inline std::unique_ptr<ProcessorState> CreateState(std::size_t offset, bool shouldCreateShortNames, SystemView system)
	{
		return std::unique_ptr<ProcessorState>(new T(offset, shouldCreateShortNames, system));
	}

	std::vector<std::pair<std::string, ParameterCreater>> ParameterCreationList =
	{
		{ ViewIndexToMap[(int)ViewTypes::Vectorscope], &CreateState<VectorScopeContent> },
		{ ViewIndexToMap[(int)ViewTypes::Oscilloscope], &CreateState<OscilloscopeContent> },
		{ ViewIndexToMap[(int)ViewTypes::Spectrum], &CreateState<SpectrumContent> }
	};

	template<typename... Args>
	std::unique_ptr<cpl::CSubView> GenerateView(ViewTypes type, Args&&... args)
	{
		switch (type)
		{
		case ViewTypes::Vectorscope: return std::make_unique<VectorScope>(args...);
		case ViewTypes::Oscilloscope: return std::make_unique<Oscilloscope>(args...);
		case ViewTypes::Spectrum: return std::make_unique<Spectrum>(args...);
		}
		CPL_RUNTIME_EXCEPTION("Unknown view generation index");
	}

	enum class RenderTypes
	{
		Software,
		openGL,
		end
	};

	enum class Utility
	{
		Freeze,
		Sync,
		IdleInBack,
		end

	};

	/// <summary>
	/// TODO: Query this at runtime
	/// </summary>
	const static std::array<int, 5> AntialisingLevels =
	{
		1,
		2,
		4,
		8,
		16
	};

	/// <summary>
	/// TODO: Query this at runtime
	/// </summary>
	const static std::vector<std::string> AntialisingStringLevels =
	{
		"1",
		"2",
		"4",
		"8",
		"16"
	};

	MainEditor::MainEditor(AudioProcessor * e, ParameterMap * parameterMap)
		: engine(e)
		, params(parameterMap)
		, AudioProcessorEditor(e)
		, CTopView(this, "Signalizer main editor")
		, rcc(this, this)
		, krenderEngine("Rendering Engine", RenderingEnginesList)
		, refreshRate(0)
		, oldRefreshRate(0)
		, unFocused(true)
		, idleInBack(false)
		, isEditorVisible(false)
		, selTab(0)
		, currentView(nullptr)
		, kioskCoords(-1, -1)
		, firstKioskMode(false)
		, hasAnyTabBeenSelected(false)
		, viewTopCoord(0)
		, kpresets(e, MainPresetName, kpresets.WithDefault)
		, kmaxHistorySize("History size")
		, tabBarTimer()
		, mouseHoversTabArea(false)
		, tabBarIsVisible(true)
	{
		// TODO: figure out why moving a viewstate causes corruption (or early deletion of moved object)
		views.reserve((std::size_t)ViewTypes::end);

		cpl::foreach_enum<ViewTypes>(
			[&](ViewTypes index)
			{
				int i = cpl::enum_cast<int>(index);
				auto & name = ViewIndexToMap[i];
				auto & localState = *params->getState(i);
				views.emplace_back(
					name,
					localState,
					[=, &localState]
					{
						return GenerateView(index, globalState, name, engine->stream, &localState);
					}
				);
			}
		);

		setOpaque(true);
		setMinimumSize(50, 50);
		setBounds(0, 0, kdefaultLength, kdefaultHeight);
		initUI();
		oglc.setComponentPaintingEnabled(false);

		nestedMouseHook.hook(this, this, true);
	}


	std::unique_ptr<StateEditor> MainEditor::createEditor()
	{
		class CContentPageDummySerialization
			: public CContentPage
			, private cpl::SafeSerializableObject
		{
			virtual cpl::SafeSerializableObject & getEditorSO() override { return *this; }
		};

		auto content = new CContentPageDummySerialization();
		content->setName(MainEditorName);
		if (auto page = content->addPage("Settings", "icons/svg/wrench.svg"))
		{
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&krefreshRate, 0);
				section->addControl(&kstableFps, 1);
				section->addControl(&kswapInterval, 0);
				section->addControl(&kvsync, 1);
				page->addSection(section, "Update");
			}
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&krenderEngine, 0);

				section->addControl(&kantialias, 1);

				page->addSection(section, "Quality");
			}
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kmaxHistorySize, 0);
				page->addSection(section, "Globals");
			}
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&krefreshState, 0);
				section->addControl(&kidle, 1);
				section->addControl(&khideTabs, 2);
				section->addControl(&kstopProcessingOnSuspend, 0);
				section->addControl(&khideWidgets, 1);
				page->addSection(section, "Globals");
			}
		}
		if (auto page = content->addPage("Colours", "icons/svg/brush.svg"))
		{
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				int numKnobsPerLine = static_cast<int>(colourControls.size() / 2);
				int rem = colourControls.size() % 2;

				for (int y = 0; y < 2; ++y)
				{
					for (int x = 0; x < numKnobsPerLine; ++x)
					{
						section->addControl(&colourControls[y * numKnobsPerLine + x], y);
					}
				}
				if (rem)
					section->addControl(&colourControls[colourControls.size() - 1], 0);
				page->addSection(section, "Colour Scheme");
			}
		}
		if (auto page = content->addPage("Presets", "icons/svg/save.svg"))
		{

			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kpresets, 0);
				page->addSection(section, "Presets");
			}
		}

		return std::unique_ptr<StateEditor>(content);
	}

	void MainEditor::focusGained(FocusChangeType cause)
	{
		if (idleInBack)
		{
			if (unFocused)
			{
				krefreshRate.bSetValue(oldRefreshRate);
			}
		}

		unFocused = false;
	}
	void MainEditor::focusLost(FocusChangeType cause)
	{
		oldRefreshRate = krefreshRate.bGetValue();
		if (idleInBack)
		{
			if (!unFocused)
			{
				krefreshRate.bSetValue(0.5f);
			}
			unFocused = true;
		}
	}

	void MainEditor::onObjectDestruction(const cpl::Utility::DestructionServer<cpl::CBaseControl>::ObjectProxy & destroyedObject)
	{
		// no-op.
	}

	void MainEditor::setTabBarVisibility(bool toggle)
	{
		if (khideTabs.bGetBoolState())
		{
			if (tabBarIsVisible == !!toggle)
				return;

			tabBarIsVisible = !!toggle;
			resized();
		}
		else
		{
			if (!tabBarIsVisible)
			{
				tabBarIsVisible = true;
				resized();
			}
		}

	}

	void MainEditor::pushEditor(StateEditor * editor)
	{
		pushEditor(std::unique_ptr<StateEditor>(editor));

	}
	void MainEditor::pushEditor(UniqueHandle<StateEditor> newEditor)
	{
		if (!newEditor.get())
			return;

		if (auto editor = getTopEditor())
		{
			editor->setVisible(false);
		}

		editorStack.push_back(std::move(newEditor));

		// beware of move construction! newEditor is now invalid!
		if (auto editor = getTopEditor())
		{
			addAndMakeVisible(editor);
			resized();
			repaint();
			if (!tabs.isOpen())
				tabs.openPanel();
		}
	}
	StateEditor * MainEditor::getTopEditor() const
	{
		return editorStack.empty() ? nullptr : editorStack.back().get();
	}

	void MainEditor::deleteEditor(MainEditor::EditorIterator i)
	{
		removeChildComponent(i->get());
		editorStack.erase(i);

		if(editorStack.empty() && tabs.isOpen())
			tabs.closePanel();
		else
		{
			if (auto editor = getTopEditor())
			{
				editor->setVisible(true);
			}
			resized();
			repaint();
		}
	}

	void MainEditor::popEditor()
	{
		if (!editorStack.empty())
		{
			return deleteEditor(editorStack.begin() + editorStack.size() - 1);
		}
	}
	void MainEditor::clearEditors()
	{
		while (!editorStack.empty())
			popEditor();
	}

	void MainEditor::setRefreshRate(int rate)
	{
		refreshRate = cpl::Math::confineTo(rate, 10, 1000);
		if (kstableFps.getValueReference().getNormalizedValue() > 0.5)
		{
			juce::HighResolutionTimer::startTimer(refreshRate);
		}
		else
		{
			juce::Timer::startTimer(refreshRate);
		}
		if (hasCurrentView())
			activeView().setApproximateRefreshRate(refreshRate);

	}
	void MainEditor::resume()
	{
		setRefreshRate(refreshRate);
	}

	void MainEditor::suspend()
	{
		juce::HighResolutionTimer::stopTimer();
		juce::Timer::stopTimer();
	}

	// these handle cases where our component is being thrown off the kiosk mode by (possibly)
	// another JUCE plugin.
	void MainEditor::componentMovedOrResized(Component& component, bool wasMoved, bool wasResized)
	{
		if (&component == activeView().getWindow())
		{
			if (activeView().getIsFullScreen() && component.isOnDesktop())
			{
				if (juce::ComponentPeer * peer = component.getPeer())
				{
					if (!peer->isKioskMode())
					{
						exitFullscreen();
						if (kkiosk.bGetValue() > 0.5)
						{
							kkiosk.bSetInternal(0.0);
						}
					}
				}
			}
		}
	}
	void MainEditor::componentParentHierarchyChanged(Component& component)
	{
		if (&component == activeView().getWindow())
		{
			if (activeView().getIsFullScreen() && component.isOnDesktop())
			{
				if (juce::ComponentPeer * peer = component.getPeer())
				{
					if (!peer->isKioskMode())
					{
						exitFullscreen();
						if (kkiosk.bGetValue() > 0.5)
						{
							kkiosk.bSetInternal(0.0);
						}
					}
				}
			}
		}
	}

	void MainEditor::valueChanged(const cpl::CBaseControl * c)
	{
		// TODO: refactor all behaviour here out to semantic functions
		// bail out early if we aren't showing anything.
		//if (!hasCurrentView())
		//	return;

		auto value = c->bGetValue();

		// freezing of displays
		if (c == &kfreeze)
		{
			engine->stream.setSuspendedState(value > 0.5);
		}
		// lower display rate if we are unfocused
		else if (c == &kidle)
		{
			idleInBack = value > 0.5 ? true : false;
		}

		else if (c == &ksettings)
		{
			if (value > 0.5)
			{
				// spawn the global setting editor
				pushEditor(this->createEditor());
				tabs.openPanel();
			}
			else
			{
				// -- remove it
				if (auto editor = getTopEditor())
				{
					if (editor->getName() == MainEditorName)
					{
						popEditor();
					}
				}
			}
		}
		else if (c == &kkiosk)
		{
			if (hasCurrentView())
			{
				if (value > 0.5)
				{
					// set a window to fullscreen.
					if (firstKioskMode)
					{
						// dont fetch the current coords if we set the view the first time..
						firstKioskMode = false;
					}
					else
					{
						if (juce::Desktop::getInstance().getKioskModeComponent() == activeView().getWindow())
							return;
						kioskCoords = activeView().getWindow()->getScreenPosition();
					}

					preFullScreenSize = getBounds().withZeroOrigin();

					removeChildComponent(activeView().getWindow());
					activeView().getWindow()->addToDesktop(juce::ComponentPeer::StyleFlags::windowAppearsOnTaskbar);

					activeView().getWindow()->setTopLeftPosition(kioskCoords.x, kioskCoords.y);
					bool useMenusAndBars = false;
					#ifdef CPL_MAC
						useMenusAndBars = true;
					#endif
					juce::Desktop::getInstance().setKioskModeComponent(activeView().getWindow(), useMenusAndBars);
					activeView().setFullScreenMode(true);
					activeView().getWindow()->setWantsKeyboardFocus(true);
					activeView().getWindow()->grabKeyboardFocus();
					// add listeners.

					activeView().getWindow()->addKeyListener(this);
					activeView().getWindow()->addComponentListener(this);

					// sets a minimal view when entering full screen
					setBounds(getBounds().withBottom(getViewTopCoordinate()));
				}
				else
				{
					exitFullscreen();
				}
			}
		}
		else if (c == &kstableFps)
		{
			if (kstableFps.getValueReference().getNormalizedValue() > 0.5)
			{
				juce::Timer::stopTimer();
				juce::HighResolutionTimer::startTimer(refreshRate);
			}
			else
			{
				juce::HighResolutionTimer::stopTimer();
				juce::Timer::startTimer(refreshRate);
			}
		}
		else if (c == &kswapInterval)
		{
			newc.swapInterval.store(
				cpl::Math::round<int>(kswapInterval.bGetValue() * kdefaultMaxSkippedFrames),
				std::memory_order_release
			);

			mtFlags.swapIntervalChanged = true;
		}
		else if (c == &kvsync)
		{
			if (kvsync.getValueReference().getNormalizedValue() > 0.5)
			{
				if (hasCurrentView())
				{
					// this is kind of stupid; the sync setting must be set after the context is created..
					struct RetrySync
					{
						RetrySync(MainEditor * h) : handle(h) {};
						MainEditor * handle;

						void operator()()
						{
							if (handle->oglc.isAttached())
								handle->oglc.setContinuousRepainting(true);
							else
								cpl::GUIUtils::FutureMainEvent(200, RetrySync(handle), handle);
						}

					};

					RetrySync(this)();
				}
			}
			else
			{
				if (hasCurrentView())
				{
					oglc.setContinuousRepainting(false);
				}
			}
		}
		// change of refresh rate
		else if (c == &krefreshRate)
		{
			refreshRate = cpl::Math::round<int>(cpl::Math::UnityScale::exp(value, 10.0, 1000.0));
			setRefreshRate(refreshRate);
		}
		// change of the rendering engine
		else if (c == &krenderEngine)
		{
			auto index = cpl::Math::distribute<RenderTypes>(value);

			switch (index)
			{
			case RenderTypes::Software:

				if (hasCurrentView() && activeView().isOpenGL())
				{
					activeView().detachFromOpenGL(oglc);
					if (oglc.isAttached())
						oglc.detach();
				}


				break;
			case RenderTypes::openGL:
				if (oglc.isAttached())
				{
					if (hasCurrentView() && !activeView().isOpenGL())
					{
						// ?? freaky
						if (cpl::CView * unknownView = dynamic_cast<cpl::CView *>(oglc.getTargetComponent()))
							unknownView->detachFromOpenGL(oglc);
						oglc.detach();
					}
				}
				if (hasCurrentView())
					activeView().attachToOpenGL(oglc);
				break;
			}
		}
		else if (c == &kantialias)
		{
			setAntialiasing();
		}
		else if (c == &krefreshState)
		{
			if (hasCurrentView())
				activeView().resetState();
		}
		else if(c == &khelp)
		{
			showAboutBox();
		}
		else if (c == &kmaxHistorySize)
		{

			struct RetryResizeCapacity
			{
				RetryResizeCapacity(MainEditor * h) : handle(h) {};
				MainEditor * handle;

				void operator()()
				{
					auto currentSampleRate = handle->engine->stream.getInfo().sampleRate.load(std::memory_order_acquire);
					if (currentSampleRate > 0)
					{
						std::int64_t value;
						std::string contents = handle->kmaxHistorySize.getInputValue();
						if (cpl::lexicalConversion(contents, value) && value >= 0)
						{
							auto msCapacity = cpl::Math::round<std::size_t>(currentSampleRate * 0.001 * value);

							handle->engine->stream.setAudioHistoryCapacity(msCapacity);

							if (contents.find_first_of("ms") == std::string::npos)
							{
								contents.append(" ms");
								handle->kmaxHistorySize.setInputValueInternal(contents);
							}

							handle->kmaxHistorySize.indicateSuccess();
						}
						else
						{
							std::string result;
							auto msCapacity = cpl::Math::round<std::size_t>(1000 * handle->engine->stream.getAudioHistoryCapacity() / handle->engine->stream.getInfo().sampleRate);
							if (cpl::lexicalConversion(msCapacity, result))
								handle->kmaxHistorySize.setInputValueInternal(result + " ms");
							else
								handle->kmaxHistorySize.setInputValueInternal("error");
							handle->kmaxHistorySize.indicateError();
						}
					}
					else
					{
						cpl::GUIUtils::FutureMainEvent(200, RetryResizeCapacity(handle), handle);
					}

				}
			};

			RetryResizeCapacity(this)();

		}
		else if (c == &kstopProcessingOnSuspend)
		{
			globalState.stopProcessingOnSuspend.store(kstopProcessingOnSuspend.bGetBoolState(), std::memory_order_release);
		}
		else if (c == &khideWidgets)
		{
			globalState.hideWidgetsOnMouseExit.store(khideWidgets.bGetBoolState(), std::memory_order_release);
		}
		else
		{
			// check if it was one of the colours
			for (unsigned i = 0; i < colourControls.size(); ++i)
			{
				if (c == &colourControls[i])
				{
					// change colour and broadcast event.
					cpl::CLookAndFeel_CPL::defaultLook().getSchemeColour(i).colour = colourControls[i].getControlColourAsColour();
					cpl::CLookAndFeel_CPL::defaultLook().updateColours();
					repaint();
				}

			}
		}
	}

	void MainEditor::setPreferredKioskCoords(juce::Point<int> preferredCoords) noexcept
	{
		firstKioskMode = true;
		kioskCoords = preferredCoords;
	}

	void MainEditor::panelOpened(cpl::CTextTabBar<> * obj)
	{
		isEditorVisible = true;
		// -- settings editor spawned the panel view
		if(!ksettings.bGetBoolState() && !getTopEditor())
		{
			pushEditor(views[selTab].getEditorDSO().getCached());
		}


		auto editorBottom = getViewTopCoordinate();

		if (getBottom() < editorBottom)
			setBounds(getBounds().withBottom(editorBottom));
		else
			resized();
		repaint();
	}

	void MainEditor::panelClosed(cpl::CTextTabBar<> * obj)
	{
		clearEditors();
		ksettings.setToggleState(false, juce::NotificationType::dontSendNotification);
		resized();
		isEditorVisible = false;
	}

	void MainEditor::setAntialiasing(int multisamplingLevel)
	{

		int sanitizedLevel = 1;

		if(multisamplingLevel == -1)
		{
			auto val = cpl::Math::confineTo(kantialias.getZeroBasedSelIndex(), 0, AntialisingLevels.size() - 1);
			sanitizedLevel = AntialisingLevels[val];
		}
		else
		{
			for(unsigned i = 0; i < AntialisingLevels.size(); ++i)
			{
				if(AntialisingLevels[i] == multisamplingLevel)
				{
					sanitizedLevel = AntialisingLevels[i];
					break;
				}
			}
		}
		// TODO: Query OpenGL for glGet GL_MAX_INTEGER_SAMPLES, to see what the maximum supported
		// multisampling level is

		if(sanitizedLevel > 0)
		{
			juce::OpenGLPixelFormat fmt;
			fmt.multisamplingLevel = sanitizedLevel;
			// true if a view exists and it is attached
			bool reattach = false;
			if(hasCurrentView())
			{
				if(activeView().isOpenGL())
				{
					activeView().detachFromOpenGL(oglc);
					reattach = true;
				}
			}
			oglc.setMultisamplingEnabled(true);
			oglc.setPixelFormat(fmt);

			if(reattach)
			{
				activeView().attachToOpenGL(oglc);
			}

		}
		else
		{
			oglc.setMultisamplingEnabled(false);
		}

	}

	int MainEditor::getRenderEngine()
	{
		return (int)cpl::Math::distribute<RenderTypes>(krenderEngine.bGetValue());
	}

	void MainEditor::suspendView(SentientViewState & svs)
	{
		if (svs.getViewDSO().hasCached())
		{
			auto view = svs.getViewDSO().getCached().get();
			view->suspend();
			if (oglc.isAttached())
			{
				view->detachFromOpenGL(oglc);
				oglc.detach();
			}
			else
			{
				// check if view is still attached to a dead context:
				if(view->isOpenGL())
				{
					// something else must have killed the context, check if it's the same

					if(view->getAttachedContext() == &oglc)
					{
						// okay, so we detach it anyway and eat the exceptions:
						view->detachFromOpenGL(oglc);
					}
				}

			}

			// serializing the main editor to a preset where the kkiosk is false, while the main editor has
			// a full screen window creates a.. interesting situation (replaced && with || to be sure).
			if (juce::Desktop::getInstance().getKioskModeComponent() == view->getWindow() || kkiosk.bGetValue() > 0.5)
			{
				exitFullscreen();
			}
			removeChildComponent(activeView().getWindow());
			view->getWindow()->removeMouseListener(this);
		}

	}

	void MainEditor::initiateView(SentientViewState & view, bool spawnNewEditor)
	{
		currentView = &view;
		addAndMakeVisible(activeView().getWindow());

		if ((RenderTypes)getRenderEngine() == RenderTypes::openGL)
		{
			// init all openGL stuff.
			if (auto oglView = dynamic_cast<cpl::COpenGLView*>(view.getViewDSO().getCached().get()))
			{
				oglView->addOpenGLEventListener(this);

				setAntialiasing();
				oglView->attachToOpenGL(oglc);
			}

		}

		if (spawnNewEditor)
		{
			pushEditor(view.getEditorDSO().getCached());
		}

		if (kkiosk.bGetBoolState())
		{
			if (firstKioskMode)
				enterFullscreenIfNeeded(kioskCoords);
			else
				enterFullscreenIfNeeded();
		}
		resized();
		activeView().getWindow()->addMouseListener(this, true);
		activeView().resume();
	}

	void MainEditor::mouseUp(const juce::MouseEvent& event)
	{
		if (hasCurrentView() && event.eventComponent == activeView().getWindow())
		{
			if (event.mods.testFlags(juce::ModifierKeys::rightButtonModifier))
			{
				kfreeze.setToggleState(!kfreeze.getToggleState(), juce::NotificationType::sendNotification);
			}
		}
		else
		{
			AudioProcessorEditor::mouseUp(event);
		}

		
	}

	void MainEditor::mouseDown(const juce::MouseEvent& event)
	{
		if (hasCurrentView() && event.eventComponent == activeView().getWindow())
		{
			if (event.mods.testFlags(juce::ModifierKeys::rightButtonModifier))
			{
				kfreeze.setToggleState(!kfreeze.getToggleState(), juce::NotificationType::sendNotification);
			}
		}
		else
		{
			AudioProcessorEditor::mouseUp(event);
		}
	}

	void MainEditor::mouseDoubleClick(const juce::MouseEvent& event)
	{
		if (hasCurrentView() && event.eventComponent == activeView().getWindow())
		{
			if (event.mods.testFlags(juce::ModifierKeys::rightButtonModifier))
			{
				kfreeze.setToggleState(!kfreeze.getToggleState(), juce::NotificationType::sendNotification);
			}
		}
		else
		{
			AudioProcessorEditor::mouseDoubleClick(event);
		}
	}


	void MainEditor::tabSelected(cpl::CTextTabBar<> * obj, int index)
	{
		index = cpl::Math::confineTo(index, 0, (int)ViewTypes::end - 1);
		hasAnyTabBeenSelected = true;
		// these lines disable the global editor if you switch view.
		//if (ksettings.getToggleState())
		//	ksettings.setToggleState(false, NotificationType::sendNotification);

		// if any editor is open currently, we have to close it and open the new.
		bool openNewEditor = false;
		// deattach old view
		if (hasCurrentView())
		{
			if (activeView().getIsFullScreen())
				setPreferredKioskCoords(activeView().getWindow()->getPosition());


			if (getTopEditor())
				openNewEditor = true;

			clearEditors();

			suspendView(*currentView);
			currentView = nullptr;
		}

		initiateView(views[index], openNewEditor);

		if (openNewEditor && ksettings.bGetBoolState())
			ksettings.bSetInternal(0.0);

		selTab = index;
	}

	template<typename Pred>
		bool MainEditor::removeAnyEditor(Pred pred)
	{
		bool alteredStack = false;
		for(std::size_t i = 0; i < editorStack.size(); ++i)
		{
			while(i < editorStack.size() && pred(editorStack[i].get()))
			{
				alteredStack = true;
				deleteEditor(editorStack.begin() + i);
			}
		}
		return alteredStack;
	}

	void MainEditor::activeTabClicked(cpl::CTextTabBar<>* obj, int index)
	{
		// TODO: spawn the clicked editor, if it doesn't exist.
		if(ksettings.bGetBoolState())
		{
			// make sure an editor is active
			if (hasCurrentView() && !std::none_of(editorStack.begin(), editorStack.end(), [](auto & edit) { return edit.get()->getName() == MainEditorName; }))
			{
				pushEditor(views[index].getEditorDSO().getCached());
			}

			removeAnyEditor([](juce::Component * e) { return e->getName() == MainEditorName; });
			ksettings.bSetInternal(0);
		}
	}


	void MainEditor::addTab(const std::string & name)
	{
		tabs.addTab(name);
	}

	void MainEditor::enterFullscreenIfNeeded(juce::Point<int> where)
	{
		// full screen set?
		if (kkiosk.bGetValue() > 0.5)
		{

			// avoid storing the current window position into kioskCoords
			// the first time we spawn the view.
			firstKioskMode = true;
			activeView().getWindow()->setTopLeftPosition(where.x, where.y);
			kkiosk.bForceEvent();
		}
	}

	void MainEditor::enterFullscreenIfNeeded()
	{
		// full screen set?
		if (hasCurrentView() && !activeView().getIsFullScreen() && kkiosk.bGetValue() > 0.5)
		{
			// avoid storing the current window position into kioskCoords
			// the first time we spawn the view.
			firstKioskMode = false;
			kkiosk.bForceEvent();
		}
	}

	void MainEditor::exitFullscreen()
	{
		juce::Desktop & instance = juce::Desktop::getInstance();
		if (hasCurrentView() && activeView().getIsFullScreen() && !this->isParentOf(activeView().getWindow()))
		{
			activeView().getWindow()->removeKeyListener(this);
			activeView().getWindow()->removeComponentListener(this);
			if (activeView().getWindow() == instance.getKioskModeComponent())
			{
				juce::Desktop::getInstance().setKioskModeComponent(nullptr);
			}

			activeView().getWindow()->setTopLeftPosition(0, 0);
			addChildComponent(activeView().getWindow());
			activeView().setFullScreenMode(false);

			if (preFullScreenSize.getWidth() > 0 && preFullScreenSize.getHeight() > 0)
			{
				// restores from minimal window
				setBounds(preFullScreenSize);
			}
			else
			{
				resized();
			}
		}

	}

	void MainEditor::serialize(cpl::CSerializer & data, cpl::Version version)
	{

		data << krefreshRate;
		data << krenderEngine;
		data << khelp;
		data << kfreeze;
		data << kidle;
		data << getBounds().withZeroOrigin();
		data << isEditorVisible;
		data << selTab;
		data << kioskCoords;
		data << hasAnyTabBeenSelected;
		data << kkiosk;

		// stuff that gets set the last.
		data << kantialias;
		data << kvsync;
		data << kswapInterval;


		for (auto & colour : colourControls)
		{
			data.getContent("Colours").getContent(colour.bGetTitle()) << colour;
		}

		for (std::size_t i = 0; i < views.size(); ++i)
		{
			data.getContent("Serialized Views").getContent(views[i].getName()) = views[i].getViewDSO().getState();
			data.getContent("Serialized Editors").getContent(views[i].getName()) = views[i].getEditorDSO().getState();
		}

		std::int64_t historySize;
		std::string contents = kmaxHistorySize.getInputValue();
		if (cpl::lexicalConversion(contents, historySize))
			data << std::max(0ll, (long long)historySize);
		else
			data << 1000ll;

		data << khideTabs;

		data << khideWidgets << kstopProcessingOnSuspend;
	}

	void MainEditor::nestedOnMouseMove(const juce::MouseEvent & e)
	{
		auto point = e.getEventRelativeTo(this);
		mouseHoversTabArea = point.y < elementSize + elementBorder;

		if (mouseHoversTabArea)
		{
			setTabBarVisibility(true);
		}
		tabBarTimer = cpl::Misc::QuickTime();
	}

	void MainEditor::nestedOnMouseExit(const juce::MouseEvent & e)
	{
		if (e.eventComponent == this)
		{
			mouseHoversTabArea = false;
			tabBarTimer = cpl::Misc::QuickTime();
		}
	}

	void MainEditor::deserialize(cpl::CSerializer & data, cpl::Version version)
	{
		cpl::CSerializer savedEditorData, savedViewData;

		if (version < cpl::Version(0, 2, 8))
		{
			// before 0.2.8, any state was stored in the view - so copy it to the editor
			savedEditorData = data.getContent("Serialized Views");
		}
		else
		{
			savedEditorData = data.getContent("Serialized Editors");
			savedViewData = data.getContent("Serialized Views");
		}
		//cpl::iCtrlPrec_t dataVal(0);
		juce::Rectangle<int> bounds;

		data >> krefreshRate;
		data >> krenderEngine;
		data >> khelp;
		data >> kfreeze;
		data >> kidle;
		data >> bounds;
		data >> isEditorVisible;
		data >> selTab;
		data >> kioskCoords;
		data >> hasAnyTabBeenSelected;
		// kind of a hack, but we don't really want to enter kiosk mode immediately.
		kkiosk.bRemoveChangeListener(this);
		data >> kkiosk;
		kkiosk.bAddChangeListener(this);

		for (auto & colour : colourControls)
		{
			auto & content = data.getContent("Colours").getContent(colour.bGetTitle());
			if (!content.isEmpty())
				content >> colour;
		}

		// sanitize bounds...
		setBounds(bounds.constrainedWithin(
			juce::Desktop::getInstance().getDisplays().getDisplayContaining(bounds.getPosition()).userArea
		).withZeroOrigin());

		// reinitiate any current views (will not be done through tab selection further down)
		for (auto & viewState : views)
		{
			{
				auto & content = savedViewData.getContent(viewState.getName());
				if (!content.isEmpty())
					viewState.getViewDSO().setState(content, content.getLocalVersion());
			}
			{
				auto & content = savedEditorData.getContent(viewState.getName());
				if (!content.isEmpty())
					viewState.getEditorDSO().setState(content, content.getLocalVersion());
			}
		}

		// will take care of opening the correct view
		if (hasAnyTabBeenSelected)
		{
			// settings this will cause setSelectedTab to
			// enter fullscreen at kioskCoords
			if (kkiosk.bGetValue() > 0.5)
				firstKioskMode = true;

			selTab = std::min(selTab, cpl::enum_cast<int32_t>(ViewTypes::end) - 1);

			if(tabs.getSelectedTab() == selTab)
			{
				// TODO: rewrite the following
				// tabSelected will NOT get called in this case, do it manually.
				tabSelected(&tabs, selTab);
			}

			tabs.setSelectedTab(selTab);
		}
		else
		{
			// if not, make sure we entered full screen if needed.
			enterFullscreenIfNeeded(kioskCoords);
		}


		// set this kind of stuff after view is initiated. note, these cause events to be fired!
		data >> kantialias;
		data >> kvsync;
		data >> kswapInterval;

		if (version >= cpl::Version(0, 2, 8))
		{
			std::int64_t historySize;
			data >> historySize;
			if (historySize > 0)
				kmaxHistorySize.setInputValue(std::to_string(historySize));

			data >> khideTabs;
		}

		if (version >= cpl::Version(0, 3, 1))
		{
			data >> khideWidgets >> kstopProcessingOnSuspend;
		}
	}

	bool MainEditor::stringToValue(const cpl::CBaseControl * ctrl, const cpl::string_ref valString, cpl::iCtrlPrec_t & val)
	{
		double newVal = 0;

		if (ctrl == &krefreshRate)
		{
			if (cpl::lexicalConversion(valString, newVal))
			{
				newVal = cpl::Math::confineTo
				(
					cpl::Math::UnityScale::Inv::exp
					(
						cpl::Math::confineTo(newVal, 10.0, 1000.0),
						10.0,
						1000.0
					),
					0.0,
					1.0
				);
				val = newVal;
				return true;
			}
		}
		else if (ctrl == &kswapInterval)
		{
			if (cpl::lexicalConversion(valString, newVal))
			{
				newVal = cpl::Math::confineTo(newVal, 0, kdefaultMaxSkippedFrames);
				val = cpl::Math::UnityScale::Inv::linear(newVal, 0.0, (double)kdefaultMaxSkippedFrames);
				return true;
			}
		}
		return false;
	}

	bool MainEditor::valueToString(const cpl::CBaseControl * ctrl, std::string & valString, cpl::iCtrlPrec_t val)
	{
		if (ctrl == &krefreshRate)
		{
			auto refreshRateVal = cpl::Math::round<int>(cpl::Math::UnityScale::exp(val, 10.0, 1000.0));
			valString = std::to_string(refreshRateVal) + " ms";
			return true;
		}
		else if (ctrl == &kswapInterval)
		{
			char buf[100];
			sprintf_s(buf, "%d frames", cpl::Math::round<int>(val * kdefaultMaxSkippedFrames));
			valString = buf;
			return true;
		}
		return false;
	}

	bool MainEditor::keyPressed(const juce::KeyPress &key, juce::Component *originatingComponent)
	{
		if (hasCurrentView() && (activeView().getWindow() == originatingComponent))
		{
			if (key.isKeyCode(key.escapeKey) && kkiosk.bGetValue() > 0.5)
			{
				exitFullscreen();

				kkiosk.bSetInternal(0.0);
				return true;
			}
		}
		return false;
	}

	MainEditor::~MainEditor()
	{
		suspendView(views[selTab]);
		notifyDestruction();
		exitFullscreen();
		juce::Timer::stopTimer();
		juce::HighResolutionTimer::stopTimer();

	}

	void MainEditor::resizeEnd()
	{
		resized();
		resume();
	}

	void MainEditor::resizeStart()
	{
		suspend();
	}

	int MainEditor::getViewTopCoordinate() const noexcept
	{
		auto editor = getTopEditor();

		if (editor)
		{
			auto maxHeight = elementSize * 5;
			auto possibleBounds = std::make_pair(getWidth() - elementBorder * 2, maxHeight);
			// content pages knows their own (dynamic) size.
			if (auto signalizerEditor = dynamic_cast<Signalizer::CContentPage *>(editor))
			{
				maxHeight = std::max(0, std::min(maxHeight, signalizerEditor->getSuggestedSize(possibleBounds).second));
			}
			return tabs.getHeight() + maxHeight + elementBorder;
		}
		else
		{
			return tabs.getBottom() + elementBorder;
		}
	}

	void MainEditor::resized()
	{
		auto const buttonSize = elementSize - elementBorder * 2;
		auto const buttonSizeW = elementSize - elementBorder * 2;
		rcc.setBounds(getWidth() - 15, getHeight() - 15, 15, 15);
		// dont resize while user is dragging
		if (rcc.isMouseButtonDown())
			return;
		// resize panel to width
		auto const top = tabBarIsVisible ? elementBorder : -elementSize;

		auto const width = getWidth();
		auto leftBorder = width - elementSize + elementBorder;
		ksettings.setBounds(1, top, buttonSizeW, buttonSize);

		kfreeze.setBounds(leftBorder, top, buttonSizeW, buttonSize);
		leftBorder -= elementSize - elementBorder;

		khelp.setBounds(leftBorder, top, buttonSizeW, buttonSize);
		leftBorder -= elementSize - elementBorder;
		kkiosk.setBounds(leftBorder, top, buttonSizeW, buttonSize);

		tabs.setBounds
		(
			ksettings.getBounds().getRight() + elementBorder,
			top,
			getWidth() - (ksettings.getWidth() + getWidth() - leftBorder + elementBorder * 3),
			elementSize - elementBorder * 2
		);


		auto editor = getTopEditor();
		if (editor)
		{
			// TODO: refactor code to merge with getViewTopCoordinate -> getEditorRect()
			auto maxHeight = elementSize * 5;
			auto possibleBounds = std::make_pair(getWidth() - elementBorder * 2, maxHeight);
			// content pages knows their own (dynamic) size.
			if (auto signalizerEditor = dynamic_cast<Signalizer::CContentPage *>(editor))
			{
				maxHeight = std::max(0, std::min(maxHeight, signalizerEditor->getSuggestedSize(possibleBounds).second));
			}
			editor->setBounds(elementBorder, tabs.getBottom(), possibleBounds.first, maxHeight);
			viewTopCoord = tabs.getHeight() + maxHeight + elementBorder;
		}
		else
		{
			if (tabBarIsVisible)
				viewTopCoord = tabs.getBottom() + elementBorder;
			else
				viewTopCoord = 0;
		}
		// full screen components resize themselves.
		if (hasCurrentView() && !activeView().getIsFullScreen())
		{
			activeView().getWindow()->setBounds(0, viewTopCoord, getWidth(), getHeight() - viewTopCoord);
		}

		//rightButtonOutlines.addRectangle(juce::Rectangle<float>(0.5f, 0.5f, getWidth() - 1.5f, editor ? editor->getBottom() : elementSize - 1.5f));
	}


	void MainEditor::timerCallback()
	{
		for (std::size_t i = 0; i < params->numSetsAndState(); ++i)
		{
			params->getSet(i)->pulseUI();
		}
		if (hasCurrentView())
		{
			auto now = cpl::Misc::QuickTime();

			if (!getTopEditor() && !mouseHoversTabArea && now - tabBarTimer > tabBarTimeout)
				setTabBarVisibility(false);

			if (idleInBack)
			{
				if (!hasKeyboardFocus(true))
					focusLost(FocusChangeType::focusChangedDirectly);
				else if (unFocused)
					focusGained(FocusChangeType::focusChangedDirectly);
			}

			if(!kvsync.bGetBoolState())
				activeView().repaintMainContent();
		}
	}

	void MainEditor::hiResTimerCallback()
	{
		if (hasCurrentView())
		{
			if (idleInBack)
			{
				const juce::MessageManagerLock mml;
				if (!hasKeyboardFocus(true))
					focusLost(FocusChangeType::focusChangedDirectly);
				else if (unFocused)
					focusGained(FocusChangeType::focusChangedDirectly);
			}

			if (!kvsync.bGetBoolState())
				activeView().repaintMainContent();
		}
	}

	void MainEditor::paint(juce::Graphics& g)
	{
		// make sure to paint everything completely opaque.
		g.setColour(cpl::GetColour(cpl::ColourEntry::Separator).withAlpha(1.0f));
		g.fillRect(getBounds().withZeroOrigin().withBottom(viewTopCoord));
		if (kkiosk.bGetValue() > 0.5)
		{
			g.setColour(cpl::GetColour(cpl::ColourEntry::Deactivated).withAlpha(1.0f));
			g.fillRect(getBounds().withZeroOrigin().withTop(viewTopCoord));
			g.setColour(cpl::GetColour(cpl::ColourEntry::AuxillaryText));
			g.drawText("View is fullscreen", getBounds().withZeroOrigin().withTop(viewTopCoord), juce::Justification::centred, true);
		}
	}


	void MainEditor::onOGLRendering(cpl::COpenGLView * view) noexcept
	{
		if (mtFlags.swapIntervalChanged.cas())
		{
			oglc.setSwapInterval(newc.swapInterval.load(std::memory_order_acquire));
			view->setSwapInterval(newc.swapInterval.load(std::memory_order_acquire));
		}
	}

	void MainEditor::onOGLContextCreation(cpl::COpenGLView * view) noexcept
	{
		mtFlags.swapIntervalChanged = true;
	}

	void MainEditor::onOGLContextDestruction(cpl::COpenGLView * view) noexcept
	{

	}

	void MainEditor::showAboutBox()
	{
		khelp.bSetInternal(1);
		using namespace cpl;
		std::string contents =
			programInfo.name + " " + programInfo.version.toString() + newl +
			"Build info: \n" + programInfo.customBuildInfo + newl +
			"Written by Janus Lynggaard Thorborg, (C) 2017" + newl +
			programInfo.name + " is free and open source (GPL v3), see more at the home page: " + newl + "www.jthorborg.com/index.html?ipage=signalizer" + newl + newl +
			"Open the readme file (contains information you must read upon first use)?";

		auto ret = Misc::MsgBox(contents, "About " + programInfo.name, Misc::MsgStyle::sYesNo | Misc::MsgIcon::iInfo);
		switch (ret) {
			case Misc::MsgButton::bYes:
				#ifdef CPL_WINDOWS
					std::system(("Notepad.exe \"" + Misc::DirectoryPath() + "/READ ME.txt\"").c_str());
				#elif defined(CPL_MAC)
					std::string cmdLine = "open \"" + Misc::DirectoryPath() + "/READ ME.txt\"";
					std::system((cmdLine).c_str());
				#else
					std::string cmdLine = "gedit \"" + Misc::DirectoryPath() + "/READ ME.txt\"";
					std::system((cmdLine).c_str());
				#endif
				break;
		}
		khelp.bSetInternal(0);
	}

	void MainEditor::initUI()
	{
		auto & lnf = cpl::CLookAndFeel_CPL::defaultLook();

		// add listeners
		krefreshRate.bAddFormatter(this);
		kfreeze.bAddChangeListener(this);
		kkiosk.bAddChangeListener(this);
		kidle.bAddChangeListener(this);
		krenderEngine.bAddChangeListener(this);
		krefreshRate.bAddChangeListener(this);
		ksettings.bAddChangeListener(this);
		kstableFps.bAddChangeListener(this);
		kvsync.bAddChangeListener(this);
		tabs.addListener(this);
		kmaxHistorySize.bAddChangeListener(this);
		kswapInterval.bAddChangeListener(this);
		kantialias.bAddChangeListener(this);
		khelp.bAddChangeListener(this);
		krefreshState.bAddChangeListener(this);
		kswapInterval.bAddFormatter(this);
		kstopProcessingOnSuspend.bAddChangeListener(this);
		khideWidgets.bAddChangeListener(this);

		// design
		kfreeze.setImage("icons/svg/freeze.svg");
		ksettings.setImage("icons/svg/gears.svg");
		khelp.setImage("icons/svg/help.svg");
		kkiosk.setImage("icons/svg/fullscreen.svg");

		//kstableFps.setSize(cpl::ControlSize::Rectangle.width, cpl::ControlSize::Rectangle.height / 2);
		//kvsync.setSize(cpl::ControlSize::Rectangle.width, cpl::ControlSize::Rectangle.height / 2);
		//krefreshState.setSize(cpl::ControlSize::Rectangle.width, cpl::ControlSize::Rectangle.height / 2);
		kstableFps.setToggleable(true);
		kvsync.setToggleable(true);
		khideTabs.setToggleable(true);
		kstopProcessingOnSuspend.setToggleable(true);
		khideWidgets.setToggleable(true);

		khideTabs.setSingleText("Auto-hide tabs");
		krefreshRate.bSetTitle("Refresh Rate");
		krefreshState.setSingleText("Reset state");
		kantialias.bSetTitle("Antialiasing");
		kidle.setSingleText("Idle in back");
		kswapInterval.bSetTitle("Swap interval");
		kstableFps.setSingleText("Stable FPS");
		kvsync.setSingleText("Vertical Sync");

		kstopProcessingOnSuspend.setSingleText("Suspend processing");
		khideWidgets.setSingleText("Hide widgets");

		// setup
		krenderEngine.setValues(RenderingEnginesList);
		kantialias.setValues(AntialisingStringLevels);

		// initiate colours
		for (unsigned i = 0; i < colourControls.size(); ++i)
		{
			auto & schemeColour = lnf.getSchemeColour(i);
			colourControls[i].setControlColour(schemeColour.colour);
			colourControls[i].bSetTitle(schemeColour.name);
			colourControls[i].bSetDescription(schemeColour.description);
			colourControls[i].bAddChangeListener(this);

		}

		// add stuff
		addAndMakeVisible(ksettings);
		addAndMakeVisible(kfreeze);
		addAndMakeVisible(khelp);
		addAndMakeVisible(kkiosk);
		addAndMakeVisible(tabs);

		tabs.setOrientation(tabs.Horizontal);
		for(auto & viewName : ViewIndexToMap)
			tabs.addTab(viewName);

		// additions
		addAndMakeVisible(rcc);
		rcc.setAlwaysOnTop(true);
		// TODO: reattach?
		//currentView = &defaultView; // note: enables callbacks on value set in this function
		addAndMakeVisible(defaultView);
		krefreshRate.bSetValue(0.12);


		// descriptions
		kstableFps.bSetDescription("Stabilize frame rate using a high precision timer.");
		kvsync.bSetDescription("Synchronizes graphic view rendering to your monitors refresh rate.");
		kantialias.bSetDescription("Set the level of hardware antialising applied.");
		krefreshRate.bSetDescription("How often the view is redrawn.");
		khelp.bSetDescription("About this program");
		kkiosk.bSetDescription("Puts the view into fullscreen mode. Press Escape to untoggle, or tab out of the view.");
		kidle.bSetDescription("If set, lowers the frame rate of the view if this plugin is not in the front.");
		ksettings.bSetDescription("Open the global settings for the plugin (presets, themes and graphics).");
		kfreeze.bSetDescription("Stops the view from updating, allowing you to examine the current point in time.");
		krefreshState.bSetDescription("Resets any processing state in the active view to the default.");
		kmaxHistorySize.bSetDescription("The maximum audio history capacity, set in the respective views. No limit, so be careful!");
		kswapInterval.bSetDescription("Determines the swap interval for the graphics context; a value of zero means the graphics will"
			"update as fast as possible, a value of 1 means it updates synced to the vertical sync, a value of N means it updates every Nth vertical frame sync.");
		khideTabs.bSetDescription("Auto-hides the top tabs and buttons when not used.");
		kstopProcessingOnSuspend.bSetDescription("If set, only the selected running view will process audio - improves performance, but views are out of sync when frozen");
		khideWidgets.bSetDescription("Hides widgets on the screen (frequency trackers, for instance) when the mouse leaves the editor");


		// initial values that should be through handlers
		// TODO: remove if changed to parameter
		kmaxHistorySize.setInputValue("1000");


		resized();
	}
};
