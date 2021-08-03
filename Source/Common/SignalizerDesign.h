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
 
	file:SignalizerDesign.h

		Interface for specific design components and other UI for Signalizer.
 
*************************************************************************************/

#ifndef SIGNALIZER_SIGNALIZERDESIGN_H
	#define SIGNALIZER_SIGNALIZERDESIGN_H

	#include <cpl/Common.h>
	#include <string>
	#include <vector>
	#include <cpl/gui/GUI.h>
	#include <cpl/gui/NewStuffAndLook.h>
	#include <map>
	#include <ctime>

	namespace Signalizer
	{


		#if _DEBUG
			typedef cpl::CGreenLineTester DummyComponent;
		#else
			typedef juce::Component DummyComponent;
		#endif

		class StateEditor
			: public DummyComponent
			, public cpl::DestructionNotifier
		{
		public:

			/// <summary>
			/// Returns an object to be used for serializing the editor state
			/// </summary>
			virtual cpl::SafeSerializableObject & getEditorSO() = 0;

		};


		class CContentPage
			: public StateEditor
			, public cpl::CIconTabBar::CTabBarListener
		{
		public:

			class SeparatorSection
			:
				public juce::Component
			{
			
			public:
				void paint(juce::Graphics & g) override
				{
					//g.setColour(juce::Colours::black);
					g.setColour(cpl::GetColour(cpl::ColourEntry::Separator));
					// antialias ends
					g.drawLine(getWidth() / 2.f, getHeight() * 0.95f, getWidth() / 2.f, getHeight() * 0.1f, 0.5f);
				}
			};

			class COrderedTabPage
			:
				public DummyComponent
			{

			public:

				void addSection(juce::Component * section, const std::string & name = "", bool acquireOwnership = true)
				{
					if (section == this)
						CPL_RUNTIME_EXCEPTION("Tab page added to self!");

					if (section)
					{
						sections.emplace_back(section, name, acquireOwnership);
						addAndMakeVisible(section);
						auto separator = new SeparatorSection();
						addAndMakeVisible(separator);
						separators.emplace_back(separator);

						arrange(false);
					}
				}

				void arrange(bool fromResized = true)
				{

					int offsetX = 5;
					int posX = 0;
					int posY = 0;
					int sHeight = 0, sWidth = 0;

					if (sections.size())
					{
						posX = 0;
						jassert(sections.size() == separators.size());
						for (unsigned int i = 0; i < sections.size(); ++i)
						{
							auto sectionBounds = sections[i].component->getBounds();
							sections[i].component->setTopLeftPosition(posX, posY);
							posX += (offsetX - 1) + sections[i].component->getWidth();
							//posY += offsetY + section.first->getY();
							sHeight = std::max(sHeight, sections[i].component->getHeight());
							// position the separators
							separators[i]->setTopLeftPosition(posX, posY);
							// maybe add one more pixel? - done!
							posX += (offsetX + 1);
						}
						sWidth = posX;
					}
					suggestedSize = std::make_pair(sWidth, sHeight);
					if (!fromResized)
					{
						setSize(suggestedSize.first, suggestedSize.second);
					}
				}

				void resized() override
				{
					arrange();
					for (auto & separator : separators)
					{
						separator->setSize(2, getHeight());
					}
				}

				virtual ~COrderedTabPage()
				{
					for (auto & section : sections)
					{
						if (section.hasOwnership)
							delete section.component;
					}
				}
				std::pair<int, int> getSuggestedSize() { return suggestedSize; }

			private:

				struct ContentEntry
				{
					ContentEntry() {}
					ContentEntry(juce::Component * c, const std::string & n, bool b) : component(c), name(n), hasOwnership(b) {}
					juce::Component * component;
					std::string name;
					bool hasOwnership;
				};
				std::vector<ContentEntry> sections;
				std::vector<std::unique_ptr<juce::Component>> separators;
				std::pair<int, int> suggestedSize;

			};

			class MatrixSection
			:
				public DummyComponent
			{
			public:

				MatrixSection()
					: suggestedWidth(), suggestedHeight(), spaceAfterLargest(true)
				{

				}

				cpl::CBaseControl * operator [](const std::string & name) const
				{
					for (unsigned y = 0; y < controls.size(); ++y)
					{
						for (unsigned x = 0; x < controls[y].size(); ++x)
						{
							if (cpl::CBaseControl * c = controls[y][x].first)
							{
								if (c->bGetTitle() == name)
									return c;
							}
						}
					}
					return nullptr;
				}

				void resized() override
				{
					arrange();
				}

				void setSpacesAfterLargestElement(bool trigger = true) noexcept
				{
					spaceAfterLargest = trigger;
				}

				void arrange(bool fromResized = true)
				{
					// the general amount of separation between elements
					const std::size_t sepY(5), sepX(5);
					// current x,y offsets
					std::size_t offX(0), offY(0);
					suggestedHeight = 0;
					suggestedWidth = 0;
					std::size_t maxHeightInPrevRow = 0;
					for (std::size_t y = 0; y < controls.size(); ++y)
					{
						maxHeightInPrevRow = 0;
						offX = 0;
						offY += sepY;

						for (std::size_t x = 0; x < controls[y].size(); ++x)
						{
							// get max height for next row
							auto bounds = controls[y][x].first->bGetView()->getBounds();
							maxHeightInPrevRow = 
								spaceAfterLargest ? std::max<std::size_t>(maxHeightInPrevRow, bounds.getHeight()) : bounds.getHeight();
							offX += sepX;

							controls[y][x].first->bGetView()->setTopLeftPosition(static_cast<int>(offX), static_cast<int>(offY));
							offX += bounds.getWidth();
						}

						offY += maxHeightInPrevRow;

						suggestedWidth = std::max(suggestedWidth, offX);
					}
					suggestedHeight += offY;
					suggestedHeight += sepY;
					suggestedWidth += sepX;
					if (!fromResized)
					{
						setSize(static_cast<int>(suggestedWidth), static_cast<int>(suggestedHeight));
					}
				}

				void addControl(cpl::CBaseControl * c, std::size_t row, bool takeOwnerShip = false)
				{
					if (!c)
						return;
					// check if it already exists

					controls.resize(std::max(row + 1, controls.size()));

					for (unsigned x = 0; x < controls[row].size(); ++x)
					{
						if (c == controls[row][x].first)
						{
							return;
						}
					}
					// okay, add it.
					controls[row].emplace_back(c, takeOwnerShip);
					addAndMakeVisible(c->bGetView());
					arrange(false);
				}

				std::pair<int, int> getSuggestedSize()
				{
					return std::make_pair(static_cast<int>(suggestedWidth), static_cast<int>(suggestedHeight));
				}
				virtual ~MatrixSection()
				{
					for (auto & row : controls)
					{
						for (auto & col : row)
						{
							// second == takeOwnerShip
							if (col.second)
								delete col.first;
						}
					}
				}

			private:

				bool spaceAfterLargest;
				std::vector<std::vector<std::pair<cpl::CBaseControl *, bool>>> controls;
				std::size_t suggestedHeight, suggestedWidth;
			};

			COrderedTabPage * addPage(const std::string name, const std::string iconName = "")
			{
				icons.addTab(iconName.size() ? iconName : name);
				nameToIndex.push_back(name);
				auto & page = pages[name];
				// first added page, select it
				if (nameToIndex.size() == 1)
					tabSelected(&icons, 0);
				return &page;
			}


			COrderedTabPage * operator [](const std::string & entry)
			{
				auto it = pages.find(entry);

				if (it != pages.end())
					return &it->second;

				return nullptr;
			}

			std::pair<int, int> getSuggestedSize()
			{
				auto ret = std::make_pair(0, 0);
				for (auto & pageCombo : pages)
				{
					auto bounds = pageCombo.second.getSuggestedSize();
					ret.first = std::max(bounds.first, ret.first);
					ret.second = std::max(bounds.second, ret.second);
				}
				return ret;
				//return selectedComponent ? selectedComponent->getSuggestedSize() : std::make_pair(0, 0);
			}

			std::pair<int, int> getSuggestedSize(std::pair<int, int> possibleBounds)
			{
				auto ret = std::make_pair(0, 0);
				for (auto & pageCombo : pages)
				{
					auto bounds = pageCombo.second.getSuggestedSize();
					ret.first = std::max(bounds.first, ret.first);
					ret.second = std::max(bounds.second, ret.second);
				}

				auto const elementSize = 25;
				auto const elementBorder = 1;

				if (ret.first > possibleBounds.first - elementSize + elementBorder)
					ret.second += 20;
				return ret;
				//return selectedComponent ? selectedComponent->getSuggestedSize() : std::make_pair(0, 0);
			}

			CContentPage()
				: selectedComponent(nullptr)
			{
				icons.setOrientation(icons.Vertical);
				icons.addListener(this);
				addAndMakeVisible(contents);
				addAndMakeVisible(icons);
			}

			void paint(juce::Graphics & g) override
			{
				g.fillAll(cpl::GetColour(cpl::ColourEntry::Activated));
				g.setColour(cpl::GetColour(cpl::ColourEntry::Separator));
				g.drawHorizontalLine(getHeight() - 1, icons.getRight(), getWidth());
				return;
			}

		protected:
			void resized() override
			{
				auto const elementSize = 25;
				auto const elementBorder = 1;

				icons.setBounds(0, 0, elementSize - elementBorder, getHeight());
				contents.setBounds(elementSize, 0, getWidth() - elementSize + elementBorder, getHeight() - elementBorder);
				if (selectedComponent)
				{
					auto suggestedSize = selectedComponent->getSuggestedSize();
					auto internalHeight = (getWidth() - elementSize + elementBorder) < suggestedSize.first ? suggestedSize.second : std::min(suggestedSize.second, getHeight());

					selectedComponent->setBounds(0, 0, suggestedSize.first, internalHeight - elementBorder);
					//selectedComponent->setBounds(0, 0, getWidth() - elementSize + elementBorder, getHeight() - elementBorder);
				}

			}

			virtual void tabSelected(BaseBarType * object, int index) override
			{
				if (object == &icons && (int)nameToIndex.size() > index)
				{

					selectedComponent = &pages[nameToIndex[index]];
					
					auto currentComponent = contents.getViewedComponent();

					contents.setViewedComponent(selectedComponent, false);

#if JUCE_MAJOR_VERSION < 4 && JUCE_MINOR_VERSION < 2
					// TODO: update this when juce actually removes the old component.

					if (currentComponent != nullptr && selectedComponent != currentComponent)
					{
						for (int i = 0; i < contents.getNumChildComponents(); ++i)
						{
							if (auto possibleContainer = contents.getChildComponent(i))
							{
								if (possibleContainer->getIndexOfChildComponent(currentComponent) != -1)
								{
									possibleContainer->removeChildComponent(currentComponent);
									break;
								}
							}
							else
							{
								break;
							}
						}
					}
#endif

					resized();
				}
			};

		private:
			COrderedTabPage * selectedComponent;
			juce::Viewport contents;
			std::map<std::string, COrderedTabPage> pages;
			std::vector<std::string> nameToIndex;
			cpl::CIconTabBar icons;
		};

		class CDefaultView
		: 
			public cpl::COpenGLView
		{
		public:
			CDefaultView()
				: lastTime(cpl::Misc::TimeCounter())
				, fractionateMoves(0)
				, COpenGLView("Signalizer default view")
			{
				std::srand(std::time(nullptr));
				setOpaque(true);
				addAndMakeVisible(bouncer);
				bouncer.setText("No view selected");
				bouncer.collision();
			}
			
			~CDefaultView()
			{
				notifyDestruction();
			}
			
			void resume() override
			{
				if(isOpenGL())
				{
					getAttachedContext()->setComponentPaintingEnabled(true);
				}
			}
			
			void suspend() override
			{
				if(isOpenGL())
				{
					// TODO: ensure old state is reversed.
					getAttachedContext()->setComponentPaintingEnabled(false);
				}
			}
			void onOpenGLRendering() override
			{
				repaintMainContent2();

				juce::OpenGLHelpers::clear(juce::Colours::black);
			}

			void onGraphicsRendering(juce::Graphics & g) override
			{
				if (!isOpenGL())
				{
					g.fillAll(juce::Colours::black);
					repaintMainContent2();
				}
			}

			void setSpeed(int msToMoveAPixel)
			{
				timerSpeed = msToMoveAPixel;
			}

			void resized() override
			{
				if (firstMove && getWidth() && getHeight())
				{
					bouncer.setTopLeftPosition(std::rand() % getWidth(), std::rand() % getHeight());
				}
			}
			void repaintMainContent2()
			{
				// calculate how much we moved

				double msSinceLastCall = cpl::Misc::TimeDifference(lastTime);
				// only repaint and such if we should actually move
				if (msSinceLastCall > timerSpeed)
				{
					lastTime = cpl::Misc::TimeCounter();
					auto preciseMoves = fractionateMoves + msSinceLastCall / timerSpeed;
					auto numMoves = static_cast<int>(floor(preciseMoves));
					fractionateMoves = preciseMoves - numMoves;
					for (int i = 0; i < numMoves; ++i)
					{
						move();
					}

					// should repaint
					cpl::GUIUtils::MainEvent(*this, [this] { bouncer.setTopLeftPosition(bounds.getPosition()); });
				}
				//bouncer.repaint();
			}
			void move()
			{
				bounds.setPosition(bounds.getX() + vx, bounds.getY() + vy);
				bounds.setSize(bouncer.getWidth(), bouncer.getHeight());
				bool collided = false;
				if (bounds.getRight() >= getWidth()) // if (x + w == getWidth()) vx = -vx;
				{
					vx = -vx;
					bounds.setX(getWidth() - bounds.getWidth());
					collided = true;
				}
				if (bounds.getBottom() >= getHeight()) // if (y + h == H) vy = -vy;
				{
					vy = -vy;
					bounds.setY(getHeight() - bounds.getHeight());
					collided = true;
				}
				if (bounds.getX() <= 0)
				{
					vx = -vx;
					bounds.setX(0);
					collided = true;
				}
				if (bounds.getY() <= 0)
				{
					vy = -vy;
					bounds.setY(0);
					collided = true;
				}
				if (collided)
					bouncer.collision();
			}
			class Bouncer
				: public juce::Component
			{
			public:

				Bouncer()
				{

					this->setComponentEffect(&glow);
				}

				void paint(juce::Graphics & g) override
				{

					g.setFont(cpl::CLookAndFeel_CPL::defaultLook().getStdFont().withHeight(fontHeight - 2));
					g.setColour(colour);
					g.drawText(text, getBounds().withZeroOrigin(), juce::Justification::centred);

				}

				void collision()
				{
					colour = juce::Colour
					(
						std::uint8_t(std::rand() % 0xFF), 
						std::uint8_t(std::rand() % 0xFF), 
						std::uint8_t(std::rand() % 0xFF), 
						std::uint8_t(0xFF)
					);
					glow.setGlowProperties(2, colour.darker());
				}

				void setText(const std::string & newText)
				{
					text = newText;
					setSize((int)cpl::CLookAndFeel_CPL::defaultLook().getStdFont().withHeight(fontHeight).getStringWidth(text), (int)std::ceil(fontHeight));

				}
				void setTextSize(float height)
				{
					fontHeight = height;
					setText(text);
				}

				juce::Colour colour;
				float fontHeight = 20;

				std::string text;

				juce::GlowEffect glow;
			};

		private:
			int vx = 1, vy = 1;
			int timerSpeed = 10;
			long long lastTime;
			double fractionateMoves;
			juce::Rectangle<int> bounds;
			Bouncer bouncer;
			bool firstMove = false;
		};




	};
#endif 
