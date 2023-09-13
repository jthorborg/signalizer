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

	file:GraphEditor.cpp

		Implementation of GraphEditor.h

*************************************************************************************/

#include "GraphEditor.h"
#include "MainEditor.h"

namespace Signalizer
{
	class Content 
		: public juce::Component
		, public juce::AsyncUpdater
		, private cpl::CBaseControl::Listener
		, private juce::TextEditor::Listener
		, private juce::KeyListener
	{
	public:

		Content(GraphEditor& parent, HostGraph& graph)
			: graph(graph), parent(parent), nameField("Name"), filterField("Filter"), isMouseDown(false), identityButton(*this)
		{
			addAndMakeVisible(nameField);
			addAndMakeVisible(filterField);

			addChildComponent(identityButton);

			nameField.bAddChangeListener(this);
			filterField.addListener(this);
			filterField.addKeyListener(this);
			model = graph.getModel();
			reflectModel();
		}

	private:

		void reflectModel()
		{
			auto& name = model.nodes[model.hostIndex].name;
			parent.setName(name);

			if (nameField.getInputValue() != name)
			{
				nameField.setInputValueInternal(name);
			}

			relayoutAllEdges();

			filterField.setVisible(!model.isAlias);
			identityButton.setVisible(model.isAlias);
		}

		static constexpr int textEditorHeight = 40;
		static constexpr int space = 5;
		static constexpr int offsetForGraph = textEditorHeight * 2 + space * 2;
		static constexpr float boundarySize = 2;

		struct ImmediateLayout
		{
			static constexpr int pinSpacing = 15;
			static constexpr int pinSize = 5;

			static juce::Rectangle<float> expandPoint(juce::Point<float> pos)
			{
				return juce::Rectangle<float>().withCentre(pos).expanded(pinSize, pinSize);
			}

			template<typename Derived>
			struct NodeBase
			{
				NodeBase(HostGraph::Model* model, std::size_t index)
					: model(model), port(-1), index(index)
				{

				}

				juce::Rectangle<float> getRect() const noexcept
				{
					return rect;
				}

				juce::Point<float> getPinPosition() const noexcept
				{
					return static_cast<const Derived*>(this)->getPinPositionFor(getCurrentPin());
				}

				juce::Rectangle<float> getPinRect() const noexcept
				{
					return expandPoint(getPinPosition());
				}

				juce::Rectangle<float> getPinRectFor(PinInt p) const noexcept
				{
					return expandPoint(static_cast<Derived*>(this)->getPinPositionFor(p));
				}

				PinInt getCurrentPin() const noexcept
				{
					return port;
				}

				const std::string& currentName() const noexcept
				{
					return model->nodes[index].name;
				}

				HostGraph::SerializedHandle getHandle() const noexcept
				{
					return model->nodes[index].node;
				}

				static int calculateHeight(int numPins) noexcept
				{
					return std::max(30, (numPins + 1) * pinSpacing);
				}

				static int pinOffset(int pinIndex) noexcept
				{
					return (pinIndex + 1) * pinSpacing;
				}

				int port;
				std::size_t index;
				HostGraph::Model* model;
				juce::Rectangle<float> rect;
			};

			struct Host : public NodeBase<Host>
			{
				using Base = NodeBase<Host>;

				Host(HostGraph::Model& model, juce::Rectangle<float> bounds)
					: Base(&model, 0)
				{
					auto width = bounds.getWidth() / 3;

					rect = juce::Rectangle<float>(
						bounds.getRight() - (space + width),
						bounds.getY() + space * 2,
						width,
						calculateHeight(Signalizer::MaxInputChannels)
					);
				}
				
				juce::Point<float> getPinPositionFor(PinInt port) const noexcept
				{
					return juce::Point<float>(rect.getX(), rect.getY() + pinOffset(port));
				}

				bool moveNextPort()
				{
					return ++port < Signalizer::MaxInputChannels;
				}
			};

			struct Node : public NodeBase<Node>
			{
				using Base = NodeBase<Node>;

				Node()
					: Base(nullptr, 0), con(-1), hasConnections(false)
				{

				}

				Node(HostGraph::Model& model, std::size_t index, float previousBottom, juce::Rectangle<float> bounds)
					: Base(&model, index), hasConnections(false), con(-1)
				{
					auto ports = model.nodes[index].portCount;
					auto height = calculateHeight(ports);

					rect = bounds
						.translated(space, 0)
						.withY(previousBottom + space * 2)
						.withSize(bounds.getWidth() / 3, height);
				}

				const HostGraph::Model::NodeView& getView() const noexcept
				{
					return model->nodes[index];
				}

				juce::Rectangle<float> getPinRect() const noexcept
				{
					return expandPoint(getPinPosition());
				}

				juce::Point<float> getPinPositionFor(PinInt port) const noexcept
				{
					return juce::Point<float>(rect.getRight(), rect.getY() + pinOffset(port));
				}

				bool moveNextPort()
				{
					con = -1;

					std::size_t pin = static_cast<std::size_t>(++port);

					if (pin >= model->nodes[index].portCount)
						return false;

					hasConnections = model->nodes[index].connectionCount > 0;

					return true;
				}

				bool moveNextConnection()
				{
					std::size_t connection = static_cast<std::size_t>(++con);

					if (connection >= model->nodes[index].connectionCount)
						return false;

					return true;
				}

				auto getCurrentConnection(const Host& host) const
				{
					checkBoundsOfConnection();

					auto connection = model->connections[model->nodes[index].connectionOffset + static_cast<std::size_t>(con)];

					return std::make_tuple(
						getPinPositionFor(connection.Source), 
						host.getPinPositionFor(connection.Destination),
						connection
					);
				}

				bool isPortConnected() const noexcept
				{
					return hasConnections;
				}

				void checkBoundsOfPort() const
				{
					CPL_RUNTIME_ASSERTION(static_cast<std::size_t>(port) < model->nodes[index].portCount);
				}

				void checkBoundsOfConnection() const
				{
					CPL_RUNTIME_ASSERTION(static_cast<std::size_t>(con) < model->nodes[index].connectionCount);
				}

				bool hasConnections;
				int con;
			};

			ImmediateLayout(HostGraph::Model& model, juce::Rectangle<float> bounds, std::string& filter) 
				: model(model), index(-1), bounds(bounds), filter(filter)
			{
			}


			Host getHost() const noexcept
			{
				return Host(model, bounds);
			}


			bool moveNextNode(Node& result)
			{
				std::size_t current;

				do
				{
					current = static_cast<std::size_t>(++index);
					if (current >= model.nodes.size())
						return false;
				} while (!filter.empty() && model.nodes[current].name.find(filter) == std::string::npos);


				auto previousBottom = result.model ? result.getRect().getBottom() : bounds.getY();

				result = Node(model, current, previousBottom, bounds);

				return true;
			}


			int index;
			juce::Rectangle<float> bounds;
			HostGraph::Model& model;
			std::string& filter;
		};

		template<typename TNode>
		void paintNode(juce::Graphics& g, TNode n, juce::Colour base, juce::Colour baseTextColour)
		{
			// draw main rect 
			juce::Rectangle<float> rect = n.getRect();

			juce::Colour rectFill = base;
			juce::Colour textColour = baseTextColour;

			if (!sourceDrag.has_value() && rect.withTrimmedRight(ImmediateLayout::pinSize).withTrimmedLeft(ImmediateLayout::pinSize).contains(lastControlPosition))
			{
				rectFill = rectFill.contrasting(isMouseDown ? 0.2f : 0.3f);
				textColour = cpl::GetColour(cpl::ColourEntry::SelectedText);
			}

			juce::Colour boundary = rectFill.contrasting(0.3f);


			// highlight?

			g.setColour(rectFill);
			g.fillRoundedRectangle(rect, 5);

			g.setColour(boundary);
			g.drawRoundedRectangle(rect, 5, boundarySize);

			g.setColour(textColour);

			g.drawText(
				n.currentName(),
				rect.reduced(5).withTrimmedLeft(5),
				juce::Justification::centredLeft
			);

			// draw pins.

			while (n.moveNextPort())
			{
				juce::Rectangle<float> pinRect = n.getPinRect();

				bool highLight = pinRect.contains(lastControlPosition) || pinRect.contains(pointOfInterest);

				g.setColour(highLight ? boundary.contrasting(isMouseDown ? 0.5f : 0.3) : boundary);

				g.fillEllipse(pinRect);
			}
		}

		auto graphRect()
		{
			return getBounds().toFloat().removeFromTop(offsetForGraph).translated(0, offsetForGraph);
		}

		template<typename Callable>
		void forEachConnection(Callable c)
		{
			ImmediateLayout layout(model, graphRect(), filter);
			ImmediateLayout::Node n;
			auto host = layout.getHost();
			std::size_t edgeCounter{};

			while (layout.moveNextNode(n))
			{
				while (n.moveNextConnection())
				{
					auto current = n.getCurrentConnection(host);
					if(!c(std::get<0>(current), std::get<1>(current), std::make_pair(n.getView(), std::get<2>(current))))
						return;
				}
			}
		}

		void relayoutAllEdges()
		{
			edges.resize(0);
			
			forEachConnection(
				[this](auto sourcePos, auto destPos, auto pair)
				{
					layoutEdge(edges.emplace_back(), sourcePos, destPos);
					return true;
				}
			);
		}

		void paint(juce::Graphics& g) override
		{
			std::call_once(focusGrab, [&] { filterField.grabKeyboardFocus(); });

			g.fillAll(cpl::GetColour(model.isAlias ? cpl::ColourEntry::Error : cpl::ColourEntry::Normal));

			if (model.isAlias)
				return;

			const auto base = cpl::GetColour(cpl::ColourEntry::Normal);
			const auto baseTextColour = cpl::GetColour(cpl::ColourEntry::ControlText);
			const auto edgeColour = juce::Colours::violet;
			ImmediateLayout layout(model, graphRect(), filter);

			ImmediateLayout::Node n;

			juce::PathStrokeType pst(3, juce::PathStrokeType::JointStyle::curved, juce::PathStrokeType::EndCapStyle::rounded);

			g.setColour(edgeColour);

			for (auto& p : edges)
			{
				g.strokePath(p, pst);
			}

			if (sourceDrag.has_value())
			{
				g.setColour(edgeColour.brighter());
				g.strokePath(sourceDrag->edge, pst);
			}

			while (layout.moveNextNode(n))
			{
				paintNode(g, n, base, baseTextColour);
			}

			paintNode(g, layout.getHost(), base, baseTextColour);
		}

		void resized() override
		{
			nameField.setBounds(getBounds().withHeight(textEditorHeight));

			// these overlap, but are never shown together.
			const auto filterSpace = nameField.getBounds().withY(nameField.getBottom() + space);
			filterField.setBounds(filterSpace);
			identityButton.setBounds(filterSpace.reduced(5, 0));

			relayoutAllEdges();
		}

		void mouseMove(const juce::MouseEvent& e) override
		{
			lastControlPosition = e.position;
			isMouseDown = false;
			repaint();
		}

		static void layoutEdge(juce::Path& edge, const juce::Point<float>& p1, const juce::Point<float>& p2)
		{
			edge.clear();
			edge.startNewSubPath(p1);

			edge.cubicTo(
				p1.x + (p2.x - p1.x) * 0.4f, p1.y,
				p2.x - (p2.x - p1.x) * 0.4f, p2.y,
				p2.x, p2.y
			);

			edge.setUsingNonZeroWinding(true);
		}

		void mouseDrag(const juce::MouseEvent& e) override
		{
			pointOfInterest = {};

			if (sourceDrag.has_value())
			{
				auto& edge = sourceDrag->edge;
				auto position = e.position;

				// Snap to close pins.
				ImmediateLayout layout(model, graphRect(), filter);
				auto n = layout.getHost();

				while (n.moveNextPort())
				{
					juce::Rectangle<float> pinRect = n.getPinRect().expanded(3);

					if (pinRect.contains(position))
					{
						pointOfInterest = position = n.getPinPosition();
						break;
					}
				}

				layoutEdge(edge, lastControlPosition, position);
			}

			repaint();
		}

		void mouseDown(const juce::MouseEvent& e) override
		{
			lastControlPosition = e.position;

			if (e.mods.isLeftButtonDown())
			{
				repaint();

				{
					ImmediateLayout layout(model, graphRect(), filter);

					ImmediateLayout::Node n;

					// Test if we're drawing a new edge
					while (layout.moveNextNode(n))
					{
						while (n.moveNextPort())
						{
							juce::Rectangle<float> pinRect = n.getPinRect();

							if (!pinRect.contains(lastControlPosition))
								continue;

							lastControlPosition = n.getPinPosition();
							sourceDrag.emplace(n.getHandle(), n.getCurrentPin());
							return;
						}
					}
				}

				// Test if we're redragging an existing edge
				forEachConnection(
					[this](auto sourcePos, auto destPos, const auto& con)
					{
						if (ImmediateLayout::expandPoint(destPos).contains(lastControlPosition))
						{
							lastControlPosition = sourcePos;
							// start dragging this edge, and disconnect the original connection
							sourceDrag.emplace(con.first.node, con.second.Source);
							layoutEdge(sourceDrag->edge, lastControlPosition, destPos);
							graph.disconnect(con.first.node, con.second);
							return false;
						}

						return true;
					}
				);
			}
		}


		void mouseUp(const juce::MouseEvent& e) override
		{
			isMouseDown = false;
			if (sourceDrag.has_value())
			{
				auto& edge = sourceDrag->edge;
				auto position = e.position;

				// Snap to close pins.				
				ImmediateLayout layout(model, graphRect(), filter);
				auto n = layout.getHost();

				while (n.moveNextPort())
				{
					juce::Rectangle<float> pinRect = n.getPinRect().expanded(3);

					if (pinRect.contains(position))
					{
						connectionRequest(n.getCurrentPin());
						break;
					}
				}
			}
			else
			{
				ImmediateLayout layout(model, graphRect(), filter);

				// test if we're clicking the host
				if (layout.getHost().getRect().contains(lastControlPosition))
				{
					// toggle all visible nodes
					std::vector<HostGraph::SerializedHandle> refs;
					ImmediateLayout::Node n;

					while (layout.moveNextNode(n))
					{
						refs.emplace_back(n.getHandle());
					}

					graph.toggleSet(refs);
				}
				else
				{
					// test if we're clicking any nodes
					ImmediateLayout::Node node;

					while (layout.moveNextNode(node))
					{
						if (node.getRect().contains(lastControlPosition))
						{
							graph.toggleSet({ node.getHandle() });
							break;
						}
					}
				}
			}

			sourceDrag.reset();
			pointOfInterest = { };
			lastControlPosition = e.position;

			repaint();
		}

		void connectionRequest(PinInt destination)
		{
			CPL_RUNTIME_ASSERTION(sourceDrag.has_value());
			graph.connect(sourceDrag->SourceNode, { sourceDrag->SourcePin, destination });
		}

		// Inherited via AsyncUpdater
		virtual void handleAsyncUpdate() override
		{
			graph.updateModel(model);
			reflectModel();
			repaint();
		}

		void assumeIdentity()
		{
			graph.assumeNonAliasedIdentity();
		}

		class Button : public cpl::CButton
		{
		public: 

			Button(Content& parent)
				: content(parent)
			{
				bSetDescription("This can temporarily happen when loading projects or pasting plugins / presets on Signalizers. When clicked, this plugin assumes a new identity.");
				setTexts("This Signalizer was duplicated, click to fix", "This Signalizer was duplicated, click to fix");
			}

			void clicked() override
			{
				content.assumeIdentity();
				cpl::CButton::clicked();
			}

			Content& content;
		};

		HostGraph& graph;
		HostGraph::Model model;
		GraphEditor& parent;
		cpl::CInputControl nameField;
		Button identityButton;
		juce::TextEditor filterField;
		juce::Point<float> lastControlPosition, pointOfInterest;
		std::vector<juce::Path> edges;
		std::string filter;
		std::once_flag focusGrab;
		
		struct DraggedEdge
		{
			DraggedEdge(const HostGraph::SerializedHandle& handle, PinInt source) : SourceNode(handle), SourcePin(source) {}

			HostGraph::SerializedHandle SourceNode;
			PinInt SourcePin;

			juce::Path edge;
		};

		std::optional<DraggedEdge> sourceDrag;
		bool isMouseDown;

		// Inherited via Listener
		virtual void onObjectDestruction(const cpl::CBaseControl::ObjectProxy& destroyedObject) override
		{
		}

		virtual void valueChanged(const cpl::CBaseControl* c) override
		{
			if (c == &nameField)
			{
				graph.setName(nameField.getInputValue());
			}
		}

		virtual void textEditorTextChanged(juce::TextEditor& te) override
		{
			if (&te == &filterField)
			{
				filter = te.getText().toStdString();
				repaint();
				relayoutAllEdges();
			}
		}

		// Inherited via KeyListener
		virtual bool keyPressed(const juce::KeyPress& key, juce::Component* originatingComponent) override
		{
			// TODO: intercept keystrokes from field to toggle on special combos
			return false;
		}
};

	GraphEditor::GraphEditor(MainEditor* editor, HostGraph& h)
		: juce::DocumentWindow("Graph editor", juce::Colours::aliceblue, juce::DocumentWindow::TitleBarButtons::allButtons)
		, editor(editor)
		, host(h)
		, content(std::make_shared<Content>(*this, h))
	{
		setUsingNativeTitleBar(true);
		setBounds(editor->getScreenBounds().withSizeKeepingCentre(400, 600));
		setVisible(true);

		setContentNonOwned(content.get(), false);
		h.addModelListener(content);
	}

	GraphEditor::~GraphEditor()
	{
		if (editor)
			editor->graphEditorDied();
	}

	void GraphEditor::mainEditorDied()
	{
		editor = nullptr;
		delete this;
	}

	void GraphEditor::closeButtonPressed()
	{
		delete this;
	}
}

