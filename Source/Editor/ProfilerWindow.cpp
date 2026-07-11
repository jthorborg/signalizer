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

#include <algorithm>
#include <cctype>
#include <memory>

#include <cpl/profiling/ProfilingModel.h>
#include <cpl/gui/controls/CValueKnobSlider.h>

#include "ProfilerWindow.h"
#include "MainEditor.h"
#include "../Common/SharedBehaviour.h"

namespace Signalizer
{
	struct EWMALaneJuceRenderer
	{
		friend struct cpl::Profiling::EWMAModel::LaneData; // Give access to operator()
		using Model = cpl::Profiling::EWMAModel;
		using Scalar = float;

		constexpr static int space = 1;

	public:

		EWMALaneJuceRenderer(
			Model::LaneData data, 
			juce::Rectangle<int> bounds, 
			juce::Point<Scalar> zoom = juce::Point<Scalar>(0, 1), 
			Scalar parentPruning = -std::numeric_limits<Scalar>::infinity(),
			std::optional<Model::Seconds> laneLength = std::nullopt // by default, the lane is "budget" long
		)
			: data(data)
			, layout(data.buildLayout(parentPruning))
			, bounds(bounds.toFloat())
			, window(zoom)
			, depthLevelsRequired(layout.maxDepthLevelsInLayout() + 1)
			, yPixelsForHeight(std::min(20, (bounds.getHeight() - space * depthLevelsRequired) / depthLevelsRequired))
			, length(laneLength.value_or(data.budget()))
		{

		}

		void paint(
			juce::Graphics& g, 
			juce::Colour spanFill, 
			juce::Colour outline, 
			std::optional<juce::Colour> text = std::nullopt,
			std::optional<juce::Colour> hot = std::nullopt)
		{
			juce::Font old = g.getCurrentFont();
			juce::Font newFont;

			newFont.setTypefaceName(juce::Font::getDefaultMonospacedFontName());
			g.setFont(newFont);

			// TODO: set up clip regions <-- still a problem
			this->g = &g;
			spanColour = spanFill;
			textColour = text.value_or(spanColour.contrasting());
			hotColour = hot.value_or(spanColour);
			outlineColour = outline;

			// Render the name of the lane as a root node
			renderNode(-1, data.getName().c_str(), Model::Seconds(0), data.duration(), std::nullopt);
			// Then all the children.
			data.visit(*this, layout);

			char buffer[1024];

			cpl::sprintfs(
				buffer,
				"FPS: %6.2f (budget: %6.2f ms)\nDelta Time: %5.2f ms",
				1.0 / data.budget().count(),
				data.budget().count() * 1000,
				data.deltaTime().count() * 1000
			);

			auto topLeft = bounds.getTopLeft();

			g.setColour(textColour);
			g.drawMultiLineText(buffer, cpl::Math::round<int>(topLeft.x), cpl::Math::round<int>(topLeft.y + 20), 400);

			// lane outlines
			float paths[2] = { 5, 5 };

			auto x = secondsToX(Model::Seconds(0));
			g.drawDashedLine(
				{ x, bounds.getY() + bounds.getHeight() * 0.5f, x, bounds.getBottom() },
				paths,
				2, // elements in paths
				2 // line thickness
			);

			x = secondsToX(data.budget());
			g.drawDashedLine(
				{ x, bounds.getY() + bounds.getHeight() * 0.5f, x, bounds.getBottom()},
				paths,
				2, // elements in paths
				2 // line thickness
			);

			g.setFont(old);
		}

	private:

		juce::Colour spanColour;
		juce::Colour outlineColour;
		juce::Colour textColour;
		juce::Colour hotColour;

		juce::Graphics* g = nullptr;

		const Model::LaneData data;
		const Model::LaneData::Layout layout;
		const juce::Rectangle<Scalar> bounds;
		const juce::Point<Scalar> window;
		const int depthLevelsRequired;
		const int yPixelsForHeight;
		const Model::Seconds length;

		Scalar secondsToX(Model::Seconds seconds) const
		{
			const auto normalized = seconds / length;
			const auto zoomed = (normalized - window.getX()) / (window.getY() - window.getX());
			return bounds.getX() + zoomed * bounds.getWidth();
		}

		void operator() (int depth, cpl::Profiling::Region::Identifier identifier, Model::Seconds start, Model::Seconds self, Model::Seconds total) const
		{
			renderNode(
				depth,
				cpl::Profiling::resolveRegion(identifier).name,
				start,
				total,
				self
			);
		}

		void renderNode(int depth, const char* name, Model::Seconds start, Model::Seconds total, std::optional<Model::Seconds> self) const
		{
			auto left = secondsToX(start);
			auto right = secondsToX(start + total);
			
			// depth + 1 since root is painted as well below
			auto bottom = bounds.getBottom() - (depth + 1) * (yPixelsForHeight + space);
			auto top = bottom - yPixelsForHeight;

			auto rect = juce::Rectangle<float>(
				left, // x
				top, // y
				right - left, // width
				static_cast<float>(yPixelsForHeight - space) // height
			);

			// Have the text clip to borders to it doesn't disappear
			auto visibleRect = rect;

			// early out if window completely culls
			if (!bounds.intersectRectangle(visibleRect))
				return;

			// no 'self' means we measure against the whole frame.
			// really only here to reuse this function drawing for the root node.
			const auto target = self ? data.duration() : data.budget();
			auto selfProportion = self ? *self / total : total / target;

			const auto finalSpanColour = spanColour.interpolatedWith(hotColour, selfProportion);

			g->setColour(finalSpanColour);
			g->fillRect(visibleRect);

			g->setColour(outlineColour);
			g->drawRect(rect);

			char buffer[2048];

			cpl::sprintfs(
				buffer,
				"%s %5.2f%%\t (%5.2f ms)",
				name,
				(total / target) * 100,
				total.count() * 1000
			);

			g->setColour(textColour);
			g->drawFittedText(buffer, visibleRect.toNearestInt(), juce::Justification::centred, 1, 1.0f);
		}
	};

	class ProfilerContent
		: public juce::Component
		, private juce::Timer
		, private cpl::ValueEntityBase::ValueEntityListener
	{
		static constexpr int kTimerFrequency = 30;
		static constexpr int kBorder = 5;
		static constexpr int kControlPaneHeight = 40;

	public:

		enum class TimeAxisModes
		{
			IndependentBudget,
			IndependentDuration,
			AlignedBudget,
			AlignedDuration
		};

		ProfilerContent(std::shared_ptr<const SharedBehaviour> behaviour)
			: behaviour(std::move(behaviour))
			, pruneRange(0.0001, 0.99)
			, pruneValue(&pruneRange, &pruneFormatter)
			, pruneControl(&pruneValue)
			, timeChoices(timeRange)
			, timeAxisValue(&timeRange, &timeChoices)
		{
			// Each lane pools 8 snapshots and producers drop (never allocate) when full,
			// so a drain rate below the fastest producer only decimates - it doesn't break anything.
			startTimerHz(kTimerFrequency);
			timeChoices.setValues({ "Budget; independent", "Duration; independent", "Budget; aligned", "Duration; aligned" });

			timeAxisControl = std::make_unique<cpl::CValueComboBox>(&timeAxisValue);

			pruneControl.bSetTitle("Prune parents");
			pruneControl.bSetDescription("Avoid showing parent profiler sections whos self-time is less than this");

			timeAxisControl->bSetTitle("Axis scaling");
			timeAxisControl->bSetDescription("Select how the lanes' time axes are scaled");

			pruneValue.setTransformedValue(0.01);
			timeAxisValue.setAsTEnum(TimeAxisModes::IndependentDuration);

			pruneControl.bSetPos(kBorder, kBorder);
			timeAxisControl->bSetPos(pruneControl.getRight() + kBorder, kBorder);

			pruneValue.addListener(this);
			timeAxisValue.addListener(this);

			addAndMakeVisible(*timeAxisControl);
			addAndMakeVisible(&pruneControl);
		}

	private:

		void valueEntityChanged(ValueEntityListener* sender, cpl::ValueEntityBase* value) override
		{
			repaint();
		}

		void timerCallback() override
		{
			// Single-drainer invariant: the lane queues are SPSC, and this timer callback
			// is the only consumer of all three lanes - the model is created, updated and
			// read exclusively on the message thread. Draining a lane from anywhere else,
			// or reading the model off-thread, breaks this.
			model.consume(*behaviour->getRenderingLane());
			model.consume(behaviour->getRealtimeLane());
			model.consume(behaviour->getAsyncDSPLane());

			repaint();
		}

		void paint(juce::Graphics& g) override
		{	
			const auto backgroundColour = cpl::GetColour(cpl::ColourEntry::Normal);
			const auto outlineColour = cpl::GetColour(cpl::ColourEntry::Separator);
			const auto foregroundColour = cpl::GetColour(cpl::ColourEntry::Auxillary);
			const auto textColour = cpl::GetColour(cpl::ColourEntry::ControlText);
			const auto hotColour = foregroundColour.interpolatedWith(cpl::GetColour(cpl::ColourEntry::Error), 0.25);

			g.fillAll(backgroundColour);

			auto localBounds = getLocalBounds()
				.withTrimmedTop(kControlPaneHeight + kBorder)
				.expanded(-kBorder, -kBorder);

			auto lanes = { 
				model.getLaneData(*behaviour->getRenderingLane().get()), 
				model.getLaneData(behaviour->getRealtimeLane()),
				model.getLaneData(behaviour->getAsyncDSPLane())
			};

			// do better than all must be alive in future.
			if (std::any_of(lanes.begin(), lanes.end(), [](auto& model) { return !model; }))
				return;

			auto totalDepthsNeeded = std::accumulate(
				lanes.begin(), 
				lanes.end(), 
				0, 
				[](int acc, const auto& l) 
				{
					return acc + l->maxDepthSeen() + 5; // +1 for levels, +1 to ensure one slot inbetween all lanes
				}
			);

			auto maxBudget = std::max_element(lanes.begin(), lanes.end(), [](const auto& a, const auto& b) { return a->budget() < b->budget(); });
			auto maxDuration = std::max_element(lanes.begin(), lanes.end(), [](const auto& a, const auto& b) { return a->duration() < b->duration(); });

			auto scalingMode = timeAxisValue.getAsTEnum<TimeAxisModes>();

			auto currentBottom = localBounds.getY();

			for (const auto& data : lanes)
			{
				std::optional<cpl::Profiling::EWMAModel::Seconds> laneLength;

				switch (scalingMode)
				{
					//case IndependentBudget: // default
					case TimeAxisModes::IndependentDuration: laneLength = data->duration(); break;
					case TimeAxisModes::AlignedBudget: laneLength = (*maxBudget)->budget(); break;
					case TimeAxisModes::AlignedDuration: laneLength = (*maxDuration)->duration(); break;
				}

				auto proportionalSpaceNeeded = (data->maxDepthSeen() + 5.f) / totalDepthsNeeded;

				auto rect = localBounds
					.toFloat()
					.withHeight(proportionalSpaceNeeded * localBounds.getHeight())
					.withY(static_cast<float>(currentBottom));

				EWMALaneJuceRenderer renderer(
					*data,  
					rect.toNearestInt(),
					viewOffsets.toFloat(), 
					static_cast<float>(pruneValue.getTransformedValue()),
					laneLength
				);

				renderer.paint(g, foregroundColour, outlineColour, std::nullopt, hotColour);

				currentBottom = rect.getBottom();
			}
		}

		// zooms view offsets
		void mouseWheelMove(const juce::MouseEvent& event, const juce::MouseWheelDetails& wheel) override
		{
			constexpr double increment = 1.2;
			auto position = event.position.getX() / (getWidth() - 1.0);

			double sign = wheel.isReversed ? 1 : -1;

			auto span = (viewOffsets.getY() - viewOffsets.getX());
			auto scale = std::pow(increment, sign * wheel.deltaY);
			auto delta = span - span * scale;

			viewOffsets.setX(viewOffsets.getX() + delta * position);
			viewOffsets.setY(viewOffsets.getY() - delta * (1 - position));

			// Increase time constant as the window gets smaller, but moderated.
			model.setTimeConstant(cpl::Profiling::EWMAModel::Seconds(std::sqrt(1 / span)));

			repaint();
		}

		// resets view offsets on left click, freezes on right click
		void mouseDoubleClick(const juce::MouseEvent& event) override
		{
			if (event.mods.isLeftButtonDown())
			{
				viewOffsets = { 0, 1 };
			}
			else
			{
				// stop timer update
				if (isTimerRunning())
					stopTimer();
				else
					startTimerHz(kTimerFrequency);
			}

			repaint();
		}

		void mouseUp(const juce::MouseEvent& e) override
		{
			if (e.mods.isLeftButtonDown())
				priorDragPosition.reset();
		}

		void mouseDown(const juce::MouseEvent& e) override
		{
			if (e.mods.isLeftButtonDown())
				priorDragPosition = e.position.getX();
		}

		// translates view offsets
		void mouseDrag(const juce::MouseEvent& event) override
		{
			if (!priorDragPosition)
				return;

			auto current = event.position.getX();

			auto delta = current - *priorDragPosition;
			auto fraction = -delta / (getWidth() - 1.0);
			auto span = viewOffsets.getY() - viewOffsets.getX();
			viewOffsets.addXY(fraction * span, fraction * span);

			priorDragPosition = current;

			repaint();
		}

		std::shared_ptr<const SharedBehaviour> behaviour;
		cpl::Profiling::EWMAModel model;
		juce::Point<double> viewOffsets {0, 1};
		std::optional<double> priorDragPosition;

		cpl::ExponentialRange<cpl::ValueT> pruneRange;
		cpl::PercentageFormatter<cpl::ValueT> pruneFormatter;
		cpl::SelfcontainedValue<> pruneValue;
		cpl::CValueKnobSlider pruneControl;

		cpl::ChoiceTransformer<cpl::ValueT> timeRange;
		cpl::ChoiceFormatter<cpl::ValueT> timeChoices;
		cpl::SelfcontainedValue<> timeAxisValue;
		std::unique_ptr<cpl::CValueComboBox> timeAxisControl;
	};

	ProfilerWindow::ProfilerWindow(MainEditor* editor, std::shared_ptr<const SharedBehaviour> behaviour)
		: juce::DocumentWindow("Profiler", juce::Colours::black, juce::DocumentWindow::TitleBarButtons::allButtons)
		, editor(editor)
	{
		setUsingNativeTitleBar(true);
		setResizable(true, false);
		setContentOwned(new ProfilerContent(behaviour), false);
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

