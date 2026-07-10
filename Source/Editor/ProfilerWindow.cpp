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

#include <cpl/profiling/ProfilingModel.h>

#include "ProfilerWindow.h"
#include "MainEditor.h"
#include "../Common/SharedBehaviour.h"

namespace Signalizer
{
	struct EWMALaneJuceRenderer
	{
		friend class cpl::Profiling::EWMAModel;
		using Model = cpl::Profiling::EWMAModel;
		using Scalar = float;

		constexpr static std::size_t space = 1;

	public:

		EWMALaneJuceRenderer(Model::LaneData data, juce::Rectangle<int> bounds, juce::Point<Scalar> zoom = juce::Point<Scalar>(0, 1))
			: data(data)
			, bounds(bounds.toFloat())
			, window(zoom)
			, numDepths(std::max(1, data.maxDepthSeen() + 1))
			, yPixelsForHeight((bounds.getHeight() - space * numDepths) / numDepths)
		{

		}

		void paint(juce::Graphics& g, juce::Colour spanFill, juce::Colour outline, std::optional<juce::Colour> text = std::nullopt)
		{
			// TODO: set up clip regions
			this->g = &g;
			spanColour = spanFill;
			textColour = text.value_or(spanColour.contrasting());
			outlineColour = outline;
			data.visit(*this);

			char buffer[2048];

			cpl::sprintfs(
				buffer,
				"%s usage: %2.2f%%\nDuration: %2.2f ms\nBudget: %2.2f ms\nDelta Time: %2.2f ms",
				data.getName().c_str(),
				data.usage() * 100,
				data.duration().count() * 1000,
				data.budget().count() * 1000,
				data.deltaTime().count() * 1000
			);

			auto topLeft = bounds.getTopLeft();

			g.setColour(textColour);
			g.drawMultiLineText(buffer, topLeft.x, topLeft.y + 20, 200);
		}

	private:

		juce::Colour spanColour;
		juce::Colour outlineColour;
		juce::Colour textColour;

		juce::Graphics* g = nullptr;

		const Model::LaneData data;
		const juce::Rectangle<Scalar> bounds;
		const juce::Point<Scalar> window;
		const int numDepths;
		const int yPixelsForHeight;

		Scalar secondsToX(Model::Seconds seconds) const
		{
			const auto normalized = seconds / data.budget();
			const auto zoomed = (normalized - window.getX()) / (window.getY() - window.getX());
			return bounds.getX() + zoomed * bounds.getWidth();
		}

		void operator() (int depth, cpl::Profiling::Region::Identifier identifier, Model::Seconds start, Model::Seconds self, Model::Seconds total) const
		{
			const auto name = cpl::Profiling::resolveRegion(identifier).name;

			auto left = secondsToX(start);
			auto right = secondsToX(start + total);
			
			auto bottom = bounds.getBottom() - depth * yPixelsForHeight;
			auto top = bottom - yPixelsForHeight;

			// TODO: early out if zoom culls
			// TODO: bound rect
			auto rect = juce::Rectangle<float>(
				left, // x
				top, // y
				right - left, // width
				yPixelsForHeight // height
			);

			g->setColour(spanColour);
			g->fillRect(rect);
			g->setColour(outlineColour);
			g->drawRect(rect);
			g->setColour(textColour);
			g->drawFittedText(name, rect.toNearestInt(), juce::Justification::centred, 1, 1.0f);
		}
	};

	class ProfilerContent
		: public juce::Component
		, private juce::Timer
	{
		static constexpr int kTimerFrequency = 30;

	public:

		ProfilerContent(std::shared_ptr<const SharedBehaviour> behaviour)
			: behaviour(std::move(behaviour))
		{
			// Each lane pools 8 snapshots and producers drop (never allocate) when full,
			// so a drain rate below the fastest producer only decimates - it doesn't break anything.
			startTimerHz(kTimerFrequency);
		}

	private:

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
			const auto foregrundColour = cpl::GetColour(cpl::ColourEntry::Auxillary);
			const auto textColour = cpl::GetColour(cpl::ColourEntry::ControlText);

			g.fillAll(backgroundColour);

			auto localBounds = getLocalBounds();
			localBounds.setHeight(localBounds.getHeight() / 3);

			for (const auto lane : { behaviour->getRenderingLane().get(), &behaviour->getRealtimeLane(), &behaviour->getAsyncDSPLane()})
			{
				auto data = model.getLaneData(*lane);

				if (!data)
					continue;

				EWMALaneJuceRenderer renderer(*data, localBounds, viewOffsets.toFloat());
				renderer.paint(g, foregrundColour, outlineColour);

				localBounds.translate(0, localBounds.getHeight());
			}
		}

		// zooms view offsets
		void mouseWheelMove(const juce::MouseEvent& event, const juce::MouseWheelDetails& wheel) override
		{
			constexpr double increment = 1.2;
			auto fraction = event.position.getX() / (getWidth() - 1.0);

			double sign = wheel.isReversed ? -1 : 1;

			auto span = (viewOffsets.getY() - viewOffsets.getX());
			span *= sign * wheel.deltaY * 0.2;

			viewOffsets.setX(viewOffsets.getX() + span * fraction);
			viewOffsets.setY(viewOffsets.getY() - span * (1 - fraction));

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

