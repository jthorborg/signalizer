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

	file:Oscilloscope.cpp

		Implementation of UI, logic and dsp for the oscilloscope.

*************************************************************************************/


#include "Oscilloscope.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include <cpl/rendering/OpenGLRasterizers.h>
#include <cpl/simd.h>
#include <cpl/LexicalConversion.h>
#include "OscilloscopeDSP.inl"
#include "StreamPreprocessing.h"

namespace Signalizer
{
	const double Oscilloscope::lowerAutoGainBounds = cpl::Math::dbToFraction(-120.0);
	const double Oscilloscope::higherAutoGainBounds = cpl::Math::dbToFraction(120.0);

    constexpr std::size_t OscilloscopeContent::LookaheadSize;
    constexpr std::size_t OscilloscopeContent::InterpolationKernelSize;
    

	Oscilloscope::Oscilloscope(
		std::shared_ptr<const SharedBehaviour>& globalBehaviour,
		std::shared_ptr<const ConcurrentConfig>& config,
		const cpl::string_ref nameId,
		std::shared_ptr<AudioStream::Output>& stream,
		std::shared_ptr<ProcessorState>& params
	)
		: COpenGLView(nameId.string())
		, globalBehaviour(globalBehaviour)
		, audioStream(stream)
		, processorSpeed(0)
		, lastFrameTick(0)
		, lastMousePos()
		, editor(nullptr)
		, state()
		, medianPos()
		, isMouseInside(false)
		, processor(std::make_shared<ProcessorShell>())
	{
		if (!(content = std::dynamic_pointer_cast<OscilloscopeContent>(params)))
		{
			CPL_RUNTIME_EXCEPTION("Cannot cast parameter set's user data to OscilloscopeContent");
		}

		processor->streamState.lock()->content = content;

		transformBuffer.resize(OscilloscopeContent::LookaheadSize);
		temporaryBuffer.resize(OscilloscopeContent::LookaheadSize);

		mtFlags.firstRun = true;
		setOpaque(true);
		processorSpeed = cpl::system::CProcessor::getMHz();
		initPanelAndControls();
		stream->addListener(processor);
	}

	void Oscilloscope::suspend()
	{
		processor->isSuspended = true;
	}

	void Oscilloscope::resume()
	{
		processor->isSuspended = false;
	}

	juce::Component * Oscilloscope::getWindow()
	{
		return this;
	}

	Oscilloscope::~Oscilloscope()
	{
		audioStream->removeListener(processor);
		notifyDestruction();
	}

	std::size_t Oscilloscope::getEffectiveChannels() const noexcept
	{
		if (state.channelMode >= OscChannels::OffsetForMono)
		{
			return state.channelMode == OscChannels::Separate ? state.numChannels : 2;
		}

		return 1;
	}

	void Oscilloscope::initPanelAndControls()
	{
		// design
		setMouseCursor(juce::MouseCursor::DraggingHandCursor);
	}

	void Oscilloscope::freeze()
	{
		state.isFrozen = true;
	}

	void Oscilloscope::unfreeze()
	{
		state.isFrozen = false;
	}

	void Oscilloscope::setLastMousePos(const juce::Point<float> position) noexcept
	{
		lastMousePos = position;
		threadedMousePos.first = position.x;
		threadedMousePos.second = position.y;
	}

	void Oscilloscope::mouseWheelMove(const juce::MouseEvent& event, const juce::MouseWheelDetails& wheel)
	{
		using V = OscilloscopeContent::ViewOffsets;

		double yp;

		if (!state.overlayChannels && state.channelMode > OscChannels::OffsetForMono)
		{
			auto heightPerScope = (getHeight() - 1.0) / getEffectiveChannels();
			yp = event.position.y;
			while (yp > heightPerScope)
				yp -= heightPerScope;

			yp = yp / heightPerScope;
		}
		else
		{
			yp = double(event.position.y) / (getHeight() - 1);
		}

		auto xp = double(event.position.x) / (getWidth() - 1);

		auto get = [&](auto i) { return content->viewOffsets[i].getTransformedValue(); };

		auto amount = wheel.deltaY;
		if (event.mods.isCtrlDown())
		{
			if (event.mods.isShiftDown())
			{
				// scale only vertically
				auto top = get(V::Top);
				auto bottom = get(V::Bottom);

				auto incY = -(top - bottom) * wheel.deltaY / 5;
				// TODO: change to pow()
				content->viewOffsets[V::Top].setTransformedValue(top + yp * incY);
				content->viewOffsets[V::Bottom].setTransformedValue(bottom - (1 - yp) * incY);
			}
			else
			{
				// TODO: fix to pow()
				content->inputGain.setNormalizedValue(content->inputGain.getNormalizedValue() + amount / 80);
			}


		}
		else if (event.mods.isShiftDown())
		{
			// scale only horizontally
			auto left = get(V::Left);
			auto right = get(V::Right);

			auto incX = -(left - right) * wheel.deltaY / 5;
			// TODO: change to pow()
			content->viewOffsets[V::Left].setTransformedValue(left + xp * incX);
			content->viewOffsets[V::Right].setTransformedValue(right - (1 - xp) * incX);
		}
		else
		{
			// zoom both graph axis

			auto left = get(V::Left);
			auto right = get(V::Right);
			auto top = get(V::Top);
			auto bottom = get(V::Bottom);

			auto incX = -(left - right) * wheel.deltaY / 5;
			auto incY = -(top - bottom) * wheel.deltaY / 5;
			// TODO: change to pow()
			content->viewOffsets[V::Left].setTransformedValue(left + xp * incX);
			content->viewOffsets[V::Right].setTransformedValue(right - (1 - xp) * incX);
			content->viewOffsets[V::Top].setTransformedValue(top + yp * incY);
			content->viewOffsets[V::Bottom].setTransformedValue(bottom - (1 - yp) * incY);

		}

	}
	void Oscilloscope::mouseDoubleClick(const juce::MouseEvent& event)
	{
		using V = OscilloscopeContent::ViewOffsets;

		if (event.mods.isLeftButtonDown())
		{
			// reset all zooming, offsets etc. when doubleclicking left
			// TODO: reset to preset
			//content->inputGain.setTransformedValue(1);
			auto & matrix = content->transform;
			matrix.getValueIndex(matrix.Position, matrix.X).setTransformedValue(0);
			matrix.getValueIndex(matrix.Position, matrix.Y).setTransformedValue(0);

			content->viewOffsets[V::Left].setNormalizedValue(0);
			content->viewOffsets[V::Right].setNormalizedValue(0);
			content->viewOffsets[V::Top].setNormalizedValue(0);
			content->viewOffsets[V::Bottom].setNormalizedValue(0);
		}
	}
	void Oscilloscope::mouseDrag(const juce::MouseEvent& event)
	{
		auto deltaDifference = event.position - lastMousePos;

		using V = OscilloscopeContent::ViewOffsets;

		auto yp = double(deltaDifference.y) / (getHeight() - 1);
		auto xp = double(deltaDifference.x) / (getWidth() - 1);

		auto get = [&](auto i) { return content->viewOffsets[i].getTransformedValue(); };

		auto left = get(V::Left);
		auto right = get(V::Right);
		auto top = get(V::Top);
		auto bottom = get(V::Bottom);

		auto addClamped = [&](auto i, double val) {
			content->viewOffsets[i].setTransformedValue(cpl::Math::confineTo(get(i) + val, 0.0, 1.0));
		};

		auto verticalFactor = state.overlayChannels ? 1 : getEffectiveChannels();

		addClamped(V::Left, xp * (left - right));
		addClamped(V::Right, xp * (left - right));
		addClamped(V::Top, verticalFactor * yp * (top - bottom));
		addClamped(V::Bottom, verticalFactor * yp * (top - bottom));

		setLastMousePos(event.position);
	}
	void Oscilloscope::mouseUp(const juce::MouseEvent& event)
	{
		// TODO: implement beginChangeGesture()
	}
	void Oscilloscope::mouseDown(const juce::MouseEvent& event)
	{
		// TODO: implement endChangeGesture()
		setLastMousePos(event.position);
	}

	void Oscilloscope::mouseMove(const juce::MouseEvent & event)
	{
		setLastMousePos(event.position);
	}

	void Oscilloscope::mouseExit(const juce::MouseEvent & e)
	{
		isMouseInside = true;
	}

	void Oscilloscope::mouseEnter(const juce::MouseEvent & e)
	{
		isMouseInside = true;
	}

	void Oscilloscope::handleFlagUpdates(Oscilloscope::StreamState& cs)
	{
		const auto windowValue = content->windowSize.getTransformedValue();

		cs.envelopeMode = cpl::enum_cast<EnvelopeModes>(content->autoGain.param.getTransformedValue());
		state.sampleInterpolation = cpl::enum_cast<SubSampleInterpolation>(content->subSampleInterpolation.param.getTransformedValue());
		state.manualGain = content->inputGain.getTransformedValue();
		state.autoGain = shared.autoGainEnvelope;
		state.antialias = content->antialias.getTransformedValue() > 0.5;
		state.diagnostics = content->diagnostics.getTransformedValue() > 0.5;
		state.primitiveSize = content->primitiveSize.getTransformedValue();
		cs.triggerMode = cpl::enum_cast<OscilloscopeContent::TriggeringMode>(content->triggerMode.param.getTransformedValue());
		state.customTrigger = content->triggerOnCustomFrequency.getNormalizedValue() > 0.5;
		state.customTriggerFrequency = content->customTriggerFrequency.getTransformedValue();
		state.colourChannelsByFrequency = content->channelColouring.param.getAsTEnum<OscilloscopeContent::ColourMode>() == OscilloscopeContent::ColourMode::SpectralEnergy;
		state.overlayChannels = content->overlayChannels.getTransformedValue() > 0.5;
		state.drawCursorTracker = content->cursorTracker.parameter.getValue() > 0.5;

		state.colourBackground = content->backgroundColour.getAsJuceColour();
		state.colourGraph = content->graphColour.getAsJuceColour();
		state.colourTracker = content->trackerColour.getAsJuceColour();

		state.timeMode = cpl::enum_cast<OscilloscopeContent::TimeMode>(content->timeMode.param.getTransformedValue());
		state.beatDivision = windowValue;
		state.dotSamples = content->dotSamples.getNormalizedValue() > 0.5;
		cs.channelMode = cpl::enum_cast<OscChannels>(content->channelConfiguration.param.getTransformedValue());

		state.triggerHysteresis = content->triggerHysteresis.parameter.getValue();
		state.triggerThreshold = content->triggerThreshold.getTransformedValue();
		state.numChannels = cs.channelData.numChannels();

		state.drawLegend = content->showLegend.getTransformedValue() > 0.5;

		cpl::foreach_enum<VO>([this](auto i) {
			state.viewOffsets[i] = content->viewOffsets[i].getTransformedValue();
		});

		for (std::size_t c = 0; c < state.numChannels; ++c)
		{
			state.colours[c] = cs.channelData.filterStates.channels[c].defaultKey = content->getColour(c).getAsJuceColour();
		}

		switch (state.timeMode)
		{
		case OscilloscopeContent::TimeMode::Beats:
			state.effectiveWindowSize = audioStream.getAudioHistorySamplerate() * (60 / (std::max(10.0, audioStream.getASyncPlayhead().getBPM()) * state.beatDivision));
			state.effectiveWindowSize = std::max(state.effectiveWindowSize, 128.0);
			break;
		case OscilloscopeContent::TimeMode::Cycles:
			// +1 to allow another sample's display when centering waveforms in the middle,
			// used by preprocessing triggers
			state.effectiveWindowSize = windowValue * triggerState.cycleSamples + 1;
			break;
		case OscilloscopeContent::TimeMode::Time:
		default:
			state.effectiveWindowSize = windowValue;
			break;
		}

		cs.triggeringProcessor->setSettings(cs.triggerMode, state.effectiveWindowSize, state.triggerThreshold, state.triggerHysteresis);
	}

	inline void Oscilloscope::ProcessorShell::onStreamAudio(AudioStream::ListenerContext& source, AudioStream::DataType** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		if (isSuspended && globalBehaviour->stopProcessingOnSuspend)
			return;

		cpl::simd::dynamic_isa_dispatch<float, AudioDispatcher>(*this, source, buffer, numChannels, numSamples);
	}

	inline void Oscilloscope::ProcessorShell::onStreamPropertiesChanged(AudioStream::ListenerContext& source, const AudioStream::AudioStreamInfo& before)
	{
		auto& access = streamState.lock();
		access->channelNames = source.getChannelNames();
		access->historyCapacity = source.getInfo().audioHistoryCapacity;
	}

	Oscilloscope::StreamState::StreamState()
		: triggeringProcessor(std::make_unique<TriggeringProcessor>())
	{

	}
	Oscilloscope::StreamState::~StreamState()
	{
	}
};
