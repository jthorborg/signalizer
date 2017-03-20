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

	file:COscilloscope.cpp

		Implementation of UI, logic and dsp for the oscilloscope.

*************************************************************************************/


#include "COscilloscope.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include "SignalizerDesign.h"
#include <cpl/rendering/OpenGLRasterizers.h>
#include <cpl/simd.h>
#include <cpl/LexicalConversion.h>
#include "VectorScopeParameters.h"
#include "COscilloscopeDSP.inl"

namespace Signalizer
{
	static std::vector<std::string> OperationalModeNames = {"Lissajous", "Polar"};
	static std::vector<std::string> EnvelopeModeNames = {"None", "RMS", "Peak Decay"};

	const double COscilloscope::lowerAutoGainBounds = cpl::Math::dbToFraction(-120.0);
	const double COscilloscope::higherAutoGainBounds = cpl::Math::dbToFraction(120.0);


	COscilloscope::COscilloscope(const std::string & nameId, AudioStream & data, ProcessorState * params)
	:
		COpenGLView(nameId),
		audioStream(data),
		processorSpeed(0),
		lastFrameTick(0),
		lastMousePos(),
		editor(nullptr),
		state(),
		filters(),
		triggerState(),
		medianPos()
	{
		if (!(content = dynamic_cast<OscilloscopeContent *>(params)))
		{
			CPL_RUNTIME_EXCEPTION("Cannot cast parameter set's user data to OscilloscopeContent");
		}

		transformBuffer.resize(OscilloscopeContent::LookaheadSize);
		temporaryBuffer.resize(OscilloscopeContent::LookaheadSize);

		mtFlags.firstRun = true;
		setOpaque(true);
		textbuf = std::unique_ptr<char>(new char[400]);
		processorSpeed = cpl::SysStats::CProcessorInfo::instance().getMHz();
		initPanelAndControls();
		listenToSource(audioStream);
	}

	void COscilloscope::suspend()
	{
	}

	void COscilloscope::resume()
	{
	}

	juce::Component * COscilloscope::getWindow()
	{
		return this;
	}

	COscilloscope::~COscilloscope()
	{
		detachFromSource();
		notifyDestruction();
	}

	void COscilloscope::initPanelAndControls()
	{
		// design
		setMouseCursor(juce::MouseCursor::DraggingHandCursor);
	}

	void COscilloscope::freeze()
	{
		state.isFrozen = true;
	}

	void COscilloscope::unfreeze()
	{
		state.isFrozen = false;
	}

	void COscilloscope::mouseWheelMove(const juce::MouseEvent& event, const juce::MouseWheelDetails& wheel)
	{
		using V = OscilloscopeContent::ViewOffsets;

		double yp;

		if (!state.overlayChannels && state.channelMode > OscChannels::OffsetForMono)
		{
			auto halfHeight = 0.5 * (getHeight() - 1);
			yp = event.position.y;
			if (yp > halfHeight)
				yp -= halfHeight;

			yp = yp / halfHeight;
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
				content->inputGain.setNormalizedValue(content->inputGain.getNormalizedValue() + amount / 20);
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
	void COscilloscope::mouseDoubleClick(const juce::MouseEvent& event)
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
	void COscilloscope::mouseDrag(const juce::MouseEvent& event)
	{
		auto & matrix = content->transform;
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

		auto verticalFactor = !state.overlayChannels && state.channelMode > OscChannels::OffsetForMono ? 2 : 1;

		addClamped(V::Left, xp * (left - right));
		addClamped(V::Right, xp * (left - right));
		addClamped(V::Top, verticalFactor * yp * (top - bottom));
		addClamped(V::Bottom, verticalFactor * yp * (top - bottom));

		lastMousePos = event.position;
	}
	void COscilloscope::mouseUp(const juce::MouseEvent& event)
	{
		// TODO: implement beginChangeGesture()
	}
	void COscilloscope::mouseDown(const juce::MouseEvent& event)
	{
		// TODO: implement endChangeGesture()
		lastMousePos = event.position;
	}

	void COscilloscope::handleFlagUpdates()
	{
		const auto windowValue = content->windowSize.getTransformedValue();

		state.envelopeMode = cpl::enum_cast<EnvelopeModes>(content->autoGain.param.getTransformedValue());
		state.sampleInterpolation = cpl::enum_cast<SubSampleInterpolation>(content->subSampleInterpolation.param.getTransformedValue());
		state.manualGain = content->inputGain.getTransformedValue();
		state.autoGain = shared.autoGainEnvelope.load(std::memory_order_relaxed);
		state.antialias = content->antialias.getTransformedValue() > 0.5;
		state.diagnostics = content->diagnostics.getTransformedValue() > 0.5;
		state.primitiveSize = content->primitiveSize.getTransformedValue();
		state.triggerMode = cpl::enum_cast<OscilloscopeContent::TriggeringMode>(content->triggerMode.param.getTransformedValue());
		state.customTrigger = content->triggerOnCustomFrequency.getNormalizedValue() > 0.5;
		state.customTriggerFrequency = content->customTriggerFrequency.getTransformedValue();
		state.colourChannelsByFrequency = content->channelColouring.param.getAsTEnum<OscilloscopeContent::ColourMode>() == OscilloscopeContent::ColourMode::SpectralEnergy;
		state.overlayChannels = content->overlayChannels.getTransformedValue() > 0.5;
		state.drawCursorTracker = content->cursorTracker.parameter.getValue() > 0.5;

		state.colourPrimary = content->primaryColour.getAsJuceColour();
		state.colourSecondary = content->secondaryColour.getAsJuceColour();
		state.colourBackground = content->backgroundColour.getAsJuceColour();
		state.colourGraph = content->graphColour.getAsJuceColour();
		state.colourTracker = content->trackerColour.getAsJuceColour();

		state.timeMode = cpl::enum_cast<OscilloscopeContent::TimeMode>(content->timeMode.param.getTransformedValue());
		state.beatDivision = windowValue;
		state.dotSamples = content->dotSamples.getNormalizedValue() > 0.5;
		state.channelMode = cpl::enum_cast<OscChannels>(content->channelConfiguration.param.getTransformedValue());

		cpl::foreach_enum<VO>([this](auto i) {
			state.viewOffsets[i] = content->viewOffsets[i].getTransformedValue();
		});

		for (std::size_t c = 0; c < channelData.channels.size(); ++c)
		{
			if (c == 0)
				channelData.channels[c].defaultKey = state.colourPrimary;
			else
				channelData.channels[c].defaultKey = state.colourSecondary;
		}

		switch (state.timeMode)
		{
		case OscilloscopeContent::TimeMode::Beats:
			state.effectiveWindowSize = audioStream.getAudioHistorySamplerate() * (60 / (std::max(10.0, audioStream.getASyncPlayhead().getBPM()) * state.beatDivision));
			state.effectiveWindowSize = std::max(state.effectiveWindowSize, 128.0);
			break;
		case OscilloscopeContent::TimeMode::Cycles:
			state.effectiveWindowSize = windowValue * triggerState.cycleSamples + 1;
			break;
		case OscilloscopeContent::TimeMode::Time:
		default:
			state.effectiveWindowSize = windowValue;
			break;
		}



		bool firstRun = false;

		if (mtFlags.firstRun.cas())
		{
			firstRun = true;
		}


		/* if (firstRun || mtFlags.initiateWindowResize)
		{

			// we will get notified asynchronously in onAsyncChangedProperties.
			if (audioStream.getAudioHistoryCapacity() && audioStream.getAudioHistorySamplerate())
			{
				// only reset this flag if there's valid data, otherwise keep checking.
				mtFlags.initiateWindowResize.cas();
				auto value = content->windowSize.getTransformedValue();
				//audioStream.setAudioHistorySize(value);
			}
		} */
	}

	inline bool COscilloscope::onAsyncAudio(const AudioStream & source, AudioStream::DataType ** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		switch (cpl::simd::max_vector_capacity<float>())
		{
		case 32:
		case 16:
		case 8:
#ifdef CPL_COMPILER_SUPPORTS_AVX
			audioProcessing<cpl::Types::v8sf>(buffer, numChannels, numSamples);
			break;
#endif
		case 4:
			audioProcessing<cpl::Types::v4sf>(buffer, numChannels, numSamples);
			break;
		default:
			audioProcessing<float>(buffer, numChannels, numSamples);
			break;
		}
		return false;
	}

};
