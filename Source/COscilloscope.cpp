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
		oldWindowSize(-1),
		detectedFreq(),
		quantizedFreq()
	{
		if (!(content = dynamic_cast<OscilloscopeContent *>(params)))
		{
			CPL_RUNTIME_EXCEPTION("Cannot cast parameter set's user data to OscilloscopeContent");
		}


		mtFlags.firstRun = true;
		setOpaque(true);
		textbuf = std::unique_ptr<char>(new char[300]);
		processorSpeed = cpl::SysStats::CProcessorInfo::instance().getMHz();
		initPanelAndControls();
		listenToSource(audioStream);
		content->getParameterSet().addRTListener(this, true);
	}

	void COscilloscope::suspend()
	{
		//TODO: possibly unsynchronized. fix to have an internal size instead
		oldWindowSize = content->windowSize.getTransformedValue();
	}

	void COscilloscope::resume()
	{
		if(oldWindowSize != -1)
		{
			//TODO: possibly unsynchronized. fix to have an internal size instead
			content->windowSize.setTransformedValue(oldWindowSize);
		}
	}

	juce::Component * COscilloscope::getWindow()
	{
		return this;
	}

	COscilloscope::~COscilloscope()
	{
		detachFromSource();
		content->getParameterSet().removeRTListener(this, true);
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
		auto amount = wheel.deltaY;
		if (event.mods.isCtrlDown())
		{
			// increase gain
			// TODO: fix to pow()
			content->inputGain.setNormalizedValue(content->inputGain.getNormalizedValue() + amount / 20);
		}
		else /* zoom graph */
		{
			auto & matrix = content->transform;
			//matrix.scale.x += amount;
			//matrix.scale.y += amount;
			auto Z = matrix.getValueIndex(matrix.Scale, matrix.Z).getTransformedValue();
			auto actualAmount = (1 + amount / 5) * Z;
			auto deltaIncrease = (actualAmount - Z) / Z;
			matrix.getValueIndex(matrix.Scale, matrix.Z).setTransformedValue(actualAmount);

			matrix.getValueIndex(matrix.Scale, matrix.X).setTransformedValue(
				matrix.getValueIndex(matrix.Scale, matrix.X).getTransformedValue() * (1 + deltaIncrease)
			);

			matrix.getValueIndex(matrix.Scale, matrix.Y).setTransformedValue(
				matrix.getValueIndex(matrix.Scale, matrix.Y).getTransformedValue() * (1 + deltaIncrease)
			);
		}

	}
	void COscilloscope::mouseDoubleClick(const juce::MouseEvent& event)
	{

		if (event.mods.isLeftButtonDown())
		{
			// reset all zooming, offsets etc. when doubleclicking left
			// TODO: reset to preset
			content->inputGain.setTransformedValue(1);
			auto & matrix = content->transform;
			matrix.getValueIndex(matrix.Position, matrix.X).setTransformedValue(0);
			matrix.getValueIndex(matrix.Position, matrix.Y).setTransformedValue(0);
		}
	}
	void COscilloscope::mouseDrag(const juce::MouseEvent& event)
	{
		auto & matrix = content->transform;
		auto factor = float(getWidth()) / getHeight();
		auto deltaDifference = event.position - lastMousePos;
		if (event.mods.isCtrlDown())
		{
			auto yvalue = matrix.getValueIndex(matrix.Rotation, matrix.Y).getTransformedValue();
			auto newValue = std::fmod(deltaDifference.x * 0.3f + yvalue, 360.f);
			while (newValue < 0)
				newValue += 360;
			matrix.getValueIndex(matrix.Rotation, matrix.Y).setTransformedValue(newValue);

			auto xvalue = matrix.getValueIndex(matrix.Rotation, matrix.X).getTransformedValue();
			newValue = std::fmod(deltaDifference.y * 0.3f + xvalue, 360.f);
			while (newValue < 0)
				newValue += 360;
			matrix.getValueIndex(matrix.Rotation, matrix.X).setTransformedValue(newValue);
		}
		else
		{

			auto xvalue = matrix.getValueIndex(matrix.Position, matrix.X).getTransformedValue();
			matrix.getValueIndex(matrix.Position, matrix.X).setTransformedValue(xvalue + deltaDifference.x / 500.f);

			auto yvalue = matrix.getValueIndex(matrix.Position, matrix.Y).getTransformedValue();
			matrix.getValueIndex(matrix.Position, matrix.Y).setTransformedValue(factor * -deltaDifference.y / 500.f + yvalue);
		}

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

	void COscilloscope::parameterChangedRT(cpl::Parameters::Handle localHandle, cpl::Parameters::Handle globalHandle, ParameterSet::BaseParameter * param)
	{
		if (param == &content->windowSize.parameter)
		{
			mtFlags.initiateWindowResize = true;
		}
	}

	void COscilloscope::handleFlagUpdates()
	{
		state.envelopeMode = cpl::enum_cast<EnvelopeModes>(content->autoGain.param.getTransformedValue());
		state.normalizeGain = state.envelopeMode != EnvelopeModes::None;
		state.envelopeCoeff = std::exp(-1.0 / (content->envelopeWindow.getNormalizedValue() * audioStream.getInfo().sampleRate));
		state.sampleInterpolation = cpl::enum_cast<SubSampleInterpolation>(content->subSampleInterpolation.param.getTransformedValue());
		state.envelopeGain = content->inputGain.getTransformedValue();
		state.antialias = content->antialias.getTransformedValue() > 0.5;
		state.diagnostics = content->diagnostics.getTransformedValue() > 0.5;
		state.primitiveSize = content->primitiveSize.getTransformedValue();

		state.colourDraw = content->drawingColour.getAsJuceColour();
		state.colourWire = content->skeletonColour.getAsJuceColour();
		state.colourBackground = content->backgroundColour.getAsJuceColour();
		state.colourGraph = content->graphColour.getAsJuceColour();


		bool firstRun = false;

		if (mtFlags.firstRun.cas())
		{
			firstRun = true;
		}

		if (firstRun || mtFlags.initiateWindowResize)
		{

			// we will get notified asynchronously in onAsyncChangedProperties.
			if (audioStream.getAudioHistoryCapacity() && audioStream.getAudioHistorySamplerate())
			{
				// only reset this flag if there's valid data, otherwise keep checking.
				mtFlags.initiateWindowResize.cas();
				auto value = content->windowSize.getTransformedValue();
				audioStream.setAudioHistorySize(value);
			}
		}
	}


	template<typename V>
		void COscilloscope::audioProcessing(typename cpl::simd::scalar_of<V>::type ** buffer, std::size_t numChannels, std::size_t numSamples)
		{
			using namespace cpl::simd;
			typedef typename scalar_of<V>::type T;
			if (numChannels != 2)
				return;

			T filterEnv[2] = { filters.envelope[0], filters.envelope[1] };

			for (std::size_t n = 0; n < numSamples; n++)
			{
				using fs = FilterStates;
				// collect squared inputs (really just a cheap abs)
				const auto lSquared = buffer[fs::Left][n] * buffer[fs::Left][n];
				const auto rSquared = buffer[fs::Right][n] * buffer[fs::Right][n];

				// average envelope
				filterEnv[fs::Left] = lSquared + state.envelopeCoeff * (filterEnv[fs::Left] - lSquared);
				filterEnv[fs::Right] = rSquared + state.envelopeCoeff * (filterEnv[fs::Right] - rSquared);



			}
			// store calculated envelope
			if (state.envelopeMode == EnvelopeModes::RMS && state.normalizeGain)
			{
				// we end up calculating envelopes even though its not possibly needed, but the overhead
				// is negligible
				double currentEnvelope = 1.0 / (2 * std::max(std::sqrt(filterEnv[0]), std::sqrt(filterEnv[1])));

				// only update filters if this mode is on.
				filters.envelope[0] = filterEnv[0];
				filters.envelope[1] = filterEnv[1];

				if (std::isnormal(currentEnvelope))
				{
					content->inputGain.getParameterView().updateFromProcessorTransformed(
						currentEnvelope,
						cpl::Parameters::UpdateFlags::All & ~cpl::Parameters::UpdateFlags::RealTimeSubSystem
					);
				}
			}

		}

	bool COscilloscope::onAsyncAudio(const AudioStream & source, AudioStream::DataType ** buffer, std::size_t numChannels, std::size_t numSamples)
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

	void COscilloscope::onAsyncChangedProperties(const AudioStream & source, const AudioStream::AudioStreamInfo & before)
	{
		mtFlags.audioWindowWasResized = true;
	}

};
