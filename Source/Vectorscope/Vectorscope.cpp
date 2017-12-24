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

	file:VectorScope.cpp

		Implementation of UI, logic and dsp for the vector scope.

*************************************************************************************/


#include "Vectorscope.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include <cpl/rendering/OpenGLRasterizers.h>
#include <cpl/simd.h>
#include <cpl/LexicalConversion.h>
#include "VectorscopeParameters.h"

namespace Signalizer
{
	enum class OperationalModes : int
	{
		Lissajous,
		Polar
	};

	const double VectorScope::lowerAutoGainBounds = cpl::Math::dbToFraction(-120.0);
	const double VectorScope::higherAutoGainBounds = cpl::Math::dbToFraction(120.0);


	VectorScope::VectorScope(const SharedBehaviour & globalBehaviour, const std::string & nameId, AudioStream & data, ProcessorState * params)
		: COpenGLView(nameId)
		, globalBehaviour(globalBehaviour)
		, audioStream(data)
		, processorSpeed(0)
		, lastFrameTick(0)
		, lastMousePos()
		, editor(nullptr)
		, state()
		, filters()
		, oldWindowSize(-1)
	{
		if (!(content = dynamic_cast<VectorScopeContent *>(params)))
		{
			CPL_RUNTIME_EXCEPTION("Cannot cast parameter set's user data to VectorScopeContent");
		}


		mtFlags.firstRun = true;
		state.secondStereoFilterSpeed = 0.25f;
		state.envelopeGain = 1;
		setOpaque(true);
		textbuf = std::unique_ptr<char>(new char[300]);
		processorSpeed = cpl::system::CProcessor::getMHz();
		initPanelAndControls();
		listenToSource(audioStream);
		content->getParameterSet().addRTListener(this, true);
	}

	void VectorScope::suspend()
	{
		//TODO: possibly unsynchronized. fix to have an internal size instead
		oldWindowSize = content->windowSize.getTransformedValue();
		state.isSuspended = true;
	}

	void VectorScope::resume()
	{
		if(oldWindowSize != -1)
		{
			//TODO: possibly unsynchronized. fix to have an internal size instead
			content->windowSize.setTransformedValue(oldWindowSize);
		}
		state.isSuspended = false;
	}

	juce::Component * VectorScope::getWindow()
	{
		return this;
	}

	VectorScope::~VectorScope()
	{
		detachFromSource();
		content->getParameterSet().removeRTListener(this, true);
		notifyDestruction();
	}

	void VectorScope::initPanelAndControls()
	{
		// design
		setMouseCursor(juce::MouseCursor::DraggingHandCursor);
	}

	void VectorScope::freeze()
	{
		state.isFrozen = true;
	}

	void VectorScope::unfreeze()
	{
		state.isFrozen = false;
	}


	void VectorScope::mouseWheelMove(const juce::MouseEvent& event, const juce::MouseWheelDetails& wheel)
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

	void VectorScope::mouseDoubleClick(const juce::MouseEvent& event)
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

	void VectorScope::mouseDrag(const juce::MouseEvent& event)
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

	void VectorScope::mouseUp(const juce::MouseEvent& event)
	{
		// TODO: implement beginChangeGesture()
	}

	void VectorScope::mouseDown(const juce::MouseEvent& event)
	{
		// TODO: implement endChangeGesture()
		lastMousePos = event.position;
	}

	void VectorScope::parameterChangedRT(cpl::Parameters::Handle localHandle, cpl::Parameters::Handle globalHandle, ParameterSet::BaseParameter * param)
	{
		if (param == &content->windowSize.parameter)
		{
			mtFlags.initiateWindowResize = true;
		}
	}

	void VectorScope::handleFlagUpdates()
	{
		state.envelopeMode = cpl::enum_cast<EnvelopeModes>(content->autoGain.param.getTransformedValue());
		state.normalizeGain = state.envelopeMode != EnvelopeModes::None;
		state.envelopeCoeff = std::exp(-1.0 / (content->envelopeWindow.getNormalizedValue() * audioStream.getInfo().sampleRate));
		state.stereoCoeff = std::exp(-1.0 / (content->stereoWindow.getNormalizedValue() * audioStream.getInfo().sampleRate));

		state.isPolar = cpl::enum_cast<OperationalModes>(content->operationalMode.param.getTransformedValue()) == OperationalModes::Polar;
		state.antialias = content->antialias.getTransformedValue() > 0.5;
		state.fadeHistory = content->fadeOlderPoints.getTransformedValue() > 0.5;
		state.fillPath = content->interconnectSamples.getTransformedValue() > 0.5;
		state.diagnostics = content->diagnostics.getTransformedValue() > 0.5;
		state.rotation = content->waveZRotation.getNormalizedValue();
		state.primitiveSize = content->primitiveSize.getTransformedValue();
		state.userGain = content->inputGain.getTransformedValue();

		state.colourDraw = content->drawingColour.getAsJuceColour();
		state.colourWire = content->skeletonColour.getAsJuceColour();
		state.colourBackground = content->backgroundColour.getAsJuceColour();
		state.colourMeter = content->meterColour.getAsJuceColour();
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


	template<typename ISA>
		void VectorScope::audioProcessing(AudioStream::DataType ** buffer, std::size_t numChannels, std::size_t numSamples)
		{
			typedef typename ISA::V V;
			using namespace cpl::simd;
			typedef typename scalar_of<V>::type T;
			if (numChannels != 2)
				return;

			T filterEnv[2] = { filters.envelope[0], filters.envelope[1] };
			T stereoPoles[2] = { state.stereoCoeff, std::pow(state.stereoCoeff, state.secondStereoFilterSpeed) };

			const T cosineRotation = (T)std::cos(M_PI * 135 / 180);
			const T sineRotation = (T)std::sin(M_PI * 135 / 180);
			const V vMatrixReal = set1<V>(cosineRotation);
			const V vMatrixImag = set1<V>(sineRotation);

			const V vDummyAngle = set1<V>((T)(M_PI * 0.25));
			const V vZero = zero<V>();

			auto const loopIncrement = elements_of<V>::value;

			// ensure a perfect multiple and no buffer overrun
			numSamples -= numSamples & (loopIncrement - 1);

			suitable_container<V> outputPhases;

			for (std::size_t n = 0; n < numSamples; n += loopIncrement)
			{
				// some of the heavy math done in vector mode..
				const V vLeft = loadu<V>(buffer[0] + n);
				const V vRight = loadu<V>(buffer[1] + n);

				// rotate 235 degrees...
				const V vX = vLeft * vMatrixReal - vRight * vMatrixImag;
				const V vY = vRight * vMatrixImag + vLeft * vMatrixReal;

				// check if both axes are zero
				const V vZeroAxes = vand((V)(vX == vZero), (V)(vY == vZero));

				// compute the phase angle and replace the zero vector elements with the dummy angle (to avoid nans, +/-infs are defined)
				const V vRadians = atan(vY / vX);
				const V vPhaseAngle = vselect(vDummyAngle, vRadians, vZeroAxes);

				outputPhases = vPhaseAngle;
				// the phase angle is discontinuous around PI, so we take the cosine
				// to avoid the discontinuity and give a small smoothing
				// for some arcane fucking reason, the phase response of this function IN THIS CONTEXT is wrong
				// testing reveals no problems, it just doesn't work right here. Uncomment if you find a fix.
#pragma message cwarn("Errornous cosine output here.")
				//outputPhases = cos(outputPhases * set1<V>(2));

				// scalar code segment for recursive IIR filters..
				for (std::size_t i = 0; i < loopIncrement; ++i)
				{
					auto const z = n + i;
					// collect squared inputs (really just a cheap abs)
					const auto lSquared = buffer[0][z] * buffer[0][z];
					const auto rSquared = buffer[1][z] * buffer[1][z];

					// average envelope
					filterEnv[0] = lSquared + state.envelopeCoeff * (filterEnv[0] - lSquared);
					filterEnv[1] = rSquared + state.envelopeCoeff * (filterEnv[1] - rSquared);

					// balance average source

					using fs = FilterStates;

					filters.balance[fs::Slow][fs::Left]  = lSquared + stereoPoles[fs::Slow] * (filters.balance[fs::Slow][fs::Left]  - lSquared);
					filters.balance[fs::Slow][fs::Right] = rSquared + stereoPoles[fs::Slow] * (filters.balance[fs::Slow][fs::Right] - rSquared);
					filters.balance[fs::Fast][fs::Left]  = lSquared + stereoPoles[fs::Fast] * (filters.balance[fs::Fast][fs::Left]  - lSquared);
					filters.balance[fs::Fast][fs::Right] = rSquared + stereoPoles[fs::Fast] * (filters.balance[fs::Fast][fs::Right] - rSquared);

					// phase averaging
					// see larger comment above.
					outputPhases[i] = cos(outputPhases[i] * consts<T>::two);
					filters.phase[0] = outputPhases[i] + stereoPoles[0] * (filters.phase[0] - outputPhases[i]);
					filters.phase[1] = outputPhases[i] + stereoPoles[1] * (filters.phase[1] - outputPhases[i]);
				}


			}
			// store calculated envelope
			if (state.envelopeMode == EnvelopeModes::RMS && state.normalizeGain)
			{
				// we end up calculating envelopes even though its not possibly needed, but the overhead
				// is negligible
				double currentEnvelope = 1.0 / (std::max(std::sqrt(filterEnv[0]), std::sqrt(filterEnv[1])));

				// only update filters if this mode is on.
				filters.envelope[0] = filterEnv[0];
				filters.envelope[1] = filterEnv[1];

				if (std::isnormal(currentEnvelope))
				{
					state.envelopeGain = currentEnvelope;
				}
			}

		}

	bool VectorScope::onAsyncAudio(const AudioStream & source, AudioStream::DataType ** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		if (state.isSuspended && globalBehaviour.stopProcessingOnSuspend.load(std::memory_order_relaxed))
			return false;

		cpl::simd::dynamic_isa_dispatch<AFloat, AudioDispatcher>(*this, buffer, numChannels, numSamples);
		return false;
	}

	void VectorScope::onAsyncChangedProperties(const AudioStream & source, const AudioStream::AudioStreamInfo & before)
	{
		mtFlags.audioWindowWasResized = true;
	}

};
