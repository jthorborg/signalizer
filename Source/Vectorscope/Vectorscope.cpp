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


	VectorScope::VectorScope(
		std::shared_ptr<const SharedBehaviour>& globalBehaviour,
		std::shared_ptr<const ConcurrentConfig>& config,
		std::shared_ptr<AudioStream::Output>& stream, 
		std::shared_ptr<VectorScopeContent>& params
	)
		: GraphicsWindow(params->getName())
		, editor(nullptr)
		, state()
		, processor(std::make_shared<Processor>(globalBehaviour))
		, config(config)
		, content(params)
		, audioStream(stream)
		, globalBehaviour(globalBehaviour)
	{
		mtFlags.firstRun = true;
		setOpaque(true);

		initPanelAndControls();
		stream->addListener(processor);
		content->getParameterSet().addRTListener(this, true);
	}

	void VectorScope::suspend()
	{
		processor->isSuspended = true;
	}

	void VectorScope::resume()
	{
		mtFlags.initiateWindowResize = true;
		processor->isSuspended = false;
	}

	juce::Component * VectorScope::getWindow()
	{
		return this;
	}

	VectorScope::~VectorScope()
	{
		audioStream->removeListener(processor);
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
			auto & tsX = content->transform;
			auto Z = tsX.getValueIndex(tsX.Scale, tsX.Z).getTransformedValue();
			auto actualAmount = (1 + amount / 5) * Z;
			auto deltaIncrease = (actualAmount - Z) / Z;
			tsX.getValueIndex(tsX.Scale, tsX.Z).setTransformedValue(actualAmount);

			tsX.getValueIndex(tsX.Scale, tsX.X).setTransformedValue(
				tsX.getValueIndex(tsX.Scale, tsX.X).getTransformedValue() * (1 + deltaIncrease)
			);

			tsX.getValueIndex(tsX.Scale, tsX.Y).setTransformedValue(
				tsX.getValueIndex(tsX.Scale, tsX.Y).getTransformedValue() * (1 + deltaIncrease)
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
			auto & tsX = content->transform;
			tsX.getValueIndex(tsX.Position, tsX.X).setTransformedValue(0);
			tsX.getValueIndex(tsX.Position, tsX.Y).setTransformedValue(0);
		}
	}

	void VectorScope::mouseDrag(const juce::MouseEvent& event)
	{
		auto & tsX = content->transform;
		auto factor = float(getWidth()) / getHeight();
		auto deltaDifference = event.position - currentMouse.getPoint();
		if (event.mods.isCtrlDown())
		{
			auto yvalue = tsX.getValueIndex(tsX.Rotation, tsX.Y).getTransformedValue();
			auto newValue = std::fmod(deltaDifference.x * 0.3f + yvalue, 360.f);
			while (newValue < 0)
				newValue += 360;
			tsX.getValueIndex(tsX.Rotation, tsX.Y).setTransformedValue(newValue);

			auto xvalue = tsX.getValueIndex(tsX.Rotation, tsX.X).getTransformedValue();
			newValue = std::fmod(deltaDifference.y * 0.3f + xvalue, 360.f);
			while (newValue < 0)
				newValue += 360;
			tsX.getValueIndex(tsX.Rotation, tsX.X).setTransformedValue(newValue);
		}
		else
		{

			auto xvalue = tsX.getValueIndex(tsX.Position, tsX.X).getTransformedValue();
			tsX.getValueIndex(tsX.Position, tsX.X).setTransformedValue(xvalue + deltaDifference.x / 500.f);

			auto yvalue = tsX.getValueIndex(tsX.Position, tsX.Y).getTransformedValue();
			tsX.getValueIndex(tsX.Position, tsX.Y).setTransformedValue(factor * -deltaDifference.y / 500.f + yvalue);
		}

		GraphicsWindow::mouseDrag(event);
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
		processor->envelopeMode = cpl::enum_cast<EnvelopeModes>(content->autoGain.param.getTransformedValue());
		processor->normalizeGain = processor->envelopeMode != EnvelopeModes::None;
		processor->envelopeCoeff = std::exp(-1.0 / (content->envelopeWindow.getNormalizedValue() * config->sampleRate));
		processor->stereoCoeff = std::exp(-1.0 / (content->stereoWindow.getNormalizedValue() * config->sampleRate));

		state.isPolar = cpl::enum_cast<OperationalModes>(content->operationalMode.param.getTransformedValue()) == OperationalModes::Polar;
		state.antialias = content->antialias.getTransformedValue() > 0.5;
		state.fadeHistory = content->fadeOlderPoints.getTransformedValue() > 0.5;
		state.fillPath = content->interconnectSamples.getTransformedValue() > 0.5;
		state.diagnostics = content->diagnostics.getTransformedValue() > 0.5;
		state.rotation = content->waveZRotation.getNormalizedValue();
		state.primitiveSize = content->primitiveSize.getTransformedValue();
		state.userGain = content->inputGain.getTransformedValue();
		state.drawLegend = content->showLegend.getTransformedValue() > 0.5;
		state.scalePolar = content->scalePolarModeToFill.getTransformedValue() > 0.5;

		state.colourWaveform = content->waveformColour.getAsJuceColour();
		state.colourWire = content->wireframeColour.getAsJuceColour();
		state.colourBackground = content->backgroundColour.getAsJuceColour();
		state.colourMeter = content->meterColour.getAsJuceColour();
		state.colourAxis = content->axisColour.getAsJuceColour();
		state.colourWidget = content->widgetColour.getAsJuceColour();


		bool firstRun = false;

		if (mtFlags.firstRun.cas())
		{
			firstRun = true;
		}

		if (firstRun || mtFlags.initiateWindowResize)
		{
			// we will get notified asynchronously in onStreamPropertiesChanged.
			if (config->historyCapacity > 0 && config->sampleRate > 0)
			{
				// only reset this flag if there's valid data, otherwise keep checking.
				mtFlags.initiateWindowResize.cas();
				auto value = content->windowSize.getTransformedValue();

				audioStream->modifyConsumerInfo(
					[&](auto & info)
					{
						info.audioHistorySize = value;
					}
				);
			}
		}
	}


	template<typename ISA>
		void VectorScope::Processor::audioProcessing(AudioStream::DataType ** buffer, std::size_t numChannels, std::size_t numSamples)
		{
			typedef typename ISA::V V;
			using namespace cpl::simd;
			typedef typename scalar_of<V>::type T;
			if (numChannels < 2)
				return;

			T filterEnv[2] { filters.envelope[0], filters.envelope[1] };
			T balance[2][2]{ { filters.balance[0][0], filters.balance[0][1] }, { filters.balance[1][0], filters.balance[1][1] } };
			T phase[2] { filters.phase[0], filters.phase[1] };

			const T stereoPoles[2] = { stereoCoeff, std::pow(stereoCoeff, secondStereoFilterSpeed) };
			const T envelope = envelopeCoeff;
			const V vMatrixReal = consts<V>::sqrt_half_two_minus;
			const V vMatrixImag = consts<V>::sqrt_half_two;

			const V vDummyAngle = consts<V>::pi_quarter;
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

				// the phase angle is discontinuous around PI, so we take the cosine
				// to avoid the discontinuity and give a small smoothing

				outputPhases = cos(vPhaseAngle * consts<V>::two);

				// scalar code segment for recursive IIR filters..
				for (std::size_t i = 0; i < loopIncrement; ++i)
				{
					auto const z = n + i;
					// collect squared inputs (really just a cheap abs)
					const auto lSquared = buffer[0][z] * buffer[0][z];
					const auto rSquared = buffer[1][z] * buffer[1][z];

					// average envelope
					filterEnv[0] = lSquared + envelope * (filterEnv[0] - lSquared);
					filterEnv[1] = rSquared + envelope * (filterEnv[1] - rSquared);

					// balance average source

					using fs = FilterStates;

					balance[fs::Slow][fs::Left]  = lSquared + stereoPoles[fs::Slow] * (balance[fs::Slow][fs::Left]  - lSquared);
					balance[fs::Slow][fs::Right] = rSquared + stereoPoles[fs::Slow] * (balance[fs::Slow][fs::Right] - rSquared);
					balance[fs::Fast][fs::Left]  = lSquared + stereoPoles[fs::Fast] * (balance[fs::Fast][fs::Left]  - lSquared);
					balance[fs::Fast][fs::Right] = rSquared + stereoPoles[fs::Fast] * (balance[fs::Fast][fs::Right] - rSquared);

					// phase averaging
					phase[fs::Slow] = outputPhases[i] + stereoPoles[fs::Slow] * (phase[fs::Slow] - outputPhases[i]);
					phase[fs::Fast] = outputPhases[i] + stereoPoles[fs::Fast] * (phase[fs::Fast] - outputPhases[i]);
				}

			}

			// store calculated envelope
			if (envelopeMode == EnvelopeModes::RMS && normalizeGain)
			{
				// we end up calculating envelopes even though its not possibly needed, but the overhead
				// is negligible
				double currentEnvelope = 1.0 / (std::max(std::sqrt(filterEnv[0]), std::sqrt(filterEnv[1])));

				// only update filters if this mode is on.
				for (std::size_t i = 0; i < 2; ++i)
				{
					filters.envelope[i] = filterEnv[i];
					filters.phase[i] = phase[i];
					for (std::size_t j = 0; j < 2; ++j)
						filters.balance[i][j] = balance[i][j];
				}
				
				if (std::isnormal(currentEnvelope))
				{
					envelopeGain = currentEnvelope;
				}
			}
		}

	void VectorScope::Processor::onStreamAudio(AudioStream::ListenerContext& source, AudioStream::DataType ** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		if (isSuspended && globalBehaviour->stopProcessingOnSuspend)
			return;

		cpl::simd::dynamic_isa_dispatch<AFloat, AudioDispatcher>(*this, buffer, numChannels, numSamples);
	}

	void VectorScope::Processor::onStreamPropertiesChanged(AudioStream::ListenerContext& ctx, const AudioStream::AudioStreamInfo & before)
	{
		audioWindowWasResized = true;
		*channelNames.lock() = ctx.getChannelNames();
	}

};
