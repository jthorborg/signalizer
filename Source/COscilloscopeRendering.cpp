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

	file:COscilloscopeRendering.cpp

		Implementation of all rendering code for the oscilloscope.

*************************************************************************************/


#include "COscilloscope.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include "SignalizerDesign.h"
#include <cpl/rendering/OpenGLRasterizers.h>
#include <cpl/simd.h>
#include <cpl/CDBMeterGraph.h>

namespace Signalizer
{
	// swapping the right channel might give an more intuitive view

	static const char * ChannelDescriptions[] = { "+L", "+R", "-L", "-R", "L", "R", "C", "S" };

	static const float quarterPISinCos = 0.707106781186547f;
	static const float circleScaleFactor = 1.1f;

	enum Textures
	{
		LPlus,
		RPlus,
		LMinus,
		RMinus,
		Left,
		Right,
		Center,
		Side
	};

	template<typename V>
	void COscilloscope::paint2DGraphics(juce::Graphics & g)
	{

		auto cStart = cpl::Misc::ClockCounter();

		if (content->diagnostics.getNormalizedValue() > 0.5)
		{
			auto fps = 1.0 / (avgFps.getAverage() / juce::Time::getHighResolutionTicksPerSecond());

			auto totalCycles = renderCycles + cpl::Misc::ClockCounter() - cStart;
			double cpuTime = (double(totalCycles) / (processorSpeed * 1000 * 1000) * 100) * fps;
			g.setColour(juce::Colours::blue);
			sprintf(textbuf.get(), "%dx%d: %.1f fps - %.1f%% cpu, deltaG = %f, deltaO = %f (rt: %.2f%% - %.2f%%), (as: %.2f%% - %.2f%%), qHZ: %.5f - HZ: %.5f - PHASE: %.5f",
				getWidth(), getHeight(), fps, cpuTime, graphicsDeltaTime(), openGLDeltaTime(),
				100 * audioStream.getPerfMeasures().rtUsage.load(std::memory_order_relaxed),
				100 * audioStream.getPerfMeasures().rtOverhead.load(std::memory_order_relaxed),
				100 * audioStream.getPerfMeasures().asyncUsage.load(std::memory_order_relaxed),
				100 * audioStream.getPerfMeasures().asyncOverhead.load(std::memory_order_relaxed),
				(double)triggerState.record.index,
				triggerState.fundamental,
				triggerState.phase);
			g.drawSingleLineText(textbuf.get(), 10, 20);

		}

		drawWireFrame<V>(g, getLocalBounds().toFloat(), state.envelopeGain);
	}

	void COscilloscope::onGraphicsRendering(juce::Graphics & g)
	{

		// do software rendering
		if(!isOpenGL())
		{
			g.fillAll(state.colourBackground.withAlpha(1.0f));
			g.setColour(state.colourBackground.withAlpha(1.0f).contrasting());
			g.drawText("Enable OpenGL in settings to use the vectorscope", getLocalBounds(), juce::Justification::centred);

			// post fps anyway
			auto tickNow = juce::Time::getHighResolutionTicks();
			avgFps.setNext(tickNow - lastFrameTick);
			lastFrameTick = tickNow;

		}


	}

	void COscilloscope::initOpenGL()
	{

	}

	void COscilloscope::closeOpenGL()
	{

	}

	void COscilloscope::onOpenGLRendering()
	{
		switch (cpl::simd::max_vector_capacity<float>())
		{
		case 32:
		case 16:
		case 8:
            #ifdef CPL_COMPILER_SUPPORTS_AVX
                vectorGLRendering<cpl::Types::v8sf>();
                break;
            #endif
		case 4:
			vectorGLRendering<cpl::Types::v4sf>();
			break;
		default:
			vectorGLRendering<float>();
			break;
		}
	}

	template<typename V>
		void COscilloscope::vectorGLRendering()
		{
			cpl::CMutex lock(bufferLock);
			auto cStart = cpl::Misc::ClockCounter();
			CPL_DEBUGCHECKGL();

			handleFlagUpdates();

			calculateFundamentalPeriod();
			calculateTriggeringOffset();

			resizeAudioStorage();


			juce::OpenGLHelpers::clear(state.colourBackground);
			{
				cpl::OpenGLRendering::COpenGLStack openGLStack;
				// set up openGL
				openGLStack.setBlender(GL_ONE, GL_ONE_MINUS_SRC_COLOR);
				openGLStack.loadIdentityMatrix();
				cpl::GraphicsND::Transform3D<GLfloat> transform(1);
				content->transform.fillTransform3D(transform);
				openGLStack.applyTransform3D(transform);
				state.antialias ? openGLStack.enable(GL_MULTISAMPLE) : openGLStack.disable(GL_MULTISAMPLE);

				// the peak filter has to run on the whole buffer each time.
				if (state.envelopeMode == EnvelopeModes::PeakDecay)
				{
					runPeakFilter<V>(audioStream.getAudioBufferViews());
				}

				openGLStack.setLineSize(static_cast<float>(oglc->getRenderingScale()) * state.primitiveSize);
				openGLStack.setPointSize(static_cast<float>(oglc->getRenderingScale()) * state.primitiveSize);

				drawWavePlot<V>(openGLStack);
				CPL_DEBUGCHECKGL();

				openGLStack.setLineSize(static_cast<float>(oglc->getRenderingScale()) * 2.0f);

				CPL_DEBUGCHECKGL();
				renderCycles = cpl::Misc::ClockCounter() - cStart;
			}

			renderGraphics(
				[&](juce::Graphics & g)
				{
					// draw graph and wireframe
					//drawWireFrame<V>(openGLStack);
					paint2DGraphics<V>(g);

				}
			);

			auto tickNow = juce::Time::getHighResolutionTicks();
			avgFps.setNext(tickNow - lastFrameTick);
			lastFrameTick = tickNow;

		}

	template<typename V>
		void COscilloscope::drawWireFrame(juce::Graphics & g, const juce::Rectangle<float> rect, const float gain)
		{
			const auto xoff = rect.getX();
			const auto yoff = rect.getY();

			auto const minSpacing = rect.getHeight() / 50;
			// quantize to multiples of 3
			auto const numLines = 2 * (std::size_t)(1.5 + 0.5 * content->pctForDivision.getNormalizedValue() * minSpacing) - 1;

			g.setColour(content->skeletonColour.getAsJuceColour());

			const auto middle = rect.getHeight() * 0.5;
			const auto offset = state.envelopeGain;
			char textBuf[200];


			for (std::size_t i = 1; i < numLines; ++i)
			{
				auto const fraction = double(i) / (numLines - 1);
				auto const coord = fraction * middle;

				auto const y = middle + coord;
				auto const dBs = 20 * std::log10(fraction / offset);
				g.drawLine(xoff, yoff + y, xoff + rect.getWidth(), yoff + y);

				sprintf_s(textBuf, "%.3f dB", dBs);

				g.drawSingleLineText(textBuf, xoff + 5, yoff + y - 10, juce::Justification::left);
				g.drawSingleLineText(textBuf, xoff + 5, yoff + rect.getHeight() - (y - 15), juce::Justification::left);

				g.drawLine(xoff, yoff + rect.getHeight() - y, xoff + rect.getWidth(), yoff + rect.getHeight() - y);
			}

			g.drawLine(xoff, yoff + middle, xoff + rect.getWidth(), yoff + middle);

			drawTimeDivisions<V>(g, rect, rect.getHeight() / numLines);
		}

	template<typename V>
	void COscilloscope::drawTimeDivisions(juce::Graphics & g, juce::Rectangle<float> rect, double horizontalFractionGranularity)
	{
		auto const minVerticalSpacing = rect.getWidth() / 75;

		auto const wantedVerticalLines = (std::size_t)(0.5 + content->pctForDivision.getNormalizedValue() * minVerticalSpacing);

		const auto xoff = rect.getX();
		const auto yoff = rect.getY();

		if (wantedVerticalLines > 0)
		{
			char textBuf[200];
			auto const windowSize = 1000 * (state.effectiveWindowSize - 1) / audioStream.getAudioHistorySamplerate();
			double msIncrease = 0;
			double roundedPower = 0;
			auto cyclesTotal = state.effectiveWindowSize / (triggerState.cycleSamples);

			if (state.timeMode == OscilloscopeContent::TimeMode::Time)
			{
				constexpr double scaleTable[] = { 1, 2, 5, 10 };
				int size = std::extent<decltype(scaleTable)>::value;

				auto getIncrement = [&](int level)
				{

					if (level < 0)
					{
						return std::pow(2, level);
					}
					else if (level >= size)
					{
						// scale by magnitude difference
						return scaleTable[level % size] * std::pow(scaleTable[size - 1], level / size);
					}
					else
					{
						return scaleTable[level];
					}
				};


				auto index = 0;
				std::size_t count = 0;

				std::size_t numLines = 0;

				for (;;)
				{
					msIncrease = getIncrement(index);
					double incz1 = getIncrement(index - 1);
					std::size_t numLinesz1 = static_cast<int>(windowSize / incz1);
					numLines = static_cast<int>(windowSize / msIncrease);

					if (numLines > wantedVerticalLines)
						index++;
					else if (numLinesz1 > wantedVerticalLines)
						break;
					else
						index--;

					count++;

					// TODO: converge problem?
					if (count > 20)
						break;
				}
			}
			else if (state.timeMode == OscilloscopeContent::TimeMode::Cycles)
			{
				auto index = 0;
				std::size_t count = 0;

				std::size_t numLines = 0;

				auto numLinesPerCycles = wantedVerticalLines / cyclesTotal;
				roundedPower = std::pow(2, std::round(std::log2(numLinesPerCycles)));

				msIncrease = 1000 * (triggerState.cycleSamples) / audioStream.getAudioHistorySamplerate();
				msIncrease /= roundedPower;

			}
			else if (state.timeMode == OscilloscopeContent::TimeMode::Beats)
			{
				roundedPower = std::pow(2, std::round(std::log2(wantedVerticalLines)));
				msIncrease = windowSize / roundedPower;
			}

	
			g.setColour(content->skeletonColour.getAsJuceColour());
			double currentMsPos = 0;

			// (first line is useless)
			currentMsPos += msIncrease;
			int i = 0;
			while (currentMsPos < windowSize)
			{
				auto const fraction = currentMsPos / windowSize;
				auto const x = xoff + fraction * rect.getWidth();

				g.drawLine(x, 0, x, rect.getHeight());
				switch (state.timeMode)
				{
					case OscilloscopeContent::TimeMode::Time: 
						sprintf_s(textBuf, "%.3f ms", currentMsPos); 
						g.drawSingleLineText(textBuf, std::round(x + 5), rect.getHeight() - 10, juce::Justification::left);
						break;
					case OscilloscopeContent::TimeMode::Cycles: 
					{
						
						sprintf_s(textBuf, "%.3f ms", currentMsPos);
						g.drawSingleLineText(textBuf, std::round(x + 5), rect.getHeight() - 25, juce::Justification::left);

						double moduloI = std::fmod(i, roundedPower) + 1;
						sprintf_s(textBuf, "%.0f/%.0f (%.2f r)",
							roundedPower > 1 ? moduloI : moduloI / roundedPower,
							std::max(1.0, roundedPower),
							(1 / roundedPower) * moduloI * cpl::simd::consts<double>::tau
						);
						g.drawSingleLineText(textBuf, std::floor(x + 5), rect.getHeight() - 10, juce::Justification::left);
						break;
					}
					case OscilloscopeContent::TimeMode::Beats:
					{
						auto multiplier = std::max(state.beatDivision, roundedPower);
						sprintf_s(textBuf, "%.3f ms", currentMsPos);
						g.drawSingleLineText(textBuf, std::round(x + 5), rect.getHeight() - 25, juce::Justification::left);

						double moduloI = std::fmod(i, roundedPower) + 1;
						sprintf_s(textBuf, "%.0f/%.0f",
							moduloI * multiplier / roundedPower, multiplier
						);
						g.drawSingleLineText(textBuf, std::floor(x + 5), rect.getHeight() - 10, juce::Justification::left);
						break;

					}
				}
				



				currentMsPos += msIncrease;

				i++;

				if (i > 100)
					break;
			}

		}


	}

	template<typename V>
		void COscilloscope::drawWavePlot(cpl::OpenGLRendering::COpenGLStack & openGLStack)
		{
			{
				cpl::OpenGLRendering::MatrixModification matrixMod;
				// and apply the gain:
				const auto gain = (GLfloat)state.envelopeGain;
				matrixMod.scale(1, gain, 1);


				auto && view = lifoStream.createProxyView();

				cpl::ssize_t
					cursor = view.cursorPosition(),
					bufferSamples = view.size(),
					roundedWindow = static_cast<cpl::ssize_t>(std::ceil(state.effectiveWindowSize)),
					quantizedCycleSamples(0);

				const auto sizeMinusOne = std::max(1.0, state.effectiveWindowSize - 1);
				const auto sampleDisplacement = 2.0 / sizeMinusOne;
				auto cycleSamples = triggerState.cycleSamples;
				auto sampleOffset = triggerState.sampleOffset;
				auto offset = -1 + -(2 * triggerState.sampleOffset / sizeMinusOne);
				cpl::ssize_t bufferOffset = 0;
				double subSampleOffset = 0;

				if (state.triggerMode != OscilloscopeContent::TriggeringMode::None)
				{
					quantizedCycleSamples = static_cast<cpl::ssize_t>(std::ceil(triggerState.cycleSamples));
					subSampleOffset = (quantizedCycleSamples - triggerState.cycleSamples) + (roundedWindow - state.effectiveWindowSize);
					offset += 2 * ((1 - subSampleOffset) / sizeMinusOne);
				}
				else
				{
					subSampleOffset = cycleSamples = sampleOffset = 0;
					offset = -1;
				}

				if(state.timeMode != OscilloscopeContent::TimeMode::Beats)
				{
					bufferOffset = roundedWindow + quantizedCycleSamples;
				}
				else
				{
					auto const realOffset = std::fmod(state.transportPosition, state.effectiveWindowSize);
					bufferOffset = static_cast<cpl::ssize_t>(std::ceil(realOffset));
					subSampleOffset = (quantizedCycleSamples - triggerState.cycleSamples) + (realOffset - bufferOffset);
				}


				cpl::OpenGLRendering::PrimitiveDrawer<1024> drawer(openGLStack, GL_LINE_STRIP);

				drawer.addColour(state.colourDraw);


				auto pointer = view.begin() + (cursor - bufferOffset);

				while (pointer < view.begin())
					pointer += bufferSamples;

				roundedWindow = std::max(2, roundedWindow);

				auto end = view.end();

				for (cpl::ssize_t i = 0; i < roundedWindow + quantizedCycleSamples; ++i)
				{
					drawer.addVertex(static_cast<GLfloat>(i * sampleDisplacement + offset), *pointer++, 0);

					if (pointer == end)
						pointer -= bufferSamples;
				}
			}

			// draw end marks (temporary? move to 2d graphics anyway.)
			cpl::OpenGLRendering::PrimitiveDrawer<1024> drawer(openGLStack, GL_LINES);

			drawer.addColour(content->skeletonColour.getAsJuceColour());
			drawer.addVertex(-1, -1, 0);
			drawer.addVertex(-1, 1, 0);
			drawer.addVertex(1, -1, 0);
			drawer.addVertex(1, 1, 0);

		}

	template<typename V>
		void COscilloscope::runPeakFilter(const AudioStream::AudioBufferAccess & audio)
		{
			// TODO: Fix to extract data per-channel
			if (state.normalizeGain && audio.getNumChannels() >= 2)
			{
				AudioStream::AudioBufferView views[2] = { audio.getView(0), audio.getView(1) };

				std::size_t numSamples = views[0].size();

				double currentEnvelope = 0;
				// since this runs in every frame, we need to scale the coefficient by how often this function runs
				// (and the amount of samples)
				double power = numSamples * (avgFps.getAverage() / juce::Time::getHighResolutionTicksPerSecond());

				double coeff = std::pow(state.envelopeCoeff, power);

				// there is a number of optimisations we can do here, mostly that we actually don't care about
				// timing, we are only interested in the current largest value in the set.
				using namespace cpl;
				using namespace cpl::simd;
				using cpl::simd::load;

				V
					vLMax = zero<V>(),
					vRMax = zero<V>(),
					vSign = consts<V>::sign_mask;

				auto const loopIncrement = elements_of<V>::value;

				auto * leftBuffer = views[0].begin();
				auto * rightBuffer = views[1].begin();

				auto stop = numSamples - (numSamples & (loopIncrement - 1));
				if (stop <= 0)
					stop = 0;

				for (std::size_t i = 0; i < stop; i += loopIncrement)
				{
					auto const vLInput = loadu<V>(leftBuffer + i);
					vLMax = max(vand(vLInput, vSign), vLMax);
					auto const vRInput = loadu<V>(rightBuffer + i);
					vRMax = max(vand(vRInput, vSign), vRMax);
				}

				// TODO: remainder?

				suitable_container<V> lmax = vLMax, rmax = vRMax;

				double highestLeft = *std::max_element(lmax.begin(), lmax.end());
				double highestRight = *std::max_element(rmax.begin(), rmax.end());

				filters.envelope[0] = std::max(filters.envelope[0] * coeff, highestLeft  * highestLeft);
				filters.envelope[1] = std::max(filters.envelope[1] * coeff, highestRight * highestRight);

				currentEnvelope = 1.0 / std::max(std::sqrt(filters.envelope[0]), std::sqrt(filters.envelope[1]));

				if (std::isnormal(currentEnvelope))
				{
					content->inputGain.getParameterView().updateFromProcessorTransformed(
						currentEnvelope,
						cpl::Parameters::UpdateFlags::All & ~cpl::Parameters::UpdateFlags::RealTimeSubSystem
					);
				}
			}
		}
};
