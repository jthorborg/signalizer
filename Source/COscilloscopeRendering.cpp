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

	static const char * ChannelDescriptions[] = { "+L", "+R", "-L", "-R", "L", "R", "C", "S"};

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
			sprintf(textbuf.get(), "%dx%d: %.1f fps - %.1f%% cpu, deltaG = %f, deltaO = %f (rt: %.2f%% - %.2f%%), (as: %.2f%% - %.2f%%), qHZ: %.5f - HZ: %.5f",
				getWidth(), getHeight(), fps, cpuTime, graphicsDeltaTime(), openGLDeltaTime(),
				100 * audioStream.getPerfMeasures().rtUsage.load(std::memory_order_relaxed),
				100 * audioStream.getPerfMeasures().rtOverhead.load(std::memory_order_relaxed),
				100 * audioStream.getPerfMeasures().asyncUsage.load(std::memory_order_relaxed),
				100 * audioStream.getPerfMeasures().asyncOverhead.load(std::memory_order_relaxed),
				quantizedFreq,
				detectedFreq);
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

			CPL_DEBUGCHECKGL();
			auto && lockedView = audioStream.getAudioBufferViews();
			handleFlagUpdates();
			auto cStart = cpl::Misc::ClockCounter();
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
					runPeakFilter<V>(lockedView);
				}

				triggerOffset = getTriggeringOffset(lockedView);

				openGLStack.setLineSize(static_cast<float>(oglc->getRenderingScale()) * state.primitiveSize);
				openGLStack.setPointSize(static_cast<float>(oglc->getRenderingScale()) * state.primitiveSize);

				drawWavePlot<V>(openGLStack, lockedView);

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

			const bool nicelyQuantize = true;
			
			if(true || nicelyQuantize)
			{
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

					g.drawSingleLineText(textBuf, xoff + 5, yoff + y - 10, juce::Justification::centredLeft);
					g.drawSingleLineText(textBuf, xoff + 5, yoff + rect.getHeight() - (y - 15), juce::Justification::centredLeft);

					g.drawLine(xoff, yoff + rect.getHeight() - y, xoff + rect.getWidth(), yoff + rect.getHeight() - y);
				}

				g.drawLine(xoff, yoff + middle, xoff + rect.getWidth(), yoff + middle);
			}

			auto const minVerticalSpacing = rect.getWidth() / 75;

			auto const wantedVerticalLines = (std::size_t)(0.5 + content->pctForDivision.getNormalizedValue() * minVerticalSpacing);

			if(wantedVerticalLines > 0)
			{

				g.setColour(content->skeletonColour.getAsJuceColour());
				char textBuf[200];


				// quantize to multiples of 3

				if (nicelyQuantize)
				{
					double scaleTable[] = { 1, 2, 5, 10 };
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

					auto const currentMS = content->windowSize.getTransformedValue() * 1000 / audioStream.getAudioHistorySamplerate();

					auto index = 0;
					std::size_t count = 0;
					double inc = 0;
					std::size_t numLines = 0;

					for (;;)
					{
						inc = getIncrement(index);
						double incz1 = getIncrement(index - 1);
						std::size_t numLinesz1 = static_cast<int>(currentMS / incz1);
						numLines = static_cast<int>(currentMS / inc);

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

					double currentPos = 0;

					// (first line is useless)
					currentPos += inc;

					while (currentPos < currentMS)
					{
						auto const fraction = currentPos / currentMS;
						auto const x = xoff + fraction * rect.getWidth();

						g.drawLine(x, 0, x, rect.getHeight());

						sprintf_s(textBuf, "%.2f ms", currentPos);

						g.drawSingleLineText(textBuf, x + 5, rect.getHeight() - 10, juce::Justification::centredLeft);

						currentPos += inc;
					}

				}
				else
				{
					for (std::size_t i = 1; i < wantedVerticalLines; ++i)
					{
						auto const fraction = double(i) / (wantedVerticalLines - 1);
						auto const x = xoff + fraction * rect.getWidth();

						auto const time = content->windowSize.getTransformedValue() * 1000 * fraction / audioStream.getAudioHistorySamplerate();

						g.drawLine(x, 0, x, rect.getHeight());

						sprintf_s(textBuf, "%.2f ms", time);

						g.drawSingleLineText(textBuf, x + 5, rect.getHeight() - 10, juce::Justification::centredLeft);
					}
				}

			}

		}


	template<typename V>
		void COscilloscope::drawWavePlot(cpl::OpenGLRendering::COpenGLStack & openGLStack, const AudioStream::AudioBufferAccess & audio)
		{

			if (audio.getNumChannels() < 1)
				return;

			cpl::OpenGLRendering::MatrixModification matrixMod;
			// and apply the gain:
			auto gain = (GLfloat)state.envelopeGain;
			matrixMod.scale(1, gain, 1);
			float sampleDisplacement = 2.0f / std::max<int>(1, static_cast<int>(audio.getNumSamples() - 1));

			cpl::OpenGLRendering::PrimitiveDrawer<1024> drawer(openGLStack, GL_LINE_STRIP);

			drawer.addColour(state.colourDraw);

			auto const offset = -1 + -(2 * triggerOffset / static_cast<int>(audio.getNumSamples() - 1));

			// TODO: glDrawArrays
			audio.iterate<1, true>
			(
				[&] (std::size_t sampleFrame, AudioStream::DataType & left)
				{
					drawer.addVertex(sampleFrame * sampleDisplacement + offset, left, 0);
				}
			);


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
