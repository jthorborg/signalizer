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


#include "Oscilloscope.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include <cpl/rendering/OpenGLRasterizers.h>
#include <cpl/simd.h>
#include <cpl/special/AxisTools.h>
#include "SampleColourEvaluators.h"
#include "OscilloscopeDSP.inl"

namespace Signalizer
{

	struct VerticalScreenSplitter
	{
		VerticalScreenSplitter(juce::Rectangle<int> clipRectangle, cpl::OpenGLRendering::COpenGLStack & stack, int amountOfPasses, bool doSeparate)
			: stack(stack)
			, width(clipRectangle.getWidth())
			, height(clipRectangle.getHeight())
			, passesDone(0)
			, numPasses(amountOfPasses)
			, separate(doSeparate)
		{
			if (doSeparate)
			{
				m.scale(1, 1.0f / amountOfPasses, 1);
				// center on upper rect, equiv. to (1 - 1 / amountOfPasses) / (1 / amountOfPasses)
				m.translate(0, amountOfPasses - 1, 0);
				stack.enable(GL_SCISSOR_TEST);
				glScissor(0, cpl::Math::round<GLint>(GLfloat((numPasses - passesDone - 1) * height) / numPasses), width, height / numPasses);

			}

		}

		void nextPass()
		{
			if (separate)
			{				
				// one unit of vertical displacement is half a clipped region, so move two to center on next rect
				m.translate(0, -2, 0);
				passesDone++;

				glScissor(0, cpl::Math::round<GLint>(GLfloat((numPasses - passesDone - 1) * height) / numPasses), width, height / numPasses);
			}

		}

		~VerticalScreenSplitter() { if (separate) stack.disable(GL_SCISSOR_TEST); }

		cpl::OpenGLRendering::MatrixModification m;
		cpl::OpenGLRendering::COpenGLStack& stack;
		GLint width, height;
		GLint passesDone, numPasses;
		bool separate;
	};

	template<typename ISA>
	void Oscilloscope::paint2DGraphics(juce::Graphics & g)
	{
		if (content->diagnostics.getNormalizedValue() > 0.5)
		{
			g.setColour(juce::Colours::blue);

			const auto perf = audioStream->getPerfMeasures();

			float averageFps, averageCpu;
			computeAverageStats(averageFps, averageCpu);

			char textbuf[1024];

			cpl::sprintfs(textbuf, "%dx%d: %.1f fps - %.1f%% cpu, deltaG = %f, deltaO = %f (rt: %.2f%% - %.2f%%), (as: %.2f%% - %.2f%%), qHZ: %.5f - HZ: %.5f - PHASE: %.5f",
				getWidth(), getHeight(), averageFps, averageCpu, graphicsDeltaTime(), openGLDeltaTime(),
				100 * perf.producerUsage,
				100 * perf.producerOverhead,
				100 * perf.consumerUsage,
				100 * perf.consumerOverhead,
				(double)triggerState.record.index,
				triggerState.fundamental,
				triggerState.sampleOffset
			);

			g.drawSingleLineText(textbuf, 10, 20);

		}

		auto bounds = getLocalBounds().toFloat();

		if (state.colourAxis.getAlpha() != 0)
		{
			auto gain = static_cast<float>(getGain());
			drawTimeDivisions<ISA>(g, bounds);

			if (!shared.overlayChannels && shared.channelMode > OscChannels::OffsetForMono)
			{
				const auto cappedChannels = getEffectiveChannels();

				auto window = bounds.withHeight(bounds.getHeight() / cappedChannels);

				for (std::size_t c = 0; c < cappedChannels; ++c)
				{
					juce::Graphics::ScopedSaveState s(g);

					g.setColour(state.colourAxis);
					g.drawLine(0, window.getY(), window.getWidth(), window.getY());
					g.reduceClipRegion(window.toType<int>());

					drawWireFrame<ISA>(g, window, gain);

					window = window.withY(window.getY() + window.getHeight());

				}

			}
			else
			{
				drawWireFrame<ISA>(g, bounds, gain);
			}
		}

		auto mouseCheck = globalBehaviour->hideWidgetsOnMouseExit ? isMouseInside.load() : true;

		if (state.drawLegend && mouseCheck)
		{
			state.legend.paint(g, state.colourWidget, state.colourBackground);
		}

		if (state.drawCursorTracker && mouseCheck && getEffectiveChannels() > 0 /* HACK */)
		{
			g.setColour(state.colourWidget);

			const auto mouseX = currentMouse.x.load(), mouseY = currentMouse.y.load();

			double estimatedSize[2] = { 180, 70 };
			double textOffset[2] = { 20, -estimatedSize[1] };

			auto xpoint = mouseX + textOffset[0];
			if (xpoint + estimatedSize[0] > getWidth())
				xpoint = mouseX - (textOffset[0] + estimatedSize[0]);

			auto ypoint = mouseY + textOffset[1] - textOffset[0];
			if (ypoint < 0)
				ypoint = mouseY + textOffset[0];

			auto fraction = (double)mouseY / (getHeight() - 1);

			juce::Point<float> verticalMouseSection{ 0.0f, static_cast<float>(getBounds().getHeight()) };

			const auto cappedChannels = getEffectiveChannels();

			const auto normalizedSpace = 1.0 / cappedChannels;

			if (!shared.overlayChannels && shared.channelMode > OscChannels::OffsetForMono)
			{
				const auto position = fraction / normalizedSpace;
				const auto rounded = static_cast<int>(position);

				verticalMouseSection = { 
					static_cast<float>(getBounds().getHeight() * rounded * normalizedSpace), 
					static_cast<float>(getBounds().getHeight() * normalizedSpace)
				};

				fraction = std::fmod(fraction, normalizedSpace);
				fraction *= cappedChannels;
			}

			g.drawLine(0, mouseY, bounds.getWidth(), mouseY);
			g.drawLine(mouseX, verticalMouseSection.getX(), mouseX, verticalMouseSection.getX() + verticalMouseSection.getY());

			fraction = cpl::Math::UnityScale::linear(fraction, state.viewOffsets[VO::Top], state.viewOffsets[VO::Bottom]);
			fraction = 2 * (0.5 - fraction) / getGain();

			juce::Rectangle<float> rect{ (float)xpoint, (float)ypoint, (float)estimatedSize[0], (float)estimatedSize[1] };

			auto rectInside = rect.withSizeKeepingCentre(estimatedSize[0] * 0.95f, estimatedSize[1] * 0.95f).toType<int>();

			// clear background
			g.setColour(state.colourBackground);
			g.fillRoundedRectangle(rect, 2);

			// reset colour
			g.setColour(state.colourWidget);
			g.drawRoundedRectangle(rect, 2, 0.7f);

			auto const horizontalFraction = cpl::Math::UnityScale::linear((double)mouseX / (getWidth() - 1), state.viewOffsets[VO::Left], state.viewOffsets[VO::Right]);
			auto samples = (state.effectiveWindowSize - 1) * horizontalFraction;

			if (state.triggerMode == OscilloscopeContent::TriggeringMode::EnvelopeHold || state.triggerMode == OscilloscopeContent::TriggeringMode::ZeroCrossing)
			{
				samples -= (state.effectiveWindowSize - 1) *  0.5;
			}

			char text[1024];

			cpl::sprintfs(text,
				"y: %+.10f\ny: %+.10f\tdB\nx: %.10f\tms\nx: %.10f\tsmps",
				fraction,
				20 * std::log10(std::abs(fraction)),
				1e3 * samples / state.sampleRate,
				samples
			);

			g.setFont(juce::Font(juce::Font::getDefaultMonospacedFontName(), cpl::TextSize::normalText * 0.9f, 0));

			g.drawFittedText(text, rectInside, juce::Justification::centredLeft, 6);
		}

	}

	void Oscilloscope::onOpenGLRendering()
	{
		cpl::simd::dynamic_isa_dispatch<float, RenderingDispatcher>(*this);
	}

	bool Oscilloscope::checkAndInformInvalidCombinations(Oscilloscope::StreamState& cs)
	{
		if (state.timeMode == OscilloscopeContent::TimeMode::Cycles && cs.triggerMode != OscilloscopeContent::TriggeringMode::Spectral)
		{
			renderGraphics(
				[&](juce::Graphics & g)
				{
					g.setColour(cpl::GetColour(cpl::ColourEntry::Error));
					g.drawFittedText("Invalid combination of time and triggering modes", getLocalBounds(), juce::Justification::centred, 5);
				}
			);

			return false;
		}

		return true;
	}

	template<typename ISA>
		void Oscilloscope::vectorGLRendering()
		{
            CPL_DEBUGCHECKGL();
            
			{
				auto cs = processor->streamState.lock();
				auto& streamState = *cs;
                handleFlagUpdates(streamState);
                
                juce::OpenGLHelpers::clear(state.colourBackground);
                
				if (state.sampleRate == 0)
					return;

                if (!checkAndInformInvalidCombinations(streamState))
                    return;
                
				auto& channelData = streamState.channelData;

				cpl::OpenGLRendering::COpenGLStack openGLStack;
				// set up openGL
				openGLStack.setBlender(GL_ONE, GL_ONE_MINUS_SRC_COLOR);
				openGLStack.loadIdentityMatrix();
				state.antialias ? openGLStack.enable(GL_MULTISAMPLE) : openGLStack.disable(GL_MULTISAMPLE);
				openGLStack.setLineSize(static_cast<float>(oglc->getRenderingScale()) * state.primitiveSize);
				openGLStack.setPointSize(static_cast<float>(oglc->getRenderingScale()) * state.primitiveSize);


				CPL_DEBUGCHECKGL();

				auto mode = streamState.channelMode;

				const auto numChannels = shared.numChannels.load();
				CPL_RUNTIME_ASSERTION((numChannels % 2) == 0);
				// should be safe to assert equality 
				CPL_RUNTIME_ASSERTION(numChannels <= channelData.filterStates.channels.size());

				const auto triggerChannel = std::min(numChannels, static_cast<std::size_t>(cs->triggeringChannel)) - 1;
				// Pre-analyse triggering and state
				switch (mode)
				{
					default: case OscChannels::Left:
						analyseAndSetupState<ISA, SampleColourEvaluator<OscChannels::Left>>({ channelData, triggerChannel / 2 }, streamState); break;
					case OscChannels::Right:
						analyseAndSetupState<ISA, SampleColourEvaluator<OscChannels::Right>>({ channelData, triggerChannel / 2 }, streamState); break;
					case OscChannels::Mid:
						analyseAndSetupState<ISA, SampleColourEvaluator<OscChannels::Mid>>({ channelData, triggerChannel / 2 }, streamState); break;
					case OscChannels::Side:
						analyseAndSetupState<ISA, SampleColourEvaluator<OscChannels::Side>>({ channelData, triggerChannel / 2 }, streamState); break;
					case OscChannels::Separate:
						analyseAndSetupState<ISA, DynamicChannelEvaluator>({ channelData, triggerChannel }, streamState); break;
					case OscChannels::MidSide:
					{
						if ((triggerChannel & 0x1) == 0)
							analyseAndSetupState<ISA, SampleColourEvaluator<OscChannels::Mid>>({ channelData, triggerChannel }, streamState);
						else
							analyseAndSetupState<ISA, SampleColourEvaluator<OscChannels::Side>>({ channelData, triggerChannel & ~0x1 }, streamState);
						break;
					}
				}

				const auto numSplits = static_cast<int>(numChannels / (mode > OscChannels::OffsetForMono ? 1 : 2));
				VerticalScreenSplitter w(getLocalBounds() * oglc->getRenderingScale(), openGLStack, numSplits, !shared.overlayChannels);

				for (std::size_t channelPair = 0; channelPair < numChannels; channelPair += 2)
				{
					EvaluatorParams params { channelData, channelPair };

					switch (mode)
					{
						default: case OscChannels::Left:
							drawWavePlot<ISA, SampleColourEvaluator<OscChannels::Left>>(openGLStack, params, streamState); break;
						case OscChannels::Right:
							drawWavePlot<ISA, SampleColourEvaluator<OscChannels::Right>>(openGLStack, params, streamState); break;
						case OscChannels::Mid:
							drawWavePlot<ISA, SampleColourEvaluator<OscChannels::Mid>>(openGLStack, params, streamState); break;
						case OscChannels::Side:
							drawWavePlot<ISA, SampleColourEvaluator<OscChannels::Side>>(openGLStack, params, streamState); break;
						case OscChannels::Separate:
						{
							drawWavePlot<ISA, DynamicChannelEvaluator>(openGLStack, params, streamState);
							w.nextPass();
							params.channelIndex++;
							drawWavePlot<ISA, DynamicChannelEvaluator>(openGLStack, params, streamState);

							break;
						}
						case OscChannels::MidSide:
						{
							drawWavePlot<ISA, SampleColourEvaluator<OscChannels::Mid>>(openGLStack, params, streamState);
							w.nextPass();
							drawWavePlot<ISA, SampleColourEvaluator<OscChannels::Side>>(openGLStack, params, streamState);
							break;
						}
					}

					w.nextPass();
				}

				CPL_DEBUGCHECKGL();
			}

			renderGraphics(
				[&](juce::Graphics & g)
				{
					// draw graph and wireframe
					paint2DGraphics<ISA>(g);
				}
			);

			postFrame();
		}

	template<typename ISA>
		void Oscilloscope::drawWireFrame(juce::Graphics & g, const juce::Rectangle<float> rect, const float gain)
		{
			const auto xoff = rect.getX();
			const auto yoff = rect.getY();

			const auto verticalDelta = state.viewOffsets[VO::Bottom] - state.viewOffsets[VO::Top];

			auto const minSpacing = (rect.getHeight() / 50) / verticalDelta;
			// quantize to multiples of 3
			auto const numLines = 2 * (std::size_t)(1.5 + 0.5 * (1 - content->pctForDivision.getNormalizedValue()) * minSpacing) - 1;

			g.setColour(state.colourAxis);

			const auto offset = gain;
			char textBuf[200];

			auto viewTransform = [&](auto pos) {
				return (pos + (state.viewOffsets[VO::Bottom] - 1)) / verticalDelta;
			};

			auto drawMarkerAt = [&](auto where, auto size) {

				auto const transformedPos = viewTransform(where);
				auto const coord = transformedPos * rect.getHeight();

				auto const waveSpace = std::abs(2 * where - 1);

				auto const y = coord;
				auto const dBs = 20 * std::log10(waveSpace / offset);

				cpl::sprintfs(textBuf, "%.3f dB", dBs);

				g.drawSingleLineText(textBuf, xoff + 5, yoff + rect.getHeight() - (y - 15), juce::Justification::left);

				g.drawLine(xoff, yoff + rect.getHeight() - y, xoff + rect.getWidth(), yoff + rect.getHeight() - y, size);
			};

			// -300 dB
			auto const zeroDBEpsilon = 0.000000000000001;

			auto
				inc = 1.0 / (numLines - 1),
				end = 1 - state.viewOffsets[VO::Top];

			auto start = cpl::Math::roundToNextMultiplier(1 - state.viewOffsets[VO::Bottom], inc);

			auto unitSpacePos = start + inc;


			while (unitSpacePos < end)
			{
				if(std::abs(unitSpacePos - 0.5) > zeroDBEpsilon)
					drawMarkerAt(unitSpacePos, 1.5f);

				unitSpacePos += inc;
			}

			if(start < 0.5 && end > 0.5)
				drawMarkerAt(0.5, 1.5f);
		}

	template<typename ISA>
	void Oscilloscope::drawTimeDivisions(juce::Graphics & g, juce::Rectangle<float> rect)
	{
		auto const horizontalDelta = (state.viewOffsets[VO::Right] - state.viewOffsets[VO::Left]);
		auto const minVerticalSpacing = (rect.getWidth() / (state.timeMode == OscilloscopeContent::TimeMode::Time ? 75 : 110)) / horizontalDelta;
		auto const wantedVerticalLines = (std::size_t)(0.5 + (1 - content->pctForDivision.getNormalizedValue()) * minVerticalSpacing);

		const auto xoff = rect.getX();
		const auto yoff = rect.getY();

		if (wantedVerticalLines > 0)
		{
			char textBuf[200];
			auto const windowSize = 1000 * (state.effectiveWindowSize - 1) / state.sampleRate;
			if (windowSize == 0)
				return;

			double msIncrease = 0;
			double roundedPower = 0;
			auto cyclesTotal = state.effectiveWindowSize / (triggerState.cycleSamples);

			if (state.timeMode == OscilloscopeContent::TimeMode::Time)
			{
				constexpr double scaleTable[] = { 1, 2, 5, 10 };
				msIncrease = cpl::special::SuitableAxisDivision<double>(scaleTable, wantedVerticalLines, windowSize);
			}
			else if (state.timeMode == OscilloscopeContent::TimeMode::Cycles)
			{
				auto numLinesPerCycles = wantedVerticalLines / cyclesTotal;
				roundedPower = std::pow(2, std::round(std::log2(numLinesPerCycles)));
				msIncrease = 1000 * (triggerState.cycleSamples) / state.sampleRate;
				msIncrease /= roundedPower;

			}
			else if (state.timeMode == OscilloscopeContent::TimeMode::Beats)
			{
				roundedPower = std::pow(2, std::round(std::log2(wantedVerticalLines)));
				msIncrease = windowSize / roundedPower;
			}


			g.setColour(state.colourAxis);

			auto transformView = [&](auto x) {
				return (x - state.viewOffsets[VO::Left]) / horizontalDelta;
			};

			auto end = state.viewOffsets[VO::Right] * windowSize;

			auto start = state.viewOffsets[VO::Left] * windowSize;
			double offset = 0;
			if (state.triggerMode == OscilloscopeContent::TriggeringMode::EnvelopeHold || state.triggerMode == OscilloscopeContent::TriggeringMode::ZeroCrossing)
			{
				offset = windowSize * 0.5;
				end -= offset;
				start -= offset;
			}

			auto multiplier = start / msIncrease;

			auto i = static_cast<int>(std::floor(multiplier)) - 1;

			auto currentMsPos = cpl::Math::roundToNextMultiplier(start, msIncrease);


			while (currentMsPos < end)
			{
				double moduloI = std::fmod(i, roundedPower) + 1;

				auto const fraction = (currentMsPos + offset) / windowSize;
				auto const x = xoff + transformView(fraction) * rect.getWidth();
				auto const samples = 1e-3 * currentMsPos * state.sampleRate;

				g.drawLine(x, yoff, x, rect.getHeight(), samples == 0 ? 1.5f : 1.0f);

				float offset = 10;

				auto textOut = [&](auto format, auto... args) {
					cpl::sprintfs(textBuf, format, args...);
					g.drawSingleLineText(textBuf, std::floor(x + 5), yoff + rect.getHeight() - offset, juce::Justification::left);
					offset += 15;
				};

				switch (state.timeMode)
				{
					case OscilloscopeContent::TimeMode::Cycles:
					{
						textOut("%.0f/%.0f (%.2f r)",
							roundedPower > 1 ? moduloI : moduloI / roundedPower,
							std::max(1.0, roundedPower),
							(1 / roundedPower) * moduloI * cpl::simd::consts<double>::tau
						);
						break;
					}
					case OscilloscopeContent::TimeMode::Beats:
					{
						textOut("%.0f/%.0f",
							(moduloI + 1), roundedPower
						);
						break;
					}
				}

				textOut("%.2f smps", samples);
				textOut("%.4f ms", currentMsPos);

				currentMsPos += msIncrease;
				i++;
			}
		}
	}

	template<typename ISA, typename Evaluator>
		void Oscilloscope::drawWavePlot(cpl::OpenGLRendering::COpenGLStack& openGLStack, const EvaluatorParams& params, Oscilloscope::StreamState& cs)
		{

			typedef cpl::OpenGLRendering::PrimitiveDrawer<1024> Renderer;

			cpl::OpenGLRendering::MatrixModification matrixMod;
			// and apply the gain:
			const auto gain = static_cast<GLfloat>(getGain());

			auto
				left = state.viewOffsets[VO::Left],
				right = state.viewOffsets[VO::Right],
				top = state.viewOffsets[VO::Top],
				bottom = state.viewOffsets[VO::Bottom];

			auto verticalDelta = bottom - top;
			auto horizontalDelta = right - left;


			cpl::ssize_t
				roundedWindow = static_cast<cpl::ssize_t>(std::ceil(state.effectiveWindowSize)),
				quantizedCycleSamples(0);

			const auto sizeMinusOne = std::max(1.0, state.effectiveWindowSize - 1);
			const auto sampleDisplacement = 1.0 / sizeMinusOne;
			cpl::ssize_t bufferOffset = 0;
			double subSampleOffset = 0, offset = 0;
			auto const pixelsPerSample = oglc->getRenderingScale() * std::abs((getWidth() - 1) / (sizeMinusOne * (horizontalDelta)));
			const auto triggerMode = state.triggerMode;

			auto interpolation = state.sampleInterpolation;
			if (pixelsPerSample < 1 && state.sampleInterpolation != SubSampleInterpolation::None)
			{
				interpolation = SubSampleInterpolation::Linear;
			}

			if (triggerMode == OscilloscopeContent::TriggeringMode::Window)
			{
				auto const realOffset = std::fmod(cs.transportPosition, state.effectiveWindowSize);
				bufferOffset = static_cast<cpl::ssize_t>(std::ceil(realOffset));
				offset = subSampleOffset = 0;
			}
			else if (triggerMode == OscilloscopeContent::TriggeringMode::EnvelopeHold || triggerMode == OscilloscopeContent::TriggeringMode::ZeroCrossing)
			{
				auto const realOffset = triggerState.sampleOffset;
				bufferOffset = static_cast<cpl::ssize_t>(std::ceil(realOffset));
				subSampleOffset = bufferOffset - realOffset;
				offset += (1 - subSampleOffset) / sizeMinusOne;
			}
			else
			{
				// TODO: FIX. This causes an extra cycle to be rendered for lanczos so it has more to eat from the edges.
				auto const cycleBuffers = interpolation == SubSampleInterpolation::Lanczos ? 2 : 1;
				// calculate fractionate offsets used for sample-space rendering
				if (state.triggerMode != OscilloscopeContent::TriggeringMode::None)
				{
					quantizedCycleSamples = static_cast<cpl::ssize_t>(std::ceil(triggerState.cycleSamples));
					subSampleOffset = cycleBuffers * (quantizedCycleSamples - triggerState.cycleSamples) + (roundedWindow - state.effectiveWindowSize);
					offset = -triggerState.sampleOffset / sizeMinusOne;
				}
				bufferOffset = roundedWindow + cycleBuffers * quantizedCycleSamples;
				offset += (1 - subSampleOffset) / sizeMinusOne;
			}


			roundedWindow = std::max<cpl::ssize_t>(2, roundedWindow);

			// modify the horizontal axis into [0, 1] instead of [-1, 1]
			matrixMod.translate(-1, 0, 0);
			matrixMod.scale(2, 1, 1);

			// apply horizontal transformation
			matrixMod.scale(1 / (horizontalDelta), 1, 1);
			matrixMod.translate(-left, 0, 0);

			// apply vertical transformation
			matrixMod.scale(1, 1.0 / verticalDelta, 0);
			matrixMod.translate(0, top + (bottom - 1), 0);
			matrixMod.scale(1, gain, 0);

			const GLfloat endCondition = static_cast<GLfloat>(roundedWindow + quantizedCycleSamples /* + 2 */);

			auto renderSampleSpace = [&](auto kernel, GLint primitive, cpl::ssize_t sampleOffset = 0)
			{
				cpl::OpenGLRendering::MatrixModification m;
				// translate triggering offset + 1
				matrixMod.translate(offset - sampleDisplacement, 0, 0);
				// scale to sample/pixels space
				matrixMod.scale(sampleDisplacement, 1, 1);

				Evaluator eval(params);
				if (!eval.isWellDefined())
					return;

				eval.startFrom(-(bufferOffset + sampleOffset), -(bufferOffset + sampleOffset));
				cpl::OpenGLRendering::PrimitiveDrawer<1024> drawer(openGLStack, primitive);
				kernel(eval, drawer);
			};

			auto dotSamples = [&] (cpl::ssize_t offset)
			{
				auto oldPointSize = openGLStack.getPointSize();

                auto normScale = 1.0 / std::sqrt(oglc->getRenderingScale());
                
				if (pixelsPerSample * normScale > 5 && state.sampleInterpolation != SubSampleInterpolation::None)
				{
                    openGLStack.setPointSize(oldPointSize * 4 * normScale);
				}

				if (state.colourChannelsByFrequency)
				{
					renderSampleSpace(
						// nested lambda auto evaluation fails on vc 2015
						[&] (Evaluator & evaluator, Renderer & drawer)
						{
							for (GLfloat i = 0; i < endCondition; i += 1)
							{
								const auto data = evaluator.evaluate();
								drawer.addColour(data.second);
								drawer.addVertex(i, data.first, 0);
								evaluator.inc();
							}
						},
						GL_POINTS,
						offset
					);
				}
				else
				{
					renderSampleSpace(
						[&] (Evaluator & evaluator, Renderer & drawer)
						{
							drawer.addColour(evaluator.getDefaultKey());

							for (GLfloat i = 0; i < endCondition; i += 1)
							{
								drawer.addVertex(i, evaluator.evaluateSampleInc(), 0);
							}
						},
						GL_POINTS,
						offset
					);
				}
				openGLStack.setPointSize(oldPointSize);
			};

			// draw dots for very zoomed displays and when there's no subsample interpolation
			if ((state.dotSamples && pixelsPerSample > 5) || state.sampleInterpolation == SubSampleInterpolation::None)
			{
				dotSamples(0);
			}

			// TODO: Add scaled rendering (getAttachedContext()->getRenderingScale())
			switch (interpolation)
			{
				case SubSampleInterpolation::Linear:
				{
					if (state.colourChannelsByFrequency)
					{
						renderSampleSpace(
							[&] (auto & evaluator, auto & drawer)
							{
								for (GLfloat i = 0; i < endCondition; i += 1)
								{
									const auto data = evaluator.evaluate();

									drawer.addColour(data.second);
									drawer.addVertex(i, data.first, 0);

									evaluator.inc();
								}
							},
							GL_LINE_STRIP
						);
					}
					else
					{
						renderSampleSpace(
							[&] (auto & evaluator, auto & drawer)
							{
								drawer.addColour(evaluator.getDefaultKey());

								for (GLfloat i = 0; i < endCondition; i += 1)
								{
									drawer.addVertex(i, evaluator.evaluateSampleInc(), 0);
								}
							},
							GL_LINE_STRIP
						);
					}

					break;
				}
				case SubSampleInterpolation::Rectangular:
				{
					if (state.colourChannelsByFrequency)
					{
						renderSampleSpace(
							[&] (auto & evaluator, auto & drawer)
							{
								typename Evaluator::ColourT oldColour = evaluator.evaluateColour();

								for (GLfloat i = 0; i < endCondition; i += 1)
								{
									const auto data = evaluator.evaluate();

									drawer.addColour(oldColour);
									drawer.addVertex(i, data.first, 0);
									drawer.addColour(data.second);
									drawer.addVertex(i + 1, data.first, 0);

									oldColour = data.second;

									evaluator.inc();
								}
							},
							GL_LINE_STRIP
						);
					}
					else
					{
						renderSampleSpace(
							[&](auto & evaluator, auto & drawer)
							{
								drawer.addColour(evaluator.getDefaultKey());
								for (GLfloat i = 0; i < endCondition; i += 1)
								{
									const auto vertex = evaluator.evaluateSampleInc();
									drawer.addVertex(i, vertex, 0);
									drawer.addVertex(i + 1, vertex, 0);
								}
							},
							GL_LINE_STRIP
						);
					}
					break;
				}
				case SubSampleInterpolation::Lanczos:
				{

					auto const KernelSize = OscilloscopeContent::InterpolationKernelSize;
					auto const KernelBufferSize = KernelSize * 2 + 1;

					double samplePos = 0;

					if (triggerMode == OscilloscopeContent::TriggeringMode::Window)
					{
						samplePos = std::fmod(cs.transportPosition, state.effectiveWindowSize) - 1;
					}
					else if (triggerMode == OscilloscopeContent::TriggeringMode::EnvelopeHold || triggerMode == OscilloscopeContent::TriggeringMode::ZeroCrossing)
					{
						samplePos = triggerState.sampleOffset;
					}
					else
					{
						// TODO: possible small graphic glitch here if the interpolation eats into the next cycle.
						// calculations different here, depending on how you interpret phase information in the frequency domain
						samplePos = triggerState.cycleSamples * 2 + state.effectiveWindowSize - triggerState.sampleOffset;
					}

					// otherwise we will have a discontinuity as the interpolation kernel moves past T = 0
					if (triggerMode == OscilloscopeContent::TriggeringMode::None || triggerMode == OscilloscopeContent::TriggeringMode::Window)
					{
						samplePos = std::ceil(samplePos);
						//if (state.timeMode == OscilloscopeContent::TimeMode::Time)
						//	samplePos += KernelSize;
					}

					// adjust for left
					double inc = horizontalDelta / (oglc->getRenderingScale() * (getWidth() - 1));
					double unitSpacePos = left;
					double samplesPerPixel = 1.0 / (pixelsPerSample);

					samplePos += -unitSpacePos / inc * samplesPerPixel;
					double currentSample = std::floor(samplePos);

					Evaluator eval(params);

					if (!eval.isWellDefined())
						return;

					eval.startFrom(-static_cast<cpl::ssize_t>(std::floor(samplePos)) - (int)KernelSize, 2 - KernelBufferSize - static_cast<cpl::ssize_t>(std::floor(samplePos)));

					AFloat kernel[KernelBufferSize];

					typename Evaluator::ColourT currentColour, nextColour;

					auto get = [&]() {
						auto data = eval.evaluate();
						eval.inc();
						currentColour = nextColour;
						nextColour = data.second;

						return data.first;
					};

					auto insert = [&] (auto val) {
						std::rotate(kernel, kernel + 1, kernel + KernelBufferSize);
						kernel[KernelBufferSize - 1] = val;
					};

					std::for_each(std::begin(kernel), std::end(kernel), [&](auto & f) { f = get(); });

					{
						cpl::OpenGLRendering::PrimitiveDrawer<1024> drawer(openGLStack, GL_LINE_STRIP);
						if (!state.colourChannelsByFrequency)
							drawer.addColour(eval.getDefaultKey());
						else
							drawer.addColour(currentColour);

						do
						{
							auto delta = currentSample - samplePos;

							while (delta > 1)
							{
								samplePos += 1;
								delta -= 1;
								insert(get());
							}

							const auto interpolatedValue = cpl::dsp::lanczosFilter<double>(kernel, KernelBufferSize, (KernelSize) + delta, KernelSize);

							if (state.colourChannelsByFrequency)
							{
								drawer.addColour(currentColour.lerp(nextColour, delta));
							}

							drawer.addVertex(unitSpacePos, interpolatedValue, 0);
							currentSample += samplesPerPixel;

							unitSpacePos += inc;

						} while (unitSpacePos < (right + inc));
					}


					break;
				}

			}

		}
};
