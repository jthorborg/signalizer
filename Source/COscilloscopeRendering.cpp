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
#include <tuple>

namespace Signalizer
{
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

			const auto verticalDelta = state.viewOffsets[VO::Bottom] - state.viewOffsets[VO::Top];

			auto const minSpacing = (rect.getHeight() / 50) / verticalDelta;
			// quantize to multiples of 3
			auto const numLines = 2 * (std::size_t)(1.5 + 0.5 * content->pctForDivision.getNormalizedValue() * minSpacing) - 1;

			g.setColour(state.colourGraph);

			const auto offset = state.envelopeGain;
			char textBuf[200];


			auto viewTransform = [&](auto pos) {
				return (pos + (state.viewOffsets[VO::Bottom] - 1)) / verticalDelta;
			};

			auto drawMarkerAt = [&](auto where) {

				auto const transformedPos = viewTransform(where);
				auto const coord = transformedPos * rect.getHeight();

				auto const waveSpace = std::abs(2 * where - 1);

				auto const y = coord;
				auto const dBs = 20 * std::log10(waveSpace / offset);

				sprintf_s(textBuf, "%.3f dB", dBs);

				g.drawSingleLineText(textBuf, xoff + 5, yoff + rect.getHeight() - (y - 15), juce::Justification::left);

				g.drawLine(xoff, yoff + rect.getHeight() - y, xoff + rect.getWidth(), yoff + rect.getHeight() - y);
			};

			// -300 dB
			auto const zeroDBEpsilon = 0.000000000000001;

			auto
				inc = 1.0 / (numLines - 1),
				end = 1 - state.viewOffsets[VO::Top];

			auto bottom = viewTransform(0);
			auto top = viewTransform(1);

			auto unitSpacePos = cpl::Math::roundToNextMultiplier(1 - state.viewOffsets[VO::Bottom], inc);


			while (unitSpacePos < end)
			{
				if(std::abs(unitSpacePos - 0.5) > zeroDBEpsilon)
					drawMarkerAt(unitSpacePos);


				unitSpacePos += inc;
			}

			drawMarkerAt(0.5);

			drawTimeDivisions<V>(g, rect, rect.getHeight() / numLines);
		}

	template<typename V>
	void COscilloscope::drawTimeDivisions(juce::Graphics & g, juce::Rectangle<float> rect, double horizontalFractionGranularity)
	{
		auto const horizontalDelta = (state.viewOffsets[VO::Right] - state.viewOffsets[VO::Left]);
		auto const minVerticalSpacing = (rect.getWidth() / (state.timeMode == OscilloscopeContent::TimeMode::Time ? 75 : 110)) / horizontalDelta;
		auto const wantedVerticalLines = (std::size_t)(0.5 + content->pctForDivision.getNormalizedValue() * minVerticalSpacing);

		const auto xoff = rect.getX();
		const auto yoff = rect.getY();

		if (wantedVerticalLines > 0)
		{
			char textBuf[200];
			auto const windowSize = 1000 * (state.effectiveWindowSize - 1) / audioStream.getAudioHistorySamplerate();
			if (windowSize == 0)
				return;

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

				int index = 0;
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

					// TODO: converge problem? - also, provide inital estimate
					if (count > 20)
						break;
				}
			}
			else if (state.timeMode == OscilloscopeContent::TimeMode::Cycles)
			{
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

	
			g.setColour(state.colourGraph);

			auto transformView = [&](auto x) {
				return (x - state.viewOffsets[VO::Left]) / horizontalDelta;
			};

			auto end = state.viewOffsets[VO::Right] * windowSize;

			auto start = state.viewOffsets[VO::Left] * windowSize;
			auto multiplier = start / msIncrease;

			auto i = static_cast<int>(std::floor(multiplier)) - 1;

			auto currentMsPos = cpl::Math::roundToNextMultiplier(start, msIncrease);

			while (currentMsPos < end)
			{
				double moduloI = std::fmod(i, roundedPower) + 1;

				auto const fraction = currentMsPos / windowSize;
				auto const x = xoff + transformView(fraction) * rect.getWidth();
				auto const samples = 1e-3 * currentMsPos * audioStream.getAudioHistorySamplerate();
				g.drawLine(x, 0, x, rect.getHeight());

				float offset = 10;

				auto textOut = [&](auto format, auto... args) {
					sprintf_s(textBuf, format, args...);
					g.drawSingleLineText(textBuf, std::floor(x + 5), rect.getHeight() - offset, juce::Justification::left);
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


	template<>
	class COscilloscope::SampleColourEvaluator<OscChannels::Left>
	{
	public:

		typedef ChannelData::AudioBuffer::ProxyView::const_iterator AudioIt;
		typedef ChannelData::ColourBuffer::ProxyView::const_iterator ColourIt;
		typedef ChannelData::AudioBuffer::ProxyView::value_type AudioT;
		typedef ChannelData::ColourBuffer::ProxyView::value_type ColourT;

		SampleColourEvaluator(COscilloscope & oscilloscope)
			: osc(oscilloscope)
			, audioView(oscilloscope.channelData.channels[0].audioData.createProxyView())
			, colourView(oscilloscope.channelData.channels[0].colourData.createProxyView())
		{

		}

		inline bool isWellDefined() const noexcept
		{
			return audioView.size() > 0 && colourView.size() > 0;
		}

		void startFrom(cpl::ssize_t audioOffset, cpl::ssize_t colourOffset)
		{
			audioPointer = audioView.begin() + audioView.cursorPosition() + audioOffset;

			while (audioPointer < audioView.begin())
				audioPointer += audioView.size();

			while (audioPointer >= audioView.end())
				audioPointer -= audioView.size(); 

			colourPointer = colourView.begin() + colourView.cursorPosition() + colourOffset;

			while (colourPointer < colourView.begin())
				colourPointer += colourView.size();

			while (colourPointer >= colourView.end())
				colourPointer -= colourView.size(); 

		}

		/*inline void move(cpl::ssize_t offset) noexcept
		{
			auto diff = distance();

			audioPointer = audioView.begin() + (diff.first + offset) % audioView.size();
			colourPointer = colourView.begin() + (diff.second + offset) % audioView.size();
		} */

		inline void inc() noexcept
		{
			audioPointer++, colourPointer++;

			if (audioPointer == audioView.end())
				audioPointer -= audioView.size();

			if (colourPointer == colourView.end())
				colourPointer -= colourView.size();
		}

		inline std::pair<std::ptrdiff_t, std::ptrdiff_t> distance() const noexcept
		{
			return std::make_pair(
				std::distance<AudioIt>(audioPointer, audioView.begin()),
				std::distance<ColourIt>(colourPointer, colourView.begin())
			);
		}

		inline std::pair<AudioT, ColourT> evaluate() const noexcept
		{
			return { *audioPointer, *colourPointer };
		}

		AudioT evaluateSample() const noexcept
		{
			return *audioPointer;
		}

		ColourT evaluateColour() const noexcept
		{
			return *colourPointer;
		}

		AudioT evaluateSampleInc() noexcept
		{
			auto ret = *audioPointer++;

			if (audioPointer == audioView.end())
				audioPointer -= audioView.size();

			return ret;
		}

		ColourT evaluateColourInc() noexcept
		{
			auto ret = *colourPointer++;

			if (colourPointer == colourView.end())
				colourPointer -= colourView.size();

			return ret;
		}

	private:

		COscilloscope & osc;
		ChannelData::AudioBuffer::ProxyView audioView;
		ChannelData::ColourBuffer::ProxyView colourView;

		AudioIt audioPointer {};
		ColourIt colourPointer {};
	};

	template<typename V>
		void COscilloscope::drawWavePlot(cpl::OpenGLRendering::COpenGLStack & openGLStack)
		{

			typedef SampleColourEvaluator<OscChannels::Left> Evaluator;
			typedef cpl::OpenGLRendering::PrimitiveDrawer<1024> Renderer;

			cpl::OpenGLRendering::MatrixModification matrixMod;
			// and apply the gain:
			const auto gain = (GLfloat)state.envelopeGain;

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


			auto interpolation = state.sampleInterpolation;
			if (pixelsPerSample < 1 && state.sampleInterpolation != SubSampleInterpolation::None)
			{
				interpolation = SubSampleInterpolation::Linear;
			}

			if (state.triggerMode == OscilloscopeContent::TriggeringMode::Window)
			{
				auto const realOffset = std::fmod(state.transportPosition, state.effectiveWindowSize);
				bufferOffset = static_cast<cpl::ssize_t>(std::ceil(realOffset));
				offset = subSampleOffset = 0;
			}
			else
			{
				// TODO: FIX. This causes an extra cycle to be rendered for lanczos so it has more to eat from the edges.
				auto const cycleBuffers = interpolation == SubSampleInterpolation::Lanczos ? 2 : 1;
				// calculate fractionate offsets used for sample-space rendering
				if (state.triggerMode == OscilloscopeContent::TriggeringMode::Spectral)
				{

					quantizedCycleSamples = static_cast<cpl::ssize_t>(std::ceil(triggerState.cycleSamples));
					subSampleOffset = cycleBuffers * (quantizedCycleSamples - triggerState.cycleSamples) + (roundedWindow - state.effectiveWindowSize);
					offset = -triggerState.sampleOffset / sizeMinusOne;
				}
				bufferOffset = roundedWindow + cycleBuffers * quantizedCycleSamples;
				offset += (1 - subSampleOffset) / sizeMinusOne;
			}
			

			roundedWindow = std::max(2, roundedWindow);

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

			const GLfloat endCondition = static_cast<GLfloat>(roundedWindow + quantizedCycleSamples + 2);

			auto renderSampleSpace = [&](auto kernel, GLint primitive, cpl::ssize_t sampleOffset = 0)
			{
				cpl::OpenGLRendering::MatrixModification m;
				// translate triggering offset + 1
				matrixMod.translate(offset - sampleDisplacement, 0, 0);
				// scale to sample/pixels space
				matrixMod.scale(sampleDisplacement, 1, 1);

				Evaluator eval(*this);

				if (!eval.isWellDefined())
					return;
				eval.startFrom(-(bufferOffset + 1 + sampleOffset), -(bufferOffset + sampleOffset));
				cpl::OpenGLRendering::PrimitiveDrawer<1024> drawer(openGLStack, primitive);
				kernel(eval, drawer);
			};

			auto dotSamples = [&] (cpl::ssize_t offset) 
			{
				auto oldPointSize = openGLStack.getPointSize();

				if (pixelsPerSample > 5 && state.sampleInterpolation != SubSampleInterpolation::None)
				{
					openGLStack.setPointSize(oldPointSize * 4);
				}

				renderSampleSpace(
					// nested lambda auto evaluation fails on vc 2015
					[&] (Evaluator & evaluator, Renderer & drawer)
					{
						drawer.addColour(state.colourPrimary);

						for (GLfloat i = 0; i < endCondition; i += 1)
						{
							drawer.addVertex(i, evaluator.evaluateSampleInc(), 0);
						}
					},
					GL_POINTS,
					offset
				);

				openGLStack.setPointSize(oldPointSize);
			};

			// draw dots for very zoomed displays and when there's no subsample interpolation
			if ((state.dotSamples && pixelsPerSample > 5) || state.sampleInterpolation == SubSampleInterpolation::None)
			{
				dotSamples(interpolation == SubSampleInterpolation::Lanczos && state.triggerMode == OscilloscopeContent::TriggeringMode::None && state.timeMode == OscilloscopeContent::TimeMode::Time
					? -(int)(OscilloscopeContent::InterpolationKernelSize) : 0);
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
								Evaluator::ColourT oldColour = evaluator.evaluateColour();

								for (GLfloat i = 0; i < endCondition; i += 1)
								{
									const auto data = evaluator.evaluate();

									drawer.addColour(oldColour);
									drawer.addVertex(i - 1, data.first, 0);
									drawer.addColour(data.second);
									drawer.addVertex(i, data.first, 0);

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
								drawer.addColour(state.colourPrimary);
								for (GLfloat i = 0; i < endCondition; i += 1)
								{
									const auto vertex = evaluator.evaluateSampleInc();
									drawer.addVertex(i - 1, vertex, 0);
									drawer.addVertex(i, vertex, 0);
								}
							},
							GL_LINE_STRIP
						);
					}
					break;
				}
				/*case SubSampleInterpolation::Lanczos:
				{

					auto const KernelSize = OscilloscopeContent::InterpolationKernelSize;
					auto const KernelBufferSize = KernelSize * 2 + 1;

					double samplePos = 0;

					if (state.triggerMode == OscilloscopeContent::TriggeringMode::Window)
					{
						samplePos = std::fmod(state.transportPosition, state.effectiveWindowSize) - 1;
					}
					else
					{
						// TODO: possible small graphic glitch here if the interpolation eats into the next cycle.
						samplePos = triggerState.cycleSamples * 2 + state.effectiveWindowSize - triggerState.sampleOffset;
					}

					// otherwise we will have a discontinuity as the interpolation kernel moves past T = 0
					if (state.triggerMode == OscilloscopeContent::TriggeringMode::None || state.triggerMode == OscilloscopeContent::TriggeringMode::Window)
					{
						samplePos = std::ceil(samplePos);
						if (state.timeMode == OscilloscopeContent::TimeMode::Time)
							samplePos += KernelSize;
					}

					// adjust for left
					double inc = horizontalDelta / (oglc->getRenderingScale() * (getWidth() - 1));
					double unitSpacePos = left;
					double samplesPerPixel = 1.0 / (pixelsPerSample);

					samplePos += -unitSpacePos / inc * samplesPerPixel;
					double currentSample = std::floor(samplePos);

					auto localPointer = (view.begin() + cursor) - static_cast<cpl::ssize_t>(std::floor(samplePos));
					auto clocalPointer = (cview.begin() + cursor + 2 - KernelBufferSize) - static_cast<cpl::ssize_t>(std::floor(samplePos));
					AFloat kernel[KernelBufferSize];

					val_typeof(*clocalPointer) currentColour, nextColour;

					auto get = [&]() {
						auto ret = *localPointer++;
						currentColour = nextColour;
						nextColour = *clocalPointer++;

						if (localPointer == end)
						{
							localPointer -= bufferSamples;
						}
	
						if(clocalPointer == cview.end())
						{
							clocalPointer -= bufferSamples;
						}

						return ret;
					};

					auto insert = [&] (auto val) {
						std::rotate(kernel, kernel + 1, kernel + KernelBufferSize);
						kernel[KernelBufferSize - 1] = val;
					};

					adjustCircular(localPointer, -((int)KernelSize) - 1);
					while (clocalPointer < cview.begin())
						clocalPointer += bufferSamples;

					std::for_each(std::begin(kernel), std::end(kernel), [&](auto & f) { f = get(); });

					{
						cpl::OpenGLRendering::PrimitiveDrawer<1024> drawer(openGLStack, GL_LINE_STRIP);
						if (!state.colourChannelsByFrequency)
							drawer.addColour(state.colourPrimary);
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
								//if(state.colourChannelsByFrequency)
								//	drawer.addColour(currentColour);
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
				}*/

			}
			
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
