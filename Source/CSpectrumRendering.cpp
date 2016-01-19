
#include "CSpectrum.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include "SignalizerDesign.h"
#include <cpl/rendering/OpenGLRasterizers.h>
#include <cpl/simd.h>
#include <cpl/dsp/filterdesign.h>
#include <array>

namespace Signalizer
{

	using namespace cpl;

	void CSpectrum::onGraphicsRendering(juce::Graphics & g)
	{
		// do software rendering
		if (!isOpenGL())
		{
			g.fillAll(state.colourBackground.withAlpha(1.0f));
			g.setColour(state.colourBackground.withAlpha(1.0f).contrasting());
			g.drawText("Enable OpenGL in settings to use the spectrum", getLocalBounds(), juce::Justification::centred);

			// post fps anyway
			auto tickNow = juce::Time::getHighResolutionTicks();
			avgFps.setNext(tickNow - lastFrameTick);
			lastFrameTick = tickNow;
		}
	}

	void CSpectrum::paint2DGraphics(juce::Graphics & g)
	{
		auto cStart = cpl::Misc::ClockCounter();

#pragma message cwarn("lel")

		// ------- draw frequency graph

		char buf[200];

		g.setColour(state.colourGrid.withMultipliedBrightness(0.5f));

		if (state.displayMode == DisplayMode::LineGraph)
		{
			auto complexScale = state.configuration == ChannelConfiguration::Complex ? 2.0f : 1.0f;
			g.setColour(state.colourGrid);
			const auto & divs = frequencyGraph.getDivisions();
			const auto & cdivs = complexFrequencyGraph.getDivisions();
			for (auto & sdiv : divs)
			{
				sprintf_s(buf, "%.2f", sdiv.frequency);
				g.drawText(buf, float(complexScale * sdiv.coord) + 5, 20, 100, 20, juce::Justification::centredLeft);

			}
			if (state.configuration == ChannelConfiguration::Complex)
			{
				auto normalizedScaleX = 1.0 / frequencyGraph.getBounds().dist();
				auto normXC = [=](double in) { return -static_cast<float>(normalizedScaleX * in * 2.0 - 1.0); };

				for (auto & sdiv : cdivs)
				{
					sprintf_s(buf, "-i*%.2f", sdiv.frequency);
					// transform back and forth from unit cartesion... should insert a TODO here.
					g.drawText(buf, getWidth() * (normXC(sdiv.coord) + 1) * 0.5 + 5, 20, 100, 20, juce::Justification::centredLeft);
				}
			}

			for (auto & dbDiv : dbGraph.getDivisions())
			{
				sprintf_s(buf, "%.2f", dbDiv.dbVal);
				g.drawText(buf, 5, float(complexScale * dbDiv.coord), 100, 20, juce::Justification::centredLeft);
			}

		}
		else
		{
			float height = getHeight();
			float baseWidth = getWidth() * 0.05f;
			
			float gradientOffset = 10.0f;

			g.setColour(state.colourGrid);
			const auto & divs = frequencyGraph.getDivisions();


			for (auto & sdiv : divs)
			{
				sprintf_s(buf, "%.2f", sdiv.frequency);
				g.drawText(buf, gradientOffset + baseWidth + 5, float(height - sdiv.coord) - 10 /* height / 2 */, 100, 20, juce::Justification::centredLeft);
			}


			// draw gradient

			juce::ColourGradient gradient;

			// fill in colours

			double gradientPos = 0.0;

			for (int i = 0; i < numSpectrumColours + 1; i++)
			{
				gradientPos += state.normalizedSpecRatios[i];
				gradient.addColour(gradientPos, state.colourSpecs[i].toJuceColour());
			}


			gradient.point1 = {gradientOffset * 0.5f, (float)getHeight() };
			gradient.point2 = {gradientOffset * 0.5f, 0.0f };

			g.setGradientFill(gradient);

			g.fillRect(0.0f, 0.0f, gradientOffset, (float)getHeight());
		}


		if (kdiagnostics.bGetValue() > 0.5)
		{
			char text[1000];
			auto fps = 1.0 / (avgFps.getAverage() / juce::Time::getHighResolutionTicksPerSecond());

			auto totalCycles = renderCycles + cpl::Misc::ClockCounter() - cStart;
			double cpuTime = (double(totalCycles) / (processorSpeed * 1000 * 1000) * 100) * fps;
			double asu = 100 * audioStream.getPerfMeasures().asyncUsage.load(std::memory_order_relaxed);
			double aso = 100 * audioStream.getPerfMeasures().asyncOverhead.load(std::memory_order_relaxed);
			g.setColour(juce::Colours::blue);
			sprintf(text, "%dx%d {%.3f, %.3f}: %.1f fps - %.1f%% cpu, deltaG = %.4f, deltaO = %.4f (rt: %.2f%% - %.2f%%, d: %llu), (as: %.2f%% - %.2f%%)",
				getWidth(), getHeight(), state.viewRect.left, state.viewRect.right,
				fps, cpuTime, graphicsDeltaTime(), openGLDeltaTime(),
				100 * audioStream.getPerfMeasures().rtUsage.load(std::memory_order_relaxed),
				100 * audioStream.getPerfMeasures().rtOverhead.load(std::memory_order_relaxed),
				audioStream.getPerfMeasures().droppedAudioFrames.load(std::memory_order_relaxed),
				asu,
				aso);
			
			g.drawSingleLineText(text, 10, 20);

		}
	}

	void CSpectrum::initOpenGL()
	{
		const int imageSize = 128;
		const float fontToPixelScale = 90 / 64.0f;

		//oglImage.resize(getWidth(), getHeight(), false);
		flags.openGLInitiation = true;
		textures.clear();

	}

	void CSpectrum::closeOpenGL()
	{
		textures.clear();
		oglImage.offload();
	}


	template<std::size_t N, std::size_t V, cpl::GraphicsND::ComponentOrder order>
		void ColourScale2(cpl::GraphicsND::UPixel<order> * __RESTRICT__ outPixel, 
			float intensity, cpl::GraphicsND::UPixel<order> colours[N + 1], float normalizedScales[N])
		{
			intensity = cpl::Math::confineTo(intensity, 0.0f, 1.0f);
			
			float accumulatedSum = 0.0f;

			for (std::size_t i = 0; i < N; ++i)
			{
				accumulatedSum += normalizedScales[i];

				if (accumulatedSum >= intensity)
				{
					std::uint16_t factor = cpl::Math::round<std::uint8_t>(0xFF * cpl::Math::UnityScale::Inv::linear<float>(intensity, accumulatedSum - normalizedScales[i], accumulatedSum));

					std::size_t 
						x1 = std::max(signed(i) - 1, 0), 
						x2 = std::min(x1 + 1, N - 1);

					for (std::size_t p = 0; p < V; ++p)
					{
						outPixel->pixel.data[p] = (((colours[x1].pixel.data[p] * (0x00FF - factor)) + 0x80) >> 8) + (((colours[x2].pixel.data[p] * factor) + 0x80) >> 8);
					}

					return;
				}
			}
		}

	void ColorScale(uint8_t * pixel, float intensity)
	{

		uint8_t red = 0, blue = 0, green = 0;
		// set blue

		if (intensity <= 0)
		{
		}
		else if (intensity < 0.16666f && intensity > 0)
		{
			blue = 6 * intensity * 0x7F;
		}
		else if (intensity < 0.3333f)
		{
			red = 6 * (intensity - 0.16666f) * 0xFF;
			blue = 0x7F - (red >> 1);

		}
		else if (intensity < 0.6666f)
		{
			red = 0xFF;
			green = 3 * (intensity - 0.3333) * 0xFF;
		}
		else if (intensity < 1)
		{
			red = 0xFF;
			green = 0xFF;
			blue = 3 * (intensity - 0.66666) * 0xFF;
		}
		else
			red = green = blue = 0xFF;
		// set green
		pixel[0] = red;
		pixel[1] = green;
		pixel[2] = blue;
		// set red

		// saturate


	}

	void CSpectrum::onOpenGLRendering()
	{
		auto cStart = cpl::Misc::ClockCounter();
		// starting from a clean slate?
		CPL_DEBUGCHECKGL();

		peakFilter.setSampleRate(fpoint(1.0 / openGLDeltaTime()));

		bool lineTransformReady = false;

		// lock the memory buffers, and do our thing.
		{
			handleFlagUpdates();
			// line graph data for ffts are rendered now.
			if (state.displayMode == DisplayMode::LineGraph)
				lineTransformReady = prepareTransform(audioStream.getAudioBufferViews());
		}
		// flags may have altered ogl state
		CPL_DEBUGCHECKGL();



		juce::OpenGLHelpers::clear(state.colourBackground);

		cpl::OpenGLEngine::COpenGLStack openGLStack;
		// set up openGL
		//openGLStack.setBlender(GL_ONE, GL_ONE_MINUS_SRC_COLOR);
		openGLStack.loadIdentityMatrix();
		CPL_DEBUGCHECKGL();

		state.antialias ? openGLStack.enable(GL_MULTISAMPLE) : openGLStack.disable(GL_MULTISAMPLE);
		CPL_DEBUGCHECKGL();

		openGLStack.setLineSize(std::max(0.001f, state.primitiveSize * 10));
		openGLStack.setPointSize(std::max(0.001f, state.primitiveSize * 10));
		CPL_DEBUGCHECKGL();

		CPL_DEBUGCHECKGL();

		/// <summary>
		/// 
		/// </summary>
		/// 
		///



		switch (state.displayMode)
		{
		case DisplayMode::LineGraph:
			// no need to lock in this case, as this display mode is exclusively switched,
			// such that only we have access to it.
			if (lineTransformReady)
			{
				doTransform();
				mapToLinearSpace();
				postProcessStdTransform();
			}
			renderLineGraph<Types::v8sf>(openGLStack); break;
		case DisplayMode::ColourSpectrum:
			// mapping and processing is already done here.
			renderColourSpectrum<Types::v8sf>(openGLStack); break;

		}




		renderCycles = cpl::Misc::ClockCounter() - cStart;
		auto tickNow = juce::Time::getHighResolutionTicks();
		avgFps.setNext(tickNow - lastFrameTick);
		lastFrameTick = tickNow;
		CPL_DEBUGCHECKGL();
		renderGraphics([&](juce::Graphics & g) { paint2DGraphics(g); });
		CPL_DEBUGCHECKGL();
	}


	template<typename V>
		void CSpectrum::renderColourSpectrum(cpl::OpenGLEngine::COpenGLStack & ogs)
		{
			CPL_DEBUGCHECKGL();

			if (!state.isFrozen)
			{


				framePixelPosition %= (getWidth());
				double bufferSmoothing = state.bufferSmoothing;
				auto approximateFrames = getApproximateStoredFrames();
				/*if (approximateFrames == 0)
					approximateFrames = framesPerUpdate;*/
				std::size_t processedFrames = 0;
				framesPerUpdate = approximateFrames + bufferSmoothing * (framesPerUpdate - approximateFrames);
				auto framesThisTime = cpl::Math::round<std::size_t>(framesPerUpdate);

				// if there's no buffer smoothing at all, we just capture every frame possible.
				// 
				bool shouldCap = state.bufferSmoothing != 0.0;

				while ((!shouldCap || (processedFrames++ < framesThisTime)) && processNextSpectrumFrame())
				{
#pragma message cwarn("Update frames per update each time inside here, but as a local variable! There may come more updates meanwhile.")
					// run the next frame through pixel filters and format it etc.

					for (int i = 0; i < getAxisPoints(); ++i)
					{
						ColourScale2<numSpectrumColours + 1, 4>(columnUpdate.data() + i, filterResults[i].magnitude, state.colourSpecs, state.normalizedSpecRatios);
						//ColorScale(&columnUpdate[i * 4], filterResults[i].magnitude);
					}
					//CPL_DEBUGCHECKGL();
					oglImage.updateSingleColumn(framePixelPosition, columnUpdate, GL_RGBA);
					//CPL_DEBUGCHECKGL();

					framePixelPosition++;
					framePixelPosition %= (getWidth());
				}
			}

			CPL_DEBUGCHECKGL();

			{
				cpl::OpenGLEngine::COpenGLImage::OpenGLImageDrawer imageDrawer(oglImage, ogs);

				imageDrawer.drawCircular((float)((double)(framePixelPosition) / (getWidth() - 1)));

			}

			CPL_DEBUGCHECKGL();

			// render grid
			{
				auto normalizedScale = 1.0 / getHeight();

				// draw vertical lines.
				const auto & lines = frequencyGraph.getLines();

				auto norm = [=](double in) { return static_cast<float>(normalizedScale * in * 2.0 - 1.0); };

				float baseWidth = 0.1f;

				float gradientOffset = 10.0f / getWidth() - 1.0f;

				OpenGLEngine::PrimitiveDrawer<128> lineDrawer(ogs, GL_LINES);

				lineDrawer.addColour(state.colourGrid.withMultipliedBrightness(0.5f));

				for (auto dline : lines)
				{
					auto line = norm(dline);
					lineDrawer.addVertex(gradientOffset, line, 0.0f);
					lineDrawer.addVertex(gradientOffset + baseWidth * 0.7f, line, 0.0f);
				}

				lineDrawer.addColour(state.colourGrid);
				const auto & divs = frequencyGraph.getDivisions();

				for (auto & sdiv : divs)
				{
					auto line = norm(sdiv.coord);
					lineDrawer.addVertex(gradientOffset, line, 0.0f);
					lineDrawer.addVertex(gradientOffset + baseWidth, line, 0.0f);
				}

			}
			CPL_DEBUGCHECKGL();
		}


	template<typename V>
	void CSpectrum::renderLineGraph(cpl::OpenGLEngine::COpenGLStack & ogs)
	{
		ogs.setBlender(GL_ONE, GL_ONE_MINUS_SRC_COLOR);
		int points = getAxisPoints() - 1;
		switch (state.configuration)
		{
		case ChannelConfiguration::MidSide:
		case ChannelConfiguration::Phase:
		case ChannelConfiguration::Separate:
		{
			OpenGLEngine::PrimitiveDrawer<256> lineDrawer(ogs, GL_LINE_STRIP);
			lineDrawer.addColour(state.colourTwo);
			for (int i = 0; i < (points + 1); ++i)
			{
				lineDrawer.addVertex((float(i) / points) * 2 - 1, filterResults[i].rightMagnitude * 2 - 1, -0.5);
			}
		}
		// (fall-through intentional)
		case ChannelConfiguration::Left:
		case ChannelConfiguration::Right:
		case ChannelConfiguration::Merge:
		case ChannelConfiguration::Side:
		case ChannelConfiguration::Complex:
		{
			OpenGLEngine::PrimitiveDrawer<256> lineDrawer(ogs, GL_LINE_STRIP);
			lineDrawer.addColour(state.colourOne);
			for (int i = 0; i < (points + 1); ++i)
			{
				lineDrawer.addVertex((float(i) / points) * 2 - 1, filterResults[i].leftMagnitude * 2 - 1, 0);
			}
		}
		default:
			break;
		}
		// render grid
		{
			OpenGLEngine::PrimitiveDrawer<128> lineDrawer(ogs, GL_LINES);

			lineDrawer.addColour(state.colourGrid.withMultipliedBrightness(0.5f));

			auto xDist = frequencyGraph.getBounds().dist();
			auto normalizedScaleX = 1.0 / xDist;
			auto normalizedScaleY = 1.0 / getHeight();
			// draw vertical lines.
			const auto & lines = frequencyGraph.getLines();
			const auto & clines = complexFrequencyGraph.getLines();
			// TODO: fix using a matrix modification instead (cheaper)
			auto normX = [=](double in) { return static_cast<float>(normalizedScaleX * in * 2.0 - 1.0); };
			auto normXC = [=](double in) { return -static_cast<float>(normalizedScaleX * in * 2.0 - 1.0); };
			//auto normXC = normX;
			auto normY = [=](double in) {  return static_cast<float>(1.0 - normalizedScaleY * in * 2.0); };



			for (auto dline : lines)
			{
				auto line = normX(dline);
				lineDrawer.addVertex(line, -1.0f, 0.0f);
				lineDrawer.addVertex(line, 1.0f, 0.0f);
			}

			for (auto dline : clines)
			{
				auto line = normXC(dline);
				lineDrawer.addVertex(line, -1.0f, 0.0f);
				lineDrawer.addVertex(line, 1.0f, 0.0f);
			}

			//m.scale(1, getHeight(), 1);
			lineDrawer.addColour(state.colourGrid);
			const auto & divs = frequencyGraph.getDivisions();
			const auto & cdivs = complexFrequencyGraph.getDivisions();

			for (auto & sdiv : divs)
			{
				auto line = normX(sdiv.coord);
				lineDrawer.addVertex(line, -1.0f, 0.0f);
				lineDrawer.addVertex(line, 1.0f, 0.0f);
			}

			for (auto & sdiv : cdivs)
			{
				auto line = normXC(sdiv.coord);
				lineDrawer.addVertex(line, -1.0f, 0.0f);
				lineDrawer.addVertex(line, 1.0f, 0.0f);
			}

			// draw horizontal lines:
			for (auto & dbDiv : dbGraph.getDivisions())
			{
				auto line = normY(dbDiv.coord);
				lineDrawer.addVertex(-1.0f, line, 0.0f);
				lineDrawer.addVertex(1.0f, line, 0.0f);
			}
		}



	}

};
