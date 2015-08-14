
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
		auto cStart = cpl::Misc::ClockCounter();

		// do software rendering
		if(!isOpenGL())
		{
			g.fillAll(state.colourBackground.withAlpha(1.0f));
			g.setColour(state.colourBackground.withAlpha(1.0f).contrasting());
			g.drawText("Enable OpenGL in settings to use the spectrum", getLocalBounds(), juce::Justification::centred);
			
			// post fps anyway
			auto tickNow = juce::Time::getHighResolutionTicks();
			avgFps.setNext(tickNow - lastFrameTick);
			lastFrameTick = tickNow;
			
		}
#pragma message cwarn("lel")

		// ------- draw frequency graph


		const auto & lines = frequencyGraph.getLines();

		char buf[200];

		g.setColour(state.colourGrid.withMultipliedBrightness(0.5f));

		if (state.displayMode == DisplayMode::LineGraph)
		{
			// draw vertical lines.
			for (auto line : lines)
			{
				g.drawLine({ (float)line, 0.0f, (float)line, (float)getHeight() }, 1.0f);
			}

			g.setColour(state.colourGrid);
			const auto & divs = frequencyGraph.getDivisions();

			for (auto & sdiv : divs)
			{
				sprintf_s(buf, "%.2f", sdiv.frequency);
				g.drawLine({ (float)sdiv.coord, 0.0f, (float)sdiv.coord, (float)getHeight() }, 1.0f);
				g.drawText(buf, (float)sdiv.coord + 5, 20, 100, 20, juce::Justification::centredLeft);
			}

			// draw horizontal lines:

			for (auto & dbDiv : dbGraph.getDivisions())
			{
				sprintf_s(buf, "%.2f", dbDiv.dbVal);
				g.drawLine({ 30, (float)dbDiv.coord, float(getWidth()), (float)dbDiv.coord }, 1.0f);
				g.drawText(buf, 5, (float)dbDiv.coord, 100, 20, juce::Justification::centredLeft);
			}

		}
		else
		{
			float height = getHeight();
			float baseWidth = getWidth() * 0.05f;
			
			float gradientOffset = 10.0f;



			// draw horizontal lines.

			for (auto line : lines)
			{
				g.drawLine({ gradientOffset, float(height - line), gradientOffset + baseWidth * 0.7f, float(height - line) }, 1.0f);
			}

			g.setColour(state.colourGrid);
			const auto & divs = frequencyGraph.getDivisions();


			for (auto & sdiv : divs)
			{
				sprintf_s(buf, "%.2f", sdiv.frequency);
				g.drawLine({ gradientOffset, float(height - sdiv.coord), gradientOffset + baseWidth, float(height - sdiv.coord) }, 1.0f);
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


		if (true || kdiagnostics.bGetValue() > 0.5)
		{
			char text[1000];
			auto fps = 1.0 / (avgFps.getAverage() / juce::Time::getHighResolutionTicksPerSecond());

			auto totalCycles = renderCycles + cpl::Misc::ClockCounter() - cStart;
			double cpuTime = (double(totalCycles) / (processorSpeed * 1000 * 1000) * 100) * fps;
			g.setColour(juce::Colours::blue);
			sprintf(text, "%dx%d: %.1f fps - %.1f%% cpu (audioThread: %.1f%% cpu, d: %d, fpu: %f), deltaG = %f, deltaO = %f, viewRect = {%f, %f}",
				getWidth(), getHeight(), fps, cpuTime, audioThreadUsage, droppedAudioFrames, framesPerUpdate, graphicsDeltaTime(), openGLDeltaTime(), state.viewRect.left, state.viewRect.right);
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

					std::size_t x1 = std::max(signed(i) - 1, 0), x2 = std::min(x1 + 1, N - 1);

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
		// starting from a clean slate?
		CPL_DEBUGCHECKGL();
		if (audioStream.empty())
			return;


		handleFlagUpdates();
		// flags may have altered ogl state
		CPL_DEBUGCHECKGL();


		auto cStart = cpl::Misc::ClockCounter();
		juce::OpenGLHelpers::clear(state.colourBackground);

		cpl::AudioBuffer * buffer;

		if (state.isFrozen)
		{
			buffer = &audioStreamCopy;
		}
		else if (isSynced)
		{
			std::vector<cpl::CMutex> locks;

			for (auto & buffer : audioStream)
				locks.emplace_back(buffer);

			for (unsigned i = 0; i < audioStream.size(); ++i)
			{
				audioStream[i].clone(audioStreamCopy[i]);

			}
			buffer = &audioStreamCopy;
		}
		else
		{
			buffer = &audioStream;

		}
		const std::size_t numSamples = (*buffer)[0].size - 1;

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
	
		prepareTransform();

		doTransform();

		mapToLinearSpace();

		switch (state.displayMode)
		{
		case DisplayMode::ColourSpectrum:
			renderColourSpectrum<Types::v8sf>(openGLStack); break;
		case DisplayMode::LineGraph:
			renderLineGraph<Types::v8sf>(openGLStack); break;
		}




		renderCycles = cpl::Misc::ClockCounter() - cStart;
		auto tickNow = juce::Time::getHighResolutionTicks();
		avgFps.setNext(tickNow - lastFrameTick);
		lastFrameTick = tickNow;
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
				int processedFrames = 0;
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
			ogs.enable(GL_TEXTURE_2D);

			cpl::OpenGLEngine::COpenGLImage::OpenGLImageDrawer imageDrawer(oglImage, ogs);
			CPL_DEBUGCHECKGL();
			/*framePixelPosition--;
			if (framePixelPosition < 0)
			{
				framePixelPosition = getWidth();
			}*/

			imageDrawer.drawCircular((float)((double)(framePixelPosition) / (getWidth() - 1)));

			/*framePixelPosition++;
			framePixelPosition %= (getWidth() + 1);*/
		}


	template<typename V>
	void CSpectrum::renderLineGraph(cpl::OpenGLEngine::COpenGLStack & ogs)
	{

		int points = getAxisPoints() - 1;
		switch (state.configuration)
		{
		case ChannelConfiguration::Phase:
		case ChannelConfiguration::Separate:
		{
			OpenGLEngine::PrimitiveDrawer<256> lineDrawer(ogs, GL_LINE_STRIP);
			lineDrawer.addColour(state.colourTwo);
			for (int i = 0; i < (points + 1); ++i)
			{
				lineDrawer.addVertex((float(i) / points) * 2 - 1, filterResults[i].rightMagnitude * 2 - 1, 0);
			}
		}
		// (fall-through intentional)
		case ChannelConfiguration::Left:
		case ChannelConfiguration::Right:
		case ChannelConfiguration::Merge:
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

	}

};
