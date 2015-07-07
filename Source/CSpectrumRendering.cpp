
#include "CSpectrum.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include "SignalizerDesign.h"
#include <cpl/rendering/OpenGLRasterizers.h>
#include <cpl/simd.h>
#include <cpl/dsp/filterdesign.h>

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
#pragma cwarn("lel")

		// ------- draw frequency graph

		g.setColour(juce::Colours::grey);
		const auto & lines = frequencyGraph.getLines();
		for (auto line : lines)
		{
			g.drawLine({ (float)line, 0.0f, (float)line, (float)getHeight() }, 1.0f);
		}

		g.setColour(juce::Colours::white);
		const auto & divs = frequencyGraph.getDivisions();

		char buf[200];

		for (auto & sdiv : divs)
		{
			sprintf_s(buf, "%.2f", sdiv.frequency);
			g.drawLine({ (float)sdiv.coord, 0.0f, (float)sdiv.coord, (float)getHeight() }, 1.0f);
			g.drawText(buf, (float)sdiv.coord + 5, 20, 100, 20, juce::Justification::centredLeft);
		}


		if (true || kdiagnostics.bGetValue() > 0.5)
		{
			char text[1000];
			auto fps = 1.0 / (avgFps.getAverage() / juce::Time::getHighResolutionTicksPerSecond());

			auto totalCycles = renderCycles + cpl::Misc::ClockCounter() - cStart;
			double cpuTime = (double(totalCycles) / (processorSpeed * 1000 * 1000) * 100) * fps;
			g.setColour(juce::Colours::blue);
			sprintf(text, "%dx%d: %.1f fps - %.1f%% cpu, deltaG = %f, deltaO = %f, viewRect = {%f, %f}",
				getWidth(), getHeight(), fps, cpuTime, graphicsDeltaTime(), openGLDeltaTime(), state.viewRect.left, state.viewRect.right);
			g.drawSingleLineText(text, 10, 20);

		}
	}

	void CSpectrum::initOpenGL()
	{
		const int imageSize = 128;
		const float fontToPixelScale = 90 / 64.0f;

		oglImage.resize(getWidth(), getHeight(), false);

		textures.clear();

	}

	void CSpectrum::closeOpenGL()
	{
		textures.clear();
		oglImage.offload();
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
		if (audioStream.empty())
			return;

		handleFlagUpdates();

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
		openGLStack.setBlender(GL_ONE, GL_ONE_MINUS_SRC_COLOR);
		openGLStack.loadIdentityMatrix();

		state.antialias ? openGLStack.enable(GL_MULTISAMPLE) : openGLStack.disable(GL_MULTISAMPLE);


		openGLStack.setLineSize(state.primitiveSize * 10);
		openGLStack.setPointSize(state.primitiveSize * 10);

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

			for (int i = 0; i < getAxisPoints(); ++i)
			{
				ColorScale(&columnUpdate[i * 3], filterResults[i].magnitude);
			}
			//CPL_DEBUGCHECKGL();
			oglImage.updateSingleColumn(framePixelPosition, columnUpdate);
			//CPL_DEBUGCHECKGL();


			cpl::OpenGLEngine::COpenGLImage::OpenGLImageDrawer imageDrawer(oglImage, ogs);



			imageDrawer.drawCircular((float)((double)(framePixelPosition) / (getWidth() - 1)));

			framePixelPosition++;
			framePixelPosition %= (getWidth() + 1);
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
