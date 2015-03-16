
#include "CVectorScope.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include "SignalizerDesign.h"
#include <cpl/rendering/OpenGLRasterizers.h>
#include <cpl/simd.h>

namespace Signalizer
{
	// swapping the right channel might give an more intuitive view

	//static const char * ChannelDescriptions[] = { "+R", "+L", "-R", "-L" };
	
	static const char * ChannelDescriptions[] = { "+L", "+R", "-L", "-R" };
	static std::vector<std::string> OperationalModeNames = {"Lissajous", "Polar"};
	
	enum class OperationalModes
	{
		Lissajous,
		Polar
		
	};
	
	enum Textures
	{
		LPlus,
		RPlus,
		LMinus,
		RMinus
	};





	void CVectorScope::paint(juce::Graphics & g)
	{

		
		auto cStart = cpl::Misc::ClockCounter();
		

		
		// do software rendering
		if(!isOpenGL())
		{
			g.fillAll(kbackgroundColour.getControlColourAsColour().withAlpha(1.0f));
			g.setColour(kbackgroundColour.getControlColourAsColour().withAlpha(1.0f).contrasting());
			g.drawText("Enable OpenGL in settings to use the vectorscope", getLocalBounds(), juce::Justification::centred);
			
			// post fps
			auto tickNow = juce::Time::getHighResolutionTicks();
			avgFps.setNext(tickNow - lastFrameTick);
			lastFrameTick = tickNow;
			
		}
		
		if (kdiagnostics.bGetValue() > 0.5)
		{
			auto fps = 1.0 / (avgFps.getAverage() / juce::Time::getHighResolutionTicksPerSecond());

			auto totalCycles = renderCycles + cpl::Misc::ClockCounter() - cStart;
			double cpuTime = (double(totalCycles) / (processorSpeed * 1000 * 1000) * 100) * fps;
			g.setColour(juce::Colours::blue);
			sprintf(textbuf.get(), "%dx%d: %.1f fps - %.1f%% cpu", 
				getWidth(), getHeight(), fps, cpuTime);
			g.drawSingleLineText(textbuf.get(), 10, 20);
			
		}
	}

	void CVectorScope::initOpenGL()
	{
		const int imageSize = 128;
		const float fontToPixelScale = 90 / 64.0f;
		
		textures.clear();

		for (int i = 0; i < ArraySize(ChannelDescriptions); ++i)
		{

			juce::Image letter(juce::Image::ARGB, imageSize, imageSize, true);
			{
				juce::Graphics g(letter);
				g.fillAll(juce::Colours::transparentBlack);
				g.setColour(Colours::white);
				g.setFont(letter.getHeight() * fontToPixelScale * 0.5);
				g.drawText(ChannelDescriptions[i], letter.getBounds().toFloat(), juce::Justification::centred, false);
			}
			textures.push_back(std::unique_ptr<juce::OpenGLTexture>((new juce::OpenGLTexture())));
			textures[i]->loadImage(letter);
		}
	}
	void CVectorScope::closeOpenGL()
	{
		textures.clear();
	}

	void CVectorScope::renderOpenGL()
	{
		if (audioStream.empty())
			return;

		auto cStart = cpl::Misc::ClockCounter();
		OpenGLHelpers::clear(kbackgroundColour.getControlColourAsColour());

		const bool fillPath = kdrawLines.bGetValue() > 0.5;
		const bool fadeHistory = kfadeOld.bGetValue() > 0.5;
		const bool antiAlias = kantiAlias.bGetValue() > 0.5;
		const float psize = static_cast<GLfloat>(kprimitiveSize.bGetValue() * 10);
		const float gain = getGain();
		const bool isPolar = kopMode.getZeroBasedSelIndex() == (int)OperationalModes::Polar;

		float sleft(0.0f), sright(0.0f);

		cpl::AudioBuffer * buffer;

		if (isFrozen)
		{

			buffer = &audioStreamCopy;
		}
		else if (isSynced)
		{
			std::vector<cpl::CMutex> locks;

			for (auto & buffer : audioStream)
				locks.emplace_back(buffer);

			for (auto i = 0; i < audioStream.size(); ++i)
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


		openGLStack.applyTransform3D(ktransform.getTransform3D());
		antiAlias ? openGLStack.enable(GL_MULTISAMPLE) : openGLStack.disable(GL_MULTISAMPLE);

		cpl::CAudioBuffer::CChannelBuffer * leftChannel, * rightChannel;
		leftChannel = &((*buffer)[0]);
		rightChannel = &((*buffer)[1]);

		openGLStack.setLineSize(psize);
		openGLStack.setPointSize(psize);

		// add a stack scope for transformations.
		if (isPolar)
		{
			using namespace cpl::simd;
			cpl::OpenGLEngine::MatrixModification matrixMod;
			// fixate the rotation to fit.
			//matrixMod.rotate(135, 0, 0, 1);
			//matrixMod.rotate(krotation.bGetValue() * 360, 0, 0, 1);
			// and apply the gain:
			matrixMod.scale(gain, gain, 1);

			typedef v4sf V;
			cpl::OpenGLEngine::PrimitiveDrawer<1024> drawer(openGLStack, fillPath ? GL_LINE_STRIP : GL_POINTS);

			drawer.addColour(kdrawingColour.getControlColourAsColour());

			static const v4sf vCosine = set1<v4sf>((float)std::cos(M_PI * 135.0 / 180));
			static const v4sf vSine = set1<v4sf>((float)std::sin(M_PI * 135.0 / 180));
			const v4sf vZero = zero<v4sf>();
			const v4sf vSign = consts<v4sf>::sign_bit;
			float sampleFade = 1.0 / numSamples;

			suitable_container<v4sf> outX, outY;

			// iterate front of buffer, then back.
			// this feels awkward, but that's kinda how circular 
			// buffers work.
			std::size_t depthCounter = 0;
			for (int section = 0; section < 2; ++section)
			{
				auto const lit = leftChannel->getIterator<4>();
				auto const rit = rightChannel->getIterator<4>();

				auto const sectionSamples = lit.sizes[section];

				for (std::size_t i = 0; i < sectionSamples; i += elements_of<v4sf>::value)
				{
					v4sf vLeft = loadu<v4sf>(lit.getIndex(section) + i);
					v4sf vRight = loadu<v4sf>(rit.getIndex(section) + i);

					// the length of the hypotenuse of the triangle, we 
					// convert the unit square to.
					auto const length = max(abs(vLeft), abs(vRight));

					// rotate our view manually (to center on Y-axis)
					v4sf vY = vLeft * vCosine - vRight * vSine;
					v4sf vX = vLeft * vSine + vRight * vCosine;

					// check for any zero elements.
					V leftIsZero = vLeft == vZero;
					V rightIsZero = vRight == vZero;
					auto mask = vand(leftIsZero, rightIsZero);
					mask = vnot(mask);
					// get the phase angle. use atan2 if you want to draw the full circle.
					auto angle = atan(vX / vY);
					// replace nan elements of angle with zero
					angle = vand(mask, angle);
					// calcuate x,y coordinates for the right triangle
					sincos(angle, &vX, &vY);

					// construct triangle.
					outX = vX * length;
					outY = vY * length;

					// draw vertices
					drawer.addVertex(outX[0], outY[0], depthCounter++ * sampleFade - 1);
					drawer.addVertex(outX[1], outY[1], depthCounter++ * sampleFade - 1);
					drawer.addVertex(outX[2], outY[2], depthCounter++ * sampleFade - 1);
					drawer.addVertex(outX[3], outY[3], depthCounter++ * sampleFade - 1);
				}
			}

		}
		else
		{
			cpl::OpenGLEngine::MatrixModification matrixMod;
			// apply the custom rotation to the waveform
			matrixMod.rotate(krotation.bGetValue() * 360, 0, 0, 1);
			// and apply the gain:
			matrixMod.scale(gain, gain, 1);

			float sampleFade = 1.0 / numSamples;

			if (!fadeHistory)
			{
				cpl::OpenGLEngine::PrimitiveDrawer<1024> drawer(openGLStack, fillPath ? GL_LINE_STRIP : GL_POINTS);

				drawer.addColour(kdrawingColour.getControlColourAsColour());
				for (std::size_t i = 0; i < numSamples; ++i)
				{
					// use glDrawArrays instead here
					sleft = leftChannel->singleCheckAccess(i);
					sright = rightChannel->singleCheckAccess(i);

					drawer.addVertex(sright, sleft, i * sampleFade - 1);
				}
			}
			else
			{
				cpl::OpenGLEngine::PrimitiveDrawer<1024> drawer(openGLStack, fillPath ? GL_LINE_STRIP : GL_POINTS);
				
				auto jcolour = kdrawingColour.getControlColourAsColour();
				float red = jcolour.getFloatRed(), blue = jcolour.getFloatBlue(),
					green = jcolour.getFloatGreen(), alpha = jcolour.getFloatGreen();

				float fade = 0;

				for (std::size_t i = 0; i < numSamples; ++i)
				{
					sleft = leftChannel->singleCheckAccess(i);
					sright = rightChannel->singleCheckAccess(i);

					fade = i * sampleFade;

					drawer.addColour(fade * red, fade * green, fade * blue, alpha);
					drawer.addVertex(sright, sleft, i * sampleFade - 1);
				
				}
			}
		}

		openGLStack.setLineSize(2.0f);

		// draw skeleton graph
		if (!isPolar)
		{
			cpl::OpenGLEngine::PrimitiveDrawer<128> drawer(openGLStack, GL_LINES);
			drawer.addColour(kskeletonColour.getControlColourAsColour());
			int nlines = 14;
			auto rel = 1.0f / nlines;

			// front vertival
			for (int i = 0; i <= nlines; ++i)
			{
				drawer.addVertex(i * rel * 2 - 1, -1.f, 0.0f);
				drawer.addVertex(i * rel * 2 - 1, 1.f, 0.0f);
			}
			// front horizontal
			for (int i = 0; i <= nlines; ++i)
			{
				drawer.addVertex(-1.0f, i * rel * 2 - 1, 0.0f);
				drawer.addVertex(1.0f, i * rel * 2 - 1, 0.0f);
			}
			// back vertical
			for (int i = 0; i <= nlines; ++i)
			{
				drawer.addVertex(i * rel * 2 - 1, -1.f, -1.0f);
				drawer.addVertex(i * rel * 2 - 1, 1.f, -1.0f);
			}
			// back horizontal
			for (int i = 0; i <= nlines; ++i)
			{
				drawer.addVertex(-1.0f, i * rel * 2 - 1, -1.0f);
				drawer.addVertex(1.0f, i * rel * 2 - 1, -1.0f);
			}
		}
		else
		{
			// draw two half circles
			auto lut = circleData.get();
			
			int numInt = lut->tableSize;
			float advance = 1.0f / (numInt - 1);
			{
				cpl::OpenGLEngine::PrimitiveDrawer<512> drawer(openGLStack, GL_LINES);
				drawer.addColour(kskeletonColour.getControlColourAsColour());

				float oldY = 0.0f;
				for (int i = 1; i < numInt; ++i)
				{
					auto fraction = advance * i;
					auto yCoordinate = lut->linearLookup(fraction);
					auto leftX = -1.0f + fraction;
					auto rightX = 1.0f - fraction;
					// left part
					drawer.addVertex(leftX - advance, oldY, 0);
					drawer.addVertex(leftX, yCoordinate, 0);
					drawer.addVertex(leftX - advance, oldY, -1);
					drawer.addVertex(leftX, yCoordinate, -1);
					
					// right part
					drawer.addVertex(rightX + advance, oldY, 0);
					drawer.addVertex(rightX, yCoordinate, 0);
					drawer.addVertex(rightX + advance, oldY, -1);
					drawer.addVertex(rightX, yCoordinate, -1);
					//drawer.addVertex(1 - fraction, yCoordinate, 0);
					
					oldY = yCoordinate;
				}
				
				// add front and back horizontal lines.
				drawer.addVertex(-1.0f, 0.0f, 0.0f);
				drawer.addVertex(1.0f, 0.0f, 0.0f);
				drawer.addVertex(-1.0f, 0.0f, -1.0f);
				drawer.addVertex(1.0f, 0.0f, -1.0f);
				
				// add critical diagonal phase lines.
				const float quarterPISinCos = 0.707106781186547f;
				drawer.addVertex(0.0f, 0.0f, 0.0f);
				drawer.addVertex(quarterPISinCos, quarterPISinCos, 0.0f);
				drawer.addVertex(0.0f, 0.0f, 0.0f);
				drawer.addVertex(-quarterPISinCos, quarterPISinCos, 0.0f);
				
				drawer.addVertex(0.0f, 0.0f, -1.0f);
				drawer.addVertex(quarterPISinCos, quarterPISinCos, -1.0f);
				drawer.addVertex(0.0f, 0.0f, -1.0f);
				drawer.addVertex(-quarterPISinCos, quarterPISinCos, -1.0f);;
			}
		}

		// Draw basic graph

		{
			cpl::OpenGLEngine::PrimitiveDrawer<12> drawer(openGLStack, GL_LINES);
			drawer.addColour(kgraphColour.getControlColourAsColour());
			// front x, y axii
			drawer.addVertex(-1.0f, 0.0f, 0.0f);
			drawer.addVertex(1.0f, 0.0f, 0.0f);
			drawer.addVertex(0.0f, 1.0f, 0.0f);
			drawer.addVertex(0.0f, isPolar ? 0.0f : -1.0f, 0.0f);
			

		}


		// draw channel rotations letters.
		{
			auto rotation = -krotation.bGetValue() * 2 * M_PI;
			const float heightToWidthFactor = float(getHeight()) / getWidth();
			using namespace cpl::simd;

			cpl::OpenGLEngine::MatrixModification m;
			// this undoes the text squashing due to variable aspect ratios.
			m.scale(heightToWidthFactor, 1.0f, 1.0f);

			// calculate coordinates using sin/cos simd pairs.
			suitable_container<v4sf> phases, xcoords, ycoords;

			// set phases (we rotate L/R etc. text around in a circle)
			phases[0] = rotation;
			phases[1] = M_PI * 0.5f + rotation;
			phases[2] = M_PI + rotation;
			phases[3] = M_PI * 1.5f + rotation;

			const float circleScaleFactor = 1.1f;

			// some registers
			v4sf 
				vsines, 
				vcosines, 
				vscale = set1<v4sf>(circleScaleFactor), 
				vadd = set1<v4sf>(1.0f - circleScaleFactor),
				vheightToWidthFactor = set1<v4sf>(1.0f / heightToWidthFactor);

			// do 8 trig functions in one go!
			cpl::sse::sincos_ps(phases, &vsines, &vcosines);

			// place the circle just outside the graph, and offset it.
			xcoords = vsines * vscale * vheightToWidthFactor + vadd;
			ycoords = vcosines * vscale + vadd;

			auto jcolour = kgraphColour.getControlColourAsColour();
			// render texture text at coordinates.
			for (int i = 0; i < 4; ++i)
			{
				cpl::OpenGLEngine::ImageDrawer text(openGLStack, *textures[i]);
				text.setColour(jcolour);
				text.drawAt({ xcoords[i], ycoords[i], 0.1f, 0.1f });
			}
		}

		renderCycles = cpl::Misc::ClockCounter() - cStart;
		auto tickNow = juce::Time::getHighResolutionTicks();
		avgFps.setNext(tickNow - lastFrameTick);
		lastFrameTick = tickNow;
	}




};
