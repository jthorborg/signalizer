
#include "CVectorScope.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include "SignalizerDesign.h"
#include <cpl/rendering/OpenGLRasterizers.h>
#include <cpl/simd.h>
#include <cpl/dsp/filterdesign.h>
namespace Signalizer
{
	// swapping the right channel might give an more intuitive view
	
	static const char * ChannelDescriptions[] = { "+L", "+R", "-L", "-R", "L", "R", "C"};
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
		RMinus,
		Left,
		Right,
		Center
	};


	void CVectorScope::paint(juce::Graphics & g)
	{

		auto cStart = cpl::Misc::ClockCounter();

		if (normalizeGain && isEditorOpen())
		{
			setGainAsFraction(envelopeGain);
		}

		// do software rendering
		if(!isOpenGL())
		{
			g.fillAll(kbackgroundColour.getControlColourAsColour().withAlpha(1.0f));
			g.setColour(kbackgroundColour.getControlColourAsColour().withAlpha(1.0f).contrasting());
			g.drawText("Enable OpenGL in settings to use the vectorscope", getLocalBounds(), juce::Justification::centred);
			
			// post fps anyway
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
		juce::OpenGLHelpers::clear(kbackgroundColour.getControlColourAsColour());
		const float quarterPISinCos = 0.707106781186547f;
		const float circleScaleFactor = 1.1f;
		const bool antiAlias = kantiAlias.bGetValue() > 0.5;
		const bool isPolar = kopMode.getZeroBasedSelIndex() == (int)OperationalModes::Polar;
		GLfloat gain = getGain();

		const float psize = static_cast<GLfloat>(kprimitiveSize.bGetValue() * 10);
		const bool fillPath = kdrawLines.bGetValue() > 0.5;
		const bool fadeHistory = kfadeOld.bGetValue() > 0.5;
		float sleft, sright;
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
			drawPolarPlot<cpl::simd::v8sf>(openGLStack, *buffer);
		}
		else // is Lissajous
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
		if(!isPolar)
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

			// some registers
			v4sf 
				vsines, 
				vcosines, 
				vscale = set1<v4sf>(circleScaleFactor), 
				vadd = set1<v4sf>(1.0f - circleScaleFactor),
				vheightToWidthFactor = set1<v4sf>(1.0f / heightToWidthFactor);

			// do 8 trig functions in one go!
			cpl::simd::sincos(phases.toType(), &vsines, &vcosines);

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
		else
		{
			const float heightToWidthFactor = float(getHeight()) / getWidth();
			using namespace cpl::simd;

			cpl::OpenGLEngine::MatrixModification m;
			// this undoes the text squashing due to variable aspect ratios.
			m.scale(heightToWidthFactor, 1.0f, 1.0f);

			auto jcolour = kgraphColour.getControlColourAsColour();
			float nadd = 1.0f - circleScaleFactor;
			float xcoord = quarterPISinCos * circleScaleFactor / heightToWidthFactor + nadd;
			float ycoord = quarterPISinCos * circleScaleFactor + nadd;
			
			// render texture text at coordinates.
			{
				cpl::OpenGLEngine::ImageDrawer text(openGLStack, *textures[Textures::Left]);
				text.setColour(jcolour);
				text.drawAt({ -xcoord + nadd, ycoord, 0.1f, 0.1f });
			}
			{
				cpl::OpenGLEngine::ImageDrawer text(openGLStack, *textures[Textures::Center]);
				text.setColour(jcolour);
				text.drawAt({ 0 + nadd * 0.5f, 1, 0.1f, 0.1f });
			}
			{
				cpl::OpenGLEngine::ImageDrawer text(openGLStack, *textures[Textures::Right]);
				text.setColour(jcolour);
				text.drawAt({ xcoord, ycoord, 0.1f, 0.1f });
			}
			
		}
		renderCycles = cpl::Misc::ClockCounter() - cStart;
		auto tickNow = juce::Time::getHighResolutionTicks();
		avgFps.setNext(tickNow - lastFrameTick);
		lastFrameTick = tickNow;
	}



	template<typename V>
		void CVectorScope::drawPolarPlot(cpl::OpenGLEngine::COpenGLStack & openGLStack, const cpl::AudioBuffer & buf)
		{
			using namespace cpl::simd;
			typedef typename scalar_of<V>::type Ty;

			cpl::OpenGLEngine::MatrixModification matrixMod;
			auto const gain = (GLfloat)getGain();
			auto const numSamples = buf[0].size;
			const bool fillPath = kdrawLines.bGetValue() > 0.5;
			const bool fadeHistory = kfadeOld.bGetValue() > 0.5;
			static const long long vectorLength = elements_of<V>::value;
			matrixMod.scale(gain, gain, 1);
			auto const colour = kdrawingColour.getControlColourAsColour();
			suitable_container<V> outX, outY, outFade;

			// simd consts
			static const Ty cosineRotation = (Ty)std::cos(M_PI * 135.0 / 180);
			static const Ty sineRotation = (Ty)std::sin(M_PI * 135.0 / 180);
			static const V vCosine = set1<V>(cosineRotation);
			static const V vSine = set1<V>(sineRotation);
			const V vZero = zero<V>();
			const V vOne = consts<V>::one;
			const V vSign = consts<V>::sign_bit;
			auto const fadePerSample = (Ty)1.0 / numSamples;
			auto const vIncrementalFade = set1<V>(fadePerSample * vectorLength);

			const float
				red = colour.getFloatRed(),
				green = colour.getFloatGreen(),
				blue = colour.getFloatBlue();

			for (int i = 0; i < vectorLength; ++i)
			{
				outFade[i] = fadePerSample * i;
			}

			V vSampleFade = outFade;


			auto const lit = buf[0].getIterator<vectorLength>();
			auto const rit = buf[1].getIterator<vectorLength>();
			if(!fadeHistory)
			{
				cpl::OpenGLEngine::PrimitiveDrawer<1024> drawer(openGLStack, fillPath ? GL_LINE_STRIP : GL_POINTS);
				drawer.addColour(red, green, blue);
				// iterate front of buffer, then back.
				for (int section = 0; section < 2; ++section)
				{

					// using long longs to safely jump out of loops with elements_of<V> > sectionSamples
					long long i = 0;

					const Ty * left = lit.getIndex(section);
					const Ty * right = rit.getIndex(section);
					const long long sectionSamples = lit.sizes[section];

					for (; i < (sectionSamples - vectorLength); i += vectorLength)
					{
						V vLeft = loadu<V>(left + i);
						V vRight = loadu<V>(right + i);

						// the length of the hypotenuse of the triangle, we 
						// convert the unit square to.
						auto const vLength = max(abs(vLeft), abs(vRight));

						// rotate our view manually (to center on Y-axis)
						V vY = vLeft * vCosine - vRight * vSine;
						V vX = vLeft * vSine + vRight * vCosine;

						// check for any zero elements.
						vLeft = (vLeft == vZero);
						vRight = (vRight == vZero);
						auto vMask = vnot(vand(vLeft, vRight));

						// get the phase angle. use atan2 if you want to draw the full circle.
						// x and y are swapped at this point, btw.
						auto vAngle = atan(vX / vY);
						// replace nan elements of angle with zero
						vAngle = vand(vMask, vAngle);
						// calcuate x,y coordinates for the right triangle
						sincos(vAngle, &vX, &vY);

						// construct triangle.
						outX = vX * vLength;
						outY = vY * vLength;

						outFade = vSampleFade - vOne;

						// draw vertices.
						for (cpl::Types::fint_t n = 0; n < vectorLength; ++n)
						{
							drawer.addVertex(outX[n], outY[n], outFade[n]);
						}

						vSampleFade += vIncrementalFade;

					}
					//continue;
					// deal with remainder, scalar route
					long long remaindingSamples = 0;
					auto currentSampleFade = outFade[vectorLength - 1];

					for (; i < sectionSamples; i++, remaindingSamples++)
					{
						Ty vLeft = left[i];
						Ty vRight = right[i];

						// the length of the hypotenuse of the triangle, we 
						// convert the unit square to.
						auto const length = std::max(std::abs(vLeft), std::abs(vRight));

						// rotate our view manually (to center on Y-axis)
						Ty vY = vLeft * cosineRotation - vRight * sineRotation;
						Ty vX = vLeft * sineRotation + vRight * cosineRotation;

						// check for any zero elements.

						// get the phase angle. use atan2 if you want to draw the full circle.
						// x and y are swapped at this point, btw.
						auto angle = std::atan(vX / vY);
						// replace nan elements of angle with zero
						angle = (vLeft == Ty(0) && vRight == Ty(0)) ? Ty(0) : angle;
						// calcuate x,y coordinates for the right triangle
						sincos(angle, &vX, &vY);

						drawer.addVertex(vX * length, vY * length, (currentSampleFade - remaindingSamples * fadePerSample));



					}
					// fractionally increase sample fade levels
					vSampleFade += set1<V>(fadePerSample * remaindingSamples);
				}
			}
			else // apply fading
			{

				const V 
					vRed = set1<V>(red), 
					vGreen = set1<V>(green),
					vBlue = set1<V>(blue);

				suitable_container<V> outRed, outGreen, outBlue;
				// iterate front of buffer, then back.

				cpl::OpenGLEngine::PrimitiveDrawer<1024> drawer(openGLStack, fillPath ? GL_LINE_STRIP : GL_POINTS);

				for (int section = 0; section < 2; ++section)
				{

					// using long longs to safely jump out of loops with elements_of<V> > sectionSamples
					long long i = 0;

					const Ty * left = lit.getIndex(section);
					const Ty * right = rit.getIndex(section);
					const long long sectionSamples = lit.sizes[section];

					for (; i < (sectionSamples - vectorLength); i += vectorLength)
					{
						V vLeft = loadu<V>(left + i);
						V vRight = loadu<V>(right + i);

						// the length of the hypotenuse of the triangle, we 
						// convert the unit square to.
						auto const vLength = max(abs(vLeft), abs(vRight));

						// rotate our view manually (to center on Y-axis)
						V vY = vLeft * vCosine - vRight * vSine;
						V vX = vLeft * vSine + vRight * vCosine;

						// check for any zero elements.
						vLeft = (vLeft == vZero);
						vRight = (vRight == vZero);
						auto vMask = vnot(vand(vLeft, vRight));

						// get the phase angle. use atan2 if you want to draw the full circle.
						// x and y are swapped at this point, btw.
						auto vAngle = atan(vX / vY);
						// replace nan elements of angle with zero
						vAngle = vand(vMask, vAngle);
						// calcuate x,y coordinates for the right triangle
						sincos(vAngle, &vX, &vY);

						// construct triangle.
						outX = vX * vLength;
						outY = vY * vLength;

						outFade = vSampleFade - vOne;

						// set colours

						outRed = vRed * vSampleFade;
						outBlue = vBlue * vSampleFade;
						outGreen = vGreen * vSampleFade;

						// draw vertices.
						for (cpl::Types::fint_t n = 0; n < vectorLength; ++n)
						{
							drawer.addColour(outRed[n], outGreen[n], outBlue[n]);
							drawer.addVertex(outX[n], outY[n], outFade[n]);
						}

						vSampleFade += vIncrementalFade;

					}
					//continue;
					// deal with remainder, scalar route
					long long remaindingSamples = 0;
					auto currentSampleFade = outFade[vectorLength - 1];

					for (; i < sectionSamples; i++, remaindingSamples++)
					{
						Ty vLeft = left[i];
						Ty vRight = right[i];

						// the length of the hypotenuse of the triangle, we 
						// convert the unit square to.
						auto const length = std::max(std::abs(vLeft), std::abs(vRight));

						// rotate our view manually (to center on Y-axis)
						Ty vY = vLeft * cosineRotation - vRight * sineRotation;
						Ty vX = vLeft * sineRotation + vRight * cosineRotation;

						// check for any zero elements.

						// get the phase angle. use atan2 if you want to draw the full circle.
						// x and y are swapped at this point, btw.
						auto angle = std::atan(vX / vY);
						// replace nan elements of angle with zero
						angle = (vLeft == Ty(0) && vRight == Ty(0)) ? Ty(0) : angle;
						// calcuate x,y coordinates for the right triangle
						sincos(angle, &vX, &vY);

						auto currentFade = (currentSampleFade - remaindingSamples * fadePerSample);
						drawer.addColour(red * (currentFade + 1), green * (currentFade + 1), blue * (currentFade + 1));
						drawer.addVertex(vX * length, vY * length, currentFade);

					}
					// fractionally increase sample fade levels
					vSampleFade += set1<V>(fadePerSample * remaindingSamples);

				}
			}
		}



};
