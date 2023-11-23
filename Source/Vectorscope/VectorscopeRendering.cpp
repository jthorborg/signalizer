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

	file:VectorScopeRendering.cpp

		Implementation of all rendering code for the vector scope.

*************************************************************************************/


#include "Vectorscope.h"
#include "VectorscopeParameters.h"
#include <cstdint>
#include <cpl/Mathext.h>
#include <cpl/rendering/OpenGLRasterizers.h>
#include <cpl/simd.h>

namespace Signalizer
{
	// swapping the right channel might give an more intuitive view

	static const char * ChannelDescriptions[] = { "+L", "+R", "-L", "-R", "L", "R", "C"};

	static const float circleScaleFactor = 1.1f;

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

	void VectorScope::paint2DGraphics(juce::Graphics & g, std::size_t numChannels)
	{
		bool paintDiag = content->diagnostics.getNormalizedValue() > 0.5;
		if (paintDiag)
		{
			g.setColour(juce::Colours::blue);

			const auto perf = audioStream->getPerfMeasures();

			char textbuf[4096];

			float averageFps, averageCpu;

			computeAverageStats(averageFps, averageCpu);

			auto scale = oglc->getRenderingScale();

			cpl::sprintfs(textbuf, "%dx%d (%.2f): %.1f fps - %.1f%% cpu, deltaT = %f (rt: %.2f%% - %.2f%%), (as: %.2f%% - %.2f%%)",
				getWidth(), getHeight(), scale, averageFps, averageCpu, openGLDeltaTime(),
				100 * perf.producerUsage,
				100 * perf.producerOverhead,
				100 * perf.consumerUsage,
				100 * perf.consumerOverhead
			);

			g.drawSingleLineText(textbuf, 10, 20);

		}

		auto mouseCheck = globalBehaviour->hideWidgetsOnMouseExit ? isMouseInside.load() : true;

		if (globalBehaviour->showLegend && mouseCheck)
		{
			state.legend.paint(g, state.colourWidget, state.colourBackground);
		}
	}

	void VectorScope::initOpenGL()
	{
		const int imageSize = 128;
		const float fontToPixelScale = 90 / 64.0f;

		textures.clear();

		for (std::size_t i = 0; i < std::extent<decltype(ChannelDescriptions)>::value; ++i)
		{

			juce::Image letter(juce::Image::ARGB, imageSize, imageSize, true);
			{
				juce::Graphics g(letter);
				g.fillAll(juce::Colours::transparentBlack);
				g.setColour(juce::Colours::white);
				g.setFont(letter.getHeight() * fontToPixelScale * 0.5);
				g.drawText(ChannelDescriptions[i], letter.getBounds().toFloat(), juce::Justification::centred, false);
			}
			textures.push_back(std::unique_ptr<juce::OpenGLTexture>((new juce::OpenGLTexture())));
			textures[i]->loadImage(letter);
		}
	}

	void VectorScope::closeOpenGL()
	{
		textures.clear();
	}

	void VectorScope::onOpenGLRendering()
	{
		cpl::simd::dynamic_isa_dispatch<float, RenderingDispatcher>(*this);
	}

	template<typename ISA>
		void VectorScope::vectorGLRendering()
		{
			std::size_t numChannels;
			CPL_DEBUGCHECKGL();
            {
                auto && lockedView = audioStream->getAudioBufferViews();
                handleFlagUpdates();

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
                    if (processor->envelopeMode == EnvelopeModes::PeakDecay)
                    {
                        runPeakFilter<ISA>(lockedView);
                    }
                    else if (processor->envelopeMode == EnvelopeModes::None)
                    {
						processor->envelopeGain = 1;
                    }

                    openGLStack.setLineSize(static_cast<float>(oglc->getRenderingScale()) * state.primitiveSize);
                    openGLStack.setPointSize(static_cast<float>(oglc->getRenderingScale()) * state.primitiveSize);

					numChannels = lockedView.getNumChannels();

                    // draw actual stereoscopic plot
                    if (numChannels >= 2)
                    {
						ColourRotation rotation(state.colourWaveform, lockedView.getNumChannels(), true);

						for (std::size_t c = 0; c < numChannels; c += 2)
						{
							if (state.isPolar)
							{
								drawPolarPlot<ISA>(openGLStack, lockedView, c, rotation);
							}
							else // is Lissajous
							{
								drawRectPlot<ISA>(openGLStack, lockedView, c, rotation);
							}
						}

                    }
                    CPL_DEBUGCHECKGL();

                    openGLStack.setLineSize(static_cast<float>(oglc->getRenderingScale()) * 2.0f);
					openGLStack.enable(GL_MULTISAMPLE);
                    // draw graph and wireframe
                    drawWireFrame<ISA>(openGLStack);
                    CPL_DEBUGCHECKGL();
                    // draw channel text(ures)
                    drawGraphText<ISA>(openGLStack, lockedView);
                    CPL_DEBUGCHECKGL();
                    // draw 2d stuff (like stereo meters)
                    drawStereoMeters<ISA>(openGLStack, lockedView);
                    CPL_DEBUGCHECKGL();
                }
            }

			renderGraphics([&](juce::Graphics & g) { paint2DGraphics(g, numChannels); });

			postFrame();
		}


	template<typename ISA>
		void VectorScope::drawGraphText(cpl::OpenGLRendering::COpenGLStack & openGLStack, const AudioStream::AudioBufferAccess & view)
		{
			openGLStack.enable(GL_TEXTURE_2D);
			openGLStack.setBlender(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			using consts = cpl::simd::consts<float>;
			// draw channel rotations letters.
			if(!state.isPolar)
			{
				auto rotation = -state.rotation * 2 * consts::pi;
				const float heightToWidthFactor = float(getHeight()) / getWidth();
				using namespace cpl::simd;

				cpl::OpenGLRendering::MatrixModification m;
				// this undoes the text squashing due to variable aspect ratios.
				m.scale(heightToWidthFactor, 1.0f, 1.0f);

				// calculate coordinates using sin/cos simd pairs.
				suitable_container<v4sf> phases, xcoords, ycoords;

				// set phases (we rotate L/R etc. text around in a circle)
				phases[0] = rotation;
				phases[1] = consts::pi * 0.5f + rotation;
				phases[2] = consts::pi + rotation;
				phases[3] = consts::pi * 1.5f + rotation;

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

				auto jcolour = state.colourAxis;
				// render texture text at coordinates.
				for (int i = 0; i < 4; ++i)
				{
					cpl::OpenGLRendering::ImageDrawer text(openGLStack, *textures[i]);
					text.setColour(jcolour);
					text.drawAt({ xcoords[i], ycoords[i], 0.1f, 0.1f });
				}
			}
			else
			{
				const float heightToWidthFactor = float(getHeight()) / getWidth();
				using namespace cpl::simd;

				cpl::OpenGLRendering::MatrixModification m;
				// this undoes the text squashing due to variable aspect ratios.
				m.scale(heightToWidthFactor, 1.0f, 1.0f);

				auto jcolour = state.colourAxis;
				float nadd = 1.0f - circleScaleFactor;
				float xcoord = consts::sqrt_half_two * circleScaleFactor / heightToWidthFactor + nadd;
				float ycoord = (state.scalePolar ? consts::half : consts::sqrt_half_two) * circleScaleFactor + nadd;

				// render texture text at coordinates.
				{
					cpl::OpenGLRendering::ImageDrawer text(openGLStack, *textures[Textures::Left]);
					text.setColour(jcolour);
					text.drawAt({ -xcoord + nadd, ycoord, 0.1f, 0.1f });
				}
				{
					cpl::OpenGLRendering::ImageDrawer text(openGLStack, *textures[Textures::Center]);
					text.setColour(jcolour);
					text.drawAt({ 0 + nadd * 0.5f, 1, 0.1f, 0.1f });
				}
				{
					cpl::OpenGLRendering::ImageDrawer text(openGLStack, *textures[Textures::Right]);
					text.setColour(jcolour);
					text.drawAt({ xcoord, ycoord, 0.1f, 0.1f });
				}
			}
		}

		struct Conditional01To11HeightTransform
		{
			Conditional01To11HeightTransform(bool doApply)
			{
				if (doApply)
				{
					m->scale(1, 2, 1);
					m->translate(0, -0.5f, 0);
				}
			}

			cpl::OpenGLRendering::MatrixModification& additional() { return m.reference(); }

			cpl::Utility::LazyStackPointer<cpl::OpenGLRendering::MatrixModification> m;
		};

	template<typename ISA>
		void VectorScope::drawWireFrame(cpl::OpenGLRendering::COpenGLStack & openGLStack)
		{
			openGLStack.setBlender(GL_ONE, GL_ONE_MINUS_SRC_COLOR);
			Conditional01To11HeightTransform m(state.scalePolar && state.isPolar);
			cpl::OpenGLRendering::PrimitiveDrawer<128> drawer(openGLStack, GL_LINES);

			// draw skeleton graph
			if (!state.isPolar)
			{
				drawer.addColour(state.colourWire);
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
				typedef typename ISA::V V;

				using namespace cpl::simd;
				using cpl::simd::abs;
				typedef typename scalar_of<V>::type Ty;
				constexpr ssize_t lanes = elements_of<V>::value;

				const auto angularResolution = static_cast<int>(4 * std::sqrt(getWidth() * getHeight()) / 75.0);

				const auto step = (consts<Ty>::pi_half) / (angularResolution);

				suitable_container<V> phases, sines, cosines;
				drawer.addColour(state.colourWire);

				Ty oldX = 1, oldY = 0;

				int i = 1;

				auto emitXY = [&](float newX, float newY)
				{
					drawer.addVertex(oldX, oldY, 0);
					drawer.addVertex(newX, newY, 0);
					drawer.addVertex(oldX, oldY, -1);
					drawer.addVertex(newX, newY, -1);

					// left part
					drawer.addVertex(-oldX, oldY, 0);
					drawer.addVertex(-newX, newY, 0);
					drawer.addVertex(-oldX, oldY, -1);
					drawer.addVertex(-newX, newY, -1);

					oldX = newX;
					oldY = newY;
				};

				for (; i < angularResolution && i + lanes < angularResolution; i += lanes)
				{
					auto baseAngle = i * step;

					for (std::size_t c = 0; c < lanes; ++c)
					{
						phases[c] = baseAngle + step * c;
					}

					V vSines, vCosines;

					sincos(phases.toType(), &vSines, &vCosines);

					sines = vSines;
					cosines = vCosines;

					for (std::size_t c = 0; c < lanes; ++c)
					{
						emitXY(cosines[c], sines[c]);
					}
				}

				// scalar loop
				for (; i < angularResolution; i++)
				{
					Ty sine, cosine;

					sincos(i * step, &sine, &cosine);
					emitXY(cosine, sine);
				}

				// connect half circles to final coordinates
				emitXY(0, 1); 

				// add front and back horizontal lines.
				drawer.addVertex(-1.0f, 0.0f, 0.0f);
				drawer.addVertex(1.0f, 0.0f, 0.0f);
				drawer.addVertex(-1.0f, 0.0f, -1.0f);
				drawer.addVertex(1.0f, 0.0f, -1.0f);

				const auto sqrtHalfTwo = consts<float>::sqrt_half_two;

				// add critical diagonal phase lines.
				drawer.addVertex(0.0f, 0.0f, 0.0f);
				drawer.addVertex(sqrtHalfTwo, sqrtHalfTwo, 0.0f);
				drawer.addVertex(0.0f, 0.0f, 0.0f);
				drawer.addVertex(-sqrtHalfTwo, sqrtHalfTwo, 0.0f);

				drawer.addVertex(0.0f, 0.0f, -1.0f);
				drawer.addVertex(sqrtHalfTwo, sqrtHalfTwo, -1.0f);
				drawer.addVertex(0.0f, 0.0f, -1.0f);
				drawer.addVertex(-sqrtHalfTwo, sqrtHalfTwo, -1.0f);
			}

			// Draw basic axis
			{
				// TODO: consider whether all rendering should use premultiplied alpha - src compositing or true transparancy
				drawer.addColour(state.colourAxis);
				// front x, y axii
				drawer.addVertex(-1.0f, 0.0f, 0.0f);
				drawer.addVertex(1.0f, 0.0f, 0.0f);
				drawer.addVertex(0.0f, 1.0f, 0.0f);
				drawer.addVertex(0.0f, state.isPolar ? 0.0f : -1.0f, 0.0f);
			}
		}

	template<typename ISA, typename ColourArray>
		void VectorScope::drawRectPlot(cpl::OpenGLRendering::COpenGLStack & openGLStack, const AudioStream::AudioBufferAccess & audio, std::size_t offset, const ColourArray& colours)
		{
			cpl::OpenGLRendering::MatrixModification matrixMod;
			// apply the custom rotation to the waveform
			matrixMod.rotate(state.rotation * 360, 0, 0, 1);
			// and apply the gain:
			const auto gain = static_cast<GLfloat>(processor->envelopeGain * state.userGain);
			matrixMod.scale(gain, gain, 1);
			float sampleFade = 1.0f / std::max<int>(1, static_cast<int>(audio.getNumSamples() - 1));

			if (!state.fadeHistory)
			{
				cpl::OpenGLRendering::PrimitiveDrawer<1024> drawer(openGLStack, state.fillPath ? GL_LINE_STRIP : GL_POINTS);

				drawer.addColour(colours[offset]);

				// TODO: glDrawArrays
				audio.iterate<2, true>
				(
					[&] (std::size_t sampleFrame, AudioStream::DataType left, AudioStream::DataType right)
					{
						drawer.addVertex(right, left, sampleFrame * sampleFade - 1);
					},
					offset
				);

			}
			else
			{
				cpl::OpenGLRendering::PrimitiveDrawer<1024> drawer(openGLStack, state.fillPath ? GL_LINE_STRIP : GL_POINTS);

				auto jcolour = colours[offset];
				float red = jcolour.getFloatRed(), blue = jcolour.getFloatBlue(),
				green = jcolour.getFloatGreen(), alpha = jcolour.getFloatGreen();

				float fade = 0;

				// TODO: glDrawArrays
				audio.iterate<2, true>
				(
					[&] (std::size_t sampleFrame, AudioStream::DataType left, AudioStream::DataType right)
					{
						fade = sampleFrame * sampleFade;
						drawer.addColour(fade * red, fade * green, fade * blue, alpha);
						drawer.addVertex(right, left, fade - 1);
					},
					offset
				);

			}


		}


	template<typename ISA, typename ColourArray>
		void VectorScope::drawPolarPlot(cpl::OpenGLRendering::COpenGLStack & openGLStack, const AudioStream::AudioBufferAccess & audio, std::size_t offset, const ColourArray& colours)
		{
			typedef typename ISA::V V;
			AudioStream::AudioBufferView views[2] = { audio.getView(0 + offset), audio.getView(1 + offset) };

			using namespace cpl::simd;
			using cpl::simd::abs;
			typedef typename scalar_of<V>::type Ty;

			Conditional01To11HeightTransform m(state.scalePolar);

			const auto gain = static_cast<GLfloat>(processor->envelopeGain * state.userGain);
			auto const numSamples = views[0].size();
			// TODO: handle all cases of potential signed overflow.
			typedef std::make_signed<std::size_t>::type ssize_t;
			constexpr ssize_t vectorLength = elements_of<V>::value;
			m.additional().scale(gain, gain, 1);
			suitable_container<V> outX, outY, outFade;

			// simd consts
			const Ty cosineRotation = consts<Ty>::sqrt_half_two_minus;
			const Ty sineRotation = consts<Ty>::sqrt_half_two;
			const V vCosine = consts<V>::sqrt_half_two_minus; // equals e^i*0.75*pi
			const V vSine = consts<V>::sqrt_half_two;
			const V vZero = zero<V>();
			const V vOne = consts<V>::one;
			const V vSign = consts<V>::sign_bit;
			auto const fadePerSample = (Ty)1.0 / numSamples;
			auto const vIncrementalFade = set1<V>(fadePerSample * vectorLength);

			const auto colour = colours[offset];

			const float
				red = colour.getFloatRed(),
				green = colour.getFloatGreen(),
				blue = colour.getFloatBlue();

			for (int i = 0; i < vectorLength; ++i)
			{
				outFade[i] = fadePerSample * i;
			}

			V vSampleFade = outFade;

			if(!state.fadeHistory)
			{
				cpl::OpenGLRendering::PrimitiveDrawer<1024> drawer(openGLStack, state.fillPath ? GL_LINE_STRIP : GL_POINTS);
				drawer.addColour(red, green, blue);
				// iterate front of buffer, then back.

				for (ssize_t section = 0; section < static_cast<ssize_t>(AudioStream::bufferIndices); ++section)
				{

					// using signed ints to safely jump out of loops with elements_of<V> > sectionSamples
					ssize_t i = 0;

					const Ty * left = views[0].getItIndex(section);
					const Ty * right = views[1].getItIndex(section);

					ssize_t sectionSamples = views[0].getItRange(section);

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
						for (cpl::ssize_t n = 0; n < vectorLength; ++n)
						{
							drawer.addVertex(outX[n], outY[n], outFade[n]);
						}

						vSampleFade += vIncrementalFade;

					}

					// deal with remainder, scalar route
					ssize_t remaindingSamples = 0;
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

				cpl::OpenGLRendering::PrimitiveDrawer<1024> drawer(openGLStack, state.fillPath ? GL_LINE_STRIP : GL_POINTS);

				for (ssize_t section = 0; section < 2; ++section)
				{

					ssize_t i = 0;

					const Ty * left = views[0].getItIndex(section);
					const Ty * right = views[1].getItIndex(section);

					ssize_t sectionSamples = views[0].getItRange(section);

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
						for (cpl::ssize_t n = 0; n < vectorLength; ++n)
						{
							drawer.addColour(outRed[n], outGreen[n], outBlue[n]);
							drawer.addVertex(outX[n], outY[n], outFade[n]);
						}

						vSampleFade += vIncrementalFade;

					}
					//continue;
					// deal with remainder, scalar route
					ssize_t remaindingSamples = 0;
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

	template<typename ISA>
		void VectorScope::drawStereoMeters(cpl::OpenGLRendering::COpenGLStack & openGLStack, const AudioStream::AudioBufferAccess & audio)
		{
			if (state.colourMeter.getBrightness() == 0)
				return;

			using namespace cpl;
			OpenGLRendering::MatrixModification m;
			m.loadIdentityMatrix();
			const float heightToWidthFactor = float(getHeight()) / getWidth();

			const float balanceX = -0.95f;
			const float balanceLength = 1.9f;
			const float stereoY = -0.85f;
			const float stereoLength = 1.7f;
			const float sideSize = 0.05f;
			const float indicatorSize = 0.05f;

			// remember, y / x
			auto& filters = processor->filters;
			float balanceQuick = std::atan(filters.balance[0][1] / filters.balance[0][0]) / (simd::consts<float>::pi * 0.5f);
			if (!std::isnormal(balanceQuick))
				balanceQuick = 0.5f;
			float balanceSlow = std::atan(filters.balance[1][1] / filters.balance[1][0]) / (simd::consts<float>::pi * 0.5f);
			if (!std::isnormal(balanceSlow))
				balanceSlow = 0.5f;

			const float stereoQuick = filters.phase[0] * 0.5f + 0.5f;
			const float stereoSlow = filters.phase[1] * 0.5f + 0.5f;

			// this undoes the squashing due to variable aspect ratios.
			OpenGLRendering::RectangleDrawer2D<> rect(openGLStack);

			// draw slow balance
			rect.setColour(state.colourMeter.withMultipliedBrightness(0.75f));
			rect.setBounds(balanceX + (balanceLength - indicatorSize) * balanceSlow, balanceX, indicatorSize, sideSize);
			rect.fill();

			// draw quick balance
			rect.setColour(state.colourMeter);
			rect.setBounds(balanceX + (balanceLength - indicatorSize * 0.25f) * balanceQuick, balanceX, indicatorSize * 0.25f, sideSize);
			rect.fill();


			// draw slow stereo
			rect.setColour(state.colourMeter.withMultipliedBrightness(0.75f));
			rect.setBounds(-balanceX, stereoY + (stereoLength - indicatorSize) * stereoSlow, sideSize * heightToWidthFactor, indicatorSize);
			rect.fill();

			// draw quick stereo
			rect.setColour(state.colourMeter);
			rect.setBounds(-balanceX, stereoY + (stereoLength - indicatorSize * 0.25f) * stereoQuick, sideSize * heightToWidthFactor, indicatorSize * 0.25f);
			rect.fill();


			// draw bounding rectangle of balance meter
			rect.setBounds(balanceX, balanceX, balanceLength, sideSize);
			rect.setColour(state.colourMeter.withMultipliedBrightness(0.5f));
			rect.renderOutline();

			// draw bounding rectangle of correlation meter
			rect.setBounds(-balanceX, stereoY, sideSize * heightToWidthFactor, stereoLength);
			rect.setColour(state.colourMeter.withMultipliedBrightness(0.5f));
			rect.renderOutline();

			auto pieceColour = state.colourBackground.contrasting().withMultipliedBrightness(state.colourMeter.getPerceivedBrightness() * 0.5f);
			rect.setColour(pieceColour);
			// draw contrasting center piece for balance
			rect.setBounds(balanceX + balanceLength * 0.5f - indicatorSize * 1/16.f, balanceX, indicatorSize * 0.125f, sideSize);
			rect.fill();

			// draw contrasting center piece for stereo
			rect.setBounds(-balanceX, stereoY + stereoLength * 0.5f - indicatorSize * 1 / 16.f, sideSize * heightToWidthFactor, indicatorSize * 0.125f);
			rect.fill();
		}


	template<typename ISA>
		void VectorScope::runPeakFilter(const AudioStream::AudioBufferAccess & audio)
		{
			typedef typename ISA::V V;

			double currentEnvelope = 1;

			if (audio.getNumChannels() >= 2)
			{
				AudioStream::AudioBufferView views[2] = { audio.getView(0), audio.getView(1) };

				std::size_t numSamples = views[0].size();

				// since this runs in every frame, we need to scale the coefficient by how often this function runs
				// (and the amount of samples)
				double power = numSamples * openGLDeltaTime();

				double coeff = std::pow(processor->envelopeCoeff.load(), power);

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

				processor->filters.envelope[0] = std::max(processor->filters.envelope[0] * coeff, highestLeft  * highestLeft);
				processor->filters.envelope[1] = std::max(processor->filters.envelope[1] * coeff, highestRight * highestRight);

				currentEnvelope = 1.0 / std::max(std::sqrt(processor->filters.envelope[0]), std::sqrt(processor->filters.envelope[1]));
			}

			if (std::isnormal(currentEnvelope))
			{
				processor->envelopeGain = currentEnvelope;
			}
		}
};
