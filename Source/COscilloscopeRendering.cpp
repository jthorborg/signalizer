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

	void COscilloscope::paint2DGraphics(juce::Graphics & g)
	{

		auto cStart = cpl::Misc::ClockCounter();

		if (content->diagnostics.getNormalizedValue() > 0.5)
		{
			auto fps = 1.0 / (avgFps.getAverage() / juce::Time::getHighResolutionTicksPerSecond());

			auto totalCycles = renderCycles + cpl::Misc::ClockCounter() - cStart;
			double cpuTime = (double(totalCycles) / (processorSpeed * 1000 * 1000) * 100) * fps;
			g.setColour(juce::Colours::blue);
			sprintf(textbuf.get(), "%dx%d: %.1f fps - %.1f%% cpu, deltaG = %f, deltaO = %f (rt: %.2f%% - %.2f%%), (as: %.2f%% - %.2f%%)",
				getWidth(), getHeight(), fps, cpuTime, graphicsDeltaTime(), openGLDeltaTime(),
				100 * audioStream.getPerfMeasures().rtUsage.load(std::memory_order_relaxed),
				100 * audioStream.getPerfMeasures().rtOverhead.load(std::memory_order_relaxed),
				100 * audioStream.getPerfMeasures().asyncUsage.load(std::memory_order_relaxed),
				100 * audioStream.getPerfMeasures().asyncOverhead.load(std::memory_order_relaxed));
			g.drawSingleLineText(textbuf.get(), 10, 20);

		}
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
	void COscilloscope::closeOpenGL()
	{
		textures.clear();
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

				openGLStack.setLineSize(static_cast<float>(oglc->getRenderingScale()) * state.primitiveSize);
				openGLStack.setPointSize(static_cast<float>(oglc->getRenderingScale()) * state.primitiveSize);

				drawWavePlot<V>(openGLStack, lockedView);

				CPL_DEBUGCHECKGL();

				openGLStack.setLineSize(static_cast<float>(oglc->getRenderingScale()) * 2.0f);


				CPL_DEBUGCHECKGL();
				// draw channel text(ures)
				drawGraphText<V>(openGLStack, lockedView);
				CPL_DEBUGCHECKGL();
				renderCycles = cpl::Misc::ClockCounter() - cStart;
			}
			renderGraphics(
				[&](juce::Graphics & g)
				{
					// draw graph and wireframe
					//drawWireFrame<V>(openGLStack);
					paint2DGraphics(g);

				}
			);

			auto tickNow = juce::Time::getHighResolutionTicks();
			avgFps.setNext(tickNow - lastFrameTick);
			lastFrameTick = tickNow;

		}


	template<typename V>
		void COscilloscope::drawGraphText(cpl::OpenGLRendering::COpenGLStack & openGLStack, const AudioStream::AudioBufferAccess & view)
		{
			openGLStack.enable(GL_TEXTURE_2D);
			openGLStack.setBlender(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			using consts = cpl::simd::consts<float>;
			// draw channel rotations letters.

			// TODO: Remove rotation
			auto rotation = 0;
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

			auto jcolour = state.colourGraph;
			// render texture text at coordinates.
			for (int i = 0; i < 4; ++i)
			{
				cpl::OpenGLRendering::ImageDrawer text(openGLStack, *textures[i]);
				text.setColour(jcolour);
				text.drawAt({ xcoords[i], ycoords[i], 0.1f, 0.1f });
			}
		}


	template<typename V>
		void COscilloscope::drawWireFrame(juce::Graphics & g, juce::Rectangle<float> rect, float gain)
		{

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

			// TODO: glDrawArrays
			audio.iterate<1, true>
			(
				[&] (std::size_t sampleFrame, AudioStream::DataType & left)
				{
					drawer.addVertex(sampleFrame * sampleDisplacement - 1, left, 0);
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
