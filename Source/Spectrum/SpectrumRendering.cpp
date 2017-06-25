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

	file:CSpectrumRendering.cpp

		Implementation of rendering code for the spectrum.

*************************************************************************************/

#include "Spectrum.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include <cpl/rendering/OpenGLRasterizers.h>
#include <cpl/simd.h>
#include <array>

namespace Signalizer
{
	static const char * Semitones[] =
	{
		"C ",
		"C#",
		"D ",
		"D#",
		"E ",
		"F ",
		"F#",
		"G ",
		"G#",
		"A ",
		"A#",
		"B "
	};

	using namespace cpl;

	/// <summary>
	///
	/// </summary>
	static std::string frequencyToSemitone(double a4Ref, double frequency)
	{
		if (!std::isnormal(frequency))
			return "nan";
		auto note = 12 * std::log2(std::abs(frequency / a4Ref)) + 49;
		auto roundedNote = cpl::Math::round<int>(note);
		auto semitoneIndex = (roundedNote - 4) % 12;
		while(semitoneIndex < 0)
			semitoneIndex += 12;
		auto octave = (roundedNote - semitoneIndex) / 12;
		auto detune = cpl::Math::round<int>(1000 * (note - roundedNote));
		char buf[100];
		sprintf_s(buf, "%s%d%+3.1fc", Semitones[semitoneIndex], octave, detune * 0.1);
		return buf;
	}

	void Spectrum::onGraphicsRendering(juce::Graphics & g)
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

	void Spectrum::paint2DGraphics(juce::Graphics & g)
	{
		auto cStart = cpl::Misc::ClockCounter();

		// ------- draw frequency graph

		char buf[200];
		bool skipText = state.colourGrid.getAlpha() == 0;

		if (state.displayMode == SpectrumContent::DisplayMode::LineGraph)
		{
			if(!skipText)
			{
				auto complexScale = state.configuration == SpectrumChannels::Complex ? 2.0f : 1.0f;
				g.setColour(state.colourGrid);

				const auto & divs = frequencyGraph.getDivisions();
				const auto & cdivs = complexFrequencyGraph.getDivisions();
				// text for frequency divisions
				for (auto & sdiv : divs)
				{
					sprintf_s(buf, "%.2f", sdiv.frequency);
					g.drawText(buf, float(complexScale * sdiv.coord) + 5, 20, 100, 20, juce::Justification::centredLeft);

				}
				// text for complex frequency divisions
				if (state.configuration == SpectrumChannels::Complex)
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
				// text for db divisions
				for (auto & dbDiv : dbGraph.getDivisions())
				{
					sprintf_s(buf, "%.2f", dbDiv.dbVal);
					g.drawText(buf, 5, float(dbDiv.coord), 100, 20, juce::Justification::centredLeft);
				}
			}
		}
		else
		{
			float height = getHeight();
			float baseWidth = getWidth() * 0.05f;

			float gradientOffset = 10.0f;

			if(!skipText)
			{
				g.setColour(state.colourGrid);
				const auto & divs = frequencyGraph.getDivisions();


				for (auto & sdiv : divs)
				{
					sprintf_s(buf, "%.2f", sdiv.frequency);
					g.drawText(buf, gradientOffset + baseWidth + 5, float(height - sdiv.coord) - 10 /* height / 2 */, 100, 20, juce::Justification::centredLeft);
				}
			}

			// draw gradient

			juce::ColourGradient gradient;

			// fill in colours

			double gradientPos = 0.0;

			for (int i = 0; i < SpectrumContent::numSpectrumColours + 1; i++)
			{
				gradientPos += state.normalizedSpecRatios[i];
				gradient.addColour(gradientPos, state.colourSpecs[i].toJuceColour());
			}


			gradient.point1 = {gradientOffset * 0.5f, (float)getHeight() };
			gradient.point2 = {gradientOffset * 0.5f, 0.0f };

			g.setGradientFill(gradient);

			g.fillRect(0.0f, 0.0f, gradientOffset, (float)getHeight());
		}

		laggedFPS = 1.0 / (avgFps.getAverage() / juce::Time::getHighResolutionTicksPerSecond());

		drawFrequencyTracking(g);
		
		if (content->diagnostics.getTransformedValue() > 0.5)
		{
			char text[1000];

			auto totalCycles = renderCycles + cpl::Misc::ClockCounter() - cStart;
			double cpuTime = (double(totalCycles) / (processorSpeed * 1000 * 1000) * 100) * laggedFPS;
			double asu = 100 * audioStream.getPerfMeasures().asyncUsage.load(std::memory_order_relaxed);
			double aso = 100 * audioStream.getPerfMeasures().asyncOverhead.load(std::memory_order_relaxed);
			g.setColour(juce::Colours::blue);
			//TODO: ensure format specifiers are correct always (%llu mostly)
			sprintf(text, "%dx%d {%.3f, %.3f}: %.1f fps - %.1f%% cpu, deltaG = %.4f, deltaO = %.4f (rt: %.2f%% - %.2f%%, d: %llu), (as: %.2f%% - %.2f%%)",
				getWidth(), getHeight(), state.viewRect.left, state.viewRect.right,
				laggedFPS, cpuTime, graphicsDeltaTime(), openGLDeltaTime(),
				100 * audioStream.getPerfMeasures().rtUsage.load(std::memory_order_relaxed),
				100 * audioStream.getPerfMeasures().rtOverhead.load(std::memory_order_relaxed),
				audioStream.getPerfMeasures().droppedAudioFrames.load(std::memory_order_relaxed),
				asu,
				aso);

			g.drawSingleLineText(text, 10, 20);

		}
	}

	void Spectrum::drawFrequencyTracking(juce::Graphics & g)
	{
		if (globalBehaviour.hideWidgetsOnMouseExit.load(std::memory_order_acquire))
		{
			if (!isMouseInside.load(std::memory_order_relaxed))
				return;
		}

		auto graphN = state.frequencyTrackingGraph;
		// for adding colour spectrums, one would need to ensure correct concurrent access to the data structures
		if (graphN == SpectrumContent::LineGraphs::None || state.displayMode == SpectrumContent::DisplayMode::ColourSpectrum)
			return;

		if(state.colourTracker.getAlpha() == 0)
			return;

		g.setColour(state.colourTracker);

		double
			mouseFrequency = 0,
			mouseDBs = 0,
			peakFrequency = 0,
			peakDBs = 0,
			mouseFraction = 0,
			mouseFractionOrth = 0,
			peakFraction = 0,
			peakFractionY = 0,
			peakX = 0,
			peakY = 0,
			mouseX = 0,
			mouseY = 0,
			mouseSlope = 0,
			adjustedScallopLoss = scallopLoss,
			peakDeviance = 0,
			peakSlope = 0,
			peakSlopeDbs = 0,
			sampleRate = audioStream.getInfo().sampleRate.load(std::memory_order_relaxed),
			maxFrequency = sampleRate / 2;

		bool frequencyIsComplex = false;

		char buf[2000];

		double viewSize = state.viewRect.dist();

		if (state.displayMode == SpectrumContent::DisplayMode::LineGraph)
		{
			mouseX = cmouse.x.load(std::memory_order_acquire);
			mouseY = cmouse.y.load(std::memory_order_acquire);

			// a possible concurrent bug
			mouseX = cpl::Math::confineTo(mouseX, 0, getAxisPoints() - 1);
			mouseY = cpl::Math::confineTo(mouseY, 0, getHeight() - 1);

			g.drawLine(static_cast<float>(mouseX), 0, static_cast<float>(mouseX), getHeight(), 1);
			g.drawLine(0, static_cast<float>(mouseY), getWidth(), static_cast<float>(mouseY), 1);

			mouseFraction = mouseX / (getWidth() - 1);
			mouseFractionOrth = mouseY / (getHeight() - 1);
			mouseSlope = 20 * std::log10(slopeMap[mouseX]);

			if (state.viewScale == SpectrumContent::ViewScaling::Logarithmic)
			{
				if (state.configuration != SpectrumChannels::Complex)
				{
					mouseFrequency = state.minLogFreq * std::pow(maxFrequency / state.minLogFreq, state.viewRect.left + viewSize * mouseFraction);
				}
				else
				{
					auto arg = state.viewRect.left + viewSize * mouseFraction;
					if(arg < 0.5)
						mouseFrequency = state.minLogFreq * std::pow(maxFrequency / state.minLogFreq, arg * 2);
					else
					{
						arg -= 0.5;
						auto power = state.minLogFreq * std::pow(maxFrequency / state.minLogFreq, 1.0 - arg * 2);
						mouseFrequency = power;
						frequencyIsComplex = true;
					}
				}
			}
			else
			{
				auto arg = state.viewRect.left + viewSize * mouseFraction;
				if (state.configuration != SpectrumChannels::Complex)
				{
					mouseFrequency = arg * maxFrequency;
				}
				else
				{
					if (arg < 0.5)
					{
						mouseFrequency = 2 * arg * maxFrequency;
					}
					else
					{
						arg -= 0.5;
						mouseFrequency = (1.0 - arg * 2) * maxFrequency;
						frequencyIsComplex = true;
					}
				}
			}
			const auto & dbs = getDBs();
			mouseDBs = cpl::Math::UnityScale::linear(mouseFractionOrth, dbs.high, dbs.low); // y coords are flipped..
		}
		else
		{
			// ?
		}

		mouseFraction = cpl::Math::confineTo(mouseFraction, 0, 1);

		// calculate nearest peak
		double const nearbyFractionToConsider = 0.03;

		auto precisionError = 0.001;
		auto interpolationError = 0.01;

		// TODO: these special cases can be handled (on a rainy day)
		if (state.configuration == SpectrumChannels::Complex || !(state.algo.load(std::memory_order_acquire) == SpectrumContent::TransformAlgorithm::FFT && graphN == SpectrumContent::LineGraphs::Transform))
		{

			if (graphN == SpectrumContent::LineGraphs::Transform)
				graphN = SpectrumContent::LineGraphs::LineMain;

			auto & results = lineGraphs[graphN].results;
			auto N = results.size();
			auto pivot = cpl::Math::round<std::size_t>(N * mouseFraction);
			auto range = cpl::Math::round<std::size_t>(N * nearbyFractionToConsider);

			auto lowerBound = range > pivot ? 0 : pivot - range;
			auto higherBound = range + pivot > N ? N : range + pivot;



			auto peak = std::max_element(results.begin() + lowerBound, results.begin() + higherBound,
				[](const UComplex & left, const UComplex & right) { return left.leftMagnitude < right.leftMagnitude; });

			// scan for continuously rising peaks at boundaries
			if (peak == results.begin() + lowerBound && lowerBound != 0)
			{
				while (true)
				{
					auto nextPeak = peak - 1;
					if (nextPeak == results.begin())
						break;
					else if (nextPeak->leftMagnitude < peak->leftMagnitude)
						break;
					else
						peak = nextPeak;
				}
			}
			else if (peak == results.begin() + (higherBound - 1))
			{
				while (true)
				{
					auto nextPeak = peak + 1;
					if (nextPeak == results.end())
						break;
					else if (nextPeak->leftMagnitude < peak->leftMagnitude)
						break;
					else
						peak = nextPeak;
				}
			}

			auto peakOffset = std::distance(results.begin(), peak);

			// interpolate using a parabolic fit
			// https://ccrma.stanford.edu/~jos/parshl/Peak_Detection_Steps_3.html

			peakFrequency = mappedFrequencies[peakOffset];


			auto offsetIsEnd = peakOffset == static_cast<ptrdiff_t>(mappedFrequencies.size() - 1);
			peakDeviance = mappedFrequencies[offsetIsEnd ? peakOffset : peakOffset + 1] - mappedFrequencies[offsetIsEnd ? peakOffset - 1 : peakOffset];

			if (state.algo.load(std::memory_order_acquire) == SpectrumContent::TransformAlgorithm::FFT)
			{
				// non-smooth interpolations suffer from peak detection losses
				if (state.binPolation != SpectrumContent::BinInterpolation::Lanczos)
					peakDeviance = std::max(peakDeviance, 0.5 * getFFTSpace<std::complex<fftType>>() / N);
			}

			peakX = peakOffset;

			peakSlope = slopeMap[peakOffset];


			peakFractionY = results[peakOffset].leftMagnitude;
			peakY = getHeight() - peakFractionY * getHeight();
			const auto & dbs = getDBs();
			peakDBs = cpl::Math::UnityScale::linear(peakFractionY, dbs.low, dbs.high);

			adjustedScallopLoss = 20 * std::log10(scallopLoss - precisionError * 0.1);

		}
		else
		{
			// search the original FFT
			auto N = getFFTSpace<std::complex<fftType>>();
			auto points = getNumFilters();
			auto lowerBound = cpl::Math::round<cpl::ssize_t>(points * (mouseFraction - nearbyFractionToConsider));
			lowerBound = cpl::Math::round<cpl::ssize_t>((N * mappedFrequencies[cpl::Math::confineTo(lowerBound, 0, points - 1)] / sampleRate));
			auto higherBound = cpl::Math::round<cpl::ssize_t>(points * (mouseFraction + nearbyFractionToConsider));
			higherBound = cpl::Math::round<cpl::ssize_t>((N * mappedFrequencies[cpl::Math::confineTo(higherBound, 0, points - 1)] / sampleRate));

			lowerBound = cpl::Math::confineTo(lowerBound, 0, N);
			higherBound = cpl::Math::confineTo(higherBound, 0, N);

			auto source = getAudioMemory<std::complex<fftType>>();

			auto peak = std::max_element(source + lowerBound, source + higherBound + 1,
				[](const std::complex<fftType> & left, const std::complex<fftType> & right) { return cpl::Math::square(left) < cpl::Math::square(right); });

			// scan for continuously rising peaks at boundaries
			if (peak == source + lowerBound && lowerBound != 0)
			{
				while (true)
				{
					auto nextPeak = peak - 1;
					if (nextPeak == source)
						break;
					else if (cpl::Math::square(*nextPeak) < cpl::Math::square(*peak))
						break;
					else
						peak = nextPeak;
				}
			}
			else if (peak == source + (higherBound - 1))
			{
				while (true)
				{
					auto nextPeak = peak + 1;
					if (nextPeak == source + N)
						break;
					else if (cpl::Math::square(*nextPeak) < cpl::Math::square(*peak))
						break;
					else
						peak = nextPeak;
				}
			}

			auto peakOffset = std::distance(source, peak);

			auto const invSize = windowScale / (getWindowSize() * 0.5);

			// interpolate using a parabolic fit
			// https://ccrma.stanford.edu/~jos/parshl/Peak_Detection_Steps_3.html
			// jos suggests doing the fit in logarithmic domain, it tends to create nans and infs we wouldn't have got otherwise -
			// explaning the various isnormal() checks
			auto alpha = 20 * std::log10(std::abs(source[peakOffset == 0 ? 0 : peakOffset - 1] * invSize));
			auto beta = 20 * std::log10(std::abs(source[peakOffset] * invSize));
			auto gamma = 20 * std::log10(std::abs(source[peakOffset == static_cast<std::ptrdiff_t>(N) ? peakOffset : peakOffset + 1] * invSize));

			auto phi = 0.5 * (alpha - gamma) / (alpha - 2 * beta + gamma);

			// global
			peakFraction = 2 * (peakOffset + (std::isnormal(phi) ? phi : 0)) / double(N);
			peakFrequency = 0.5 * peakFraction * sampleRate;
			// translate to local

			peakDBs = beta - 0.25*(alpha - gamma) * phi;
			if (!std::isnormal(peakDBs))
				peakDBs = 20 * std::log10(std::abs(source[peakOffset]) / (N * 0.5));

			peakX = frequencyGraph.fractionToCoordTransformed(peakFraction);

			peakSlope = slopeMap[cpl::Math::confineTo(cpl::Math::round<std::size_t>(peakX), 0, getNumFilters() - 1)];

			peakDBs += 20 * std::log10(peakSlope);
			const auto & dbs = getDBs();
			peakY = cpl::Math::UnityScale::Inv::linear(peakDBs, dbs.low, dbs.high);

			// notice we only adjust the resulting peak value with the slope if we are inside a non-processed transform
			// (like the raw fft), because it hasn't itself been 'sloped' yet
			peakY = getHeight() - peakY * getHeight();


			// deviance firstly considers bin width in frequency, scaled by bin position (precision gets better as
			// frequency increases) and scaled by assumed precision interpolation
			auto normalizedDeviation = (1.0 - peakFraction) / N;
			peakDeviance = 2 * interpolationError * normalizedDeviation * sampleRate + precisionError;
			adjustedScallopLoss = 1.0 - ((1.0 - scallopLoss) * interpolationError + normalizedDeviation); // subtract error
			adjustedScallopLoss = 20 * std::log10(adjustedScallopLoss - precisionError * 0.2);
		}

		// draw a line to the peak from the mouse
		if (std::isnormal(peakX) && std::isnormal(peakY))
		{

			g.drawLine((float)mouseX, (float)mouseY, (float)peakX, (float)peakY, 1.5f);

			// this mechanism, although technically delayed by a frame, avoids recalculating
			// the costly scalloping loss each frame.
			auto truncatedPeak = (int)peakX;
			if (truncatedPeak != lastPeak)
			{
				scallopLoss = getScallopingLossAtCoordinate(truncatedPeak);
				lastPeak = truncatedPeak;
			}

		}

		peakSlopeDbs = 20 * std::log10(peakSlope);

		// adjust for complex things
		bool peakIsComplex = peakFrequency > sampleRate * 0.5;
		if (peakIsComplex)
			peakFrequency = sampleRate - peakFrequency;

		// TODO: calculate at runtime, at some point.

		double estimatedSize[2] = { static_cast<double>(peakIsComplex || frequencyIsComplex ? 170 : 150), 135 };
		double textOffset[2] = { 20, -estimatedSize[1] };

		if (peakDBs > 1000)
			peakDBs = std::numeric_limits<double>::infinity();
		else if(peakDBs < -1000)
			peakDBs = -std::numeric_limits<double>::infinity();

		peakState.clearNonNormals();
		peakState.design(content->trackerSmoothing.parameter.getValue() * 1000, laggedFPS);

		peakState.process(peakDBs, peakFrequency);

		// TODO: Translate these back
		peakDBs = peakState.getPeakDBs();
		peakFrequency = peakState.getFrequency();

		auto reference = content->referenceTuning.getTransformedValue();
		std::string mouseNote = frequencyToSemitone(reference, mouseFrequency);
		std::string freqNote = frequencyToSemitone(reference, peakFrequency);


		// TODO: use a monospace font for this part
		// also: is printf-style really more readable than C++ formatting..
		sprintf_s(buf,
			u8"+x:  %s%11.5f Hz\n"
			u8"+x:  %s\n"
			u8"+y:  %+9.5f dB\n"
			u8"+/:  %+7.3f dB\n"
			u8"\u039Bx:  %s%11.5f Hz\n"
			u8"\u039Bx:  %s\n"
			u8"\u039B~:  %6.3f Hz\u03C3\n"
			u8"\u039By:  %+9.5f dB\n"
			u8"\u039B/:  %+7.3f dB\n"
			u8"\u039BSL: +%6.4f dB\u03C3 ",
			frequencyIsComplex ? "-i*" : "", mouseFrequency,
			mouseNote.c_str(),
			mouseDBs,
			mouseSlope,
			peakIsComplex ? "-i*" : "", peakFrequency,
			freqNote.c_str(),
			peakDeviance,
			peakDBs,
			peakSlopeDbs,
			-adjustedScallopLoss
		);

		// render text rectangle
		auto xpoint = mouseX + textOffset[0] ;
		if (xpoint + estimatedSize[0] > getWidth())
			xpoint = mouseX - (textOffset[0] + estimatedSize[0]);

		auto ypoint = mouseY + textOffset[1] - textOffset[0];
		if (ypoint  < 0)
			ypoint = mouseY + textOffset[0];

		juce::Rectangle<float> rect { (float)xpoint, (float)ypoint, (float)estimatedSize[0], (float)estimatedSize[1] };

		auto rectInside = rect.withSizeKeepingCentre(estimatedSize[0] * 0.95f, estimatedSize[1] * 0.95f).toType<int>();

		// clear background
		g.setColour(state.colourBackground);
		g.fillRoundedRectangle(rect, 2);


		// reset colour
		g.setColour(state.colourTracker);
		g.drawRoundedRectangle(rect, 2, 0.7f);

		g.setFont(juce::Font(juce::Font::getDefaultMonospacedFontName(), cpl::TextSize::normalText * 0.9f, 0));

		g.drawFittedText(juce::CharPointer_UTF8(buf), rectInside, juce::Justification::centredLeft, 6);

	}

	void Spectrum::initOpenGL()
	{
		//oglImage.resize(getWidth(), getHeight(), false);
		flags.openGLInitiation = true;
		textures.clear();

	}

	void Spectrum::closeOpenGL()
	{
		textures.clear();
		oglImage.offload();
	}


	template<std::size_t N, std::size_t V, cpl::GraphicsND::ComponentOrder order>
		void ColourScale2(cpl::GraphicsND::UPixel<order> * CPL_RESTRICT outPixel,
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

    void Spectrum::onOpenGLRendering()
    {
		cpl::simd::dynamic_isa_dispatch<float, RenderingDispatcher>(*this);
    }

    template<typename ISA>
    void Spectrum::vectorGLRendering()
	{
		auto cStart = cpl::Misc::ClockCounter();
        {

            cpl::CFastMutex audioLock;
            // starting from a clean slate?
            CPL_DEBUGCHECKGL();
            juce::OpenGLHelpers::clear(state.colourBackground);


            for (std::size_t i = 0; i < SpectrumContent::LineGraphs::LineEnd; ++i)
                lineGraphs[i].filter.setSampleRate(fpoint(1.0 / openGLDeltaTime()));

            bool lineTransformReady = false;

            // lock the memory buffers, and do our thing.
            {
                handleFlagUpdates();
                // line graph data for ffts are rendered now.
                if (state.displayMode == SpectrumContent::DisplayMode::LineGraph)
                {
                    audioLock.acquire(audioResource);
                    lineTransformReady = prepareTransform(audioStream.getAudioBufferViews());
                }

            }
            // flags may have altered ogl state
            CPL_DEBUGCHECKGL();

            cpl::OpenGLRendering::COpenGLStack openGLStack;
            // set up openGL
            //openGLStack.setBlender(GL_ONE, GL_ONE_MINUS_SRC_COLOR);
            openGLStack.loadIdentityMatrix();
            CPL_DEBUGCHECKGL();
            CPL_DEBUGCHECKGL();

            openGLStack.setLineSize(static_cast<float>(oglc->getRenderingScale()));
            openGLStack.setPointSize(static_cast<float>(oglc->getRenderingScale()));

            CPL_DEBUGCHECKGL();

            switch (state.displayMode)
            {
            case SpectrumContent::DisplayMode::LineGraph:
                // no need to lock in this case, as this display mode is exclusively switched,
                // such that only we have access to it.
                if (lineTransformReady)
                {
                    audioLock.acquire(audioResource);
                    doTransform();
                    mapToLinearSpace();
                    postProcessStdTransform();
                }
                renderLineGraph<ISA>(openGLStack); break;
            case SpectrumContent::DisplayMode::ColourSpectrum:
                // mapping and processing is already done here.
                renderColourSpectrum<ISA>(openGLStack); break;

            }
            
        }
		CPL_DEBUGCHECKGL();
		renderGraphics([&](juce::Graphics & g) { paint2DGraphics(g); });
		CPL_DEBUGCHECKGL();
        renderCycles = cpl::Misc::ClockCounter() - cStart;
        auto tickNow = juce::Time::getHighResolutionTicks();
        avgFps.setNext(tickNow - lastFrameTick);
        lastFrameTick = tickNow;
	}


	template<typename ISA>
		void Spectrum::renderColourSpectrum(cpl::OpenGLRendering::COpenGLStack & ogs)
		{
			CPL_DEBUGCHECKGL();
			auto pW = oglImage.getWidth();
			if (!pW)
				return;

			if (!state.isFrozen)
			{
				framePixelPosition %= pW;
				auto approximateFrames = getApproximateStoredFrames();
				/*if (approximateFrames == 0)
					approximateFrames = framesPerUpdate;*/
				std::size_t processedFrames = 0;
				framesPerUpdate = approximateFrames + content->frameUpdateSmoothing.getTransformedValue() * (framesPerUpdate - approximateFrames);
				auto framesThisTime = cpl::Math::round<std::size_t>(framesPerUpdate);

				// if there's no buffer smoothing at all, we just capture every frame possible.
				//
				bool shouldCap = content->frameUpdateSmoothing.getTransformedValue() != 0.0;

				while ((!shouldCap || (processedFrames++ < framesThisTime)) && processNextSpectrumFrame())
				{
#pragma message cwarn("Update frames per update each time inside here, but as a local variable! There may come more updates meanwhile.")
					// run the next frame through pixel filters and format it etc.


					for (int i = 0; i < getAxisPoints(); ++i)
					{
//#define SIGNALIZER_VISUALDEBUGTEST
#ifdef SIGNALIZER_VISUALDEBUGTEST
						if (framePixelPosition & 1 && i & 1)
						{
							columnUpdate[i] = { 0xff, 0xFF, 0xff, 0xff };
						}
						else
						{
							columnUpdate[i] = { 0x00, 0x00, 0x00, 0x00 };
						}
#else
						ColourScale2<SpectrumContent::numSpectrumColours + 1, 4>(
							columnUpdate.data() + i,
							lineGraphs[SpectrumContent::LineGraphs::LineMain].results[i].magnitude,
							state.colourSpecs,
							state.normalizedSpecRatios
						);
#endif
					}
					//CPL_DEBUGCHECKGL();
					oglImage.updateSingleColumn(framePixelPosition, columnUpdate, GL_RGBA);
					//CPL_DEBUGCHECKGL();

					framePixelPosition++;
					framePixelPosition %= pW;
				}
			}

			CPL_DEBUGCHECKGL();

			{
				cpl::OpenGLRendering::COpenGLImage::OpenGLImageDrawer imageDrawer(oglImage, ogs);

				imageDrawer.drawCircular((float)((double)(framePixelPosition) / (pW - 1)));
				//imageDrawer.drawCircular((float)content->frameUpdateSmoothing.getNormalizedValue());
			}

			CPL_DEBUGCHECKGL();
			ogs.setLineSize(std::max(0.001f, static_cast<float>(oglc->getRenderingScale())));
			// render grid
			{
				auto normalizedScale = 1.0 / getHeight();

				// draw vertical lines.
				const auto & lines = frequencyGraph.getLines();

				auto norm = [=](double in) { return static_cast<float>(normalizedScale * in * 2.0 - 1.0); };

				float baseWidth = 0.1f;

				float gradientOffset = 10.0f / getWidth() - 1.0f;

				OpenGLRendering::PrimitiveDrawer<128> lineDrawer(ogs, GL_LINES);

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


	template<typename ISA>
	void Spectrum::renderLineGraph(cpl::OpenGLRendering::COpenGLStack & ogs)
	{
		int points = getAxisPoints() - 1;
		// render the flood fill with alpha
		ogs.setBlender(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		ogs.setLineSize(static_cast<float>(oglc->getRenderingScale()));

		OpenGLRendering::MatrixModification m;
		m.translate(-1, -1, 0);
		m.scale(static_cast<GLfloat>(1.0 / (points * 0.5)), 2, 1);

		// removes most of the weird black lines on flood fills.
		ogs.disable(GL_MULTISAMPLE);

		if (state.alphaFloodFill != 0.0f)
		{
			// Flood fill
			auto dbs = getDBs();
			const GLfloat endPoint = dbs.high > dbs.low ? 0 : 1;

			for (int k = SpectrumContent::LineGraphs::LineEnd - 1; k >= 0; --k)
			{
				switch (state.configuration)
				{
				case SpectrumChannels::MidSide:
				case SpectrumChannels::Phase:
				case SpectrumChannels::Separate:
				{
					OpenGLRendering::PrimitiveDrawer<512> lineDrawer(ogs, GL_LINES);
					lineDrawer.addColour(state.colourTwo[k].withAlpha(state.alphaFloodFill));
					for (int i = 0; i < (points + 1); ++i)
					{
						lineDrawer.addVertex(i, lineGraphs[k].results[i].rightMagnitude, -0.5);
						lineDrawer.addVertex(i, endPoint, -0.5);
					}
				}
				// (fall-through intentional)
				case SpectrumChannels::Left:
				case SpectrumChannels::Right:
				case SpectrumChannels::Merge:
				case SpectrumChannels::Side:
				case SpectrumChannels::Complex:
				{
					OpenGLRendering::PrimitiveDrawer<512> lineDrawer(ogs, GL_LINES);
					lineDrawer.addColour(state.colourOne[k].withAlpha(state.alphaFloodFill));
					for (int i = 0; i < (points + 1); ++i)
					{
						lineDrawer.addVertex(i, lineGraphs[k].results[i].leftMagnitude, 0);
						lineDrawer.addVertex(i, endPoint, 0);
					}
				}
				default:
					break;
				}
			}
		}

		state.antialias ? ogs.enable(GL_MULTISAMPLE) : ogs.disable(GL_MULTISAMPLE);

		// render the line graphs
		ogs.setBlender(GL_ONE, GL_ONE_MINUS_SRC_COLOR);
		ogs.setLineSize(std::max(0.001f, static_cast<float>(oglc->getRenderingScale() * state.primitiveSize)));
		// draw back to front
		for (int k = SpectrumContent::LineGraphs::LineEnd - 1; k >= 0; --k)
		{
			switch (state.configuration)
			{
			case SpectrumChannels::MidSide:
			case SpectrumChannels::Phase:
			case SpectrumChannels::Separate:
			{
				OpenGLRendering::PrimitiveDrawer<256> lineDrawer(ogs, GL_LINE_STRIP);
				lineDrawer.addColour(state.colourTwo[k]);
				for (int i = 0; i < (points + 1); ++i)
				{
					lineDrawer.addVertex(i, lineGraphs[k].results[i].rightMagnitude, -0.5);
				}
			}
			// (fall-through intentional)
			case SpectrumChannels::Left:
			case SpectrumChannels::Right:
			case SpectrumChannels::Merge:
			case SpectrumChannels::Side:
			case SpectrumChannels::Complex:
			{
				OpenGLRendering::PrimitiveDrawer<256> lineDrawer(ogs, GL_LINE_STRIP);
				lineDrawer.addColour(state.colourOne[k]);
				for (int i = 0; i < (points + 1); ++i)
				{
					lineDrawer.addVertex(i, lineGraphs[k].results[i].leftMagnitude, 0);
				}
			}
			default:
				break;
			}
		}

		ogs.setBlender(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		//ogs.setBlender(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		ogs.setLineSize(std::max(0.001f, static_cast<float>(oglc->getRenderingScale())));

		// render grid
		if (state.colourGrid.getAlpha() != 0)
		{
			auto normalizedScaleY = 1.0 / getHeight();
			// TODO: fix using a matrix modification instead (cheaper)
			auto normY = [=](double in) {  return static_cast<float>(normalizedScaleY * in * 2.0); };

			{
				const float cscale = state.configuration == SpectrumChannels::Complex ? 2 : 1;
				const float width = getWidth();

				OpenGLRendering::PrimitiveDrawer<128> lineDrawer(ogs, GL_LINES);
				lineDrawer.addColour(state.colourGrid.withMultipliedBrightness(0.5f));

				// draw vertical lines.
				const auto & lines = frequencyGraph.getLines();
				const auto & clines = complexFrequencyGraph.getLines();

				for (auto dline : lines)
				{
					auto line = cscale * static_cast<float>(dline);
					lineDrawer.addVertex(line, -1.0f, 0.0f);
					lineDrawer.addVertex(line, 1.0f, 0.0f);
				}

				for (auto dline : clines)
				{
					auto line = static_cast<float>(width - cscale *  dline);
					lineDrawer.addVertex(line, -1.0f, 0.0f);
					lineDrawer.addVertex(line, 1.0f, 0.0f);
				}

				//m.scale(1, getHeight(), 1);
				lineDrawer.addColour(state.colourGrid);
				const auto & divs = frequencyGraph.getDivisions();
				const auto & cdivs = complexFrequencyGraph.getDivisions();

				for (auto & sdiv : divs)
				{
					auto line = cscale * static_cast<float>(sdiv.coord);
					lineDrawer.addVertex(line, -1.0f, 0.0f);
					lineDrawer.addVertex(line, 1.0f, 0.0f);
				}

				for (auto & sdiv : cdivs)
				{
					auto line = static_cast<float>(width - cscale * sdiv.coord);
					lineDrawer.addVertex(line, -1.0f, 0.0f);
					lineDrawer.addVertex(line, 1.0f, 0.0f);
				}
			}


			m.loadIdentityMatrix();

			OpenGLRendering::PrimitiveDrawer<128> lineDrawer(ogs, GL_LINES);

			lineDrawer.addColour(state.colourGrid.withMultipliedBrightness(0.5f));

			// draw horizontal lines:
			for (auto & dbDiv : dbGraph.getDivisions())
			{
				auto line = 1 - (float)dbDiv.fraction * 2;
				lineDrawer.addVertex(-1.0f, line, 0.0f);
				lineDrawer.addVertex(1.0f, line, 0.0f);
			}
		}



	}

};
