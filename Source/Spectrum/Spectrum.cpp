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

	file:Spectrum.cpp

		Implementation of UI and logic for the spectrum view.

*************************************************************************************/

#include "Spectrum.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include <cpl/rendering/OpenGLRasterizers.h>
#include <cpl/simd.h>
#include <cpl/LexicalConversion.h>
#include <array>

namespace Signalizer
{

	Spectrum::Spectrum(const SharedBehaviour& globalBehaviour, const std::string & nameId, AudioStream & stream, ProcessorState * processorState)
		: COpenGLView(nameId)
		, globalBehaviour(globalBehaviour)
		, audioStream(stream)
		, processorSpeed(0)
		, lastFrameTick(0)
		, lastMousePos()
		, state()
		, framePixelPosition()
		, frequencyGraph({ 0, 1 }, { 0, 1 }, 1, 10)
		, complexFrequencyGraph({ 0, 1 }, { 0, 1 }, 1, 10)
		, flags()
		, droppedAudioFrames()
		, audioThreadUsage()
		, relayWidth()
		, relayHeight()
		, cmouse()
		, lastPeak()
		, scallopLoss()
		, oldWindowSize(-1)
		, framesPerUpdate()
		, laggedFPS()
		, isMouseInside(false)
	{
		setOpaque(true);
		if (!(content = dynamic_cast<SpectrumContent *>(processorState)))
		{
			CPL_RUNTIME_EXCEPTION("Cannot cast parameter set's user data to SpectrumContent");
		}

		content->getParameterSet().addRTListener(this, true);

        processorSpeed = cpl::system::CProcessor::getMHz();
		initPanelAndControls();
		flags.firstChange = true;

		state.viewRect = { 0.0, 1.0 }; // default full-view
		state.sampleRate = 0;
		state.newWindowSize.store(cpl::Math::round<std::size_t>(content->windowSize.getTransformedValue()), std::memory_order_release);

		oldViewRect = state.viewRect;
		oglImage.setFillColour(juce::Colours::black);
		listenToSource(stream);

		state.minLogFreq = 10;

		//setWindowSize(200);


		state.iAuxMode = true;
		state.antialias = true;
		state.primitiveSize = 0.1f;
		sfbuf.sampleBufferSize = 200;
		resetStaticViewAssumptions();
	}


	void Spectrum::suspend()
	{
		state.isSuspended = true;
		oldWindowSize = content->windowSize.getTransformedValue();
	}

	void Spectrum::resume()
	{
		state.isSuspended = false;
		if (oldWindowSize != -1)
		{
			//TODO: possibly unsynchronized. fix to have an internal size instead
			content->windowSize.setTransformedValue(oldWindowSize);
		}

		flags.openGLInitiation = true;
	}

	Spectrum::~Spectrum()
	{
		content->getParameterSet().removeRTListener(this, true);
		detachFromSource();
#pragma message cwarn("Fix this as well.")
		SFrameBuffer::FrameVector * frame;
		while (sfbuf.frameQueue.popElement(frame))
			delete frame;

		notifyDestruction();
	}


	void Spectrum::setDBs(double low, double high, bool updateControls)
	{
		content->lowDbs.setTransformedValue(low);
		content->highDbs.setTransformedValue(high);
		/*low = cpl::Math::confineTo(low, kMinDbs, kMaxDbs);
		high = cpl::Math::confineTo(high, kMinDbs, kMaxDbs);

		// ensure we always have a minimum of 3 dBs of range
		// except in the case when high < kMinDbs. too lazy
		if (low > (high - minDBRange))
		{
			low = high - minDBRange;
		}

		if (state.dynRange.low != low || state.dynRange.high != high)
		{
			if (updateControls)
			{
				klowDbs.bSetValue(cpl::Math::UnityScale::Inv::linear<float>(low, kMinDbs, kMaxDbs));
				khighDbs.bSetValue(cpl::Math::UnityScale::Inv::linear<float>(high, kMinDbs, kMaxDbs));
			}

			state.dynRange.low = low;
			state.dynRange.high = high;

			flags.dynamicRangeChange = true;
		} */


	}


	Spectrum::DBRange Spectrum::getDBs() const noexcept
	{
		return{ content->lowDbs.getTransformedValue(), content->highDbs.getTransformedValue()};
	}

	void Spectrum::initPanelAndControls()
	{
		// preliminary initialization - this will update all controls to match audio properties.
		// it may seem like a hack, but it's well-defined and avoids code duplications.
		onAsyncChangedProperties(audioStream, audioStream.getInfo());
		setMouseCursor(juce::MouseCursor::DraggingHandCursor);
	}

	void Spectrum::serialize(cpl::CSerializer::Archiver & archive, cpl::Version version)
	{

	}

	void Spectrum::deserialize(cpl::CSerializer::Builder & builder, cpl::Version version)
	{
	}

	void Spectrum::resized()
	{
		flags.resized = true;
	}

	void Spectrum::freeze()
	{
		state.isFrozen = true;
	}

	void Spectrum::unfreeze()
	{
		state.isFrozen = false;
	}


	void Spectrum::mouseWheelMove(const juce::MouseEvent& event, const juce::MouseWheelDetails& wheel)
	{
		double newFreqPos(0), newDBPos(0);

		switch (state.displayMode)
		{
			case SpectrumContent::DisplayMode::ColourSpectrum:
				newFreqPos = double(getAxisPoints() - event.position.y) / getAxisPoints();
				newDBPos = 0.5;
				break;
			case SpectrumContent::DisplayMode::LineGraph:
				newFreqPos = double(event.position.x) / getAxisPoints();
				newDBPos = double(getHeight() - event.position.y) / getHeight();
				break;
		}

		// shift down equals modification of dbs instead.
		if (!event.mods.isShiftDown())
		{
			// delta difference, scales
			auto left = content->viewLeft.getTransformedValue();
			auto right = content->viewRight.getTransformedValue();

			auto delta = left - right;
			auto inc = -delta * wheel.deltaY / 5;
			// TODO: change to pow()
			content->viewLeft.setTransformedValue(left + newFreqPos * inc);
			content->viewRight.setTransformedValue(right - (1 - newFreqPos) * inc);

			flags.viewChanged = true;
		}
		else
		{
			auto dbs = getDBs();
			auto delta = dbs.high - dbs.low;
			//if (std::abs(delta) <= minDBRange)
			//	return;
#ifdef CPL_MAC
			// OS X internally is extremely inconsistent between drivers, mouses trackpads and what not
			// best solution seems to be just to consider both axi and get some weird results once in a while
			auto inc = delta * (wheel.deltaY + wheel.deltaX) / 5;
#else
			auto inc = delta * wheel.deltaY / 5;
#endif
			dbs.low += newDBPos * inc;
			dbs.high -= (1 - newDBPos) * inc;

			setDBs(dbs.low, dbs.high, true);

			flags.dynamicRangeChange = true;
		}

	}

	void Spectrum::calculateSpectrumColourRatios()
	{
#pragma message cwarn("Exclude colours that are zero.")
		double acc = 0.0;

		std::array<double, SpectrumContent::numSpectrumColours> vals;

		for (std::size_t i = 0; i < vals.size(); ++i)
		{
			vals[i] = std::max(0.0001, content->specRatios[i].getNormalizedValue());
			acc += vals[i];
		}
		// to avoid accumulating sum >= 1.0f
		acc += std::numeric_limits<float>::epsilon();

		state.normalizedSpecRatios[0] = 0;
		for (std::size_t i = 0; i < vals.size(); ++i)
		{
			state.normalizedSpecRatios[i + 1] = static_cast<float>(vals[i] / acc);
		}

	}

	void Spectrum::mouseMove(const juce::MouseEvent & event)
	{
		cmouse.x.store(event.position.x, std::memory_order_release);
		cmouse.y.store(event.position.y, std::memory_order_release);

		flags.mouseMove = true;
	}

	void Spectrum::mouseDoubleClick(const juce::MouseEvent& event)
	{
		if (event.mods.testFlags(juce::ModifierKeys::leftButtonModifier))
		{
			content->viewLeft.setNormalizedValue(0); content->viewRight.setNormalizedValue(0);
			flags.viewChanged = true;
		}
	}

	void Spectrum::mouseDrag(const juce::MouseEvent& event)
	{
		if (event.mods.isLeftButtonDown())
		{
			auto mouseDelta = event.position - lastMousePos;
			auto left = content->viewLeft.getTransformedValue();
			auto right = content->viewRight.getTransformedValue();
			auto freqDelta = left - right;
			auto freqInc = (freqDelta * (state.displayMode == SpectrumContent::DisplayMode::LineGraph ? mouseDelta.x : mouseDelta.y)) / getAxisPoints();

			auto dbs = getDBs();

			auto dynDelta = dbs.high - dbs.low;
			auto dynInc = (dynDelta * (state.displayMode == SpectrumContent::DisplayMode::LineGraph ? mouseDelta.y / getHeight() : mouseDelta.x / getWidth()));

			dbs.high += dynInc;
			dbs.low += dynInc;

			setDBs(dbs.low, dbs.high, true);

			content->viewLeft.setTransformedValue(left + freqInc);
			content->viewRight.setTransformedValue(right + freqInc);

			lastMousePos = event.position;
			flags.dynamicRangeChange = true;
			flags.viewChanged = true;
		}

		cmouse.x.store(event.position.x, std::memory_order_release);
		cmouse.y.store(event.position.y, std::memory_order_release);

		flags.mouseMove = true;
	}

	void Spectrum::mouseUp(const juce::MouseEvent& event)
	{

	}

	void Spectrum::mouseDown(const juce::MouseEvent& event)
	{
		lastMousePos = event.position;
	}

	void Spectrum::mouseExit(const juce::MouseEvent & e)
	{
		isMouseInside.store(false, std::memory_order_relaxed);
	}

	void Spectrum::mouseEnter(const juce::MouseEvent & e)
	{
		isMouseInside.store(true, std::memory_order_relaxed);
	}

	void Spectrum::parameterChangedRT(cpl::Parameters::Handle localHandle, cpl::Parameters::Handle globalHandle, ParameterSet::BaseParameter * param)
	{
		using namespace cpl;
		// TODO: create parameter indices and turn into switch statement
		if (param == &content->windowSize.parameter)
		{
			state.newWindowSize.store(cpl::Math::round<std::size_t>(content->windowSize.getTransformedValue()), std::memory_order_release);
			flags.initiateWindowResize = true;
		}
		if (param == &content->viewScaling.param.parameter || param == &content->viewLeft.parameter || param == &content->viewRight.parameter)
		{
			flags.viewChanged = true;
		}
		else if (param == &content->displayMode.param.parameter)
		{
			flags.displayModeChange = true;
		}
		else if (param == &content->viewScaling.param.parameter)
		{
			flags.viewChanged = true;
		}
		else if (param == &content->lowDbs.parameter || param == &content->highDbs.parameter)
		{
			flags.dynamicRangeChange = true;
		}
		else if (param == &content->pctForDivision.parameter)
		{
			flags.frequencyGraphChange = true;
			flags.dynamicRangeChange = true;
		}
		else if (param == &content->algorithm.param.parameter)
		{
			flags.resetStateBuffers = true;
		}
		else if (param == &content->dspWin.alpha || param == &content->dspWin.beta || param == &content->dspWin.symmetry || param == &content->dspWin.type)
		{
			flags.windowKernelChange = true;
		}
		else if (param == &content->freeQ.parameter)
		{
			cresonator.setFreeQ(content->freeQ.getTransformedValue() > 0.5);
			flags.windowKernelChange = true;
		}
		else if (param == &content->spectrumStretching.parameter)
		{
			// TODO: only do when state.displayMode == colourspectrum? incurs sync issues
			flags.resized = true;
		}
		// TODO: consider group flag?
		else if (param == &content->slope.base || param == &content->slope.pivot || param == &content->slope.slope)
		{
			flags.slopeMapChanged = true;
		}
	}


	std::size_t Spectrum::getWindowSize() const noexcept
	{
		return state.windowSize;
	}

	void Spectrum::handleFlagUpdates()
	{
		cpl::CMutex audioLock;
		if (flags.internalFlagHandlerRunning)
			CPL_RUNTIME_EXCEPTION("Function is NOT reentrant!");

		flags.internalFlagHandlerRunning = true;
		bool firstRun = false;
		bool remapResonator = false;
		bool remapFrequencies = false;
		bool glImageHasBeenResized = false;


		if (flags.firstChange.cas())
		{
			flags.initiateWindowResize = true;
			flags.audioWindowWasResized = true;
			flags.displayModeChange = true;
			flags.audioStreamChanged = true;
			firstRun = true;
		}


		state.algo.store(content->algorithm.param.getAsTEnum<SpectrumContent::TransformAlgorithm>(), std::memory_order_release);
		state.frequencyTrackingGraph = cpl::enum_cast<SpectrumContent::LineGraphs>(content->frequencyTracker.param.getTransformedValue() + SpectrumContent::LineGraphs::None);
		state.dspWindow.store(content->dspWin.getWindowType(), std::memory_order_release);
		state.binPolation = content->binInterpolation.param.getAsTEnum<SpectrumContent::BinInterpolation>();
		state.colourGrid = content->gridColour.getAsJuceColour();
		state.colourBackground = content->backgroundColour.getAsJuceColour();
		state.colourTracker = content->trackerColour.getAsJuceColour();
		state.viewRect = { content->viewLeft.getTransformedValue(), content->viewRight.getTransformedValue() };

		for (std::size_t i = 0; i < SpectrumContent::LineEnd; ++i)
		{
			state.colourOne[i] = content->lines[i].colourOne.getAsJuceColour();
			state.colourTwo[i] = content->lines[i].colourTwo.getAsJuceColour();
			lineGraphs[i].filter.setDecayAsFraction(content->lines[i].decay.getTransformedValue(), 0.1);
		}

		if (state.algo.load(std::memory_order_relaxed) != SpectrumContent::TransformAlgorithm::FFT)
		{
			state.colourSpecs[0] = state.colourBackground;

			for (std::size_t i = 0; i < SpectrumContent::numSpectrumColours; ++i)
			{
				state.colourSpecs[i + 1] = content->specColours[i].getAsJuceColour();
			}

			calculateSpectrumColourRatios();
		}


		state.primitiveSize = content->primitiveSize.getTransformedValue();
		state.alphaFloodFill = content->floodFillAlpha.getTransformedValue();

		auto newconf = content->channelConfiguration.param.getAsTEnum<SpectrumChannels>();

		if (newconf != state.configuration)
		{
			state.configuration = newconf;
			if (newconf != SpectrumChannels::Complex)
			{
				complexFrequencyGraph.clear();
			}
			flags.viewChanged = true;
		}

		if (flags.audioStreamChanged.cas())
		{
			audioLock.acquire(audioResource);
			state.sampleRate.store(static_cast<float>(audioStream.getAudioHistorySamplerate()), std::memory_order_release);
			flags.viewChanged = true;
		}

		// TODO: on numFilters change (and resizing of buffers), lock the working/audio buffers so that async processing doesn't corrupt anything.
		float sampleRate = getSampleRate();

		if (flags.displayModeChange.cas())
		{
			// ensures any concurrent processing modes gets to finish.
			audioLock.acquire(audioResource);
			state.displayMode = cpl::enum_cast<SpectrumContent::DisplayMode>(content->displayMode.param.getTransformedValue());
			flags.resized = true;
			flags.resetStateBuffers = true;
		}

		std::size_t axisPoints = state.displayMode == SpectrumContent::DisplayMode::LineGraph ? getWidth() : getHeight();

		if (axisPoints != state.axisPoints)
		{
			audioLock.acquire(audioResource);
			flags.resized = true;
			state.axisPoints = state.numFilters = axisPoints;
		}

		const std::size_t numFilters = getNumFilters();

		// TODO: insert messagemanagerlock or rework
		auto divLimitParam = content->pctForDivision.getTransformedValue();
		auto const divLimit = 5 + (state.configuration == SpectrumChannels::Complex ? 0.25 : 1) * (numFilters * 0.02 + 0.5 * (numFilters * divLimitParam));
		auto const divLimitY = 5 + 0.6 * (getHeight() * divLimitParam);

		oglImage.setFillColour(state.colourBackground);

		if (flags.initiateWindowResize)
		{
			// we will get notified asynchronously in onAsyncChangedProperties.
			if (audioStream.getAudioHistoryCapacity() && audioStream.getAudioHistorySamplerate())
			{
				// only reset this flag if there's valid data, otherwise keep checking.
				flags.initiateWindowResize.cas();
				audioStream.setAudioHistorySize(state.newWindowSize.load(std::memory_order_acquire));
			}


		}
		if (flags.audioWindowWasResized.cas())
		{
			audioLock.acquire(audioResource);
			// TODO: possible difference between parameter and audiostream?

			auto current = audioStream.getAudioHistorySize();

			state.windowSize = getValidWindowSize(current);
			cresonator.setWindowSize(8, getWindowSize());
			remapResonator = true;
			flags.audioMemoryResize = true;
		}

		if (flags.audioMemoryResize.cas())
		{
			audioLock.acquire(audioResource);
			const auto bufSize = cpl::Math::nextPow2Inc(state.windowSize);
			// some cases it is nice to have an extra entry (see handling of
			// separating real and imaginary transforms)
			audioMemory.resize((bufSize + 1) * sizeof(std::complex<double>));
			windowKernel.resize(bufSize);
			flags.windowKernelChange = true;
		}

		if (flags.openGLInitiation.cas())
		{
			// will re-load image if necessary
			oglImage.resize(getWidth(), getHeight(), true);
			glImageHasBeenResized = true;
		}
		if (flags.resized.cas())
		{
			audioLock.acquire(audioResource);

			for (std::size_t i = 0; i < SpectrumContent::LineGraphs::LineEnd; ++i)
			{
				lineGraphs[i].resize(numFilters); lineGraphs[i].zero();
			}

			slopeMap.resize(numFilters);
			workingMemory.resize(numFilters * 2 * sizeof(std::complex<double>));

			columnUpdate.resize(getHeight());
			// avoid doing it twice.
			if (!glImageHasBeenResized)
			{
				oglImage.resize(std::max<std::size_t>(1, cpl::Math::round<std::size_t>(getWidth() / content->spectrumStretching.getTransformedValue())), getHeight(), true);
				glImageHasBeenResized = true;
			}

			flags.frequencyGraphChange = true;

			flags.viewChanged = true;
			remapFrequencies = true;
			flags.dynamicRangeChange = true;
			flags.slopeMapChanged = true;
		}

		if (flags.dynamicRangeChange.cas())
		{
			dbGraph.setBounds({ 0.0, (double) getHeight() });
			dbGraph.setDivisionLimit(0.4 * divLimitY);
			auto dynRange = getDBs();
			dbGraph.setLowerDbs(dynRange.low);
			dbGraph.setUpperDbs(dynRange.high);
			dbGraph.compileDivisions();
		}


		if (flags.viewChanged.cas())
		{
			state.viewScale = cpl::enum_cast<SpectrumContent::ViewScaling>(content->viewScaling.param.getTransformedValue());

			if (state.configuration != SpectrumChannels::Complex)
			{
				frequencyGraph.setBounds({ 0.0, (double)axisPoints });
				frequencyGraph.setView({ state.viewRect.left * axisPoints, state.viewRect.right * axisPoints });
				frequencyGraph.setMaxFrequency(sampleRate / 2);
				frequencyGraph.setScaling(state.viewScale == SpectrumContent::ViewScaling::Linear ? frequencyGraph.Linear : frequencyGraph.Logarithmic);
			}
			else
			{

				frequencyGraph.setBounds({ 0.0, axisPoints * 0.5 });
				frequencyGraph.setView({ state.viewRect.left * axisPoints, state.viewRect.right * axisPoints });
				//frequencyGraph.setView({ state.viewRect.left * axisPoints * 0.5, (state.viewRect.right - 0.5) * axisPoints * 0.5});
				frequencyGraph.setMaxFrequency(sampleRate / 2);
				frequencyGraph.setScaling(state.viewScale == SpectrumContent::ViewScaling::Linear ? frequencyGraph.Linear : frequencyGraph.Logarithmic);

				complexFrequencyGraph.setBounds({ 0.0, axisPoints * 0.5 });
				complexFrequencyGraph.setView({ (1 - state.viewRect.right) * axisPoints, (1 - state.viewRect.left) * axisPoints });

				//complexFrequencyGraph.setBounds({ axisPoints * 0.5, axisPoints * 1.0 });
				//complexFrequencyGraph.setView({ (1 - state.viewRect.right) * axisPoints, (1 - (state.viewRect.right - state.viewRect.left) * 0.5 - 0.5) * axisPoints });



				complexFrequencyGraph.setMaxFrequency(sampleRate / 2);
				complexFrequencyGraph.setScaling(state.viewScale == SpectrumContent::ViewScaling::Linear ? frequencyGraph.Linear : frequencyGraph.Logarithmic);
			}

			remapFrequencies = true;
			flags.frequencyGraphChange = true;

			if(state.displayMode == SpectrumContent::DisplayMode::ColourSpectrum && oldViewRect != state.viewRect)
				oglImage.freeLinearVerticalTranslation(oldViewRect, state.viewRect);

			oldViewRect = state.viewRect;

			audioLock.acquire(audioResource);
			for (std::size_t i = 0; i < SpectrumContent::LineGraphs::LineEnd; ++i)
				lineGraphs[i].zero();

			resetStaticViewAssumptions();
		}

		if (remapFrequencies)
		{
			audioLock.acquire(audioResource);
			mappedFrequencies.resize(numFilters);

			double viewSize = state.viewRect.dist();


			switch (state.viewScale)
			{
				case SpectrumContent::ViewScaling::Linear:
				{
					double halfSampleRate = sampleRate * 0.5;
					double complexFactor = state.configuration == SpectrumChannels::Complex ? 2.0 : 1.0;
					double freqPerPixel = halfSampleRate / (numFilters - 1);

					for (std::size_t i = 0; i < numFilters; ++i)
					{
						mappedFrequencies[i] = static_cast<float>(complexFactor * state.viewRect.left * halfSampleRate + complexFactor * viewSize * i * freqPerPixel);
					}

					break;
				}
				case SpectrumContent::ViewScaling::Logarithmic:
				{
					double sampleSize = (numFilters - 1);

					double minFreq = state.minLogFreq;

					double end = sampleRate / 2;
					if (state.configuration != SpectrumChannels::Complex)
					{
						for (std::size_t i = 0; i < numFilters; ++i)
						{
							mappedFrequencies[i] = static_cast<float>(minFreq * std::pow(end / minFreq, state.viewRect.left + viewSize * (i / sampleSize)));
						}

					}
					else
					{
						for (std::size_t i = 0; i < numFilters; ++i)
						{
							auto arg = state.viewRect.left + viewSize * i / sampleSize;
							if (arg < 0.5)
							{
								mappedFrequencies[i] = static_cast<float>(minFreq * std::pow(end / minFreq, arg * 2));
							}
							else
							{
								arg -= 0.5;
								auto power = minFreq * std::pow(end / minFreq, 1.0 - arg * 2);
								mappedFrequencies[i] = static_cast<float>(end + (end - power));
							}
						}

					}
					break;
				}
			}
			remapResonator = true;
			flags.slopeMapChanged = true;
		}

		if (flags.slopeMapChanged.cas())
		{
			cpl::PowerSlopeValue::PowerFunction slopeFunction = content->slope.derive();

			for (std::size_t i = 0; i < numFilters; ++i)
			{
				slopeMap[i] = slopeFunction.b * std::pow(mappedFrequencies[i], slopeFunction.a);
			}
		}

		if (flags.windowKernelChange.cas())
		{
			windowScale = content->dspWin.generateWindow<fftType>(windowKernel, getWindowSize());
			remapResonator = true;
		}

		if (remapResonator)
		{
			audioLock.acquire(audioResource);
			auto window = content->dspWin.getWindowType();
			cresonator.mapSystemHz(mappedFrequencies, mappedFrequencies.size(), cpl::dsp::windowCoefficients<fpoint>(window).second, sampleRate);
			flags.frequencyGraphChange = true;
			relayWidth = getWidth();
			relayHeight = getHeight();
		}

		if (flags.frequencyGraphChange.cas())
		{
			frequencyGraph.setDivisionLimit(divLimit);
			frequencyGraph.compileGraph();
			if (state.configuration == SpectrumChannels::Complex)
			{
				complexFrequencyGraph.setDivisionLimit(divLimit);
				complexFrequencyGraph.compileGraph();
			}
		}

		if (flags.resetStateBuffers.cas())
		{
			audioLock.acquire(audioResource);
			cresonator.resetState();
			for (std::size_t i = 0; i < SpectrumContent::LineGraphs::LineEnd; ++i)
				lineGraphs[i].zero();
			std::memset(audioMemory.data(), 0, audioMemory.size() /* * sizeof(char) */);
			std::memset(workingMemory.data(), 0, workingMemory.size() /* * sizeof(char) */);
		}

		// reset all flags through value-initialization
		flags.internalFlagHandlerRunning = false;
	}

	std::size_t Spectrum::getValidWindowSize(std::size_t in) const noexcept
	{
		std::size_t n = std::min(audioStream.getAudioHistoryCapacity(), in);
		return n;
	}

	void Spectrum::setWindowSize(std::size_t size)
	{
		state.newWindowSize.store(getValidWindowSize(size), std::memory_order_release);
		flags.initiateWindowResize = true;
	}


	void Spectrum::onAsyncChangedProperties(const AudioStream & source, const AudioStream::AudioStreamInfo & before)
	{
		flags.audioStreamChanged = true;
		flags.audioWindowWasResized = true;
	}

	void Spectrum::resetStaticViewAssumptions()
	{
		lastPeak = -1;
	}

};
