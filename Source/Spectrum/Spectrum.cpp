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

	Spectrum::Spectrum(
		std::shared_ptr<const SharedBehaviour>& globalBehaviour,
		std::shared_ptr<const ConcurrentConfig>& config,
		std::shared_ptr<AudioStream::Output>& stream,
		std::shared_ptr<SpectrumContent> params
	)
		: GraphicsWindow(params->getName())
		, globalBehaviour(globalBehaviour)
		, audioStream(stream)
		, state()
		, framePixelPosition()
		, frequencyGraph({ 0, 1 }, { 0, 1 }, 1, 10)
		, complexFrequencyGraph({ 0, 1 }, { 0, 1 }, 1, 10)
		, flags()
		, lastPeak()
		, scallopLoss()
		, oldWindowSize(-1)
		, framesPerUpdate()
		, processor(std::make_shared<ProcessorShell>(globalBehaviour))
		, config(config)
	{
		setOpaque(true);
		if (!(content = std::dynamic_pointer_cast<SpectrumContent>(params)))
		{
			CPL_RUNTIME_EXCEPTION("Cannot cast parameter set's user data to SpectrumContent");
		}

		content->getParameterSet().addRTListener(this, true);

		initPanelAndControls();

		// "first change" flags
		flags.initiateWindowResize = true;
		flags.audioWindowWasResized = true;
		flags.displayModeChange = true;
		flags.audioStreamChanged = true;

		state.viewRect = { 0.0, 1.0 }; // default full-view
		state.sampleRate = 0;
		state.newWindowSize = cpl::Math::round<std::size_t>(content->windowSize.getTransformedValue());

		oldViewRect = state.viewRect;
		oglImage.setFillColour(juce::Colours::black);

		state.minLogFreq = 10;

		{
			auto access = processor->streamState.lock();
			access->constant.setStorage(10, 16, state.transformSize);
		}
		audioStream->addListener(processor);

		state.antialias = true;
		state.primitiveSize = 0.1f;
		resetStaticViewAssumptions();
	}


	void Spectrum::suspend()
	{
		processor->isSuspended = true;
		oldWindowSize = content->windowSize.getTransformedValue();
	}

	void Spectrum::resume()
	{
		processor->isSuspended = false;
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
	
		audioStream->removeListener(processor);
		notifyDestruction();
	}


	void Spectrum::setDBs(double low, double high, bool updateControls)
	{
		content->lowDbs.setTransformedValue(low);
		content->highDbs.setTransformedValue(high);
	}


	Spectrum::DBRange Spectrum::getDBs() const noexcept
	{
		return{ content->lowDbs.getTransformedValue(), content->highDbs.getTransformedValue()};
	}

	void Spectrum::initPanelAndControls()
	{
		setMouseCursor(juce::MouseCursor::DraggingHandCursor);
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

	void Spectrum::calculateSpectrumColourRatios(Constant& constant)
	{
		double acc = 0.0;

		std::array<double, SpectrumContent::numSpectrumColours> vals;

		for (std::size_t i = 0; i < vals.size(); ++i)
		{
			vals[i] = std::max(0.0001, content->specRatios[i].getNormalizedValue());
			acc += vals[i];
		}
		// to avoid accumulating sum >= 1.0f
		acc += std::numeric_limits<float>::epsilon();

		constant.normalizedSpecRatios[0] = 0;
		for (std::size_t i = 0; i < vals.size(); ++i)
		{
			constant.normalizedSpecRatios[i + 1] = static_cast<float>(vals[i] / acc);
		}

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
			// TODO: use delta?
			auto mouseDelta = event.position - this->currentMouse.getPoint();
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

			flags.dynamicRangeChange = true;
			flags.viewChanged = true;
		}

		flags.mouseMove = true;

		GraphicsWindow::mouseDrag(event);
	}


	void Spectrum::parameterChangedRT(cpl::Parameters::Handle localHandle, cpl::Parameters::Handle globalHandle, ParameterSet::BaseParameter * param)
	{
		using namespace cpl;
		// TODO: create parameter indices and turn into switch statement
		if (param == &content->windowSize.parameter)
		{
			state.newWindowSize  = cpl::Math::round<std::size_t>(content->windowSize.getTransformedValue());
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

	void Spectrum::handleFlagUpdates(StreamState& stream)
	{
		bool remapResonator = false;
		bool remapFrequencies = false;
		bool glImageHasBeenResized = false;
		bool calculateLegend = false;
		// did audio stream change since last sync?
		if (state.audioStreamChanged.consumeChanges(stream.audioStreamChangeVersion))
		{
			flags.audioStreamChanged = true;
			flags.audioWindowWasResized = true;
			calculateLegend = true;
		}

		stream.constant.sampleBufferSize = getBlobSamples();
		
		stream.constant.algo = state.algo = content->algorithm.param.getAsTEnum<SpectrumContent::TransformAlgorithm>();
		state.frequencyTrackingGraph = cpl::enum_cast<SpectrumContent::LineGraphs>(content->frequencyTracker.param.getTransformedValue() + SpectrumContent::LineGraphs::None);
		stream.constant.dspWindow = content->dspWin.getWindowType();
		stream.constant.binPolation = state.binPolation = content->binInterpolation.param.getAsTEnum<SpectrumContent::BinInterpolation>();
		state.colourGrid = content->gridColour.getAsJuceColour();
		state.colourBackground = content->backgroundColour.getAsJuceColour();
		state.colourWidget = content->widgetColour.getAsJuceColour();
		state.viewRect = { content->viewLeft.getTransformedValue(), content->viewRight.getTransformedValue() };
		state.drawLegend = content->showLegend.getNormalizedValue() > 0.5f;
		stream.constant.lowDBs = content->lowDbs.getTransformedValue();	
		stream.constant.clipDB = content->lowDbs.getTransformer().transform(0);
		stream.constant.highDBs = content->highDbs.getTransformedValue();

		const auto pairs = stream.pairs.size();

		for (std::size_t i = 0; i < SpectrumContent::LineEnd; ++i)
		{
			state.colourOne[i] = ColourRotation(content->lines[i].colourOne.getAsJuceColour(), pairs, false);
			state.colourTwo[i] = ColourRotation(content->lines[i].colourTwo.getAsJuceColour(), pairs, false);

			double unitFrameTime;
			if (state.displayMode == SpectrumContent::DisplayMode::ColourSpectrum)
				unitFrameTime = content->blobSize.getTransformedValue() / 1000;
			else
				unitFrameTime = openGLDeltaTime();
			stream.constant.filter[i].setSampleRate(fpoint(1.0 / unitFrameTime));
			stream.constant.filter[i].setDecayAsFraction(content->lines[i].decay.getTransformedValue(), 0.1);
		}

		if (state.displayMode == SpectrumContent::DisplayMode::ColourSpectrum)
		{
			stream.constant.colourSpecs[0] = ColourRotation(state.colourBackground, pairs, false);

			for (std::size_t i = 0; i < SpectrumContent::numSpectrumColours; ++i)
			{
				stream.constant.colourSpecs[i + 1] = ColourRotation(content->specColours[i].getAsJuceColour(), pairs, false);
			}

			calculateSpectrumColourRatios(stream.constant);
		}


		state.primitiveSize = content->primitiveSize.getTransformedValue();
		state.alphaFloodFill = content->floodFillAlpha.getTransformedValue();

		auto newconf = content->channelConfiguration.param.getAsTEnum<SpectrumChannels>();

		if (newconf != state.configuration)
		{
			stream.constant.configuration = state.configuration = newconf;
			if (newconf != SpectrumChannels::Complex)
			{
				complexFrequencyGraph.clear();
			}
			flags.viewChanged = true;
			calculateLegend = true;
		}

		if (flags.audioStreamChanged.cas())
		{
			// TODO: Globals
			state.sampleRate = stream.streamLocalSampleRate;
			stream.constant.sampleRate = state.sampleRate;
			flags.viewChanged = true;
		}

		// TODO: on numFilters change (and resizing of buffers), lock the working/audio buffers so that async processing doesn't corrupt anything.
		float sampleRate = getSampleRate();

		if (flags.displayModeChange.cas())
		{
			stream.constant.displayMode = state.displayMode = cpl::enum_cast<SpectrumContent::DisplayMode>(content->displayMode.param.getTransformedValue());
			flags.resized = true;
			flags.resetStateBuffers = true;
			calculateLegend = true;
		}

		std::size_t axisPoints = state.displayMode == SpectrumContent::DisplayMode::LineGraph ? getWidth() : getHeight();

		if (axisPoints != state.axisPoints)
		{
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
			// TODO: globals
			// we will get notified asynchronously in onStreamPropertiesChanged.
			// TODO: cache them there?
			if (config->historyCapacity > 0 && config->sampleRate > 0)
			{
				// only reset this flag if there's valid data, otherwise keep checking.
				flags.initiateWindowResize.cas();

				audioStream->modifyConsumerInfo(
					[&](auto& info)
					{
						info.audioHistorySize = state.newWindowSize;
					}
				);
			}

		}
		if (flags.audioWindowWasResized.cas())
		{
			// TODO: possible difference between parameter and audiostream?
			state.windowSize = getValidWindowSize(config->historySize);
			remapResonator = true;
			flags.audioMemoryResize = true;
		}

		const auto bufSize = cpl::Math::nextPow2Inc(state.windowSize);

		stream.constant.setStorage(state.axisPoints, state.windowSize, state.transformSize);

		if (flags.audioMemoryResize.cas())
		{
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
				frequencyGraph.setMaxFrequency(sampleRate / 2);
				frequencyGraph.setScaling(state.viewScale == SpectrumContent::ViewScaling::Linear ? frequencyGraph.Linear : frequencyGraph.Logarithmic);

				complexFrequencyGraph.setBounds({ 0.0, axisPoints * 0.5 });
				complexFrequencyGraph.setView({ (1 - state.viewRect.right) * axisPoints, (1 - state.viewRect.left) * axisPoints });

				complexFrequencyGraph.setMaxFrequency(sampleRate / 2);
				complexFrequencyGraph.setScaling(state.viewScale == SpectrumContent::ViewScaling::Linear ? frequencyGraph.Linear : frequencyGraph.Logarithmic);
			}

			remapFrequencies = true;
			flags.frequencyGraphChange = true;

			if(state.displayMode == SpectrumContent::DisplayMode::ColourSpectrum && oldViewRect != state.viewRect)
				oglImage.freeLinearVerticalTranslation(oldViewRect, state.viewRect);

			oldViewRect = state.viewRect;

			for (auto& pair : stream.pairs)
				pair.clearLineGraphStates();

			resetStaticViewAssumptions();
		}

		if (remapFrequencies)
		{
			stream.constant.remapFrequencies(state.viewRect, state.viewScale, state.minLogFreq);
			remapResonator = true;
			flags.slopeMapChanged = true;
		}

		if (flags.slopeMapChanged.cas())
		{
			stream.constant.generateSlopeMap(content->slope.derive());
		}

		if (flags.windowKernelChange.cas())
		{
			remapResonator = true;
			stream.constant.regenerateWindowKernel(content->dspWin);
			state.windowScale = stream.constant.windowKernelScale;
		}

		if (remapResonator)
		{
			auto window = content->dspWin.getWindowType();
			stream.constant.remapResonator(content->freeQ.getTransformedValue() > 0.5, cpl::dsp::windowCoefficients<fpoint>(window).second);
			flags.frequencyGraphChange = true;
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
			for(auto& pair : stream.pairs)
				pair.clearAudioState();
		}

		if (calculateLegend)
			recalculateLegend(stream);
	}


	// TODO: Get rid of this pair
	std::size_t Spectrum::getValidWindowSize(std::size_t in) const noexcept
	{
		std::size_t n = std::min(config->historyCapacity.load(), in);
		return n;
	}

	void Spectrum::setWindowSize(std::size_t size)
	{
		state.newWindowSize = getValidWindowSize(size);
		flags.initiateWindowResize = true;
	}

	void Spectrum::resetStaticViewAssumptions()
	{
		lastPeak = -1;
	}

	void Spectrum::recalculateLegend(StreamState& cs)
	{
		const auto legendOffsetX = state.displayMode == SpectrumContent::DisplayMode::ColourSpectrum ? gradientWidth : 0;
		state.legend.reset({ legendOffsetX + 10, 10 });

		const auto numPairs = cs.pairs.size();

		auto staticColour = [&](int index, int subline, std::size_t channel)
		{
			return (index == 1 ? state.colourTwo : state.colourOne)[subline][channel];
		};

		auto c = [&](int index, std::size_t channel) -> LegendCache::ColourItem
		{
			if (state.displayMode == SpectrumContent::DisplayMode::LineGraph)
				return std::make_pair(staticColour(index, 0, channel), staticColour(index, 1, channel));

			return cs.constant.generateSpectrogramGradient(channel);
		};

		switch (cs.constant.configuration)
		{
		default:

		case SpectrumChannels::Left:
			for (std::size_t p = 0; p < numPairs; ++p)
				state.legend.addLine(cs.channelNames[p * 2], c(0, p));
			break;
		case SpectrumChannels::Right:
			for (std::size_t p = 0; p < numPairs; ++p)
				state.legend.addLine(cs.channelNames[p * 2 + 1], c(1, p));
			break;
		case SpectrumChannels::Mid:
			for (std::size_t p = 0; p < numPairs; ++p)
				state.legend.addLine(cs.channelNames[p * 2] + " + " + cs.channelNames[p * 2 + 1], c(0, p));
			break;
		case SpectrumChannels::Side:
			for (std::size_t p = 0; p < numPairs; ++p)
				state.legend.addLine(cs.channelNames[p * 2] + " - " + cs.channelNames[p * 2 + 1], c(1, p));
			break;
		case SpectrumChannels::Separate:
			for (std::size_t p = 0; p < numPairs; ++p)
			{
				state.legend.addLine(cs.channelNames[p * 2], c(0, p));
				state.legend.addLine(cs.channelNames[p * 2 + 1], c(1, p));
			}
			break;
		case SpectrumChannels::MidSide:
		{
			for (std::size_t p = 0; p < numPairs; ++p)
			{
				state.legend.addLine(cs.channelNames[p * 2] + " + " + cs.channelNames[p * 2 + 1], c(0, p));
				state.legend.addLine(cs.channelNames[p * 2] + " - " + cs.channelNames[p * 2 + 1], c(1, p));
			}
			break;
		}
		case SpectrumChannels::Phase:
		{
			for (std::size_t p = 0; p < numPairs; ++p)
			{
				state.legend.addLine("|" + cs.channelNames[p * 2] + "| + |" + cs.channelNames[p * 2 + 1] + "|", c(0, p));
				state.legend.addLine(cs.channelNames[p * 2] + " / " + cs.channelNames[p * 2 + 1], c(1, p));
			}
			break;
		}
		case SpectrumChannels::Complex:
		{
			for (std::size_t p = 0; p < numPairs; ++p)
			{
				state.legend.addLine(cs.channelNames[p * 2] + " + i*" + cs.channelNames[p * 2 + 1], c(0, p));
			}
			break;
		}
		}
	}
};
