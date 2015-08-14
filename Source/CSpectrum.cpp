
#include "CSpectrum.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include "SignalizerDesign.h"
#include <cpl/rendering/OpenGLRasterizers.h>
#include <cpl/simd.h>
#include <cpl/LexicalConversion.h>
#include <array>

namespace Signalizer
{

	static const std::vector<std::string> ViewScaleNames = { "Linear", "Logarithmic" };
	static const std::vector<std::string> AlgorithmNames = { "FFT", "Resonator" };
	static const std::vector<std::string> ChannelNames = { "Left", "Right", "Merge", "Phase", "Separate" };
	static const std::vector<std::string> DisplayModeNames = { "Line graph", "Colour spectrum" };
	static const std::vector<std::string> BinInterpolationNames = { "None", "Linear", "Lanczos" };
	// the minimum level of dbs to display
	static const double kMinDbs = -24 * 16;
	// the maximum level of dbs to display
	static const double kMaxDbs = 24 * 4;

	const double CSpectrum::minDBRange = 3.0;

	std::unique_ptr<juce::Component> CSpectrum::createEditor()
	{
		auto content = new Signalizer::CContentPage();
		if (auto page = content->addPage("Settings", "icons/svg/gear.svg"))
		{
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kviewScaling, 0);
				section->addControl(&kalgorithm, 1);
				section->addControl(&kchannelConfiguration, 0);
				section->addControl(&kbinInterpolation, 1);
				section->addControl(&kdspWindow, 0);
				section->addControl(&kdisplayMode, 1);
				page->addSection(section);
			}
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&klowDbs, 1);
				section->addControl(&khighDbs, 0);
				section->addControl(&kblobSize, 0);
				section->addControl(&kwindowSize, 1);
				section->addControl(&kdecayRate, 0);
				section->addControl(&kpctForDivision, 1);
				page->addSection(section);
			}
		}
		if (auto page = content->addPage("Rendering", "icons/svg/painting.svg"))
		{
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kline1Colour, 0);
				section->addControl(&kline2Colour, 1);
				section->addControl(&kgridColour, 0);
				section->addControl(&kbackgroundColour, 1);
				page->addSection(section);
			}

			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				for (int i = 0; i < ArraySize(kspecColours); ++i)
				{
					section->addControl(&kspecColours[i], 0);
					section->addControl(&kspecRatios[i], 1);
				}
				page->addSection(section);
			}


		}
		if (auto page = content->addPage("Utility", "icons/svg/wrench.svg"))
		{
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kframeUpdateSmoothing, 0);
				page->addSection(section);
			}
			page->addSection(&presetManager, "Preset", false);
		}
		editor = content;
		editor->addComponentListener(this);
		state.isEditorOpen = editor ? true : false;
		return std::unique_ptr<juce::Component>(content);

	}

	CSpectrum::CSpectrum(cpl::AudioBuffer & data)
		:
		audioStream(data),
		processorSpeed(0),
		audioStreamCopy(2),
		lastFrameTick(0),
		lastMousePos(),
		editor(nullptr),
		state(),
		framePixelPosition(),
		isSuspended(),
		frequencyGraph({ 0, 1 }, { 0, 1 }, 1, 10),
		flags(),
		audioFifo(20, 1000),
		droppedAudioFrames(),
		audioThreadUsage(),
		relayWidth(), relayHeight(), audioThreadIsInitiated(false),
		presetManager(this, "spectrum"),
		newWindowSize()
	{
		setOpaque(true);
		processorSpeed = juce::SystemStats::getCpuSpeedInMegaherz();
		initPanelAndControls();

		// see note in CKnobSlider::load()
		kbackgroundColour.bForceEvent();

		state.viewRect = { 0.0, 1.0 }; // default full-view
		oldViewRect = state.viewRect;
		state.displayMode = DisplayMode::ColourSpectrum;
		oglImage.setFillColour(juce::Colours::black);
		listenToSource(data[0]);

		state.minLogFreq = 10;
		peakFilter.setSampleRate((float)getSampleRate());
		setWindowSize(200);
		flags.firstChange = true;

		state.iAuxMode = true;

		sfbuf.sampleBufferSize = 200;
	

	}

	void CSpectrum::componentBeingDeleted(Component & component)
	{
		if (&component == editor)
		{
			editor = nullptr;
			state.isEditorOpen = false;
		}
	}



	void CSpectrum::suspend()
	{
		isSuspended = true;
	}

	void CSpectrum::resume()
	{
		isSuspended = false;
		flags.openGLInitiation = true;
	}

	CSpectrum::~CSpectrum()
	{
		// the audio thread is created inside this flag.
		// and that flag is set by this thread. 
		if (!flags.firstChange)
		{
			// so the thread has been created; wait for it to enter function space.
			cpl::Misc::WaitOnCondition(10000, 
				[&]() 
				{
					return audioThreadIsInitiated.load(std::memory_order_relaxed);
				}
			);

			// signal the audio thread and join it:
			audioFifo.releaseConsumer();

			// weird error-checking
			if (audioThread.get_id() != std::thread::id())
			{
				// try to join it. if it isn't joinable at this point, something very bad has happened.
				if (audioThread.joinable())
					audioThread.join();
				else
					CPL_RUNTIME_EXCEPTION("audio thread crashed.");
			}
		}



#pragma message cwarn("Fix this as well.")
		SFrameBuffer::FrameVector * frame;
		while (sfbuf.frameQueue.popElement(frame))
			delete frame;


		notifyDestruction();
		if (editor)
			editor->removeComponentListener(this);
	}


	void CSpectrum::setDBs(double low, double high, bool updateControls)
	{
		low = cpl::Math::confineTo(low, kMinDbs, kMaxDbs);
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
		}


	}


	CSpectrum::DBRange CSpectrum::getDBs() const noexcept
	{
		return state.dynRange;
	}
	void CSpectrum::initPanelAndControls()
	{
		// ------ listeners --------
		kviewScaling.bAddPassiveChangeListener(this);
		kalgorithm.bAddPassiveChangeListener(this);
		kchannelConfiguration.bAddPassiveChangeListener(this);
		kdisplayMode.bAddPassiveChangeListener(this);
		klowDbs.bAddChangeListener(this);
		khighDbs.bAddChangeListener(this);
		kdecayRate.bAddPassiveChangeListener(this);
		kwindowSize.bAddPassiveChangeListener(this);
		kline1Colour.bAddPassiveChangeListener(this);
		kline2Colour.bAddPassiveChangeListener(this);
		kgridColour.bAddPassiveChangeListener(this);
		kpctForDivision.bAddPassiveChangeListener(this);
		kdspWindow.bAddPassiveChangeListener(this);
		kblobSize.bAddPassiveChangeListener(this);
		kbackgroundColour.bAddPassiveChangeListener(this);
		kframeUpdateSmoothing.bAddPassiveChangeListener(this);
		kbinInterpolation.bAddPassiveChangeListener(this);
		int count = 0;

		for (int i = 0; i < ArraySize(kspecColours); ++i)
		{
			kspecColours[i].bAddPassiveChangeListener(this);
			kspecRatios[i].bAddPassiveChangeListener(this);
			kspecColours[i].bSetTitle("Spectrum " + std::to_string(i + 1));
			kspecRatios[i].bSetTitle("Gradient ratio " + std::to_string(i + 1));
			kspecRatios[i].bSetDescription("How large a part of the gradient this colour occupies (the 4 values are normalized).");
			kspecColours[i].bSetDescription("The four colours (together with the background colour) represents the linear interpolating colour function for the intensity found in the graph, for the colour spectrum, composing a gradient.");
		}

		// ------- formatters -------
		kwindowSize.bAddFormatter(this);
		kdecayRate.bAddFormatter(this);
		klowDbs.bAddFormatter(this);
		khighDbs.bAddFormatter(this);
		kblobSize.bAddFormatter(this);
		kframeUpdateSmoothing.bAddFormatter(this);
		// ------ titles -----------
		kviewScaling.bSetTitle("Graph scale");
		kalgorithm.bSetTitle("Transform algorithm");
		kchannelConfiguration.bSetTitle("Channel configuration");
		kdisplayMode.bSetTitle("Display mode");
		kdspWindow.bSetTitle("Window kernel");
		kframeUpdateSmoothing.bSetTitle("Upd. smoothing");
		kbinInterpolation.bSetTitle("Bin interpolation");
		klowDbs.bSetTitle("Lower limit");
		khighDbs.bSetTitle("Upper limit");
		kwindowSize.bSetTitle("Window size");
		kdecayRate.bSetTitle("Filter decay rate");

		kline1Colour.bSetTitle("Graph 1 colour");
		kline2Colour.bSetTitle("Graph 2 colour");
		kgridColour.bSetTitle("Grid colour");
		kbackgroundColour.bSetTitle("Background colour");

		kpctForDivision.bSetTitle("Grid div. space");
		kblobSize.bSetTitle("Update speed");

		// ------ content ----------
		kviewScaling.setValues(ViewScaleNames);
		kalgorithm.setValues(AlgorithmNames);
		kchannelConfiguration.setValues(ChannelNames);
		kdisplayMode.setValues(DisplayModeNames);
		kbinInterpolation.setValues(BinInterpolationNames);
		std::vector<std::string> windows;

		for (int i = 0; (cpl::dsp::WindowTypes) i != cpl::dsp::WindowTypes::End; ++i)
		{
			windows.push_back(cpl::dsp::Windows::stringFromEnum((cpl::dsp::WindowTypes)i));
		}

		kdspWindow.setValues(windows);

		//kviewScaling.setZeroBasedIndex(0); kalgorithm.setZeroBasedIndex(0); 


		// ------ descriptions -----
		kviewScaling.bSetDescription("Set the scale of the frequency-axis of the coordinate system.");
		kalgorithm.bSetDescription("Select the algorithm used for transforming the incoming audio data.");
		kchannelConfiguration.bSetDescription("Select how the audio channels are interpreted.");
		kdisplayMode.bSetDescription("Select how the information is displayed; line graphs are updated each frame while the colour spectrum maintains the previous history.");
		kdspWindow.bSetDescription("The window function describes a kernel applied to the input signal that alters the spectral leakage, through controlling the ratio between main lobe width and side-lobes.");
		kbinInterpolation.bSetDescription("Choice of interpolation for transform algorithms that produce a discrete set of values instead of an continuous function.");

		klowDbs.bSetDescription("The lower limit of the displayed dynamic range.");
		khighDbs.bSetDescription("The upper limit of the displayed dynamic range");
		kwindowSize.bSetDescription("The window size of the audio data, affects time/frequency resolution.");
		kdecayRate.bSetDescription("Allows the filters to decay more slowly, but still reacting to peaks.");
		kline1Colour.bSetDescription("The colour of the first graph.");
		kline2Colour.bSetDescription("The colour of the second graph.");
		kgridColour.bSetDescription("The colour of the dB/frequency grid.");
		kbackgroundColour.bSetDescription("The colour of the background.");
		kpctForDivision.bSetDescription("The minimum amount of free space that triggers a recursed frequency grid division; smaller values draw more frequency divisions.");
		kblobSize.bSetDescription("Controls how much audio data a horizontal unit represents; effectively controls the update rate of the colour spectrum.");
		kframeUpdateSmoothing.bSetDescription("Reduces jitter in spectrum updates at the (possible) expense of higher graphical latency.");
	

		setMouseCursor(juce::MouseCursor::DraggingHandCursor);
	}

	void CSpectrum::save(cpl::CSerializer::Archiver & archive, long long int version)
	{
		archive << kviewScaling;
		archive << kalgorithm;
		archive << kchannelConfiguration;
		archive << kdisplayMode;
		archive << kdspWindow;
		archive << khighDbs;
		archive << klowDbs;
		archive << kdecayRate;
		archive << kwindowSize;
		archive << kpctForDivision;
		archive << kline1Colour;
		archive << kline2Colour;
		archive << kgridColour;
		archive << kblobSize;
		archive << kbackgroundColour;
		archive << kframeUpdateSmoothing;

		for (int i = 0; i < numSpectrumColours; ++i)
		{
			archive << kspecColours[i];
			archive << kspecRatios[i];
		}

	}

	void CSpectrum::load(cpl::CSerializer::Builder & builder, long long int version)
	{
		try
		{
			builder >> kviewScaling;
			builder >> kalgorithm;
			builder >> kchannelConfiguration;
			builder >> kdisplayMode;
			builder >> kdspWindow;
			// set high first, so low isn't capped
			builder >> khighDbs;
			builder >> klowDbs;
			builder >> kdecayRate;
			builder >> kwindowSize;
			builder >> kpctForDivision;
			builder >> kline1Colour;
			builder >> kline2Colour;
			builder >> kgridColour;
			builder >> kblobSize;
			builder >> kbackgroundColour;
			builder >> kframeUpdateSmoothing;

			for (int i = 0; i < numSpectrumColours; ++i)
			{
				builder >> kspecColours[i];
				builder >> kspecRatios[i];
			}

		}
		catch (std::exception & e)
		{
			CPL_NOOP;
		}
	}

	void CSpectrum::resized()
	{
		flags.resized = true;
	}

	void CSpectrum::freeze()
	{
		state.isFrozen = true;
		std::vector<cpl::CMutex> locks;
		// lock streams firstly if we are synced.
		if (isSynced)
		{
			for (auto & buffer : audioStream)
			{
				locks.emplace_back(buffer);
			}
		}

		for (unsigned i = 0; i < audioStream.size(); ++i)
		{
			audioStream[i].clone(audioStreamCopy[i]);
		}
	}

	void CSpectrum::unfreeze()
	{
		state.isFrozen = false;
	}


	void CSpectrum::mouseWheelMove(const MouseEvent& event, const MouseWheelDetails& wheel)
	{
		double newFreqPos(0), newDBPos(0);

		switch (state.displayMode)
		{
			case DisplayMode::ColourSpectrum:
				newFreqPos = double(getAxisPoints() - event.position.y) / getAxisPoints();
				newDBPos = 0.5;
				break;
			case DisplayMode::LineGraph:
				newFreqPos = double(event.position.x) / getAxisPoints();
				newDBPos = double(getHeight() - event.position.y) / getHeight();
				break;
		}

		// shift down equals modification of dbs instead.
		if (!event.mods.isShiftDown())
		{
			// delta difference, scales
			volatile auto delta = state.viewRect.right - state.viewRect.left;
			volatile auto inc = delta * wheel.deltaY / 5;
			state.viewRect.left += newFreqPos * inc;
			state.viewRect.right -= (1 - newFreqPos) * inc;

			state.viewRect.left = std::max(0.0, state.viewRect.left);
			state.viewRect.right = std::min(state.viewRect.right, 1.0);

			flags.viewChanged = true;
		}
		else
		{
			auto dbs = getDBs();
			volatile auto delta = dbs.high - dbs.low;
			if (std::abs(delta) <= minDBRange)
				return;
			volatile auto inc = delta * wheel.deltaY / 5;
			dbs.low += newDBPos * inc;
			dbs.high -= (1 - newDBPos) * inc;

			setDBs(dbs.low, dbs.high, true);

			flags.dynamicRangeChange = true;
		}

	}

	void Signalizer::CSpectrum::calculateSpectrumColourRatios()
	{
		double acc = 0.0;

		std::array<double, numSpectrumColours> vals;

		for (int i = 0; i < vals.size(); ++i)
		{
			vals[i] = std::max(0.0001, kspecRatios[i].bGetValue());
			acc += vals[i];
		}
		// to avoid accumulating sum >= 1.0f
		acc += std::numeric_limits<float>::epsilon();

		state.normalizedSpecRatios[0] = 0;
		for (int i = 0; i < vals.size(); ++i)
		{
			state.normalizedSpecRatios[i + 1] = static_cast<float>(vals[i] / acc);
		}

	}

	void CSpectrum::mouseDoubleClick(const MouseEvent& event)
	{
		state.viewRect = { 0.0, 1.0 };
		flags.viewChanged = true;
	}
	void CSpectrum::mouseDrag(const MouseEvent& event)
	{
		if (event.mods.isLeftButtonDown())
		{
			auto mouseDelta = event.position - lastMousePos;
			
			auto freqDelta = state.viewRect.left - state.viewRect.right;
			auto freqInc = (freqDelta * (state.displayMode == DisplayMode::LineGraph ? mouseDelta.x : mouseDelta.y)) / getAxisPoints();

			auto dbs = getDBs();

			auto dynDelta = dbs.high - dbs.low;
			auto dynInc = (dynDelta * (state.displayMode == DisplayMode::LineGraph ? mouseDelta.y / getHeight() : mouseDelta.x / getWidth()));

			dbs.high += dynInc;
			dbs.low += dynInc;

			setDBs(dbs.low, dbs.high, true);

			state.viewRect.left += freqInc;
			state.viewRect.right += freqInc;

			state.viewRect.left = std::max(0.0, state.viewRect.left);
			state.viewRect.right = std::min(state.viewRect.right, 1.0);


			lastMousePos = event.position;
			flags.dynamicRangeChange = true;
			flags.viewChanged = true;
		}
	}
	void CSpectrum::mouseUp(const MouseEvent& event)
	{

	}
	void CSpectrum::mouseDown(const MouseEvent& event)
	{
		lastMousePos = event.position;
	}

	void CSpectrum::valueChanged(const cpl::CBaseControl * ctrl)
	{
		using namespace cpl;

		if (ctrl == &kdecayRate)
		{
			peakFilter.setDecayAsFraction(ctrl->bGetValue(), 0.01);
		}
		else if (ctrl == &kviewScaling)
		{
			state.viewScale = kviewScaling.getZeroBasedSelIndex<ViewScaling>();
			flags.viewChanged = true;
		}
		else if (ctrl == &kdisplayMode)
		{
			newDisplayMode = kdisplayMode.getZeroBasedSelIndex<DisplayMode>();
			flags.resized = true;
		}
		else if (ctrl == &kchannelConfiguration)
		{
			state.configuration = kchannelConfiguration.getZeroBasedSelIndex<ChannelConfiguration>();
		}
		else if (ctrl == &kwindowSize)
		{
			setWindowSize(static_cast<std::size_t>(ctrl->bGetValue() * audioStream[0].maxSize()));
		}
		else if (ctrl == &kpctForDivision)
		{
			flags.frequencyGraphChange = true;
			flags.dynamicRangeChange = true;
		}
		else if (ctrl == &kalgorithm)
		{
			state.algo = kalgorithm.getZeroBasedSelIndex<TransformAlgorithm>();
		}
		else if (ctrl == &kline1Colour)
		{
			state.colourOne = kline1Colour.getControlColourAsColour();
		}
		else if (ctrl == &kline2Colour)
		{
			state.colourTwo = kline2Colour.getControlColourAsColour();
		}
		else if (ctrl == &kgridColour)
		{
			state.colourGrid = kgridColour.getControlColourAsColour();
		}
		else if (ctrl == &kbackgroundColour)
		{
			state.colourBackground = kbackgroundColour.getControlColourAsColour();
			state.colourSpecs[0] = state.colourBackground;
		}
		else if (ctrl == &kdspWindow)
		{
			state.dspWindow = kdspWindow.getZeroBasedSelIndex<dsp::WindowTypes>();
			flags.windowKernelChange = true;
		}
		else if (ctrl == &kblobSize)
		{
			state.audioBlobSizeMs = Math::round<std::size_t>(Math::UnityScale::exp(kblobSize.bGetValue(), 1.0, 1000.0));
		}
		else if (ctrl == &kframeUpdateSmoothing)
		{
			state.bufferSmoothing = Math::UnityScale::linear(kframeUpdateSmoothing.bGetValue(), 0.0, 0.996);
		}
		else if (ctrl == &kbinInterpolation)
		{
			state.binPolation = kbinInterpolation.getZeroBasedSelIndex<BinInterpolation>();
		}
		else
		{
			for (int i = 0; i < numSpectrumColours; ++i)
			{
				if (ctrl == (kspecColours + i))
				{
					state.colourSpecs[i + 1] = kspecColours[i].getControlColourAsColour();
					break;
				}
				else if (ctrl == (kspecRatios + i))
				{
					// one ratio change affects the normalization, so recalculate all:
					calculateSpectrumColourRatios();
					break;
				}
			}
		}
	}
	
	bool CSpectrum::valueChanged(cpl::CBaseControl * ctrl)
	{
		if (ctrl == &khighDbs || ctrl == &klowDbs)
		{
			auto low = klowDbs.bGetValue();
			auto high = khighDbs.bGetValue();
			low = cpl::Math::UnityScale::linear<float>(low, kMinDbs, kMaxDbs);
			high = cpl::Math::UnityScale::linear<float>(high, kMinDbs, kMaxDbs);
			setDBs(low, high, false);
			auto const & actualDbs = getDBs();
			// need to repaint low dbs, in case it got changed:
			auto newVal = cpl::Math::UnityScale::Inv::linear(actualDbs.low, kMinDbs, kMaxDbs);
			if (newVal != klowDbs.bGetValue())
				klowDbs.bSetInternal(cpl::Math::UnityScale::Inv::linear(actualDbs.low, kMinDbs, kMaxDbs));

		}

		return false;
	}

	std::size_t CSpectrum::getWindowSize() const noexcept
	{
		return state.windowSize;
	}

	void CSpectrum::handleFlagUpdates()
	{

		if (flags.internalFlagHandlerRunning)
			CPL_RUNTIME_EXCEPTION("Function is NOT reentrant!");

		flags.internalFlagHandlerRunning = true;

		bool remapResonator = false;
		bool remapFrequencies = false;
		bool glImageHasBeenResized = false;
		state.displayMode = newDisplayMode;


		auto const sampleRate = getSampleRate();
		auto const axisPoints = state.displayMode == DisplayMode::LineGraph ? getWidth() : getHeight();

		if (axisPoints != state.axisPoints)
		{
			flags.resized = true;
			state.axisPoints = state.numFilters = axisPoints;
		}

		auto const numFilters = getNumFilters();

		auto const divLimit = 5 + numFilters * 0.02 + 0.5 * (numFilters * kpctForDivision.bGetValue());
		oglImage.setFillColour(state.colourBackground);

		if (flags.audioWindowResize.cas())
		{
			std::atomic_thread_fence(std::memory_order_acquire);

			state.windowSize = newWindowSize;
			for (auto & buf : audioStream)
			{
				buf.setSize(state.windowSize);
			}


			cresonator.setWindowSize(8, getWindowSize());


			//std::memset(memory.data(), 0, memory.size());
			//setBandwidth();
			remapResonator = true;
			flags.audioMemoryResize = true;

		}
		if (flags.audioMemoryResize.cas())
		{
			const auto bufSize = cpl::Math::nextPow2Inc(state.windowSize);
			audioMemory.resize(bufSize * sizeof(std::complex<double>));
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
			filterResults.resize(numFilters);
			filterStates.resize(numFilters);
			workingMemory.resize(numFilters * 2 * sizeof(std::complex<double>));

			columnUpdate.resize(getHeight());
			if (!glImageHasBeenResized)
			{
				oglImage.resize(getWidth(), getHeight(), true);
				glImageHasBeenResized = true;
			}

			flags.frequencyGraphChange = true;

			flags.viewChanged = true;
			remapFrequencies = true;
			flags.dynamicRangeChange = true;
		}
		if (flags.dynamicRangeChange.cas())
		{
			dbGraph.setBounds({ 0.0, (double) getHeight() });
			dbGraph.setDivisionLimit(0.14 * divLimit);
			auto dynRange = getDBs();
			dbGraph.setLowerDbs(dynRange.low);
			dbGraph.setUpperDbs(dynRange.high);
			dbGraph.compileDivisions();
		}


		if (flags.viewChanged.cas())
		{
			frequencyGraph.setBounds({ 0.0, (double)getAxisPoints() });
			frequencyGraph.setView({ state.viewRect.left * axisPoints, state.viewRect.right * axisPoints });
			frequencyGraph.setMaxFrequency(sampleRate / 2);
			frequencyGraph.setScaling(state.viewScale == ViewScaling::Linear ? frequencyGraph.Linear : frequencyGraph.Logarithmic);

			remapFrequencies = true;
			flags.frequencyGraphChange = true;

			if(state.displayMode == DisplayMode::ColourSpectrum && oldViewRect != state.viewRect)
				oglImage.freeLinearVerticalTranslation(oldViewRect, state.viewRect);

			oldViewRect = state.viewRect;

		}

		if (remapFrequencies)
		{
			mappedFrequencies.resize(numFilters);

			double start, stop;
			double viewSize = state.viewRect.dist();


			switch (state.viewScale)
			{
				case ViewScaling::Linear:
				{
					double halfSampleRate = sampleRate * 0.5;
					double freqPerPixel = halfSampleRate / (numFilters - 1);

					for (int i = 0; i < numFilters; ++i)
					{
						mappedFrequencies[i] = static_cast<float>(state.viewRect.left * halfSampleRate + viewSize * i * freqPerPixel);
					}
					break;
				}
				case ViewScaling::Logarithmic:
				{
					double sampleSize = (numFilters - 1);

					double minFreq = state.minLogFreq;

					double end = sampleRate / 2;

					for (int i = 0; i < numFilters; ++i)
					{
						mappedFrequencies[i] = static_cast<float>(minFreq * std::pow(end / minFreq, state.viewRect.left + viewSize * (i / sampleSize)));
					}
					break;
				}
			}
			remapResonator = true;
			
		}

		if (remapResonator)
		{
			cresonator.mapSystemHz(mappedFrequencies, mappedFrequencies.size(), sampleRate);
			flags.frequencyGraphChange = true;
			relayWidth = getWidth();
			relayHeight = getHeight();
		}

		if (flags.frequencyGraphChange.cas())
		{
			frequencyGraph.setDivisionLimit(divLimit);
			frequencyGraph.compileGraph();
		}

		if (flags.windowKernelChange.cas())
		{
			computeWindowKernel();
		}

		if (flags.resetStateBuffers.cas())
		{
			cresonator.resetState();
			std::memset(filterStates.data(), 0, filterStates.size() * sizeof(UComplexFilter<fpoint>));
			std::memset(filterResults.data(), 0, filterResults.size() * sizeof(UComplexFilter<fpoint>));
			std::memset(audioMemory.data(), 0, audioMemory.size() /* * sizeof(char) */);
			std::memset(workingMemory.data(), 0, workingMemory.size() /* * sizeof(char) */);
		}

		if (flags.firstChange.cas())
		{
			audioThread = std::thread(&CSpectrum::audioConsumerThread, this);
			framesPerUpdate = getOptimalFramesPerUpdate();
		}

		// reset all flags through value-initialization
		flags.internalFlagHandlerRunning = false;
	}

	void CSpectrum::setWindowSize(std::size_t size)
	{
		std::size_t n = std::min(audioStream[0].maxSize(), size);
		n -= (n & 0x7); // must be a multiple of 8, due to vectorization
		if (n < 16)
			n = 16;


		newWindowSize = n;
		flags.audioWindowResize = true;


	}

	bool CSpectrum::valueToString(const cpl::CBaseControl * ctrl, std::string & buffer, cpl::iCtrlPrec_t value)
	{
		using namespace cpl;
		char buf[200];
		if (ctrl == &klowDbs || ctrl == &khighDbs)
		{
			auto val = Math::UnityScale::linear<double>(value, kMinDbs, kMaxDbs);
			sprintf_s(buf, "%.2f dBs", val);
			buffer = buf;
			return true;
		}
		else if (ctrl == &kdecayRate)
		{
			sprintf_s(buf, "%.2f dBs/10ms", 20 * std::log10(value));
			buffer = buf;
			return true;
		}
		else if (ctrl == &kwindowSize)
		{
			sprintf_s(buf, "%d", (int)(audioStream[0].maxSize() * value));
			buffer = buf;
			return true;
		}
		else if (ctrl == &kblobSize)
		{
			sprintf_s(buf, "%d ms", Math::round<int>(Math::UnityScale::exp(value, 1.0, 1000.0)));
			buffer = buf;
			return true;
		}
		else if (ctrl == &kframeUpdateSmoothing)
		{
			sprintf_s(buf, "%.2f", Math::UnityScale::linear(kframeUpdateSmoothing.bGetValue(), 0.0, 0.996));
			buffer = buf;
			return true;
		}
		return false;
	}


	bool CSpectrum::stringToValue(const cpl::CBaseControl * ctrl, const std::string & buffer, cpl::iCtrlPrec_t & value)
	{
		using namespace cpl;
		double newVal(0);

		if (ctrl == &klowDbs || ctrl == &khighDbs)
		{
			if (cpl::lexicalConversion(buffer, newVal))
			{
				value = cpl::Math::UnityScale::Inv::linear(cpl::Math::confineTo(newVal, kMinDbs, kMaxDbs), kMinDbs, kMaxDbs);
				return true;
			}
		}
		else if (ctrl == &kdecayRate)
		{
			if (cpl::lexicalConversion(buffer, newVal))
			{
				value = cpl::Math::confineTo(std::pow(10, newVal / 20), 0.0, 1.0);
				return true;
			}
			
		}
		else if (ctrl == &kwindowSize)
		{
			if (cpl::lexicalConversion(buffer, newVal))
			{
				value = cpl::Math::confineTo(newVal / audioStream[0].maxSize(), 0.0, 1.0);
				return true;
			}
		}
		else if (ctrl == &kblobSize)
		{
			if (cpl::lexicalConversion(buffer, newVal))
			{
				value = Math::confineTo(Math::UnityScale::Inv::exp(newVal, 1.0, 1000.0), 0.0, 1.0);
				return true;
			}
		}
		else if (ctrl == &kframeUpdateSmoothing)
		{
			if (cpl::lexicalConversion(buffer, newVal))
			{
				value = Math::UnityScale::Inv::linear(Math::confineTo(value, 0.0, 1.0), 0.0, 0.996);
				return true;
			}
		}
		return false;
	}
	void CSpectrum::onObjectDestruction(const cpl::CBaseControl::ObjectProxy & destroyedObject)
	{
		// hmmm.....
	}


	template<typename V>
		void CSpectrum::audioProcessing(float ** buffer, std::size_t numChannels, std::size_t numSamples)
		{
			using namespace cpl::simd;

			if (numChannels != 2)
				return;




		}



};
