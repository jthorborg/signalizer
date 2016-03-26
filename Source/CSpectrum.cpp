
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
	static const std::vector<std::string> ChannelNames = { "Left", "Right", "Mid/Merge", "Side", "Phase", "Separate", "Mid/Side", "Complex" };
	static const std::vector<std::string> DisplayModeNames = { "Line graph", "Colour spectrum" };
	static const std::vector<std::string> BinInterpolationNames = { "None", "Linear", "Lanczos" };

	// the minimum level of dbs to display
	const double CSpectrum::kMinDbs = -24 * 16;
	// the maximum level of dbs to display
	const double CSpectrum::kMaxDbs = 24 * 4;
	const double StretchMax = 20;
	const double CSpectrum::minDBRange = 3.0;
	const double CSpectrum::primitiveMaxSize = 10;

	std::unique_ptr<juce::Component> CSpectrum::createEditor()
	{
		auto content = new Signalizer::CContentPage();
		if (auto page = content->addPage("Settings", "icons/svg/gear.svg"))
		{
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kviewScaling, 0);
				section->addControl(&kchannelConfiguration, 0);
				section->addControl(&kdisplayMode, 1);
				section->addControl(&kfrequencyTracker, 1);
				page->addSection(section);
			}
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&klowDbs, 1);
				section->addControl(&khighDbs, 0);
				section->addControl(&kblobSize, 0);
				section->addControl(&kwindowSize, 1);
				section->addControl(&kpctForDivision, 0);
				section->addControl(&kspectrumStretching, 1);
				page->addSection(section);
			}
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
				{
					section->addControl(&klines[i].decay, i & 1);
				}
				page->addSection(section);
			}
		}
		if (auto page = content->addPage("Algorithm", "icons/svg/formulae.svg"))
		{
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kalgorithm, 0);
				section->addControl(&kbinInterpolation, 1);
				page->addSection(section);
			}
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kdspWin, 0);
				page->addSection(section);
			}
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kfreeQ, 0);
				page->addSection(section);
			}
		}
		if (auto page = content->addPage("Rendering", "icons/svg/brush.svg"))
		{
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kgridColour, 0);
				section->addControl(&kbackgroundColour, 1);
				for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
				{
					section->addControl(&klines[i].colourOne, 0);
					section->addControl(&klines[i].colourTwo, 1);
				}
				page->addSection(section);
			}

			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				for (std::size_t i = 0; i < CPL_ARRAYSIZE(kspecColours); ++i)
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
				section->addControl(&kprimitiveSize, 1);
				page->addSection(section);
			}

			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kfloodFillAlpha, 0);
				section->addControl(&kdiagnostics, 1);
				page->addSection(section);
			}
			page->addSection(&presetManager, "Preset", false);
		}
		editor = content;
		editor->addComponentListener(this);
		state.isEditorOpen = editor ? true : false;
		return std::unique_ptr<juce::Component>(content);

	}

	CSpectrum::CSpectrum(AudioStream & stream)
	:
		COpenGLView("Spectrum view"),
		audioStream(stream),
		processorSpeed(0),
		lastFrameTick(0),
		lastMousePos(),
		editor(nullptr),
		state(),
		framePixelPosition(),
		isSuspended(),
		frequencyGraph({ 0, 1 }, { 0, 1 }, 1, 10),
		complexFrequencyGraph({ 0, 1 }, { 0, 1 }, 1, 10),
		flags(),
		droppedAudioFrames(),
		audioThreadUsage(),
		relayWidth(), relayHeight(), 
		presetManager(this, "spectrum"),
		newc(),
		cmouse(),
		lastPeak(),
		scallopLoss(),
		oldWindowSize(-1)
	{
		setOpaque(true);
		processorSpeed = juce::SystemStats::getCpuSpeedInMegaherz();
		initPanelAndControls();
		flags.firstChange = true;
		// see note in CKnobSlider::load()
		kbackgroundColour.bForceEvent();

		state.viewRect = { 0.0, 1.0 }; // default full-view
#pragma message cwarn("Following variable mustn't be zero. Feels weird to set it here.")
		// TODO: fix
		state.audioBlobSizeMs = 50;
		
		oldViewRect = state.viewRect;
		oglImage.setFillColour(juce::Colours::black);
		listenToSource(stream);

		state.minLogFreq = 10;

		newc.stretch = 1;

		//setWindowSize(200);


		state.iAuxMode = true;
		state.antialias = true;
		state.primitiveSize = 0.1f;
		sfbuf.sampleBufferSize = 200;
		resetStaticViewAssumptions();
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
		oldWindowSize = kwindowSize.bGetValue();
	}

	void CSpectrum::resume()
	{
		isSuspended = false;
		// fix for other views altering our window size.
			if(oldWindowSize != -1)
			{
				kwindowSize.bSetValue(oldWindowSize);
				kwindowSize.bForceEvent(); // in case they (still) are the same
			}

		flags.openGLInitiation = true;
	}

	CSpectrum::~CSpectrum()
	{
		detachFromSource();
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
		// preliminary initialization - this will update all controls to match audio properties.
		// it may seem like a hack, but it's well-defined and avoids code duplications.
		onAsyncChangedProperties(audioStream, audioStream.getInfo());

		// ------ listeners --------
		kviewScaling.bAddPassiveChangeListener(this);
		kalgorithm.bAddPassiveChangeListener(this);
		kchannelConfiguration.bAddPassiveChangeListener(this);
		kdisplayMode.bAddPassiveChangeListener(this);
		kfrequencyTracker.bAddPassiveChangeListener(this);
		klowDbs.bAddChangeListener(this);
		khighDbs.bAddChangeListener(this);
		kdspWin.bAddPassiveChangeListener(this);
		kwindowSize.bAddPassiveChangeListener(this);
		kspectrumStretching.bAddPassiveChangeListener(this);
		kgridColour.bAddPassiveChangeListener(this);
		kpctForDivision.bAddPassiveChangeListener(this);
		kblobSize.bAddPassiveChangeListener(this);
		kbackgroundColour.bAddPassiveChangeListener(this);
		kframeUpdateSmoothing.bAddPassiveChangeListener(this);
		kbinInterpolation.bAddPassiveChangeListener(this);
		kfreeQ.bAddPassiveChangeListener(this);
		kprimitiveSize.bAddPassiveChangeListener(this);
		kfloodFillAlpha.bAddPassiveChangeListener(this);

		for (int i = 0; i < CPL_ARRAYSIZE(kspecColours); ++i)
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
		klowDbs.bAddFormatter(this);
		khighDbs.bAddFormatter(this);
		kblobSize.bAddFormatter(this);
		kframeUpdateSmoothing.bAddFormatter(this);
		kspectrumStretching.bAddFormatter(this);
		kprimitiveSize.bAddFormatter(this);
		// ------ titles -----------
		kviewScaling.bSetTitle("Graph scale");
		kalgorithm.bSetTitle("Transform algorithm");
		kchannelConfiguration.bSetTitle("Channel configuration");
		kdisplayMode.bSetTitle("Display mode");
		kfrequencyTracker.bSetTitle("Frequency tracking");
		kframeUpdateSmoothing.bSetTitle("Upd. smoothing");
		kbinInterpolation.bSetTitle("Bin interpolation");
		klowDbs.bSetTitle("Lower limit");
		khighDbs.bSetTitle("Upper limit");
		kwindowSize.bSetTitle("Window size");
		kfreeQ.bSetTitle("Unbound Q");
		kdiagnostics.bSetTitle("Diagnostics");
		kdiagnostics.setToggleable(true);
		kfreeQ.setToggleable(true);
		kspectrumStretching.bSetTitle("Spectrum stretch");
		kprimitiveSize.bSetTitle("Primitive size");
		kfloodFillAlpha.bSetTitle("Flood fill %");
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


		std::vector<std::string> frequencyTrackingOptions;

		frequencyTrackingOptions.push_back("None");
		frequencyTrackingOptions.push_back("Transform");
		for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
		{
			if(i == LineGraphs::LineMain)
				frequencyTrackingOptions.push_back("Main graph");
			else
				frequencyTrackingOptions.push_back("Aux graph " + std::to_string(i));
		}

		kfrequencyTracker.setValues(frequencyTrackingOptions);

		//kviewScaling.setZeroBasedIndex(0); kalgorithm.setZeroBasedIndex(0); 


		// ------ descriptions -----
		kviewScaling.bSetDescription("Set the scale of the frequency-axis of the coordinate system.");
		kalgorithm.bSetDescription("Select the algorithm used for transforming the incoming audio data.");
		kchannelConfiguration.bSetDescription("Select how the audio channels are interpreted.");
		kdisplayMode.bSetDescription("Select how the information is displayed; line graphs are updated each frame while the colour spectrum maintains the previous history.");
		kbinInterpolation.bSetDescription("Choice of interpolation for transform algorithms that produce a discrete set of values instead of an continuous function.");
		kdiagnostics.bSetDescription("Toggle diagnostic information in top-left corner.");
		klowDbs.bSetDescription("The lower limit of the displayed dynamic range.");
		khighDbs.bSetDescription("The upper limit of the displayed dynamic range");
		kwindowSize.bSetDescription("The window size of the audio data, affects time/frequency resolution.");
		kgridColour.bSetDescription("The colour of the dB/frequency grid.");
		kbackgroundColour.bSetDescription("The colour of the background.");
		kpctForDivision.bSetDescription("The minimum amount of free space that triggers a recursed frequency grid division; smaller values draw more frequency divisions.");
		kblobSize.bSetDescription("Controls how much audio data a horizontal unit represents; effectively controls the update rate of the colour spectrum.");
		kframeUpdateSmoothing.bSetDescription("Reduces jitter in spectrum updates at the (possible) expense of higher graphical latency.");
		kfreeQ.bSetDescription("Frees the quality factor from being bounded by the window size for transforms that support it. "
			"Although it (possibly) makes response time slower, it also makes the time/frequency resolution exact, and is a choice for analyzing static material.");
		kspectrumStretching.bSetDescription("Stretches the spectrum horizontally, emulating a faster update rate (useful for transforms which are not continuous).");
		kfrequencyTracker.bSetDescription("Specifies which pair of graphs that is evaluated for nearby peak estimations.");
		kprimitiveSize.bSetDescription("The size of the rendered primitives (eg. lines or points).");
		kfloodFillAlpha.bSetDescription("For line graphs, add a flood fill of the same colour for each line with the following alpha %");
		klines[LineGraphs::LineMain].colourOne.bSetDescription("The colour of the first channel of the main graph.");
		klines[LineGraphs::LineMain].colourTwo.bSetDescription("The colour of the second channel of the main graph.");

		klines[LineGraphs::LineMain].colourOne.bSetTitle("Graph 1 colour");
		klines[LineGraphs::LineMain].colourTwo.bSetTitle("Graph 2 colour");
		klines[LineGraphs::LineMain].colourOne.bAddPassiveChangeListener(this);
		klines[LineGraphs::LineMain].colourTwo.bAddPassiveChangeListener(this);


		klines[LineGraphs::LineMain].decay.bSetTitle("Main decay");
		klines[LineGraphs::LineMain].decay.bSetDescription("Decay rate of the main graph channels; allows the graph to decay more slowly, but still reacting to peaks.");
		klines[LineGraphs::LineMain].decay.bAddPassiveChangeListener(this);
		klines[LineGraphs::LineMain].decay.bAddFormatter(this);

		for (std::size_t i = LineGraphs::LineMain + 1; i < LineGraphs::LineEnd; ++i)
		{
			auto graphNumber = std::to_string(i);
			klines[i].colourOne.bSetTitle("Aux 1 colour");
			klines[i].colourTwo.bSetTitle("Aux 2 colour");
			klines[i].colourOne.bAddPassiveChangeListener(this);
			klines[i].colourTwo.bAddPassiveChangeListener(this);
			klines[i].decay.bAddPassiveChangeListener(this);
			klines[i].decay.bAddFormatter(this);
			klines[i].decay.bSetTitle("Aux " + graphNumber + " decay");
			klines[i].decay.bSetDescription("Decay rate of auxillary graph " + graphNumber + " channels; allows the graph to decay more slowly, but still reacting to peaks.");
			klines[i].colourOne.bSetDescription("The colour of the first channel of auxillary graph " + graphNumber + ".");
			klines[i].colourTwo.bSetDescription("The colour of the second channel of auxillary graph " + graphNumber + ".");
		}

		setMouseCursor(juce::MouseCursor::DraggingHandCursor);
	}

	void CSpectrum::serialize(cpl::CSerializer::Archiver & archive, cpl::Version version)
	{
		archive << kviewScaling;
		archive << kalgorithm;
		archive << kchannelConfiguration;
		archive << kdisplayMode;
		archive << khighDbs;
		archive << klowDbs;
		archive << kwindowSize;
		archive << kpctForDivision;

		for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
		{
			archive << klines[i].colourOne;
			archive << klines[i].colourTwo;
			archive << klines[i].decay;
		}

		archive << kgridColour;
		archive << kblobSize;
		archive << kbackgroundColour;
		archive << kframeUpdateSmoothing;

		for (int i = 0; i < numSpectrumColours; ++i)
		{
			archive << kspecColours[i];
			archive << kspecRatios[i];
		}

		archive << kbinInterpolation;
		archive << state.viewRect;
		archive << kdspWin;
		archive << kfreeQ;
		archive << kspectrumStretching;
		archive << kfrequencyTracker;
		archive << kprimitiveSize;
		archive << kfloodFillAlpha;
	}

	void CSpectrum::deserialize(cpl::CSerializer::Builder & builder, cpl::Version version)
	{

		builder >> kviewScaling;
		builder >> kalgorithm;
		builder >> kchannelConfiguration;
		builder >> kdisplayMode;
		// set high first, so low isn't capped
		builder >> khighDbs;
		builder >> klowDbs;
		builder >> kwindowSize;
		builder >> kpctForDivision;

		for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
		{
			builder >> klines[i].colourOne;
			builder >> klines[i].colourTwo;
			builder >> klines[i].decay;
		}

		builder >> kgridColour;
		builder >> kblobSize;
		builder >> kbackgroundColour;
		builder >> kframeUpdateSmoothing;

		for (int i = 0; i < numSpectrumColours; ++i)
		{
			builder >> kspecColours[i];
			builder >> kspecRatios[i];
		}
		builder >> kbinInterpolation;
		builder >> state.viewRect;
		builder >> kdspWin;
		builder >> kfreeQ;
		builder >> kspectrumStretching;
		builder >> kfrequencyTracker;
		builder >> kprimitiveSize;
		builder >> kfloodFillAlpha;
	}

	void CSpectrum::resized()
	{
		flags.resized = true;
	}

	void CSpectrum::freeze()
	{
		state.isFrozen = true;
		std::vector<cpl::CMutex> locks;
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
#pragma message cwarn("Exclude colours that are zero.")
		double acc = 0.0;

		std::array<double, numSpectrumColours> vals;

		for (std::size_t i = 0; i < vals.size(); ++i)
		{
			vals[i] = std::max(0.0001, kspecRatios[i].bGetValue());
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

	void CSpectrum::mouseMove(const MouseEvent & event)
	{
		cmouse.x.store(event.position.x, std::memory_order_release);
		cmouse.y.store(event.position.y, std::memory_order_release);

		flags.mouseMove = true;
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


		if (ctrl == &kviewScaling)
		{
			state.viewScale = kviewScaling.getZeroBasedSelIndex<ViewScaling>();
			flags.viewChanged = true;
		}
		else if (ctrl == &kdisplayMode)
		{
			newc.displayMode.store(kdisplayMode.getZeroBasedSelIndex<DisplayMode>(), std::memory_order_release);
			setTransformOptions();
			flags.displayModeChange = true;
		}
		else if (ctrl == &kchannelConfiguration)
		{
			newc.configuration.store(kchannelConfiguration.getZeroBasedSelIndex<ChannelConfiguration>(), std::memory_order_release);
			flags.viewChanged = true;
		}
		else if (ctrl == &kwindowSize)
		{
			struct RetryResizer
			{
				RetryResizer(CSpectrum * h) : handle(h) {};
				CSpectrum * handle;
					
				void operator()()
				{
					auto currentCapacity = handle->audioStream.getAudioHistoryCapacity();
					if (currentCapacity > 0)
					{
						handle->setWindowSize(cpl::Math::round<std::size_t>(handle->kwindowSize.bGetValue() * currentCapacity));
					}
					else
					{
						GUIUtils::FutureMainEvent(200, RetryResizer(handle), handle);
					}
						
				}
			};
			RetryResizer(this)();
		}
		else if (ctrl == &kpctForDivision)
		{
			flags.frequencyGraphChange = true;
			flags.dynamicRangeChange = true;
			newc.divLimit = kpctForDivision.bGetValue();
		}
		else if (ctrl == &kalgorithm)
		{
			state.algo = kalgorithm.getZeroBasedSelIndex<TransformAlgorithm>();
			flags.resetStateBuffers = true;
			setTransformOptions();
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
		else if (ctrl == &kdspWin)
		{
			state.dspWindow = kdspWin.getParams().wType.load(std::memory_order_acquire);
			flags.windowKernelChange = true;
		}
		else if (ctrl == &kblobSize)
		{
			state.audioBlobSizeMs = Math::UnityScale::exp(kblobSize.bGetValue(), 0.5, 1000.0);
		}
		else if (ctrl == &kframeUpdateSmoothing)
		{
			state.bufferSmoothing = Math::UnityScale::linear(kframeUpdateSmoothing.bGetValue(), 0.0, 0.996);
		}
		else if (ctrl == &kbinInterpolation)
		{
			state.binPolation = kbinInterpolation.getZeroBasedSelIndex<BinInterpolation>();
		}
		else if (ctrl == &kfreeQ)
		{
			cresonator.setFreeQ(kfreeQ.bGetBoolState());
			flags.windowKernelChange = true;
		}
		else if (ctrl == &kspectrumStretching)
		{
			newc.stretch.store(cpl::Math::UnityScale::linear(ctrl->bGetValue(), 1.0, StretchMax), std::memory_order_release);
			if(state.displayMode == DisplayMode::ColourSpectrum)
				flags.resized = true;
		}
		else if (ctrl == &kfrequencyTracker)
		{
			newc.frequencyTrackingGraph.store(kfrequencyTracker.getZeroBasedSelIndex() + LineGraphs::None, std::memory_order_release);
		}
		else if (ctrl == &kprimitiveSize)
		{
			newc.primitiveSize.store(static_cast<float>(ctrl->bGetValue()), std::memory_order_release);
		}
		else if (ctrl == &kfloodFillAlpha)
		{
			newc.alphaFloodFill.store(static_cast<float>(ctrl->bGetValue()), std::memory_order_release);
		}
		else
		{
			for (int i = 0; i < numSpectrumColours; ++i)
			{
				if (ctrl == (kspecColours + i))
				{
					state.colourSpecs[i + 1] = kspecColours[i].getControlColourAsColour();
					return;
				}
				else if (ctrl == (kspecRatios + i))
				{
					// one ratio change affects the normalization, so recalculate all:
					calculateSpectrumColourRatios();
					return;
				}
			}

			for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
			{
				if (ctrl == &klines[i].colourOne)
				{
					state.colourOne[i] = klines[i].colourOne.getControlColourAsColour(); break;
				}
				else if	(ctrl == &klines[i].colourTwo)
				{
					state.colourTwo[i] = klines[i].colourTwo.getControlColourAsColour(); break;
				}
				else if (ctrl == &klines[i].decay)
				{
					lineGraphs[i].filter.setDecayAsFraction(ctrl->bGetValue(), 0.1); break;
				}
			}
		}
	}
	
	void CSpectrum::setTransformOptions()
	{
		if (state.algo == TransformAlgorithm::FFT)
		{
			// enable all windows
			for (std::size_t i = 0; i < (size_t)cpl::dsp::WindowTypes::End; i++)
			{
				kdspWin.getWindowList().setEnabledStateFor(i, true);
			}
		}
		else if(state.algo == TransformAlgorithm::RSNT)
		{
			// disable the windows unsupported by resonating algorithms
			for (std::size_t i = 0; i < (size_t)cpl::dsp::WindowTypes::End; i++)
			{
				if(cpl::dsp::windowHasFiniteDFT((cpl::dsp::WindowTypes)i))
					kdspWin.getWindowList().setEnabledStateFor(i, true);
				else
					kdspWin.getWindowList().setEnabledStateFor(i, false);
			}
		}

		if (newc.displayMode.load(std::memory_order_acquire) == DisplayMode::ColourSpectrum)
		{
			// disable all multichannel configurations
			for (std::size_t i = 0; i < (size_t)ChannelConfiguration::End; i++)
			{
				if (i > (size_t)ChannelConfiguration::OffsetForMono)
				{
					kchannelConfiguration.setEnabledStateFor(i, false);
				}
			}
		}
		else
		{
			// enable them.
			for (std::size_t i = 0; i < (size_t)ChannelConfiguration::End; i++)
			{
				kchannelConfiguration.setEnabledStateFor(i, true);
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
			{
				klowDbs.bSetInternal(cpl::Math::UnityScale::Inv::linear(actualDbs.low, kMinDbs, kMaxDbs));
				klowDbs.bRedraw();
			}


		}

		return false;
	}

	std::size_t CSpectrum::getWindowSize() const noexcept
	{
		return state.windowSize;
	}

	void CSpectrum::handleFlagUpdates()
	{
		cpl::CMutex audioLock;
		if (flags.internalFlagHandlerRunning)
			CPL_RUNTIME_EXCEPTION("Function is NOT reentrant!");

		flags.internalFlagHandlerRunning = true;
		bool firstRun = false;
		bool remapResonator = false;
		bool remapFrequencies = false;
		bool glImageHasBeenResized = false;

		// TODO: on numFilters change (and resizing of buffers), lock the working/audio buffers so that async processing doesn't corrupt anything.
		auto const sampleRate = getSampleRate();
		state.primitiveSize = newc.primitiveSize.load(std::memory_order_relaxed);
		state.alphaFloodFill = newc.alphaFloodFill.load(std::memory_order_relaxed);
		auto newconf = newc.configuration.load(std::memory_order_acquire);
		if (newconf != state.configuration)
		{
			state.configuration = newconf;
			if (newconf != ChannelConfiguration::Complex)
			{
				complexFrequencyGraph.clear();
			}
		}

		if (flags.displayModeChange.cas())
		{
			// ensures any concurrent processing modes gets to finish.
			audioLock.acquire(audioResource);
			state.displayMode = newc.displayMode.load(std::memory_order_acquire);
			flags.resized = true;
			flags.resetStateBuffers = true;
		}

		std::size_t axisPoints = state.displayMode == DisplayMode::LineGraph ? getWidth() : getHeight();

		if (axisPoints != state.axisPoints)
		{
			audioLock.acquire(audioResource);
			flags.resized = true;
			state.axisPoints = state.numFilters = axisPoints;
		}

		const std::size_t numFilters = getNumFilters();

		// TODO: insert messagemanagerlock or rework
		auto const divLimit = 5 + (state.configuration == ChannelConfiguration::Complex ? 0.25 : 1) * (numFilters * 0.02 + 0.5 * (numFilters * newc.divLimit.load(std::memory_order_relaxed)));
		auto const divLimitY = 5 + 0.6 * (getHeight() * newc.divLimit.load(std::memory_order_relaxed));

		oglImage.setFillColour(state.colourBackground);

		if (flags.firstChange.cas())
		{
			//framesPerUpdate = getOptimalFramesPerUpdate();
			flags.audioWindowWasResized = true;
			firstRun = true;
		}

		if (flags.initiateWindowResize.cas())
		{
			// we will get notified asynchronously in onAsyncChangedProperties.
			audioStream.setAudioHistorySize(newc.windowSize.load(std::memory_order_acquire));

		}
		if (flags.audioWindowWasResized.cas())
		{
			audioLock.acquire(audioResource);
			// TODO: rework this shit.

			auto current = audioStream.getAudioHistorySize();
			auto capacity = audioStream.getAudioHistoryCapacity();

			if (!firstRun)
			{
				cpl::GUIUtils::MainEvent(*this, 
					[=] 
					{
						if (capacity == 0)
							kwindowSize.bSetInternal(0);
						else
							kwindowSize.bSetInternal(double(audioStream.getAudioHistorySize()) / capacity);
						kwindowSize.bRedraw();
					}
				);

			}


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

			for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
			{
				lineGraphs[i].resize(numFilters); lineGraphs[i].zero();
			}


			workingMemory.resize(numFilters * 2 * sizeof(std::complex<double>));

			columnUpdate.resize(getHeight());
			// avoid doing it twice.
			if (!glImageHasBeenResized)
			{
				oglImage.resize(std::max<std::size_t>(1, cpl::Math::round<std::size_t>(getWidth() / newc.stretch.load(std::memory_order_acquire))), getHeight(), true);
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
			dbGraph.setDivisionLimit(0.4 * divLimitY);
			auto dynRange = getDBs();
			dbGraph.setLowerDbs(dynRange.low);
			dbGraph.setUpperDbs(dynRange.high);
			dbGraph.compileDivisions();
		}


		if (flags.viewChanged.cas())
		{
			if (state.configuration != ChannelConfiguration::Complex)
			{
				frequencyGraph.setBounds({ 0.0, (double)axisPoints });
				frequencyGraph.setView({ state.viewRect.left * axisPoints, state.viewRect.right * axisPoints });
				frequencyGraph.setMaxFrequency(sampleRate / 2);
				frequencyGraph.setScaling(state.viewScale == ViewScaling::Linear ? frequencyGraph.Linear : frequencyGraph.Logarithmic);
			}
			else
			{

				frequencyGraph.setBounds({ 0.0, axisPoints * 0.5 });
				frequencyGraph.setView({ state.viewRect.left * axisPoints, state.viewRect.right * axisPoints });
				//frequencyGraph.setView({ state.viewRect.left * axisPoints * 0.5, (state.viewRect.right - 0.5) * axisPoints * 0.5});
				frequencyGraph.setMaxFrequency(sampleRate / 2);
				frequencyGraph.setScaling(state.viewScale == ViewScaling::Linear ? frequencyGraph.Linear : frequencyGraph.Logarithmic);

				complexFrequencyGraph.setBounds({ 0.0, axisPoints * 0.5 });
				complexFrequencyGraph.setView({ (1 - state.viewRect.right) * axisPoints, (1 - state.viewRect.left) * axisPoints });
				
				//complexFrequencyGraph.setBounds({ axisPoints * 0.5, axisPoints * 1.0 });
				//complexFrequencyGraph.setView({ (1 - state.viewRect.right) * axisPoints, (1 - (state.viewRect.right - state.viewRect.left) * 0.5 - 0.5) * axisPoints });

				
				
				complexFrequencyGraph.setMaxFrequency(sampleRate / 2);
				complexFrequencyGraph.setScaling(state.viewScale == ViewScaling::Linear ? frequencyGraph.Linear : frequencyGraph.Logarithmic);
			}

			remapFrequencies = true;
			flags.frequencyGraphChange = true;

			if(state.displayMode == DisplayMode::ColourSpectrum && oldViewRect != state.viewRect)
				oglImage.freeLinearVerticalTranslation(oldViewRect, state.viewRect);

			oldViewRect = state.viewRect;

			audioLock.acquire(audioResource);
			for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
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
				case ViewScaling::Linear:
				{
					double halfSampleRate = sampleRate * 0.5;
					double complexFactor = state.configuration == ChannelConfiguration::Complex ? 2.0 : 1.0;
					double freqPerPixel = halfSampleRate / (numFilters - 1);

					for (std::size_t i = 0; i < numFilters; ++i)
					{
						mappedFrequencies[i] = static_cast<float>(complexFactor * state.viewRect.left * halfSampleRate + complexFactor * viewSize * i * freqPerPixel);
					}

					break;
				}
				case ViewScaling::Logarithmic:
				{
					double sampleSize = (numFilters - 1);

					double minFreq = state.minLogFreq;

					double end = sampleRate / 2;
					if (state.configuration != ChannelConfiguration::Complex)
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
			
		}
		if (flags.windowKernelChange.cas())
		{
			windowScale = kdspWin.generateWindow<fftType>(windowKernel, getWindowSize());
			remapResonator = true;
		}

		if (remapResonator)
		{
			audioLock.acquire(audioResource);
			auto window = kdspWin.getParams().wType.load(std::memory_order_acquire);
			cresonator.mapSystemHz(mappedFrequencies, mappedFrequencies.size(), cpl::dsp::windowCoefficients<fpoint>(window).second, sampleRate);
			flags.frequencyGraphChange = true;
			relayWidth = getWidth();
			relayHeight = getHeight();
		}

		if (flags.frequencyGraphChange.cas())
		{
			frequencyGraph.setDivisionLimit(divLimit);
			frequencyGraph.compileGraph();
			if (state.configuration == ChannelConfiguration::Complex)
			{
				complexFrequencyGraph.setDivisionLimit(divLimit);
				complexFrequencyGraph.compileGraph();
			}
		}

		if (flags.resetStateBuffers.cas())
		{
			audioLock.acquire(audioResource);
			cresonator.resetState();
			for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
				lineGraphs[i].zero();
			std::memset(audioMemory.data(), 0, audioMemory.size() /* * sizeof(char) */);
			std::memset(workingMemory.data(), 0, workingMemory.size() /* * sizeof(char) */);
		}

		// reset all flags through value-initialization
		flags.internalFlagHandlerRunning = false;
	}

	std::size_t CSpectrum::getValidWindowSize(std::size_t in) const noexcept
	{
		std::size_t n = std::min(audioStream.getAudioHistoryCapacity(), in);
		n -= (n & 0x7); // must be a multiple of 8, due to vectorization
		if (n < 16)
			n = 0;
		return n;
	}

	void CSpectrum::setWindowSize(std::size_t size)
	{
		newc.windowSize.store(getValidWindowSize(size), std::memory_order_release);
		flags.initiateWindowResize = true;
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
		else if (ctrl == &kwindowSize)
		{
			auto bufLength = cpl::Math::round<int>(value * audioStream.getAudioHistoryCapacity());
			sprintf(buf, "%d smps", bufLength);
			buffer = buf;
			return true;
		}
		else if (ctrl == &kblobSize)
		{
			sprintf_s(buf, "%.2f ms", Math::UnityScale::exp(value, 0.50, 1000.0));
			buffer = buf;
			return true;
		}
		else if (ctrl == &kframeUpdateSmoothing)
		{
			sprintf_s(buf, "%.2f", Math::UnityScale::linear(kframeUpdateSmoothing.bGetValue(), 0.0, 0.996));
			buffer = buf;
			return true;
		}
		else if (ctrl == &kspectrumStretching)
		{
			sprintf_s(buf, "%.2fx", cpl::Math::UnityScale::linear(ctrl->bGetValue(), 1.0, StretchMax));
			buffer = buf;
			return true;
		}
		else if (ctrl == &kprimitiveSize)
		{
			sprintf(buf, "%.2f pts", ctrl->bGetValue() * primitiveMaxSize);
			buffer = buf;
			return true;
		}
		else
		{
			for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
			{
				if (ctrl == &klines[i].decay)
				{
					sprintf_s(buf, "%.2f dB/s", 20 * std::log10(value));
					buffer = buf;
					return true;
				}
			}
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
		else if (ctrl == &kwindowSize)
		{
			if (buffer.find("ms") != std::string::npos)
			{
				if (cpl::lexicalConversion(buffer, newVal))
				{
					value = cpl::Math::confineTo(newVal / (1000 * audioStream.getAudioHistoryCapacity() / audioStream.getInfo().sampleRate), 0.0, 1.0);
					return true;
				}
			}
			else
			{
				if (cpl::lexicalConversion(buffer, newVal))
				{
					value = cpl::Math::confineTo(newVal / audioStream.getAudioHistoryCapacity(), 0.0, 1.0);
					return true;
				}
			}
		}
		else if (ctrl == &kblobSize)
		{
			if (cpl::lexicalConversion(buffer, newVal))
			{
				value = Math::confineTo(Math::UnityScale::Inv::exp(newVal, 0.5, 1000.0), 0.0, 1.0);
				return true;
			}
		}
		else if (ctrl == &kframeUpdateSmoothing)
		{
			if (cpl::lexicalConversion(buffer, newVal))
			{
				value = Math::UnityScale::Inv::linear(Math::confineTo(newVal, 0.0, 1.0), 0.0, 0.996);
				return true;
			}
		}
		else if (ctrl == &kspectrumStretching)
		{
			if (cpl::lexicalConversion(buffer, newVal))
			{
				value = Math::UnityScale::Inv::linear(Math::confineTo(newVal, 1.0, StretchMax), 1.0, StretchMax);
				return true;
			}
		}
		else if (ctrl == &kprimitiveSize)
		{
			if (cpl::lexicalConversion(buffer, newVal))
			{
				value = cpl::Math::confineTo(newVal / primitiveMaxSize, 0.0, 1.0);
				return true;
			}
		}
		else
		{
			for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
			{
				if (ctrl == &klines[i].decay)
				{
					if (cpl::lexicalConversion(buffer, newVal))
					{
						value = cpl::Math::confineTo(std::pow(10, newVal / 20), 0.0, 1.0);
						return true;
					}
				}
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

	void CSpectrum::onAsyncChangedProperties(const AudioStream & source, const AudioStream::AudioStreamInfo & before)
	{
		flags.audioWindowWasResized = true;
	}

	void CSpectrum::resetStaticViewAssumptions()
	{
		lastPeak = -1;
	}

};
