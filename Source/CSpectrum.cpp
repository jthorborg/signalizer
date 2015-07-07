
#include "CSpectrum.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include "SignalizerDesign.h"
#include <cpl/rendering/OpenGLRasterizers.h>
#include <cpl/simd.h>
#include <cpl/LexicalConversion.h>

namespace Signalizer
{

	static const std::vector<std::string> ViewScaleNames = { "Linear", "Logarithmic" };
	static const std::vector<std::string> AlgorithmNames = { "FFT", "Resonator" };
	static const std::vector<std::string> ChannelNames = { "Left", "Right", "Merge", "Phase", "Separate" };
	static const std::vector<std::string> DisplayModeNames = { "Line graph", "Colour spectrum" };
	// the minimum level of dbs to display
	static const double kMinDbs = -24 * 16;
	// the maximum level of dbs to display
	static const double kMaxDbs = 24 * 4;

	std::unique_ptr<juce::Component> CSpectrum::createEditor()
	{
		auto content = new Signalizer::CContentPage();
		if (auto page = content->addPage("Settings", "icons/svg/wrench.svg"))
		{
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kviewScaling, 0);
				section->addControl(&kalgorithm, 1);
				section->addControl(&kchannelConfiguration, 0);
				section->addControl(&kdspWindow, 0);
				section->addControl(&kwindowSize, 1);
				section->addControl(&kdisplayMode, 1);
				section->addControl(&klowDbs, 1);
				section->addControl(&khighDbs, 0);

				section->addControl(&kdecayRate, 0);
				section->addControl(&kpctForDivision, 0);
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
				
				page->addSection(section);
			}
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
		frequencyGraph({0, 1}, {0, 1}, 1, 10),
		flags()
	{
		setOpaque(true);
		processorSpeed = juce::SystemStats::getCpuSpeedInMegaherz();
		initPanelAndControls();

		state.viewRect = { 0.0, 1.0 }; // default full-view

		state.displayMode = DisplayMode::ColourSpectrum;

		listenToSource(data[0]);

		state.minLogFreq = 10;
		peakFilter.setSampleRate((float)getSampleRate());
		setWindowSize(200);
		flags.firstChange = true;

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
	}

	CSpectrum::~CSpectrum()
	{

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
		if (low > (high - 3))
		{
			low = high - 3;
		}
		state.dynRange.low = low;
		state.dynRange.high = high;
		if (updateControls)
		{
			klowDbs.bSetValue(cpl::Math::UnityScale::Inv::linear<float>(low, kMinDbs, kMaxDbs));
			khighDbs.bSetValue(cpl::Math::UnityScale::Inv::linear<float>(high, kMinDbs, kMaxDbs));
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

		kwindowSize.bAddFormatter(this);
		kdecayRate.bAddFormatter(this);
		klowDbs.bAddFormatter(this);
		khighDbs.bAddFormatter(this);

		// ------ titles -----------
		kviewScaling.bSetTitle("Graph scale");
		kalgorithm.bSetTitle("Transform algorithm");
		kchannelConfiguration.bSetTitle("Channel configuration");
		kdisplayMode.bSetTitle("Display mode");
		kdspWindow.bSetTitle("Window kernel");

		klowDbs.bSetTitle("Lower limit");
		khighDbs.bSetTitle("Upper limit");
		kwindowSize.bSetTitle("Window size");
		kdecayRate.bSetTitle("Filter decay rate");

		kline1Colour.bSetTitle("Graph 1 colour");
		kline2Colour.bSetTitle("Graph 2 colour");
		kgridColour.bSetTitle("Grid colour");
		kpctForDivision.bSetTitle("Grid div. space");


		// ------ content ----------
		kviewScaling.setValues(ViewScaleNames);
		kalgorithm.setValues(AlgorithmNames);
		kchannelConfiguration.setValues(ChannelNames);
		kdisplayMode.setValues(DisplayModeNames);

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


		klowDbs.bSetDescription("The lower limit of the displayed dynamic range.");
		khighDbs.bSetDescription("The upper limit of the displayed dynamic range");
		kwindowSize.bSetDescription("The window size of the audio data, affects time/frequency resolution.");
		kdecayRate.bSetDescription("Allows the filters to decay more slowly, but still reacting to peaks.");
		kline1Colour.bSetDescription("The colour of the first graph.");
		kline2Colour.bSetDescription("The colour of the second graph.");
		kgridColour.bSetDescription("The colour of the dB/frequency grid.");
		kpctForDivision.bSetDescription("The minimum amount of free space that triggers a recursed frequency grid division; smaller values draw more frequency divisions.");

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
	}

	void CSpectrum::load(cpl::CSerializer::Builder & builder, long long int version)
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
		double newPos(0);

		switch (state.displayMode)
		{
			case DisplayMode::ColourSpectrum:
				newPos = double(getAxisPoints() - event.position.y) / getAxisPoints();
				break;
			case DisplayMode::LineGraph:
				newPos = double(event.position.x) / getAxisPoints();
				break;
		}


		// delta difference, scales
		volatile auto delta = state.viewRect.left - state.viewRect.right;
		volatile auto inc = delta * wheel.deltaY / 5;
		state.viewRect.left += newPos * inc;
		state.viewRect.right -= (1 - newPos) * inc;

		state.viewRect.left = std::max(0.0, state.viewRect.left);
		state.viewRect.right = std::min(state.viewRect.right, 1.0);

		flags.viewChanged = true;
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
			
			auto delta = state.viewRect.left - state.viewRect.right;
			auto inc = (delta * (state.displayMode == DisplayMode::LineGraph ? mouseDelta.x : mouseDelta.y)) / getAxisPoints();

			state.viewRect.left += inc;
			state.viewRect.right += inc;

			state.viewRect.left = std::max(0.0, state.viewRect.left);
			state.viewRect.right = std::min(state.viewRect.right, 1.0);
			lastMousePos = event.position;
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
			state.displayMode = kdisplayMode.getZeroBasedSelIndex<DisplayMode>();
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
			state.gridColour = kgridColour.getControlColourAsColour();
		}
		else if (ctrl == &kdspWindow)
		{
			state.dspWindow = kdspWindow.getZeroBasedSelIndex<cpl::dsp::WindowTypes>();
			flags.windowKernelChange = true;
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
				klowDbs.bSetValue(cpl::Math::UnityScale::Inv::linear(actualDbs.low, kMinDbs, kMaxDbs), true);

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
		bool graphShouldBeCompiled = false;

		auto const numFilters = getNumFilters();
		auto const sampleRate = getSampleRate();
		auto const axisPoints = getAxisPoints();

		if (flags.audioWindowResize)
		{

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
		if (flags.audioMemoryResize)
		{
			const auto bufSize = cpl::Math::nextPow2(state.windowSize);
			audioMemory.resize(bufSize * sizeof(std::complex<double>));
			windowKernel.resize(bufSize);
			flags.windowKernelChange = true;
		}

		if (flags.resized)
		{
			filterResults.resize(numFilters);
			filterStates.resize(numFilters);
			workingMemory.resize(numFilters * 2 * sizeof(std::complex<double>));

			oglImage.resize(getWidth(), getHeight(), true);
			columnUpdate.resize(getHeight() * 4);

			graphShouldBeCompiled = true;

			flags.viewChanged = true;
			remapFrequencies = true;
		}

		if (flags.viewChanged)
		{
			frequencyGraph.setBounds({ 0.0, (double)getAxisPoints() });
			frequencyGraph.setView({ state.viewRect.left * axisPoints, state.viewRect.right * axisPoints });
			frequencyGraph.setMaxFrequency(sampleRate / 2);
			frequencyGraph.setScaling(state.viewScale == ViewScaling::Linear ? frequencyGraph.Linear : frequencyGraph.Logarithmic);

			remapFrequencies = true;
			graphShouldBeCompiled = true;
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
				}
			}
			remapResonator = true;
			
		}

		if (remapResonator)
		{
			cresonator.mapSystemHz(mappedFrequencies, mappedFrequencies.size(), sampleRate);
			graphShouldBeCompiled = true;
		}

		if (graphShouldBeCompiled || flags.frequencyGraphChange)
		{
			frequencyGraph.setDivisionLimit(5 + getNumFilters() * 0.05 + 0.95 * (getNumFilters() * kpctForDivision.bGetValue()));
			frequencyGraph.compileGraph();
		}

		if (flags.windowKernelChange)
		{
			computeWindowKernel();
		}

		// reset all flags through value-initialization
		flags = Flags();

	}

	void CSpectrum::setWindowSize(std::size_t size)
	{
		std::size_t n = std::min(audioStream[0].maxSize(), size);
		n -= (n & 0x7); // must be a multiple of 8, due to vectorization
		if (n < 16)
			n = 16;


		state.windowSize = n;
		flags.audioWindowResize = true;


	}

	bool CSpectrum::valueToString(const cpl::CBaseControl * ctrl, std::string & buffer, cpl::iCtrlPrec_t value)
	{
		char buf[200];
		if (ctrl == &klowDbs || ctrl == &khighDbs)
		{
			auto val = cpl::Math::UnityScale::linear<double>(value, kMinDbs, kMaxDbs);
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
		
		return false;
	}


	bool CSpectrum::stringToValue(const cpl::CBaseControl * ctrl, const std::string & buffer, cpl::iCtrlPrec_t & value)
	{
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
