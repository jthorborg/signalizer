#include "CSpectrum.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Graphics.h>
#include <cpl/Mathext.h>
#include <cpl/ffts.h>
#include <cpl/dsp.h>
#include <cpl/dsp/CPeakFilter.h>
#include <cpl/resampling.h>
#include <cpl/Utility.h>
#include <vector>

const bool showPhase = true;



namespace Signalizer
{
	// swapping the right channel might give an more intuitive view
	#define SWAP_RIGHT

	// the minimum level of dbs to display
	const double kMinDbs = -24 * 16;
	// the maximum level of dbs to display
	const double kMaxDbs = 24 * 4;

	static const juce::SoftwareImageType sit;

	struct CSerializedVectorData
	{
		union
		{
			struct
			{
				float fillLines, antiAlias, fadeHistory, sync, rgbHistory;
			};
			float params[5];
		} renderButtons;
		float viewType;
		bool panelClosed;
		float kWindow, kRotation, kGain, kGraph, yOffset, xOffset, kcolor;

	};

	static CSerializedVectorData DefaultSettings
	{
		{
			{ 
				1.0f,	// draw lines
				0.0f,	// no anti alias
				1.0f,	// fade history
				0.f,	// dont synchronize
				0.f		// dont fade rgb
			}			
		},
		0.0f,			// Spectral view
		false,			// show panel
		0.100f,			// 215 ms buffer
		135.0f/360.0f,	// 225 degrees rotation
		0.5f,			// no gain
		0.45f,			// medium graph visibility
		0.0f,			// no vertical offset
		0.0f,			// no horizontal offset
		0.481001705		// purple color

	};

	bool CSpectrum::restore(const void * data, std::size_t size)
	{
		const CSerializedVectorData * se = reinterpret_cast<const CSerializedVectorData * >(data);
		if (se && size == sizeof(CSerializedVectorData))
		{

			kwindow.bSetValue(se->kWindow);
			kgain.bSetValue(se->kGain);
			kgraph.bSetValue(se->kGraph);
			kviewType->bSetValue(se->viewType);
			kstdColor.bSetValue(se->kcolor);
			yoffset = se->yOffset;
			xoffset = se->xOffset;
			panel.bSetValue(se->panelClosed ? 1.0f : 0.0f);
			kdecay.bSetValue(0.25);
			setDBs(-120, 6, true);
			return true;
		}
		else
		{
			return restore(&DefaultSettings, sizeof(CSerializedVectorData));
		}
		return false;

	}

	bool CSpectrum::serialize(juce::MemoryBlock & data)
	{
		CSerializedVectorData se;
		se.kWindow = kwindow.bGetValue();
		se.kGain = kgain.bGetValue();
		se.xOffset = xoffset;
		se.yOffset = yoffset;
		se.kGraph = kgraph.bGetValue();
		se.kcolor = kstdColor.bGetValue();
		se.panelClosed = panel.bGetValue() > 0.5 ? true : false;
		se.viewType = kviewType->bGetValue();
		//for (unsigned i = 0; i < renderButtons->getNumButtons(); ++i)
		//	se.renderButtons.params[i] = renderButtons->getButton(i).bGetValue();

		data.append(&se, sizeof(se));
		
		return false;
	}

	CSpectrum::CSpectrum(AudioBuffer & data)
		:
		audioData(data),
		windowSize(1000),
		kwindow("Window Size"),
		kgain("Gain"),
		kgraph("Graph intensity"),
		kstdColor("Std. Colour"),
		kdecay("Decay Amount", cpl::CKnobEx::type::db),
		klowDbs("Lower Range"), khighDbs("Upper Range"),
		kaux("Auxiliary"), kaux2("Aux2"),
		processorSpeed(0),
		renderWindow(*this),
		lastFrameTick(0),
		kviewType(nullptr),
		kscaleButtons(nullptr),
		koptions(nullptr),
		firstResize(true),
		isFrozen(false),
		isXLog(false),
		isYLog(false),
		ydrag(0),
		xdrag(0),
		debug(false),
		expViewScale(0),
		mouseListener(*this),
		numChannels(2),
		transformer(data[0].sampleRate, transformer.vectorized),
		resonator(),
		qSettings(Quality::Bound),
		imagePtr(0),
		viewType(ViewTypes::SpectralView),
		sdft()
	{



		textbuf = std::unique_ptr<char>(new char[300]);

		if (e)
			e->onViewConstruction(this);

		processorSpeed = juce::SystemStats::getCpuSpeedInMegaherz();
		addAndMakeVisible(renderWindow);
		flt.setDecayAsDbs(-60);
		initPanelAndControls();
		listenToSource(audioData[0]);
	}

	void CSpectrum::buttonSelected(cpl::CBaseControl * c, int index)
	{
		if (c == kscaleButtons)
		{
			switch (index)
			{
			case ScaleButtons::X:
				isXLog = true;
				expViewScale = (float)isXLog;
				break;
			case ScaleButtons::Y:
				isYLog = true;
				break;
			}

		}
		else if (c == koptions)
		{
			auto flags = transformer.getFlags();
			switch (index)
			{
			case OptionButtons::threaded:
				transformer.setFlags(flags | transformer.threaded);
				break;
			case OptionButtons::parallel:
				transformer.setFlags(flags | transformer.accelerated);
				break;
			}

		}
		else if (c == kalgorithm)
		{
			algorithmType = cpl::distribute<Algorithm>(kalgorithm->bGetValue());
		}
		else if (c == kviewType)
		{
			viewType = cpl::distribute<ViewTypes>(kviewType->bGetValue());
			viewChanged();
		}

	}



	void CSpectrum::buttonDeselected(cpl::CBaseControl * c, int index)
	{
		if (c == kscaleButtons)
		{
			switch (index)
			{
			case ScaleButtons::X:
				isXLog = false;
				expViewScale = (float)isXLog;
				break;
			case ScaleButtons::Y:
				isYLog = false;
				break;
			}

		}
		else if (c == koptions)
		{
			auto flags = transformer.getFlags();
			switch (index)
			{
			case OptionButtons::threaded:
				transformer.setFlags(flags & ~transformer.threaded);
				break;
			case OptionButtons::parallel:
				transformer.setFlags(flags & ~transformer.accelerated);
				break;
			}
		}
	}

	void CSpectrum::initPanelAndControls()
	{
		// lel

		kaux.bSetValue(0.12);


		addAndMakeVisible(panel);
		kwindow.bSetListener(this);
		kgain.bSetListener(this);
		kdecay.bSetListener(this);
		klowDbs.bSetListener(this);
		khighDbs.bSetListener(this);
		kaux.bSetListener(this);
		kaux2.bSetListener(this);
		renderWindow.addMouseListener(&mouseListener, false);

		panel.setName("Spectrum properties...");
		panel.setOrientation(cpl::CControlContainer::Orientation::bottom);
		panel.setBounds(0, getHeight() - 100, getWidth(), 100);
		panel.bSetListener(this);
		panel.setNested(false);
		// create render buttons
		kviewType = new cpl::CButtonGroup<cpl::CRenderButton>
		(
			{ "Spectral View", "Spectrogram" },
			this,
			cpl::CButtonGroup<cpl::CRenderButton>::Behaviour::radio
		);

		// create option buttons
		koptions = new cpl::CButtonGroup<cpl::CRenderButton>
		(
			{ "Fill", "Threaded", "Parallel" },
			this,
			cpl::CButtonGroup<cpl::CRenderButton>::Behaviour::polyToggle
		);
		// create channel configuration
		kchannelConf = new cpl::CButtonGroup<cpl::CRenderButton>
		(
			{ "Left", "Right", "Merge", "Phase Diff", "Seperate"},
			nullptr,
			cpl::CButtonGroup<cpl::CRenderButton>::Behaviour::radio
		);
		// create algorithm buttons
		kalgorithm = new cpl::CButtonGroup<cpl::CRenderButton>
		(
			{ "FFT", "Resonator", "Minimum Q DFT" },
			this
		);
		kalgorithm->setText("Algorithm");
		kalgorithm->bSetValue(1);
		kchannelConf->setText("Channels");
		kviewType->setText("View Types");
		koptions->setText("Options...");
		// add rendering group
		auto group = new cpl::CControlGroup();
		group->setText("Rendering");
		group->addControl(kchannelConf, true);
		group->addControl(&kaux);
		group->addControl(&kaux2);
		group->addControl(kalgorithm, true);
		group->addControl(kviewType, true);
		group->addControl(koptions, true);
		panel.addControl(group, true);

		group = new cpl::CControlGroup();
		group->setText("Utility");
		group->addControl(&kdecay);
		group->addControl(&klowDbs);
		group->addControl(&khighDbs);
		group->addControl(&kwindow);
		panel.addControl(group, true);

		group = new cpl::CControlGroup();
		group->setText("Graph");

		kscaleButtons = new cpl::CButtonGroup<cpl::CRenderButton>(
			{ "X: Linear", "Y: Linear" }, 
			this,
			cpl::CButtonGroup<cpl::CRenderButton>::Behaviour::polyToggle
		);
		kscaleButtons->getButton(ScaleButtons::X).setToggledText("X: Logarithmic");
		kscaleButtons->getButton(ScaleButtons::Y).setToggledText("Y: Logarithmic");
		group->addControl(kscaleButtons, true);
		group->addControl(&kgraph);
		group->addControl(&kstdColor);
		panel.addControl(group, true);

		panel.resizeAccordingly();
		renderWindow.setMouseCursor(juce::MouseCursor::DraggingHandCursor);
	}

	void CSpectrum::rearrange(int width, int height, bool callback)
	{
		if (!width || !height)
			return;
		if (firstResize)
		{
			//initPanelAndControls();
			firstResize = false;
		}
		// new code

		if(callback)
			panel.setBounds(0, getHeight() - panel.getHeight(), getWidth(), panel.getHeight());

		displaySize.setBounds(0, 0, getWidth(), getHeight() - panel.getHeight());
		renderWindow.setBounds(displaySize);
		lowestCoord = 0;
		highestCoord = viewType == ViewTypes::SpectralView ? displaySize.getWidth() : displaySize.getHeight();

		viewChanged();

		debug = true;

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
		pixel[0] = blue;
		pixel[1] = green;
		pixel[2] = red;
		// set red

		// saturate


	}
	void CSpectrum::viewChanged()
	{
		numFilters = viewType == ViewTypes::SpectralView ? displaySize.getWidth() : displaySize.getHeight();
		reallocBuffers();
		if (viewType == ViewTypes::SpectroGram && spectrumImg.getBounds() != displaySize)
		{
			spectrumImg = Image(Image::PixelFormat::RGB, displaySize.getWidth(), displaySize.getHeight(), true, sit);
			//spectrumGraphics = std::unique_ptr<juce::Graphics>(new juce::Graphics(spectrumImg));
		}
		imagePtr = 0;


		rezoomed();
	}

	void CSpectrum::reallocBuffers()
	{
		auto const channelConfiguration = cpl::distribute<ChannelConfiguration>(kchannelConf->bGetValue());
		if (channelConfiguration > ChannelConfiguration::Merge)
			numChannels = 2;
		else
			numChannels = 1;
		mappedFrequencies.resize(numFilters);
		filterResults.resize(numFilters);
		filterStates.resize(numFilters);
		workspace.resize(numFilters * 2 * sizeof(std::complex<double>));
	}









	void CSpectrum::cbMouseWheelMove(const MouseEvent& event, const MouseWheelDetails& wheel)
	{
		float newPos;

		switch (viewType)
		{
			case ViewTypes::SpectroGram:
				newPos = float(getHeight() - event.y) / getHeight();
				break;
			case ViewTypes::SpectralView:
				newPos = float(event.x) / getWidth();
				break;
		}


		/* djd */
		auto delta = highestCoord - lowestCoord;
		auto inc = delta * wheel.deltaY / 5;
		lowestCoord += newPos * inc;
		highestCoord -= (1 - newPos) * inc;

		rezoomed();
	}
	void CSpectrum::rezoomed()
	{
		mapFrequencies();

	}

	void CSpectrum::cbMouseDoubleClick(const MouseEvent& event)
	{
		lowestCoord = xoffset = yoffset = xdrag = ydrag = 0;
		highestCoord = numFilters;
		rezoomed();
	}
	void CSpectrum::cbMouseDrag(const MouseEvent& event)
	{
		switch (kviewType->getToggledIndex())
		{
		case ViewTypes::SpectralView:
			xdrag = event.getDistanceFromDragStartX() / 500.0f;
			break;
		case ViewTypes::SpectroGram:
			xdrag = event.getDistanceFromDragStartY() / 500.0f;
			break;
		}

	}
	void CSpectrum::cbMouseUp(const MouseEvent& event)
	{
		if (event.mods.testFlags(ModifierKeys::leftButtonModifier))
		{
			xoffset += xdrag;
			yoffset += ydrag;
			xdrag = ydrag = 0.0f;
		}
		if (event.mods.testFlags(ModifierKeys::rightButtonModifier))
		{
			unfreeze();
		}
	}
	void CSpectrum::cbMouseDown(const MouseEvent& event)
	{
		if (event.mods.testFlags(ModifierKeys::rightButtonModifier))
		{
			freeze();
		}
	}







	void CSpectrum::renderNormal(juce::Graphics & g)
	{

		auto cStart = cpl::Misc::ClockCounter();
		auto fps = 1.0 / (avgFps.getAverage() / juce::Time::getHighResolutionTicksPerSecond());
		flt.setSampleRate(fps);



		auto gain = getGain();
		auto height = getHeight() - panel.getHeight();


		const bool frameLag = false;


		prepareTransform();
		auto pStop = cpl::Misc::ClockCounter();
		doTransform();
		auto cCalcstop = cpl::Misc::ClockCounter();

		mapToLinearSpace();
		auto cMapstop = cpl::Misc::ClockCounter();
		g.fillAll(juce::Colours::black);
		g.setColour(juce::Colours::green);


		auto const color1 = juce::Colours::purple;
		auto const color2 = juce::Colours::green;

		const bool fill = (koptions->getButton(OptionButtons::drawCurve).bGetValue() > 0.5);


		auto const w = getWidth();
		auto const channelConfiguration = cpl::distribute<ChannelConfiguration>(kchannelConf->bGetValue());
		auto const mode = cpl::distribute<ViewTypes>(kviewType->bGetValue());
		auto const bottom = height + panel.getHeight();
		switch (mode)
		{
		case ViewTypes::SpectroGram:
		{
			g.setOpacity(1.0f);

			auto const ilen = spectrumImg.getBounds().getWidth();
			auto const iheight = spectrumImg.getBounds().getHeight();
			
			juce::Image::BitmapData pixels(spectrumImg, juce::Image::BitmapData::ReadWriteMode::writeOnly);

			for (int i = 0; i < iheight; ++i)
			{

				auto ptr = pixels.getPixelPointer(imagePtr, iheight - (i + 1));
				/*val_typeof(*ptr) value = 0;
				
				if (filterResults[i].magnitude > 0)
					value = filterResults[i].magnitude * 0xFF;
				*ptr = value;*/

				//if (filterResults[i].magnitude > 0)
					ColorScale(ptr, filterResults[i].magnitude);
				
				//*ptr = value;
			}

			imagePtr++;
			imagePtr %= ilen;
			// draw first half
			g.drawImage(spectrumImg, 0, 0, ilen - imagePtr, iheight, imagePtr, 0, ilen - imagePtr, iheight);
			// draw second half
			g.drawImage(spectrumImg, ilen - imagePtr, 0, imagePtr, iheight, 0, 0, imagePtr, iheight);

			break;
		}
		case ViewTypes::SpectralView:
		{
			switch (channelConfiguration)
			{
			case ChannelConfiguration::Left:
			case ChannelConfiguration::Right:
			case ChannelConfiguration::Merge:
				g.setColour(color1);
				if (fill)
				{
					for (int x = 0; x < w; x++)
					{
						auto const mag = filterResults[x].magnitude;
						g.drawVerticalLine(x, height - height * mag, bottom);
					}
				}
				else
				{
					float old = filterResults[0].magnitude;
					for (int x = 0; x < w; x++)
					{
						auto const mag = height - height * filterResults[x].magnitude;
						g.drawLine(x - 1, mag, x, old);
						//g.drawVerticalLine(x, std::min(old, mag), std::max(old, mag));
						old = mag;
					}
				}
				break;
			case ChannelConfiguration::Separate:
				if (fill)
				{
					g.setColour(color1);
					for (int x = 0; x < w; x++)
					{
						auto const mag = filterResults[x].leftMagnitude;
						g.drawVerticalLine(x, height - height * mag, bottom);
					}
					g.setColour(color2);
					for (int x = 0; x < w; x++)
					{
						auto const mag = filterResults[x].rightMagnitude;
						g.drawVerticalLine(x, height - height * mag, bottom);
					}
				}
				else
				{
					g.setColour(color1);
					float old = filterResults[0].magnitude;
					for (int x = 0; x < w; x++)
					{
						auto const mag = height - height * filterResults[x].leftMagnitude;
						g.drawLine(x - 1, mag, x, old);
						//g.drawVerticalLine(x, std::min(old, mag), std::max(old, mag));
						old = mag;
					}
					g.setColour(color2);
					old = filterResults[0].magnitude;
					for (int x = 0; x < w; x++)
					{
						auto const mag = height - height * filterResults[x].rightMagnitude;
						g.drawLine(x - 1, mag, x, old);
						//g.drawVerticalLine(x, std::min(old, mag), std::max(old, mag));
						old = mag;
					}
				}
				break;
			case ChannelConfiguration::Phase:
				cpl::Graphics::UPixel c1 = color1.getPixelARGB().getARGB();
				cpl::Graphics::UPixel c2 = color2.getPixelARGB().getARGB();
				cpl::Graphics::UPixel t0(0);
				if (fill)
				{
					for (int x = 0; x < w; x++)
					{

						auto const phase = filterResults[x].phase;
						t0 = c1 * phase;
						t0 += c2 * (1 - phase);
						//std::uint32_t newColor = static_cast<std::uint32_t>(intc1 * phase + intc2 * (1 - phase));
						//g.setColour(juce::Colour(newColor));
						g.setColour(juce::Colour(t0.p));
						auto const mag = height - height * filterResults[x].magnitude;
						g.drawVerticalLine(x, mag, bottom);
					}
				}
				else
				{
					float old = filterResults[0].magnitude;
					g.setColour(color1);
					for (int x = 0; x < w; x++)
					{
						auto const mag = height - height * filterResults[x].magnitude;
						g.drawLine(x - 1, mag, x, old);
						//g.drawVerticalLine(x, std::min(old, mag), std::max(old, mag));
						old = mag;
					}
					old = filterResults[0].phase;
					g.setColour(color2);
					for (int x = 0; x < w; x++)
					{
						auto const mag = height * filterResults[x].phase * filterResults[x].magnitude;
						g.drawLine(x - 1, mag, x, old);
						//g.drawVerticalLine(x, std::min(old, mag), std::max(old, mag));
						old = mag;
					}

				}
				break;

			}
			break;
		}
		}


		auto cyclesNow = cpl::Misc::ClockCounter();
		auto totalCycles = renderCycles + cyclesNow - cStart;
		auto preCycles = pStop - cStart;
		auto calcCycles = cCalcstop - pStop;
		auto mapCycles = cMapstop - cCalcstop;
		double cpuTime = (double(totalCycles) / (processorSpeed * 1000 * 1000) * 100) * fps;
		g.setColour(juce::Colours::blue);
		sprintf(textbuf.get(), "%dx%d -> (%f, %f): %.1f fps - %.1f%% cpu, pre: %.1f, calc: %.1f, map: %.1f",
			displaySize.getWidth(), displaySize.getHeight(), lowestCoord, highestCoord, fps, cpuTime, 
			100 * double(preCycles) / totalCycles, 
			100 * double(calcCycles) / totalCycles, 
			100 * double(mapCycles) / totalCycles);
		g.drawSingleLineText(textbuf.get(), 10, 20);
			

	}

	cpl::Graphics::RGBPixel CSpectrum::getSTDColor()
	{
		cpl::Graphics::RGBPixel ret;
		ret.setColour(this->kstdColor.bGetValue());
		return ret;
	}

	void CSpectrum::renderOpenGL()
	{
		auto cStart = cpl::Misc::ClockCounter();

		renderCycles = cpl::Misc::ClockCounter() - cStart;

		auto tickNow = juce::Time::getHighResolutionTicks();
		avgFps.setNext(tickNow - lastFrameTick);
		lastFrameTick = tickNow;
	}

	void CSpectrum::paint(Graphics & g)
	{
		g.setColour(juce::Colours::black);
		g.fillAll();
	}

	void CSpectrum::freeze()
	{
		isFrozen = true;
		/*for (unsigned i = 0; i < ownData.size(); ++i)
		{
			audioData[i].clone(ownData[i]);
		}*/
	}

	void CSpectrum::unfreeze()
	{
		isFrozen = false;
	}

	double CSpectrum::getGain()
	{
		float val = kgain.bGetValue();
		if (val > 0.5)
		{
			return 1 + (val * 2 - 1) * 9;

		}
		else if (val < 0.5)
		{
			return val * 2;
		}

		return 1;

	}
	int CSpectrum::getWindowSize()
	{
		std::size_t n = cpl::Math::round<int>(kwindow.bGetValue() * audioData[0].maxSize());
		n -= (n & 0x7); // must be a multiple of 7, due to vectorization
		if (n < 16)
			n = 16;
		return n;
	}

	bool CSpectrum::setFullScreenMode(bool toggle)
	{
		panel.bSetValue(1);
		isFullScreen = toggle;
		//resized();
		return true;
	}
	CSpectrum::~CSpectrum()
	{

		if (e)
			e->onViewDestruction(this);

	}

	/*void g()
	{
		vector<complex<double>> signal(size);
		for (int i = 0; i < size; ++i)
			signal[i] = complex<double>(left[i], right[i]);

		fft(signal, size);

		complex<double> leftBin = signal[bin] + conj(signal[size - bin]) / 2;
		complex<double> rightBin = -(signal[bin] - conj(signal[size - bin]) / 2;

	}*/



	bool CSpectrum::valueChanged(cpl::CBaseControl * ctrl)
	{
		

		char buf[200];
		if (ctrl == &kgain)
		{
			sprintf(buf, "%.2f dB (%d%%)", cpl::Math::fractionToDB(getGain()), (int)cpl::Misc::Round(getGain() * 100));
			kgain.bSetText(buf);
			return true;
		}
		else if (ctrl == &kaux)
		{

			auto const qDBs = cpl::Math::UnityScale::exp<double>(kaux.bGetValue(), 48, 0.1);

			sprintf(buf, "-%.2f dBs", qDBs);
			kaux.bSetText(buf);
			sdft.setQ(qDBs);
			if (mappedFrequencies.size())
				sdft.mapSystemHz(mappedFrequencies, mappedFrequencies.size(), audioData[0].sampleRate);
			return true;
		}
		else if (ctrl == &kaux2)
		{
			sprintf(buf, "%.5f rads", kaux2.bGetValue() * 2 * M_PI);
			kaux2.bSetText(buf);
			sdft.setVectorDist(kaux2.bGetValue() * 2 * M_PI);
			if (mappedFrequencies.size())
				sdft.mapSystemHz(mappedFrequencies, mappedFrequencies.size(), audioData[0].sampleRate);
			return true;
		}
		else if (ctrl == &panel)
		{
			// panel switched, resize our display
			rearrange(getWidth(), getHeight(), false);

		}
		else if (ctrl == &klowDbs || ctrl == &khighDbs)
		{
			auto low = klowDbs.bGetValue();
			auto high = khighDbs.bGetValue();
			low = cpl::Math::UnityScale::linear<float>(low, kMinDbs, kMaxDbs);
			high = cpl::Math::UnityScale::linear<float>(high, kMinDbs, kMaxDbs);
			setDBs(low, high, false);
			auto const & actualDbs = getDBs();
			sprintf_s(buf, "%.2f dBs", actualDbs.low);
			klowDbs.bSetText(buf);
			sprintf_s(buf, "%.2f dBs", actualDbs.high);
			khighDbs.bSetText(buf);
			// need to repaint low dbs, in case it got changed:

			klowDbs.bSetInternal(cpl::Math::UnityScale::Inv::linear(actualDbs.low, kMinDbs, kMaxDbs));
			klowDbs.bRedraw();
			return true;
		}
		else if (ctrl == &kdecay)
		{
			// make a similar range to an exponential function here, which also maps 0
			auto x = ctrl->bGetValue();
			flt.setDecayAsFraction(x, 0.01f);
			sprintf_s(buf, "%.2f dBs/10ms", 20 * log10(x));
			ctrl->bSetText(buf);
			return true;
		}
		else if (ctrl == &kwindow)
		{

			setWindowSize(cpl::Math::round<int>(kwindow.bGetValue() * audioData[0].maxSize()));

			sprintf(buf, "%d smps", getWindowSize());
			ctrl->bSetText(buf);


			return true;
		}
		return false;
	}

	void CSpectrum::resized()
	{
		rearrange(getWidth(), getHeight());
	}
	void CSpectrum::repaintMainContent()
	{

		renderWindow.repaint();
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
		dbs.low = low;
		dbs.high = high;
		if (updateControls)
		{
			klowDbs.bSetValue(cpl::Math::UnityScale::Inv::linear<float>(low, kMinDbs, kMaxDbs));
			khighDbs.bSetValue(cpl::Math::UnityScale::Inv::linear<float>(high, kMinDbs, kMaxDbs));
		}
	}
	void CSpectrum::setDBs(CSpectrum::DBRange & newDbs, bool updateControls)
	{
		setDBs(newDbs.low, newDbs.high);
	}
	CSpectrum::DBRange CSpectrum::getDBs()
	{
		return dbs;
	}


};
