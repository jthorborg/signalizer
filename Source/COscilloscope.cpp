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
 
	file:COscilloscope.cpp

		Implementation of COscilloscope.h
 
*************************************************************************************/

#include "COscilloscope.h"
#include <cstdint>

namespace Signalizer
{
	static NativeImageType sit;

	struct ScanResult
	{
		size_t offset, direction;

	};

	template<class filter>
	ScanResult scanNextZeroCrossing(filter f, int maxSamples, int offset = 0)
	{
		ScanResult ret = { 0, 0 };
		float fnew = 0.0f, fold = 0.0f;

		fnew = f(offset);

		if (fnew > 0)
		{
			do
			{
				fold = fnew;

				fnew = f(++offset);


			} while (fnew > 0.f && offset < maxSamples);

			ret.direction = 1;
		}
		else if (fnew < 0)
		{
			do
			{
				fold = fnew;

				fnew = f(++offset);


			} while (fnew < 0.f && offset < maxSamples);

			ret.direction = -1;
		}
		else
		{
			ret.direction = 0;

		}


		ret.offset = offset;
		return ret;

	}



	COscilloscope::COscilloscope(AudioBuffer & data)
		: audioData(data), pst(1), waveForm(Image::PixelFormat::RGB, 800, 600, false, sit),
		oldPhaseDir(0), windowSize(1000),
		kchannels("Channel Drawing", "mono (left)|merge|seperate|split", ""),
		ksync("Synchronization", "lowpass zero's|free running", ""),
		kpeakMode("Point Sampling", "highest|average", ""),
		kwindow("Window Size", cpl::CKnobEx::type::ms),
		kcutoff("LP Adjust"), lpDelta(1), waveFormGraphics(nullptr), kgraphIntensity("Graph Intensity", cpl::CKnobEx::type::pct),
		kgain("Gain"), processorSpeed(0)
	{
		kchannels.bSetListener(this);
		ksync.bSetListener(this);
		kwindow.bSetListener(this);
		kcutoff.bSetListener(this);
		kgain.bSetListener(this);
		kgain.bSetValue(0.5);
		kcutoff.bSetValue(0.5);
		kchannels.bSetValue(0.26);
		kgraphIntensity.bSetValue(0.3);
		ksync.bSetValue(1);
		kpeakMode.bSetValue(0);
		kwindow.bSetValue(1000);
		kcutoff.bForceEvent();
		kpeakMode.bForceEvent();
		kwindow.bForceEvent();
		kchannels.bForceEvent();
		ksync.bForceEvent();
		kwindow.bForceEvent();
		addAndMakeVisible(kgain);
		addAndMakeVisible(kcutoff);
		addAndMakeVisible(kchannels);
		addAndMakeVisible(ksync);
		addAndMakeVisible(kwindow);
		addAndMakeVisible(kgraphIntensity);
		addAndMakeVisible(kpeakMode);
		
		for(int channel = 0; channel < audioData.size(); channel++)
		{
			kcolors.push_back(
							  new cpl::CColorKnob("Graph " + std::to_string(channel + 1) + " color",
												  cpl::CColorKnob::type::SingleChannel,
												  2 - channel % 2)
							  );
			kcolors[channel]->setColor(0xFF);
			//addAndMakeVisible(kcolors[channel]);
			
		}
		
		
		resized();
		filter.setResonance(0.1);
		
		if(e)
			e->onViewConstruction(this);
		
		processorSpeed = juce::SystemStats::getCpuSpeedInMegaherz();
	}
	void COscilloscope::repaintMainContent()
	{

		if (isOpenGL())
			;
		else
			repaint(0, 0, getWidth(), displaySize.getHeight());

	}
	struct CSerializedOscData
	{
		float kWindow, kCutoff, kChannels, kSync, kPeakMode, kGraphIntensity, kGain;
		
	};
	
	
	bool COscilloscope::restore(const void * data, std::size_t size)
	{
		const CSerializedOscData * se = reinterpret_cast<const CSerializedOscData * >(data);
		if(se && size == sizeof(CSerializedOscData))
		{
			kwindow.bSetValue(se->kWindow);
			kcutoff.bSetValue(se->kCutoff);
			kchannels.bSetValue(se->kChannels);
			ksync.bSetValue(se->kSync);
			kpeakMode.bSetValue(se->kPeakMode);
			kgraphIntensity.bSetValue(se->kGraphIntensity);
			kgain.bSetValue(se->kGain);
			updateVariables();
			return true;
			
		}

		return false;
		
	}
	
	bool COscilloscope::serialize(juce::MemoryBlock & data)
	{
		CSerializedOscData se;
		se.kWindow = kwindow.bGetValue();
		se.kCutoff = kcutoff.bGetValue();
		se.kChannels = kchannels.bGetValue();
		se.kGain = kgain.bGetValue();
		se.kGraphIntensity = kgraphIntensity.bGetValue();
		se.kPeakMode = kpeakMode.bGetValue();
		se.kSync = ksync.bGetValue();
		se.kWindow = kwindow.bGetValue();
		
		data.append(&se, sizeof(se));
		
		return true;
	}
	
	bool COscilloscope::setFullScreenMode(bool toggle)
	{
		kwindow.bSetVisible(!toggle);
		kcutoff.bSetVisible(!toggle);
		kchannels.bSetVisible(!toggle);
		ksync.bSetVisible(!toggle);
		kpeakMode.bSetVisible(!toggle);
		kgraphIntensity.bSetVisible(!toggle);
		kgain.bSetVisible(!toggle);
		isFullScreen = toggle;
		resized();
		return true;
	}
	COscilloscope::~COscilloscope()
	{
		
		if(e)
			e->onViewDestruction(this);
		
		for(auto & colorKnob : kcolors)
			delete colorKnob;
		if(waveFormGraphics)
			delete waveFormGraphics;
	}
	int COscilloscope::getWindowSize()
	{

		int dispWidth = displaySize.getWidth() - displaySize.getX();
		auto minMs = (dispWidth / audioData[0].sampleRate) * 1000;

		return std::ceil(minMs + kwindow.bGetValue() * (1000 - minMs));

	}
	
	double COscilloscope::getGain()
	{
		float val = kgain.bGetValue();
		if(val > 0.5)
		{
			return 1 + (val * 2 - 1) * 9;
			
		}
		else if(val < 0.5)
		{
			return val * 2;
		}
		
		return 1;
		
	}
	
	inline float fractionToDB(float val)
	{
		return (20.f * log10f (val));
	}

	bool COscilloscope::valueChanged(cpl::CBaseControl * ctrl)
	{
		char buf[200];
		if(ctrl == &kgain)
		{
			sprintf(buf, "%.2f dB (%d%%)", fractionToDB(getGain()), (int)cpl::Misc::Round(getGain() * 100));
			kgain.bSetText(buf);
			return true;
		}
		else if (ctrl == &kchannels)
		{

		}
		else if (ctrl == &ksync)
		{

		}
		else if (ctrl == &kwindow)
		{
			int bufLength = getWindowSize();
			for (auto & buf : audioData)
			{
				buf.setLength(bufLength);
			}

			sprintf(buf, "%d ms", bufLength);
			kwindow.bSetText(buf);

			updateVariables();
			kcutoff.bForceEvent();
			return true;
		}
		else if (ctrl == &kcutoff)
		{
			auto val = ctrl->bGetValue() * 2 - 1;


			if (val > 0)
			{
				lpDelta = 1 + val * 4;
				updateVariables();
				sprintf(buf, "(%.1f) x %.1f", lpDelta, cutoffFrequency);
			}
			else if (val < 0)
			{
				lpDelta = (val + 1);
				updateVariables();
				sprintf(buf, "(%.1f) x %.1f", lpDelta, cutoffFrequency);
			}
			else
			{
				lpDelta = 1;
				updateVariables();
				sprintf(buf, "(1) x %.1f", cutoffFrequency);
			}

			ctrl->bSetText(buf);
			return true;

		}
		return false;
	}

	void COscilloscope::resized()
	{
		resize(getWidth(), getHeight());
	}

	void COscilloscope::updateVariables()
	{
		if (getWidth() > 0)
		{
			auto width = displaySize.getWidth() - displaySize.getX();
			windowSize = getWindowSize();
			sizeFactor = audioData[0].size / width;
			cutoffFrequency = (1000.0 / windowSize) * 40 * width / 800 * lpDelta;
			filter.setCutoff(cutoffFrequency);
		}

	}
	void COscilloscope::resize(int width, int height)
	{
		if (!width || !height)
			return;


		// place knobs

		int knobTop = getHeight() - kchannels.getHeight();
		int knobWidth = kchannels.getWidth();
		int leftMargin = 20;
		int counter = 0;
		
		int maxDisplaySize = std::min((size_t) getWidth(), audioData[0].size);
		if(!isFullScreen)
		{
			kchannels.bSetPos(leftMargin + counter++ * knobWidth, knobTop);
			ksync.bSetPos(leftMargin + counter++ * knobWidth, knobTop);
			kwindow.bSetPos(leftMargin + counter++ * knobWidth, knobTop);
			kpeakMode.bSetPos(leftMargin + counter++ * knobWidth, knobTop);
			kcutoff.bSetPos(leftMargin + counter++ * knobWidth, knobTop);
			kgraphIntensity.bSetPos(leftMargin + counter++ * knobWidth, knobTop);
			kgain.bSetPos(leftMargin + counter++ * knobWidth, knobTop);
			for(auto & colorKnob : kcolors)
				colorKnob->bSetPos(leftMargin + counter++ * knobWidth, knobTop);
			// resize image
			displaySize.setBounds(40, 0, maxDisplaySize, knobTop);
			waveForm = Image(Image::PixelFormat::RGB, maxDisplaySize + leftMargin, getHeight(), false, sit);
			if (waveFormGraphics)
			{
				delete waveFormGraphics;
				waveFormGraphics = nullptr;
			}
			waveFormGraphics = new Graphics(waveForm);
		}
		else
		{
			displaySize.setBounds(40, 0, maxDisplaySize, getHeight());
			waveForm = Image(Image::PixelFormat::RGB, maxDisplaySize + leftMargin, getHeight(), false, sit);
			if (waveFormGraphics)
			{
				delete waveFormGraphics;
				waveFormGraphics = nullptr;
			}
			waveFormGraphics = new Graphics(waveForm);
		}
		updateVariables();
	}
	void COscilloscope::paint(Graphics & g)
	{

		auto clockStart = cpl::Misc::ClockCounter();

		g.fillAll(Colours::black);
		g.setColour(Colours::green);
		waveFormGraphics->fillAll(Colours::black);
		auto val = cpl::Misc::Round(kchannels.bGetValue() * 3);
		switch(val)
		{
			case 0: // mono (left) mode
			{
				paintGraph(*waveFormGraphics);
				paintChannelBitmap(g, 0);
				break;
			}
			case 1: // merge mode
			{
				paintGraph(*waveFormGraphics);
				paintMergedBitmap(g);
				break;
			}
			case 2: // seperate mode
			{
				paintGraph(*waveFormGraphics);
				auto numChannels = audioData.size();

				for(int channel = 0; channel < numChannels; channel ++)
				{
					paintChannelBitmap(g, channel, 2 - channel);
				}
				break;
			}
			case 3: // split
			{
				auto oldSize = displaySize;
				
				auto numChannels = audioData.size();
				
				auto newHeight = displaySize.getHeight() / numChannels;
				displaySize.setHeight(newHeight);
				auto oldTop = displaySize.getY();
				for(int channel = 0; channel < numChannels; channel ++)
				{

					displaySize.setY(oldTop + channel * newHeight);
					
					
					paintGraph(*waveFormGraphics);
					paintChannelBitmap(g, channel, 1);
					
				}

				
				displaySize = oldSize;
			}
		}
		auto renderStop = cpl::Misc::ClockCounter();
		g.drawImageAt(waveForm, 0, 0);
		
		
		if(!isFullScreen)
		{
			g.setFont(18);
			g.setColour(Colours::blue);
			auto cyclesNow = cpl::Misc::ClockCounter();
			auto renderCycles = renderStop - clockStart;
			auto transferCycles = cyclesNow - renderStop;
			auto totalCycles = cyclesNow - clockStart;
			auto processorTime = (totalCycles * (1000.0 / refreshRate)) / (processorSpeed * 1000 * 1000);
			auto renderTime = (double)renderCycles/totalCycles;
			auto transferTime = (double)transferCycles / totalCycles;
			sprintf_s(buf,
					  "Diag: %dX%d;  Cpu time: %.1f%% (render: %.1f%%, transfer: %.1f%%, sample:pixel factor %u:1",
					 getWidth(), getHeight(), processorTime * 100, renderTime * 100, transferTime * 100, sizeFactor);
			g.drawSingleLineText(buf, 20, 40);
		}
	}

	void COscilloscope::paintGraph(Graphics & g)
	{
		unsigned char colorIntensity = (unsigned char)kgraphIntensity.bGetValue() * 0xFF;
		if(!colorIntensity)
			return;
		auto top = displaySize.getY();	
		auto height = displaySize.getHeight();
		auto middle = top + height / 2;

		auto left = displaySize.getX();
		auto width = displaySize.getWidth();
		auto displayWidth = width - left;
		auto bottom = top + height;
		
		char stringbuf[200];
		g.setColour(Colour(colorIntensity , colorIntensity, colorIntensity));
		g.drawHorizontalLine(middle, 0, getWidth());
		g.drawLine(left, top, left, bottom, 2);
		
		// draw offset imposed by zero crossing syncs

		g.setFont(cpl::TextSize::smallerText);
		

		
		int numLevelMarkers = height / 2 / 25;

		

		
		float diff = (float)height / 2 / numLevelMarkers;
		auto gain = getGain();
		for(int i = 0; i < numLevelMarkers; i++)
		{
			auto yCoord = top + i * diff;
			g.drawLine(left, yCoord, width, yCoord);
			sprintf(stringbuf, "% .1f dB", fractionToDB(((float)(numLevelMarkers - i) / numLevelMarkers) * 1/gain) );
			g.drawSingleLineText(stringbuf, 0, yCoord );
			
			g.drawLine(left, top + height - i * diff, width, top + height - i * diff);
		}
		
		float pxlPerMs = (float)displayWidth / windowSize;
		float timeDivisor = 0.01;

		int maxLines = displayWidth / 38;
		int actualNumLines;
		while ((actualNumLines = (displayWidth / (timeDivisor * pxlPerMs))) > maxLines)
			timeDivisor *= 10;

		/*

		float divisor = 1000;
		int maxLines = (float) width / 25;
		while (msPerPxl * divisor > maxLines)
			divisor /= 10;
		maxLines = divisor * msPerPxl;
		int numPerfectTimeDivisions = width / (msPerPxl * divisor);//windowsSize;
		*/
		for(int i = 0; i < actualNumLines; ++i)
		{
			auto xCoord = left + i * timeDivisor * pxlPerMs;
			g.drawLine(xCoord, top, xCoord, bottom);
			sprintf(stringbuf, "%.1f ms", timeDivisor * i);
			g.drawSingleLineText(stringbuf, xCoord + 5, middle + 5 + cpl::TextSize::smallerText);
		}

		sprintf(stringbuf, "0 = +%u // %f", offset, pxlPerMs);

		g.drawSingleLineText(stringbuf, left + 10, top + 10);
	}


	
	void COscilloscope::paintMergedBitmap(Graphics & g)
	{
		auto left = displaySize.getX();
		auto height = displaySize.getHeight();
		auto top = displaySize.getY();
		auto middle = height / 2;
		auto gain = getGain();
		auto bottom = top + height;
		
		auto width = displaySize.getWidth();
		uint8_t colorIntensity = (uint8_t)this->kcolors[0]->bGetValue() * 0xFF;
		if (bottom > top)
		{
			//	current point	old point	temporary buffer
			float hPoint = 0, oldHPoint = 0, temp = 0;
			
			Image::BitmapData data(waveForm, Image::BitmapData::ReadWriteMode::writeOnly);
			
			
			const bool peakmode = kpeakMode.bGetValue() > 0.5;
			/*
			 calculate the next zero crossing in phase with previous crossings.
			 this is evaluated through a 4-pole lowpass filter, set to cutoff frequency
			 that maps to around 10 cycles per 100 pixels displaywidth
			 */
			offset = 1;
			//auto startpoint = filter.process(engine->audioBuffer[0][0]);
			
			int maxSamples = width * sizeFactor / 2;
			
			if (ksync.bGetValue() < 0.5)
			{
				ScanResult sr{ offset, 0 };
				do
				{
					sr = scanNextZeroCrossing
					(
					 [&](int q) { return filter.process(audioData[0][q]); },
					 maxSamples,
					 sr.offset
					 );
				} while (sr.offset < maxSamples && sr.direction &&  sr.direction != oldPhaseDir);
				
				oldPhaseDir = sr.direction;
				offset = sr.offset;
			}
			
			for (auto x = 0; x < width - 1 - left; ++x)
			{
				// store / reset points
				temp = hPoint = 0;
				for(unsigned channel = 0; channel < audioData.size(); channel++)
				{
					if (peakmode)
					{
						for (unsigned z = 0; z < sizeFactor; z++)
						{
							temp += audioData[channel][x * sizeFactor + z + offset];
						}
						hPoint += temp / sizeFactor;
					}
					else
					{
						// interpolate next point (preserve peaks)
						for (unsigned z = 0; z < sizeFactor; z++)
						{
							temp = audioData[channel][x * sizeFactor + z + offset];
							// kinda dogdy, but this is exactly the same as a std::fabs() function
							// we reinterpret a float as r-value integer, where we set the sign bit to zero, and
							// attach a r-value reference (reinterpretation again) float to it.
							// this is around a 12x optimization compared to std::fabs, even in release.
							if ((const float &&)(*((int *)&temp) & ~0x80000000) >(const float &&)(*((int *)&hPoint) & ~0x80000000))
								hPoint = temp;
						}
					}
				}
				float fHighestPoint = std::max(hPoint, oldHPoint) * gain;
				float fLowestPoint = std::min(hPoint, oldHPoint) * gain;
				
				int start = fLowestPoint * middle + middle;
				int stop = fHighestPoint * middle + middle;
				//int pos_bottom_left = std::min(hPointRight, oldHPointRight) * middle + middle;
				//int pos_top_left = std::max(hPointRight, oldHPointRight) * middle + middle;
				
				int direction = oldHPoint > hPoint;
				
				if (stop <= top)
					stop = top;
				if (stop >= bottom)
					stop = bottom - 1;
				if (start <= top)
					start = top;
				if (start >= bottom)
					start = bottom - 1;
				// can be zero -> division by zero later
				float difference = stop - start;
				if (difference == 0.0f)
					difference = 1.f;
				int index = 0;
				std::uint8_t * p1, *p2;
				
				int xCoord = x + left;
				do
				{
					// this might need some work
					p1 = data.getPixelPointer(xCoord + !direction, start + index) + 1;
					/*
					 here we add to the existing pixel a percentage times the factor of the remaining bits divided by two
					 this has the effect of never overflowing the byte, but still add more color.
					 
					 */
					*p1 = *p1 + ((254 - *p1) >> 1) *  (index + 1) / difference; // draw a straight line from bottom to top
					p1 = data.getPixelPointer(xCoord + direction, start + index) + 1;
					*p1 = *p1 + ((254 - *p1) >> 1)*  (difference - index) / difference; // draw a straight line from bottom to top
					index++;
				} while (stop > (start + index));
				oldHPoint = hPoint;
			}
			
			
		}// bottom > top
	}
	

	void COscilloscope::paintChannelBitmap(Graphics & g, int channel, uint32_t color)
	{
		auto left = displaySize.getX();
		auto height = displaySize.getHeight();
		auto middle = height / 2;
		auto top = displaySize.getY();
		auto bottom = top + height;
		uint8_t colorIntensity = (uint8_t)this->kcolors[channel]->getColor();
		auto width = displaySize.getWidth();
		auto gain = getGain();
		uint32_t shift = color % 2;
		if (bottom > top)
		{
			//	current point	old point	temporary buffer
			float hPoint = 0, oldHPoint = 0, temp = 0;

			Image::BitmapData data(waveForm, Image::BitmapData::ReadWriteMode::writeOnly);


			const bool peakmode = kpeakMode.bGetValue() > 0.5;
			/*
			calculate the next zero crossing in phase with previous crossings.
			this is evaluated through a 4-pole lowpass filter, set to cutoff frequency
			that maps to around 10 cycles per 100 pixels displaywidth
			*/
			offset = 1;
			//auto startpoint = filter.process(engine->audioBuffer[0][0]);

			int maxSamples = width * sizeFactor / 2;

			if (ksync.bGetValue() < 0.5)
			{
				ScanResult sr{ offset, 0 };
				do
				{
					sr = scanNextZeroCrossing
						(
						[&](int q) { return filter.process(audioData[0][q]); },
						maxSamples,
						sr.offset
						);
				} while (sr.offset < maxSamples && sr.direction &&  sr.direction != oldPhaseDir);

				oldPhaseDir = sr.direction;
				offset = sr.offset;
			}

			for (auto x = 0; x < width - 1 - left; ++x)
			{
				// store / reset points
				temp = hPoint = 0;

				if (peakmode)
				{
					for (int z = 0; z < sizeFactor; z++)
					{
						temp += audioData[channel][x * sizeFactor + z + offset];
					}
					hPoint = temp / sizeFactor;
				}
				else
				{
					// interpolate next point (preserve peaks)
					for (int z = 0; z < sizeFactor; z++)
					{
						temp = audioData[channel][x * sizeFactor + z + offset];
						// kinda dogdy, but this is exactly the same as a std::fabs() function
						// we reinterpret a float as r-value integer, where we set the sign bit to zero, and 
						// attach a r-value reference (reinterpretation again) float to it.
						// this is around a 12x optimization compared to std::fabs, even in release.
						if ((const float &&)(*((int *)&temp) & ~0x80000000) >(const float &&)(*((int *)&hPoint) & ~0x80000000))
							hPoint = temp;
					}
				}
				float fHighestPoint = std::max(hPoint, oldHPoint) * gain;
				float fLowestPoint = std::min(hPoint, oldHPoint) * gain;

				int start = top + fLowestPoint * middle + middle;
				int stop = top + fHighestPoint * middle + middle;
				//int pos_bottom_left = std::min(hPointRight, oldHPointRight) * middle + middle;
				//int pos_top_left = std::max(hPointRight, oldHPointRight) * middle + middle;

				int direction = oldHPoint > hPoint;

				if (stop <= top)
					stop = top;
				if (stop >= bottom)
					stop = bottom - 1;
				if (start <= top)
					start = top;
				if (start >= bottom)
					start = bottom - 1;
				// can be zero -> division by zero later
				float difference = stop - start;
				if (difference == 0.0f)
					difference = 1.f;
				int index = 0;
				std::uint8_t * p1;
				
				int xCoord = x + left;
				do
				{
					// this might need some work
					p1 = data.getPixelPointer(xCoord + !direction, start + index) + color;
					/*
						here we add to the existing pixel a percentage times the factor of the remaining bits divided by two
						this has the effect of never overflowing the byte, but still add more color.
						
					*/
					*p1 = *p1 + ((254 - *p1) >> 1) *  (index + 1) / difference; // draw a straight line from bottom to top
					p1 = data.getPixelPointer(xCoord + direction, start + index) + color;
					*p1 = *p1 + ((254 - *p1) >> 1)*  (difference - index) / difference; // draw a straight line from bottom to top
					index++;
				} while (stop > (start + index));
				oldHPoint = hPoint;
			}


		 }// bottom > top
	}
	

	void COscilloscope::paintBitmap(Graphics & g)
	{
		g.setColour(Colours::green);
		auto left = displaySize.getX();
		auto height = displaySize.getHeight();
		auto middle = height / 2;
		auto top = 0;
		auto bottom = top + height;
		
		auto width = displaySize.getWidth();
		
		if (bottom > top)
		{
			float hPointLeft = 0, oldHPointLeft = 0, tempLeft = 0;
			float hPointRight = 0, oldHPointRight = 0, tempRight = 0;
			
			
			
			Image::BitmapData data(waveForm, Image::BitmapData::ReadWriteMode::writeOnly);
			
			
			const bool peakmode = kpeakMode.bGetValue() > 0.5;
			/*
			 calculate the next zero crossing in phase with previous crossings.
			 this is evaluated through a 4-pole lowpass filter, set to cutoff frequency
			 that maps to around 10 cycles per 100 pixels displaywidth
			 */
			offset = 1;
			//auto startpoint = filter.process(engine->audioBuffer[0][0]);
			
			int maxSamples = width * sizeFactor / 2;
			
			if (ksync.bGetValue() < 0.5)
			{
				ScanResult sr{ offset, 0 };
				do
				{
					sr = scanNextZeroCrossing
					(
					 [&](int q) { return filter.process(audioData[0][q]); },
					 maxSamples,
					 sr.offset
					 );
				} while (sr.offset < maxSamples && sr.direction &&  sr.direction != oldPhaseDir);
				
				oldPhaseDir = sr.direction;
				offset = sr.offset;
			}
			
			for (auto x = left; x < width - 1; ++x)
			{
				// store / reset points
				tempRight = tempLeft = hPointLeft = hPointRight = 0;
				
				if (peakmode)
				{
					for (int z = 0; z < sizeFactor; z++)
					{
						tempLeft += audioData[0][x * sizeFactor + z + offset];
					}
					hPointLeft = tempLeft / sizeFactor;
				}
				else
				{
					// interpolate next point (preserve peaks)
					for (int z = 0; z < sizeFactor; z++)
					{
						tempLeft = audioData[0][x * sizeFactor + z + offset];
						tempRight = audioData[1][x * sizeFactor + z + offset];
						// kinda dogdy, but this is exactly the same as a std::fabs() function
						// we reinterpret a float as r-value integer, where we set the sign bit to zero, and
						// attach a r-value reference (reinterpretation again) float to it.
						// this is around a 12x optimization compared to std::fabs, even in release.
						if ((const float &&)(*((int *)&tempLeft) & ~0x80000000) >(const float &&)(*((int *)&hPointLeft) & ~0x80000000))
							hPointLeft = tempLeft;
						if ((const float &&)(*((int *)&tempRight) & ~0x80000000) > (const float &&)(*((int *)&hPointRight) & ~0x80000000))
							hPointRight = tempRight;
					}
				}
				float fHighestPoint = std::max(hPointLeft, oldHPointLeft);
				float fLowestPoint = std::min(hPointLeft, oldHPointLeft);
				
				int start = fLowestPoint * middle + middle;
				int stop = fHighestPoint * middle + middle;
				//int pos_bottom_left = std::min(hPointRight, oldHPointRight) * middle + middle;
				//int pos_top_left = std::max(hPointRight, oldHPointRight) * middle + middle;
				
				int direction = oldHPointLeft > hPointLeft;
				
				if (stop <= top)
					stop = top;
				if (stop >= bottom)
					stop = bottom - 1;
				if (start <= top)
					start = top;
				if (start >= bottom)
					start = bottom - 1;
				// can be zero -> division by zero later
				float difference = stop - start;
				if (difference == 0.0f)
					difference = 1.f;
				int index = 0;
				std::uint8_t * p1, *p2;
				do
				{
					// this might need some work
					p1 = data.getPixelPointer(x + !direction, start + index) + 1;
					/*
					 here we add to the existing pixel a percentage times the factor of the remaining bits divided by two
					 this has the effect of never overflowing the byte, but still add more color.
					 
					 */
					*p1 = *p1 + ((254 - *p1) >> 1) *  (index + 1) / difference; // draw a straight line from bottom to top
					p1 = data.getPixelPointer(x + direction, start + index) + 1;
					*p1 = *p1 + ((254 - *p1) >> 1)*  (difference - index) / difference; // draw a straight line from bottom to top
					index++;
				} while (stop > (start + index));
				oldHPointLeft = hPointLeft;
			}
			
			
			
			g.drawImageAt(waveForm, 0, 0);
			g.drawSingleLineText(std::to_string(cutoffFrequency), 100, 200);
		}// bottom > top
	}
};