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
 
	file:COscilloscope.h
		
		Interface for the oscilloscope view.
 
*************************************************************************************/

#ifndef COSCILLOSCOPE_H
#define COSCILLOSCOPE_H

#include <cpl/Common.h>
#include <cpl/CLowPass.h>
#include <cpl/CAudioBuffer.h>
#include <cpl/GraphicComponents.h>
#include <vector>

namespace Signalizer
{
	using namespace juce;
	class Engine;

	class COscilloscope : public cpl::SubView, public cpl::CCtrlListener
	{
	public:
		COscilloscope(AudioBuffer & data);

		virtual ~COscilloscope() noexcept;
		void resize(int width, int height);
		void paint(Graphics & g) override;
		void repaintMainContent() override;
		bool serialize(juce::MemoryBlock & data) override;
		bool restore(const void * data, std::size_t size) override;
		bool setFullScreenMode(bool toggle) override;
	private:
		
		char buf[300];

		unsigned long long processorSpeed; // clocks / sec
		
		Image waveForm;
		Graphics * waveFormGraphics;
		float lpDelta;
		// ???
		::cpl::CLowPass filter;
		int oldPhaseDir;
		int windowSize;
		unsigned sizeFactor, offset;
		juce::Rectangle<int> displaySize;
		PathStrokeType pst;
		float cutoffFrequency;
		AudioBuffer & audioData;
		std::vector<unsigned int> graphColors;
		bool valueChanged(cpl::CBaseControl *) override;
		void resized() override;
		void updateVariables();
		double getGain();
		int getWindowSize();
		void paintChannelBitmap(Graphics & g, int channel, uint32_t color = 1);
		void paintMergedBitmap(Graphics & g);
		void paintBitmap(Graphics & g);
		
		void paintGraph(Graphics & g);

		cpl::CValueKnobEx kchannels, ksync, kpeakMode;
		cpl::CKnobEx kwindow, kcutoff, kgraphIntensity, kgain;
		
		std::vector<cpl::CColorKnob *> kcolors;

	};





};

#endif