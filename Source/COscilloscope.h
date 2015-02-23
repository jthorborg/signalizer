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