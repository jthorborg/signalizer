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
 
	file:PluginProcessor.h

		Defines the processing interface of the plugin.
 
*************************************************************************************/
#ifndef PLUGINPROCESSOR_H_INCLUDED
#define PLUGINPROCESSOR_H_INCLUDED

#include <cpl/Common.h>
#include <cpl/CAudioStream.h>
#include <cpl/CSerializer.h>
#include <cpl/gui/CViews.h>
#include "CommonSignalizer.h"
#include <cpl/gui/widgets/CPresetWidget.h>


namespace Signalizer
{
	class MainEditor;
};

class SignalizerAudioProcessor  : public juce::AudioProcessor, cpl::CView::EventListener
{
	friend class SignalizerAudioProcessorEditor;
	friend class Signalizer::MainEditor;

public:
    //==============================================================================
    SignalizerAudioProcessor();
    ~SignalizerAudioProcessor() noexcept;

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock);
    void releaseResources();

    void processBlock (AudioSampleBuffer& buffer, MidiBuffer& midiMessages);

    //==============================================================================
    AudioProcessorEditor* createEditor();
    bool hasEditor() const;

    //==============================================================================
    const String getName() const;

    int getNumParameters();

    float getParameter (int index);
    void setParameter (int index, float newValue);

    const String getParameterName (int index);
    const String getParameterText (int index);

    const String getInputChannelName (int channelIndex) const;
    const String getOutputChannelName (int channelIndex) const;
    bool isInputChannelStereoPair (int index) const;
    bool isOutputChannelStereoPair (int index) const;

    bool acceptsMidi() const;
    bool producesMidi() const;
    bool silenceInProducesSilenceOut() const;
    double getTailLengthSeconds() const;

    //==============================================================================
    int getNumPrograms();
    int getCurrentProgram();
    void setCurrentProgram (int index);
    const String getProgramName (int index);
    void changeProgramName (int index, const String& newName);

    //==============================================================================
    void getStateInformation (MemoryBlock& destData);
    void setStateInformation (const void* data, int sizeInBytes);
	
	void onServerDestruction(cpl::DestructionNotifier * v) override;
	
private:
    //==============================================================================
	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SignalizerAudioProcessor)
	cpl::CPresetWidget::SerializerType serializedData;
	Signalizer::MainEditor * editor;
	Signalizer::AudioStream stream;
	bool hasDefaultPresetBeenLoaded;
	int nChannels;
};

#endif 
