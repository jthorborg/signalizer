/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic startup code for a Juce application.

  ==============================================================================
*/

#ifndef PLUGINPROCESSOR_H_INCLUDED
#define PLUGINPROCESSOR_H_INCLUDED
//#define JUCE_ENABLE_REPAINT_DEBUGGING 1
#include "../JuceLibraryCode/JuceHeader.h"

#include <cpl/CAudioBuffer.h>
#include <cpl/GraphicComponents.h>
#include <cpl/lib/AlignedAllocator.h>
#include <cpl/CSerializer.h>

//==============================================================================
/**
*/

class SignalizerAudioProcessorEditor;
namespace Signalizer
{
	class MainEditor;
}


class SignalizerAudioProcessor  : public AudioProcessor, cpl::View::EventListener
{
	friend class SignalizerAudioProcessorEditor;
	friend class Signalizer::MainEditor;
public:
    //==============================================================================
    SignalizerAudioProcessor();
    ~SignalizerAudioProcessor();

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
	
	void onViewConstruction(cpl::View * view) override;
	void onViewDestruction(cpl::View * view) override;
	
private:
    //==============================================================================
	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SignalizerAudioProcessor)
	cpl::CSerializer serializedData;
	Signalizer::MainEditor * editor;
	cpl::AudioBuffer audioBuffer;
};

#endif  // PLUGINPROCESSOR_H_INCLUDED
