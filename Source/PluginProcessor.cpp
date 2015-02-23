/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic startup code for a Juce application.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "MainEditor.h"
#include <cpl/CPresetManager.h>

//==============================================================================
SignalizerAudioProcessor::SignalizerAudioProcessor()
: audioBuffer(2), editor(nullptr)
{
	juce::File result;
	cpl::CPresetManager::instance().loadDefaultPreset(serializedData, result);

}

SignalizerAudioProcessor::~SignalizerAudioProcessor()
{
}

void SignalizerAudioProcessor::onViewConstruction(cpl::View * view)
{
	
	
}

void SignalizerAudioProcessor::onViewDestruction(cpl::View * view)
{
	serializedData.clear();
	view->save(serializedData, 1);

	editor = nullptr;
}

//==============================================================================
void SignalizerAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
	for (auto & buf : audioBuffer)
		buf.setSampleRate(getSampleRate());
    // Use this method as the place to do any pre-playback
    // initialisation that you need..
}

void SignalizerAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

void SignalizerAudioProcessor::processBlock (AudioSampleBuffer& buffer, MidiBuffer& midiMessages)
{
    // This is the place where you'd normally do the guts of your plugin's
    // audio processing...
	std::vector<cpl::CFastMutex> locks;
	auto buffers = buffer.getArrayOfWritePointers();

	std::size_t numSamples = buffer.getNumSamples();

	for (unsigned i = 0; i < audioBuffer.size(); ++i)
	{
		audioBuffer[i].raiseAudioEvent(buffers, audioBuffer.size(),  numSamples);
		locks.emplace_back(audioBuffer[i]);

	}


	for (unsigned i = 0; i < numSamples; ++i)
	{
		for (unsigned channel = 0; channel < audioBuffer.size(); ++channel)
		{
			audioBuffer[channel].setNextSample(buffers[channel][i]);

			// ..do something to the data...
		}
	}
    // In case we have more outputs than inputs, we'll clear any output
    // channels that didn't contain input data, (because these aren't
    // guaranteed to be empty - they may contain garbage).
    for (int i = getNumInputChannels(); i < getNumOutputChannels(); ++i)
    {
        buffer.clear (i, 0, buffer.getNumSamples());
    }
}

//==============================================================================
bool SignalizerAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}


AudioProcessorEditor* SignalizerAudioProcessor::createEditor()
{
	editor = new Signalizer::MainEditor(this);
	editor->addEventListener(this);
	if (!serializedData.isEmpty())
		editor->load(serializedData, 1);
	return editor;
}

//==============================================================================
void SignalizerAudioProcessor::getStateInformation (MemoryBlock& destData)
{
	if (editor)
	{
		serializedData.clear();
		editor->save(serializedData, 1);
	}
	if (!serializedData.isEmpty())
	{
		auto compiledData = serializedData.compile();
		destData.append(compiledData.getBlock(), compiledData.getSize());
	}
}

void SignalizerAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
	try
	{
		serializedData.clear();
		serializedData.build(cpl::CSerializer::WeakContentWrapper(data, sizeInBytes));
		if (editor)
		{
			editor->load(serializedData, 1);
		}
	}
	catch (const std::exception & e)
	{
		cpl::Misc::MsgBox(std::string("Error loading state information:\n") + e.what(), "Signalizer");
	}

}


//==============================================================================
const String SignalizerAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

int SignalizerAudioProcessor::getNumParameters()
{
    return 0;
}

float SignalizerAudioProcessor::getParameter (int index)
{
    return 0.0f;
}

void SignalizerAudioProcessor::setParameter (int index, float newValue)
{
}

const String SignalizerAudioProcessor::getParameterName (int index)
{
    return String::empty;
}

const String SignalizerAudioProcessor::getParameterText (int index)
{
    return String::empty;
}

const String SignalizerAudioProcessor::getInputChannelName (int channelIndex) const
{
    return String (channelIndex + 1);
}

const String SignalizerAudioProcessor::getOutputChannelName (int channelIndex) const
{
    return String (channelIndex + 1);
}

bool SignalizerAudioProcessor::isInputChannelStereoPair (int index) const
{
    return true;
}

bool SignalizerAudioProcessor::isOutputChannelStereoPair (int index) const
{
    return true;
}

bool SignalizerAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool SignalizerAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool SignalizerAudioProcessor::silenceInProducesSilenceOut() const
{
    return false;
}

double SignalizerAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int SignalizerAudioProcessor::getNumPrograms()
{
    return 1;
}

int SignalizerAudioProcessor::getCurrentProgram()
{
    return 0;
}

void SignalizerAudioProcessor::setCurrentProgram (int index)
{
}

const String SignalizerAudioProcessor::getProgramName (int index)
{
    return String::empty;
}

void SignalizerAudioProcessor::changeProgramName (int index, const String& newName)
{
}



//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new SignalizerAudioProcessor();
}
