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
 
	file:PluginProcessor.cpp

		Implementation of pluginprocessor.cpp
 
*************************************************************************************/

#include "PluginProcessor.h"
#include "MainEditor.h"
#include <cpl/CPresetManager.h>
#include <cpl/Protected.h>

//==============================================================================
SignalizerAudioProcessor::SignalizerAudioProcessor()
: 
	editor(nullptr), 
	hasDefaultPresetBeenLoaded(false),
	stream(16, true),
	nChannels(2),
	serializedData("HostState")
{
	serializedData.getArchiver().setMasterVersion(cpl::programInfo.version);
}

SignalizerAudioProcessor::~SignalizerAudioProcessor() noexcept
{
}

void SignalizerAudioProcessor::onServerDestruction(cpl::DestructionNotifier * v)
{
	if (v == editor)
	{
		serializedData.clear();
		serializedData.getArchiver().setMasterVersion(cpl::programInfo.version);
		editor->serializeObject(serializedData.getArchiver(), serializedData.getArchiver().getMasterVersion());

		editor = nullptr;
	}
}

//==============================================================================
void SignalizerAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
	Signalizer::AudioStream::AudioStreamInfo info = stream.getInfo();

	info.anticipatedChannels = nChannels;
	info.anticipatedSize = samplesPerBlock;
	info.callAsyncListeners = true;
	info.callRTListeners = true;
	info.isFrozen = false;
	info.isSuspended = false;
	info.sampleRate = sampleRate;
	info.storeAudioHistory = true;

	stream.initializeInfo(info);
}

void SignalizerAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

void SignalizerAudioProcessor::processBlock (AudioSampleBuffer& buffer, MidiBuffer& midiMessages)
{

	if (nChannels != buffer.getNumChannels())
	{
		// woah, what?
		CPL_BREAKIFDEBUGGED();
	}

	// stream will take it from here.
	stream.processIncomingRTAudio(buffer.getArrayOfWritePointers(), buffer.getNumChannels(), buffer.getNumSamples());

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
	//cpl::Misc::MsgBox("Attach debugger");

	editor = new Signalizer::MainEditor(this);
	editor->addEventListener(this);
	if (!serializedData.isEmpty())
		editor->deserializeObject(serializedData.getBuilder(), serializedData.getBuilder().getMasterVersion());

	return editor;
}

//==============================================================================
void SignalizerAudioProcessor::getStateInformation (MemoryBlock& destData)
{
	if (editor)
	{
		serializedData.clear();
		serializedData.getArchiver().setMasterVersion(cpl::programInfo.version);
		editor->serializeObject(serializedData.getArchiver(), cpl::programInfo.version);
	}
	if (!serializedData.isEmpty())
	{
		auto compiledData = serializedData.compile(true);
		destData.append(compiledData.getBlock(), compiledData.getSize());
	}
}

void SignalizerAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
	try
	{
		serializedData.clear();
		serializedData.build(cpl::WeakContentWrapper(data, sizeInBytes));
		// TODO: can make version check here.
		if (editor)
		{
			editor->deserializeObject(serializedData.getBuilder(), serializedData.getBuilder().getMasterVersion());
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
    return cpl::programInfo.name;
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
    return true;
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
