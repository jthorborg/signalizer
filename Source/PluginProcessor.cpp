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
#include <cpl/Mathext.h>
#include <cpl/infrastructure/Values/Values.h>

namespace Signalizer
{

	std::unique_ptr<ParameterSet> CreateParameters(int offset, bool shouldCreateShortNames, ParameterSet::AutomatedProcessor & processor)
	{


		auto content = std::make_unique<Content>();
		auto ret = std::make_unique<ParameterSet>("Vectorscope", "VS.", processor, offset);

		ret->registerSingleParameter(content->autoGain.generateUpdateRegistrator());
		ret->registerSingleParameter(content->operationalMode.generateUpdateRegistrator());
		ret->registerSingleParameter(content->envelopeWindow.generateUpdateRegistrator());
		ret->registerSingleParameter(content->stereoWindow.generateUpdateRegistrator());
		ret->registerSingleParameter(content->inputGain.generateUpdateRegistrator());
		ret->registerSingleParameter(content->windowSize.generateUpdateRegistrator());
		ret->registerSingleParameter(content->waveZRotation.generateUpdateRegistrator());
		ret->registerSingleParameter(content->antialias.generateUpdateRegistrator());
		ret->registerSingleParameter(content->fadeOlderPoints.generateUpdateRegistrator());
		ret->registerSingleParameter(content->interconnectSamples.generateUpdateRegistrator());
		ret->registerSingleParameter(content->diagnostics.generateUpdateRegistrator());

		ret->registerParameterBundle(&content->drawingColour, "Draw");
		ret->registerParameterBundle(&content->graphColour, "Graph");
		ret->registerParameterBundle(&content->backgroundColour, "BckG");
		ret->registerParameterBundle(&content->skeletonColour, "Sklt");
		ret->registerParameterBundle(&content->primitiveSize, "PSize");
		ret->registerParameterBundle(&content->transform, "SPC");
		ret->setUserData(content.release());

		ret->seal();

		return std::move(ret);
	}

	//==============================================================================
	AudioProcessor::AudioProcessor()
		: editor(nullptr)
		, hasDefaultPresetBeenLoaded(false)
		, stream(16, true)
		, nChannels(2)
		, serializedData("HostState")
	{
		serializedData.getArchiver().setMasterVersion(cpl::programInfo.version);

		parameterSet = CreateParameters(0, false, *this);
	}

	void AudioProcessor::automatedTransmitChangeMessage(int parameter, ParameterSet::FrameworkType value)
	{
		sendParamChangeMessageToListeners(parameter, value);
	}

	void AudioProcessor::automatedBeginChangeGesture(int parameter)
	{
		beginParameterChangeGesture(parameter);
	}

	void AudioProcessor::automatedEndChangeGesture(int parameter)
	{
		endParameterChangeGesture(parameter);
	}

	AudioProcessor::~AudioProcessor() noexcept
	{
	}

	void AudioProcessor::onServerDestruction(cpl::DestructionNotifier * v)
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
	void AudioProcessor::prepareToPlay(double sampleRate, int samplesPerBlock)
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

	void AudioProcessor::releaseResources()
	{
		// When playback stops, you can use this as an opportunity to free up any
		// spare memory, etc.
	}

	void AudioProcessor::processBlock(AudioSampleBuffer& buffer, MidiBuffer& midiMessages)
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
			buffer.clear(i, 0, buffer.getNumSamples());
		}
	}

	//==============================================================================
	bool AudioProcessor::hasEditor() const
	{
		return true; // (change this to false if you choose to not supply an editor)
	}


	AudioProcessorEditor* AudioProcessor::createEditor()
	{
		//cpl::Misc::MsgBox("Attach debugger");

		editor = new Signalizer::MainEditor(this, parameterSet.get());
		editor->addEventListener(this);
		if (!serializedData.isEmpty())
			editor->deserializeObject(serializedData.getBuilder(), serializedData.getBuilder().getMasterVersion());

		return editor;
	}

	//==============================================================================
	void AudioProcessor::getStateInformation(MemoryBlock& destData)
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

	void AudioProcessor::setStateInformation(const void* data, int sizeInBytes)
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
	const String AudioProcessor::getName() const
	{
		return cpl::programInfo.name;
	}

	int AudioProcessor::getNumParameters()
	{
		return parameterSet->size();
	}

	float AudioProcessor::getParameter(int index)
	{
		if (auto param = parameterSet->findParameter(index))
		{
			return param->getValueNormalized<float>();
		}

		return 0;
	}

	void AudioProcessor::setParameter(int index, float newValue)
	{
		if (auto param = parameterSet->findParameter(index))
		{
			return param->updateFromHostNormalized(newValue);
		}
	}

	const String AudioProcessor::getParameterName(int index)
	{
		if (auto param = parameterSet->findParameter(index))
		{
			return param->getExportedName();
		}
		
		return "<no such parameter>";
	}

	const String AudioProcessor::getParameterText(int index)
	{
		if (auto param = parameterSet->findParameter(index))
		{
			return param->getDisplayText();
		}

		return "<no such parameter>";
	}

	const String AudioProcessor::getInputChannelName(int channelIndex) const
	{
		return String(channelIndex + 1);
	}

	const String AudioProcessor::getOutputChannelName(int channelIndex) const
	{
		return String(channelIndex + 1);
	}

	bool AudioProcessor::isInputChannelStereoPair(int index) const
	{
		return true;
	}

	bool AudioProcessor::isOutputChannelStereoPair(int index) const
	{
		return true;
	}

	bool AudioProcessor::acceptsMidi() const
	{
#if JucePlugin_WantsMidiInput
		return true;
#else
		return false;
#endif
	}

	bool AudioProcessor::producesMidi() const
	{
#if JucePlugin_ProducesMidiOutput
		return true;
#else
		return false;
#endif
	}

	bool AudioProcessor::silenceInProducesSilenceOut() const
	{
		return true;
	}

	double AudioProcessor::getTailLengthSeconds() const
	{
		return 0.0;
	}

	int AudioProcessor::getNumPrograms()
	{
		return 1;
	}

	int AudioProcessor::getCurrentProgram()
	{
		return 0;
	}

	void AudioProcessor::setCurrentProgram(int index)
	{
	}

	const String AudioProcessor::getProgramName(int index)
	{
		return String::empty;
	}

	void AudioProcessor::changeProgramName(int index, const String& newName)
	{
	}

};


//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new Signalizer::AudioProcessor();
}
