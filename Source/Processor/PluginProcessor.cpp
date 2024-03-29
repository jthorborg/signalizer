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
#include "../Editor/MainEditor.h"
#include <cpl/CPresetManager.h>
#include <cpl/Protected.h>
#include <cpl/Mathext.h>
#include <cpl/infrastructure/values/Values.h>
#include <array>

namespace Signalizer
{
	extern std::vector<std::pair<std::string, ContentCreater>> ContentCreationList;
	extern std::string MainPresetName;
	extern std::string DefaultPresetName;

	constexpr int supportedChannels = 2;

	AudioProcessor::AudioProcessor()
		: AudioProcessor(AudioStream::create(true, 16))
	{

	}

	std::shared_ptr<const ConcurrentConfig> AudioProcessor::getConcurrentConfig()
	{
		return std::const_pointer_cast<const ConcurrentConfig>(config);
	}

	AudioProcessor::AudioProcessor(AudioStream::IO&& io)
		: config(std::make_shared<ConcurrentConfig>())
		, realtimeInput(std::move(std::get<0>(io)))
		, realtimeOutput(std::get<1>(io))
		, graph(std::make_shared<HostGraph>(std::get<1>(io)))
		, lastRecordedInputCount(1)
		, dsoEditor(
			[this] { return std::make_unique<MainEditor>(this, &this->parameterMap); },
			[](MainEditor & editor, cpl::CSerializer & sz, cpl::Version v) { editor.serializeObject(sz, v); },
			[](MainEditor & editor, cpl::CSerializer & sz, cpl::Version v) { editor.deserializeObject(sz, v); }
		)
	{

		SystemView view { getConcurrentConfig(), *this};

		for (std::size_t i = 0; i < ContentCreationList.size(); ++i)
		{
			parameterMap.insert({
				ContentCreationList[i].first,
				ContentCreationList[i].second(parameterMap.numParams(), view)
			});
		}

		juce::File location;

		// load the default preset
		try
		{
			SerializerType serializer(MainPresetName);

			cpl::CPresetManager::instance().loadPreset(
				cpl::CPresetManager::instance().getPresetDirectory() + "default." + MainPresetName + "." + cpl::programInfo.programAbbr,
				serializer,
				location
			);

			if (!serializer.isEmpty())
			{
				deserialize(serializer.getBuilder(), serializer.getBuilder().getLocalVersion());
			}
		}
		catch (std::exception & e)
		{
			cpl::Misc::MsgBox(std::string("Error reading state information from default preset:\n") + e.what(), cpl::programInfo.name);
		}

		// initialize audio stream with some default values, fixes a bug with the time knobs that rely on a valid sample rate being set.
		// TODO: convert them to be invariant.

		realtimeInput.initializeInfo(
			[&](AudioStream::ProducerInfo& info) 
			{
				info.channels = supportedChannels;
				info.anticipatedSize = 512;
				info.sampleRate = 48000;
			}
		);
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
		notifyDestruction();
	}

	//==============================================================================
	void AudioProcessor::prepareToPlay(double sampleRate, int samplesPerBlock)
	{
		lastRecordedInputCount = getNumInputChannels();
		lastRecordedBufferSize = samplesPerBlock;
		surrogateArray.resize(samplesPerBlock);

		realtimeInput.initializeInfo(
			[&](AudioStream::ProducerInfo& info)
			{
				info.channels = supportedChannels;
				info.anticipatedSize = samplesPerBlock;
				info.sampleRate = sampleRate;
			}
		);

		if (!hasAnyLayoutBeenApplied)
		{
			hasAnyLayoutBeenApplied = true;
			graph->applyDefaultLayoutFromRuntime();
		}
	}

	void AudioProcessor::releaseResources()
	{
	}

	void AudioProcessor::processBlock(juce::AudioSampleBuffer& buffer, juce::MidiBuffer& midiMessages)
	{
		int bufferSize = buffer.getNumSamples();

		if (!NONTERMINAL_ASSUMPTION(lastRecordedInputCount == getNumInputChannels()))
			return;

		if (!NONTERMINAL_ASSUMPTION(getNumInputChannels() <= buffer.getNumChannels()))
			return;

		if (!NONTERMINAL_ASSUMPTION(bufferSize <= lastRecordedBufferSize))
			return;

		if (realtimeInput.isAnyoneListening())
		{
			// TODO: Fix this when the mix graph listener supports dynamically changing channels
			std::array<const float*, supportedChannels> inputs;
			auto readPointers = buffer.getArrayOfReadPointers();

			const auto available = std::min(getNumInputChannels(), supportedChannels);

			int i = 0;
			for (; i < available; ++i)
			{
				inputs[i] = readPointers[i];
			}

			for (; i < supportedChannels; ++i)
			{
				inputs[i] = surrogateArray.data();
			}

			if (auto ph = getPlayHead())
				realtimeInput.processIncomingRTAudio(inputs.data(), supportedChannels, buffer.getNumSamples(), *ph);
			else
				realtimeInput.processIncomingRTAudio(inputs.data(), supportedChannels, buffer.getNumSamples(), AudioStream::Playhead::empty());
		}

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


	juce::AudioProcessorEditor* AudioProcessor::createEditor()
	{
		std::lock_guard<std::mutex> lock(editorCreationMutex);
		return dsoEditor.getUnique().acquire();
	}

	//==============================================================================
	void AudioProcessor::getStateInformation(juce::MemoryBlock& destData)
	{
		cpl::CPresetWidget::SerializerType serializer("HostState");
		auto& archive = serializer.getArchiver();
		archive.setMasterVersion(cpl::programInfo.version);

		archive << this;
		archive["host-graph"] << *graph;

		if (!serializer.isEmpty())
		{
			auto compiledData = serializer.compile(true);
			destData.append(compiledData.getBlock(), compiledData.getSize());
		}
	}

	void AudioProcessor::setStateInformation(const void* data, int sizeInBytes)
	{
		SerializerType serializer("HostState");

		try
		{
			serializer.build(cpl::WeakContentWrapper(data, sizeInBytes));
		}
		catch (const std::exception & e)
		{
			cpl::Misc::MsgBox(std::string("Error reading state information:\n") + e.what(), "Signalizer");
		}

		if (serializer.isEmpty())
			return;

		try
		{
			auto& builder = serializer.getBuilder();
			auto version = serializer.getBuilder().getLocalVersion();
			deserialize(builder, version);

			auto& serializedHostGraph = builder["host-graph"];

			if (!serializedHostGraph.isEmpty())
			{
				serializedHostGraph >> *graph;
				hasAnyLayoutBeenApplied |= graph->getDidEverDeserializeTopology();
			}

		}
		catch (const std::exception & e)
		{
			cpl::Misc::MsgBox(std::string("Error serializing state information:\n") + e.what(), "Signalizer");
		}
	}

	void AudioProcessor::deserialize(cpl::CSerializer & serializer, cpl::Version version)
	{
		auto & editorContent = version < cpl::Version(0, 2, 8) ? serializer : serializer.getContent("Editor");

		if (!editorContent.isEmpty())
		{
			std::lock_guard<std::mutex> lock(editorCreationMutex);

			if (dsoEditor.hasCached() && !juce::MessageManager::getInstance()->isThisTheMessageThread())
			{
				// writing to this value with acquire/release ensures we see any changes concurrently here
				std::atomic_int editorSerializationState {0};

				cpl::GUIUtils::MainEvent(
					*this,
					[&]
					{
						// interesting: host closed our window after calling setStateInformation but still while waiting for it to end..
						if (!dsoEditor.hasCached())
						{
							editorSerializationState.store(2, std::memory_order_release);
							return;
						}

						// if the object dies after this, it's a framework implementation issues (deleting windows not on the main thread??)
						dsoEditor.setState(editorContent, editorContent.getLocalVersion());
						editorSerializationState.store(1, std::memory_order_release);
					}
				);

				int exitState = 0;
				// maybe use future
				if (!cpl::Misc::WaitOnCondition(2000, [&] { return !(exitState = editorSerializationState.load(std::memory_order_acquire)); }))
					return;

				// editor death, dso should have content - no new window can exist due to lock
				if (exitState == 2)
				{
					dsoEditor.setState(editorContent, editorContent.getLocalVersion());
				}
				// if not - content should already be stored concurrently
			}
			else
			{
				// state is protected by either being on the main thread or having the window creation lock
				dsoEditor.setState(editorContent, editorContent.getLocalVersion());
			}
		}

		auto& parameterStates = serializer.getContent("Parameters");
		if (!parameterStates.isEmpty())
		{
			for (std::size_t i = 0; i < parameterMap.numSetsAndState(); ++i)
			{
				auto & serializedParameterState = parameterStates.getContent(parameterMap.getSet(i)->getName());
				if (!serializedParameterState.isEmpty())
					serializedParameterState >> *parameterMap.getState(i);
			}
		}

		auto& engineState = serializer.getContent("Engine");
		if (!engineState.isEmpty() && engineState.getLocalVersion() >= cpl::programInfo.version)
		{
			engineState >> config->historyCapacity;
		}

	}

	void AudioProcessor::serialize(cpl::CSerializer & serializer, cpl::Version version)
	{
		std::lock_guard<std::mutex> lock(editorCreationMutex);

		auto & editorContent = serializer.getContent("Editor");

		if (dsoEditor.hasCached() && !juce::MessageManager::getInstance()->isThisTheMessageThread())
		{
			// writing to this value with acquire/release ensures we see any changes concurrently here
			std::atomic_int editorSerializationState = {0};

			cpl::GUIUtils::MainEvent(
				*this,
				[&]
				{
					// interesting: host closed our window after calling getStateInformation but still while waiting for it to end..
					if (!dsoEditor.hasCached())
					{
						editorSerializationState.store(2, std::memory_order_release);
						return;
					}

					// if the object dies after this, it's a framework implementation issues (deleting windows not on the main thread??)
					editorContent = dsoEditor.getState();
					editorSerializationState.store(1, std::memory_order_release);
				}
			);

			int exitState = 0;
			// maybe use future
			if (!cpl::Misc::WaitOnCondition(2000, [&] { return !(exitState = editorSerializationState.load(std::memory_order_acquire)); }))
				return;

			// editor death, dso should have content - no new window can exist due to lock
			if (exitState == 2)
			{
				editorContent = dsoEditor.getState();
			}
			// if not - content should already be stored concurrently
		}
		else
		{
			// state is protected by either being on the main thread or having the window creation lock
			editorContent = dsoEditor.getState();

		}

		auto& parameterState = serializer.getContent("Parameters");
		parameterState.clear();
		parameterState.setMasterVersion(cpl::programInfo.version);

		for (std::size_t i = 0; i < parameterMap.numSetsAndState(); ++i)
		{
			parameterState.getContent(parameterMap.getSet(i)->getName()) << *parameterMap.getState(i);
		}

		auto& engineState = serializer.getContent("Engine");
		engineState.clear();
		engineState.setMasterVersion(cpl::programInfo.version);

		engineState << config->historyCapacity;
	}

	//==============================================================================
	const juce::String AudioProcessor::getName() const
	{
		return cpl::programInfo.name;
	}

	int AudioProcessor::getNumParameters()
	{
		return static_cast<int>(parameterMap.numParams());
	}

	float AudioProcessor::getParameter(int index)
	{
		return parameterMap.findParameter(index)->getValueNormalized<float>();

	}

	void AudioProcessor::setParameter(int index, float newValue)
	{
		return parameterMap.findParameter(index)->updateFromHostNormalized(newValue);
	}

	const juce::String AudioProcessor::getParameterName(int index)
	{
		return parameterMap.findParameter(index)->getExportedName();
	}

	const juce::String AudioProcessor::getParameterText(int index)
	{
		return parameterMap.findParameter(index)->getDisplayText();
	}

	const juce::String AudioProcessor::getInputChannelName(int channelIndex) const
	{
		return juce::String(channelIndex + 1);
	}

	const juce::String AudioProcessor::getOutputChannelName(int channelIndex) const
	{
		return juce::String(channelIndex + 1);
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

	const juce::String AudioProcessor::getProgramName(int index)
	{
		return juce::String::empty;
	}

	void AudioProcessor::changeProgramName(int index, const juce::String& newName)
	{
	}

};


//==============================================================================
// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new Signalizer::AudioProcessor();
}
