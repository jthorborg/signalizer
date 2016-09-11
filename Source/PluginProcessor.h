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
#include "MainEditor.h"

namespace Signalizer
{
	class MainEditor;

	class AudioProcessor final
		: public juce::AudioProcessor
		, public cpl::SafeSerializableObject
		, cpl::DestructionNotifier
		, ParameterSet::AutomatedProcessor
	{
		friend class MainEditor;

	public:

		typedef cpl::CPresetWidget::SerializerType SerializerType;

		//==============================================================================
		AudioProcessor();
		~AudioProcessor() noexcept;

		//==============================================================================
		void prepareToPlay(double sampleRate, int samplesPerBlock) override;
		void releaseResources() override ;

		void processBlock(AudioSampleBuffer& buffer, MidiBuffer& midiMessages) override;

		//==============================================================================
		AudioProcessorEditor* createEditor() override;
		bool hasEditor() const override;

		//==============================================================================
		const String getName() const override;

		int getNumParameters() override;

		float getParameter(int index) override;
		void setParameter(int index, float newValue) override;

		const String getParameterName(int index) override;
		const String getParameterText(int index) override;

		const String getInputChannelName(int channelIndex) const override;
		const String getOutputChannelName(int channelIndex) const override;
		bool isInputChannelStereoPair(int index) const override;
		bool isOutputChannelStereoPair(int index) const override;

		bool acceptsMidi() const override;
		bool producesMidi() const override;
		bool silenceInProducesSilenceOut() const override;
		double getTailLengthSeconds() const override;

		//==============================================================================
		int getNumPrograms() override;
		int getCurrentProgram() override;
		void setCurrentProgram(int index) override;
		const String getProgramName(int index) override;
		void changeProgramName(int index, const String& newName) override;

		//==============================================================================
		void getStateInformation(MemoryBlock& destData) override;
		void setStateInformation(const void* data, int sizeInBytes) override;

		void deserialize(cpl::CSerializer & se, cpl::Version version) override;
		void serialize(cpl::CSerializer & se, cpl::Version version) override;

	private:

		virtual void automatedTransmitChangeMessage(int parameter, ParameterSet::FrameworkType value) override;
		virtual void automatedBeginChangeGesture(int parameter) override;
		virtual void automatedEndChangeGesture(int parameter) override;

		//==============================================================================
		JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioProcessor)
		Signalizer::AudioStream stream;
		int nChannels;
		ParameterMap parameterMap;
		DecoupledStateObject<MainEditor> dsoEditor;
		std::mutex editorCreationMutex;
	};
};
#endif 
