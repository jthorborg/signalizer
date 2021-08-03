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

	file:VectorScope.h

		Interface for the Oscilloscope view

*************************************************************************************/

#ifndef SIGNALIZER_OSCILLOSCOPE_H
	#define SIGNALIZER_OSCILLOSCOPE_H

	#include "Signalizer.h"
	#include <cpl/Utility.h>
	#include <memory>
	#include <cpl/simd.h>
	#include "OscilloscopeParameters.h"
	#include <cpl/dsp/SmoothedParameterState.h>
	#include <utility>
	#include "ChannelData.h"
	#include <cpl/gui/CViews.h>

	namespace cpl
	{
		namespace OpenGLRendering
		{
			class COpenGLStack;
		};
	};

	namespace Signalizer
	{

		struct EvaluatorParams
		{
			ChannelData& data;

			std::size_t
				channelIndex = 0,
				colourIndex = 0;
		};

		class TriggeringProcessor;

		class Oscilloscope final
			: public cpl::COpenGLView
			, private AudioStream::Listener
		{
		public:

			friend class TriggeringProcessor;

			static const double higherAutoGainBounds;
			static const double lowerAutoGainBounds;

			Oscilloscope(const SharedBehaviour & globalBehaviour, const std::string & nameId, AudioStream & data, ProcessorState * params);
			virtual ~Oscilloscope();

		protected:
			
			// Component overrides
			void onGraphicsRendering(juce::Graphics & g) override;

			// OpenGLRender overrides
			void onOpenGLRendering() override;
			void initOpenGL() override;
			void closeOpenGL() override;
			// View overrides
			juce::Component * getWindow() override;
			void suspend() override;
			void resume() override;
			void freeze() override;
			void unfreeze() override;

			// cbasecontrol overrides
			bool isEditorOpen() const;

			void mouseWheelMove(const juce::MouseEvent& event, const juce::MouseWheelDetails& wheel) override;
			void mouseDoubleClick(const juce::MouseEvent& event) override;
			void mouseDrag(const juce::MouseEvent& event) override;
			void mouseUp(const juce::MouseEvent& event) override;
			void mouseDown(const juce::MouseEvent& event) override;
			void mouseMove(const juce::MouseEvent& event) override;
			void mouseExit(const juce::MouseEvent & e) override;
			void mouseEnter(const juce::MouseEvent & e) override;

			bool onAsyncAudio(const AudioStream & source, AudioStream::DataType ** buffer, std::size_t numChannels, std::size_t numSamples) override;
			//void onAsyncChangedProperties(const AudioStream & source, const AudioStream::AudioStreamInfo & before) override;

			/// <summary>
			/// Handles all set flags in mtFlags.
			/// Make sure the audioStream is locked while doing this.
			/// </summary>
			virtual void handleFlagUpdates();

		private:
			void deserialize(cpl::CSerializer::Builder & builder, cpl::Version version) override {};
			void serialize(cpl::CSerializer::Archiver & archive, cpl::Version version) override {};

			struct RenderingDispatcher
			{
				template<typename ISA> static void dispatch(Oscilloscope & o) { o.vectorGLRendering<ISA>(); }
			};

			struct AudioDispatcher
			{
				template<typename ISA> static void dispatch(Oscilloscope & o, AFloat ** buffer, std::size_t numChannels, std::size_t numSamples) 
				{ 
					o.audioEntryPoint<ISA>(buffer, numChannels, numSamples);
				}
			};

			std::size_t getEffectiveChannels() const noexcept;

			template<typename ISA>
				void vectorGLRendering();

			// vector-accelerated drawing, rendering and processing
			template<typename ISA, typename Eval>
				void drawWavePlot(cpl::OpenGLRendering::COpenGLStack &, const EvaluatorParams& params);

			template<typename ISA>
				void drawWireFrame(juce::Graphics & g, juce::Rectangle<float> rect, float gain);

			template<typename ISA>
				void drawTimeDivisions(juce::Graphics & g, juce::Rectangle<float> rect);

			template<typename ISA, typename Eval>
				void calculateFundamentalPeriod(const EvaluatorParams& params);

			template<typename ISA, typename Eval>
				void calculateTriggeringOffset(const EvaluatorParams& params);

			void resizeAudioStorage();

			template<typename ISA>
				void runPeakFilter();

			template<typename ISA, typename Eval>
				void analyseAndSetupState(const EvaluatorParams& params);

			template<typename ISA>
			void preAnalyseAudio(AFloat ** buffer, std::size_t numChannels, std::size_t numSamples);

			template<typename ISA, class Analyzer>
				void executeSamplingWindows(AFloat ** buffer, std::size_t numChannels, std::size_t numSamples);

			template<typename ISA>
				void audioEntryPoint(AFloat ** buffer, std::size_t numChannels, std::size_t numSamples);

			template<typename ISA>
				void audioProcessing(AFloat ** buffer, std::size_t numChannels, std::size_t numSamples, ChannelData::Buffer &);

			template<typename ISA>
				void paint2DGraphics(juce::Graphics & g);

			inline AFloat& smoothEnvelopeState(std::size_t i) { return channelData.filterStates.channels[i].envelope; }

			bool checkAndInformInvalidCombinations();

			void initPanelAndControls();

			/// <summary>
			/// Returns the combined gain value for the current frame.
			/// </summary>
			double getGain();

			void setLastMousePos(const juce::Point<float> position) noexcept;

			struct Flags
			{
				cpl::ABoolFlag
					firstRun,
					/// <summary>
					/// Set this to resize the audio windows (like, when the audio window size (as in fft size) is changed.
					/// The argument to the resizing is the state.newWindowSize
					/// </summary>
					initiateWindowResize,
					/// <summary>
					/// Set this if the audio buffer window size was changed from somewhere else.
					/// </summary>
					audioWindowWasResized;
			} mtFlags;

			// contains frame-updated non-atomic structures
			struct StateOptions
			{
				bool isFrozen, antialias, diagnostics, dotSamples, customTrigger, overlayChannels, colourChannelsByFrequency, drawCursorTracker, isSuspended, drawLegend;
				float primitiveSize;

				double effectiveWindowSize;
				double windowTimeOffset;
				double beatDivision;
				double customTriggerFrequency;
				double triggerHysteresis;
				double triggerThreshold;
				double autoGain, manualGain;

				double viewOffsets[4];
				std::int64_t transportPosition;
				std::size_t numChannels;
				juce::Colour colourBackground, colourGraph, colours[OscilloscopeContent::NumColourChannels], colourSecondary, colourTracker;

				EnvelopeModes envelopeMode;
				SubSampleInterpolation sampleInterpolation;
				OscilloscopeContent::TriggeringMode triggerMode;
				OscilloscopeContent::TimeMode timeMode;
				OscChannels channelMode;

				std::vector<std::string> channelNames;

			} state;

			struct SharedStateOptions
			{
				std::atomic<double> 
					autoGainEnvelope;
			} shared;

			using VO = OscilloscopeContent::ViewOffsets;
			cpl::CBoxFilter<double, 60> avgFps;
			juce::MouseCursor displayCursor;
			OscilloscopeContent * content;
			AudioStream & audioStream;
			//cpl::AudioBuffer audioStreamCopy;
			juce::Component * editor;
			std::vector<std::string> channelNames;
			std::unique_ptr<char> textbuf;
			unsigned long long processorSpeed; // clocks / sec
			juce::Point<float> lastMousePos;
			long long lastFrameTick, renderCycles;
			std::atomic_bool isMouseInside;
			/// <summary>
			/// Updates are not guaranteed to be in order
			/// </summary>
			std::pair<std::atomic<float>, std::atomic<float>> threadedMousePos;
			cpl::aligned_vector<std::complex<double>, 32> transformBuffer;
			cpl::aligned_vector<double, 16> temporaryBuffer;
			const SharedBehaviour & globalBehaviour;

			struct BinRecord
			{
				/// <summary>
				/// The fundamental frequency, quantized to the originating DFT bin
				/// </summary>
				std::size_t index;

				double 
					/// <summary>
					/// Value of this bin
					/// </summary>
					value, 
					/// <summary>
					/// The offset quantized bin creating the fundamental
					/// </summary>
					offset;

				double omega() const noexcept { return index + offset; }
			};

			struct TriggerData
			{ 
				std::unique_ptr<TriggeringProcessor> triggeringProcessor;
				/// <summary>
				/// The fundamental frequency (in hertz) in the selected window offset in time.
				/// </summary>
				double fundamental;
				BinRecord record{};
				/// <summary>
				/// The amount of samples per fundamental period
				/// </summary>
				double cycleSamples;
				/// <summary>
				/// The phase offset for the detected fundamental at T = 0 - state.effectiveWindowSize
				/// </summary>
				double phase;
				/// <summary>
				/// The sample offset for the detected fundamental at T = 0 - state.effectiveWindowSize
				/// </summary>
				double sampleOffset;
			} triggerState;


			struct MedianData
			{
				// must be a power of two
				static const std::size_t FilterSize = 8;
				BinRecord record{};
			};


			std::size_t medianPos;
			std::array<MedianData, MedianData::FilterSize> medianTriggerFilter;

			cpl::CMutex::Lockable bufferLock;
			ChannelData channelData;

			class DefaultKey;

			class SampleColourEvaluatorBase
			{
			public:

				typedef ChannelData::AudioBuffer::ProxyView::const_iterator AudioIt;
				typedef ChannelData::ColourBuffer::ProxyView::const_iterator ColourIt;
				typedef ChannelData::AudioBuffer::ProxyView::value_type AudioT;
				typedef ChannelData::ColourBuffer::ProxyView::value_type ColourT;

			};

			template<OscChannels channelConfiguration, std::size_t ColourIndice>
				class SampleColourEvaluator;

			class DynamicChannelEvaluator;

			template<std::size_t ChannelIndex, std::size_t ColourIndice>
				class SimpleChannelEvaluator;

			template<std::size_t ChannelIndex, std::size_t ColourIndice, typename BinaryFunction>
				class MidSideEvaluatorBase;
		};

	};

#endif
