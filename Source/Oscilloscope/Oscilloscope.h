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
	#include <cpl/ffts.h>

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
				channelIndex = 0;
		};

		class TriggeringProcessor;

		class Oscilloscope final
			: public GraphicsWindow
			, private AudioStream::Listener
		{
		public:

			friend class TriggeringProcessor;

			static const double higherAutoGainBounds;
			static const double lowerAutoGainBounds;

			Oscilloscope(
				std::shared_ptr<const SharedBehaviour>& globalBehaviour,
				std::shared_ptr<const ConcurrentConfig>& config,
				std::shared_ptr<AudioStream::Output>& stream,
				std::shared_ptr<OscilloscopeContent> params
			);

			virtual ~Oscilloscope();

		protected:
			
			// OpenGLRender overrides
			void onOpenGLRendering() override;
			// View overrides
			juce::Component * getWindow() override;
			void suspend() override;
			void resume() override;

			void mouseWheelMove(const juce::MouseEvent& event, const juce::MouseWheelDetails& wheel) override;
			void mouseDoubleClick(const juce::MouseEvent& event) override;
			void mouseDrag(const juce::MouseEvent& event) override;

		private:

			struct Flags
			{
				cpl::ABoolFlag
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

			// TODO: scope to gfx?
			// contains frame-updated non-atomic structures
			struct StateOptions
			{
				bool antialias, diagnostics, dotSamples, customTrigger, colourChannelsByFrequency, drawCursorTracker, drawLegend;
				float primitiveSize;

				OscilloscopeContent::TriggeringMode triggerMode;

				double effectiveWindowSize;
				double windowTimeOffset;
				double beatDivision;
				double customTriggerFrequency;
				double triggerHysteresis;
				double triggerThreshold;
				double autoGain, manualGain;
				double sampleRate;

				double viewOffsets[4];
				juce::Colour colourBackground, colourAxis, colourSecondary, colourWidget;

				SubSampleInterpolation sampleInterpolation;
				OscilloscopeContent::TimeMode timeMode;
				ChangeVersion::Listener audioStreamChanged;
				LegendCache legend;

			} state {};

			struct SharedStateOptions
			{
				// TODO: accessed from gfx and main thread
				cpl::relaxed_atomic<std::size_t> numChannels;
				cpl::relaxed_atomic<OscChannels> channelMode;
				cpl::relaxed_atomic<bool> overlayChannels;

			} shared {};

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

				/// <summary>
				/// The fundamental frequency (in hertz) in the selected window offset in time.
				/// </summary>
				double fundamental;
				BinRecord record;
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
			} triggerState{};

			struct StreamState
			{
				friend struct AudioDispatcher;
			public:

				StreamState();
				~StreamState();

				ChannelData channelData;
				// TODO: Should go back to be in place
				std::unique_ptr<TriggeringProcessor> triggeringProcessor;
				std::shared_ptr<OscilloscopeContent> content;
				std::vector<std::string> channelNames;
				std::int64_t historyCapacity;
				std::int64_t transportPosition;
				double bpm {};
				double sampleRate{};
				double envelopeGain{};
				EnvelopeModes envelopeMode;
				OscChannels channelMode;

				OscilloscopeContent::TriggeringMode triggerMode;
				ChangeVersion audioStreamChangeVersion;

				template<typename ISA>
				void preAnalyseAudio(AudioStream::ListenerContext& ctx, AFloat** buffer, std::size_t numChannels, std::size_t numSamples);

				template<typename ISA, class Analyzer>
				void executeSamplingWindows(AudioStream::ListenerContext& ctx, AFloat** buffer, std::size_t numChannels, std::size_t numSamples);

				template<typename ISA>
				void audioProcessing(
					const AudioStream::Info& info,
					const AudioStream::Playhead& playhead,
					const AudioStream::DataType* const* buffer,
					const std::size_t numChannels,
					const std::size_t numSamples,
					ChannelData::Buffer&
				);

				template<typename ISA>
				void audioEntryPoint(AudioStream::ListenerContext& source, AudioStream::DataType** buffer, std::size_t numChannels, std::size_t numSamples);
			};

			void deserialize(cpl::CSerializer::Builder & builder, cpl::Version version) override {};
			void serialize(cpl::CSerializer::Archiver & archive, cpl::Version version) override {};

			struct RenderingDispatcher
			{
				template<typename ISA> static void dispatch(Oscilloscope & o) { o.vectorGLRendering<ISA>(); }
			};


			std::size_t getEffectiveChannels() const noexcept;

			template<typename ISA>
				void vectorGLRendering();

			// vector-accelerated drawing, rendering and processing
			template<typename ISA, typename Eval>
				void drawWavePlot(cpl::OpenGLRendering::COpenGLStack &, const EvaluatorParams& params, StreamState& cs);

			template<typename ISA>
				void drawWireFrame(juce::Graphics & g, juce::Rectangle<float> rect, float gain);

			template<typename ISA>
				void drawTimeDivisions(juce::Graphics & g, juce::Rectangle<float> rect);

			template<typename ISA, typename Eval>
				void calculateFundamentalPeriod(const EvaluatorParams& params);

			template<typename ISA, typename Eval>
				void calculateTriggeringOffset(const EvaluatorParams& params);

			template<typename ISA>
				void runPeakFilter(ChannelData& data);

			template<typename ISA, typename Eval>
				void analyseAndSetupState(const EvaluatorParams& params, Oscilloscope::StreamState& cs);

			template<typename ISA>
				void paint2DGraphics(juce::Graphics & g);

			bool checkAndInformInvalidCombinations(Oscilloscope::StreamState&);

			void initPanelAndControls();

			/// <summary>
			/// Returns the combined gain value for the current frame.
			/// </summary>
			double getGain();

			void handleFlagUpdates(StreamState&);
			void recalculateLegend(Oscilloscope::StreamState& cs, ColourRotation primaryRotation, ColourRotation secondaryRotation);

			struct ProcessorShell : public AudioStream::Listener
			{
				std::shared_ptr<const SharedBehaviour> globalBehaviour;
				CriticalSection<StreamState> streamState;
				cpl::relaxed_atomic<bool> isSuspended;

				void onStreamAudio(AudioStream::ListenerContext& source, AudioStream::DataType** buffer, std::size_t numChannels, std::size_t numSamples) override;
				void onStreamPropertiesChanged(AudioStream::ListenerContext& source, const AudioStream::AudioStreamInfo& before) override;

				ProcessorShell(std::shared_ptr<const SharedBehaviour>& behaviour)
					: globalBehaviour(behaviour)
				{
				}
			};

			struct AudioDispatcher
			{
				template<typename ISA> static void dispatch(ProcessorShell& shell, AudioStream::ListenerContext& source, AudioStream::DataType** buffer, std::size_t numChannels, std::size_t numSamples)
				{
					shell.streamState.lock()->audioEntryPoint<ISA>(source, buffer, numChannels, numSamples);
				}
			};

			struct MedianData
			{
				// must be a power of two
				static const std::size_t FilterSize = 8;
				BinRecord record{};
			};

			using VO = OscilloscopeContent::ViewOffsets;
			using ChannelFloat = ChannelData::AudioBuffer::ProxyView::value_type;
			std::shared_ptr<OscilloscopeContent> content;
			std::shared_ptr<AudioStream::Output> audioStream;

			cpl::aligned_vector<std::complex<ChannelFloat>, 32> transformBuffer, triggerWork;
			std::shared_ptr<const SharedBehaviour> globalBehaviour;
			std::size_t medianPos;
			std::array<MedianData, MedianData::FilterSize> medianTriggerFilter;
			cpl::dsp::UniFFT<std::complex<ChannelFloat>> triggerFFT;

			std::shared_ptr<ProcessorShell> processor;

			class DefaultKey;

			class SampleColourEvaluatorBase
			{
			public:
				typedef ChannelData::AudioBuffer::ProxyView::const_iterator AudioIt;
				typedef ChannelData::ColourBuffer::ProxyView::const_iterator ColourIt;
				typedef ChannelFloat AudioT;
				typedef ChannelData::ColourBuffer::ProxyView::value_type ColourT;
			};

			template<OscChannels ChannelConfiguration>
				class SampleColourEvaluator;

			class DynamicChannelEvaluator;

			template<std::size_t ChannelOffset>
				class SimpleChannelEvaluator;

			template<std::size_t ChannelOffset, typename BinaryFunction>
				class MidSideEvaluatorBase;
		};

	};

#endif
