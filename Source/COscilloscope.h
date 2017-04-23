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

	file:CVectorScope.h

		Interface for the Oscilloscope view

*************************************************************************************/

#ifndef SIGNALIZER_COSCILLOSCOPE_H
	#define SIGNALIZER_COSCILLOSCOPE_H

	#include "CommonSignalizer.h"
	#include <cpl/Utility.h>
	#include <cpl/gui/controls/Controls.h>
	#include <cpl/gui/widgets/Widgets.h>
	#include <memory>
	#include <cpl/simd.h>
	#include "OscilloscopeParameters.h"
	#include <cpl/dsp/LinkwitzRileyNetwork.h>
	#include <cpl/dsp/SmoothedParameterState.h>
	#include <utility>
	#include "SharedBehaviour.h"
	#include "StreamPreprocessing.h"

	namespace cpl
	{
		namespace OpenGLRendering
		{
			class COpenGLStack;
		};
	};

	namespace Signalizer
	{


		class COscilloscope final
			: public cpl::COpenGLView
			, private AudioStream::Listener
		{
		public:

			static const double higherAutoGainBounds;
			static const double lowerAutoGainBounds;

			COscilloscope(const SharedBehaviour & globalBehaviour, const std::string & nameId, AudioStream & data, ProcessorState * params);
			virtual ~COscilloscope();

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
				template<typename ISA> static void dispatch(COscilloscope & o) { o.vectorGLRendering<ISA>(); }
			};

			struct AudioDispatcher
			{
				template<typename ISA> static void dispatch(COscilloscope & o, AFloat ** buffer, std::size_t numChannels, std::size_t numSamples) 
				{ 

					if (numChannels > 1)
					{
						AFloat * localPointers[2];
						localPointers[0] = buffer[0];
						localPointers[1] = buffer[1];
						o.preprocessAudio<ISA>(localPointers, 2, numSamples);
						o.audioProcessing<ISA>(localPointers, 2, numSamples);
					}
					else
					{
						AFloat * localBuffer = buffer[0];
						o.preprocessAudio<ISA>(&localBuffer, 1, numSamples);
						o.audioProcessing<ISA>(&localBuffer, 1, numSamples);
					}

				}
			};


			template<typename ISA>
				void vectorGLRendering();

			// vector-accelerated drawing, rendering and processing
			template<typename ISA, typename Eval>
				void drawWavePlot(cpl::OpenGLRendering::COpenGLStack &);

			template<typename ISA>
				void drawWireFrame(juce::Graphics & g, juce::Rectangle<float> rect, float gain);

			template<typename ISA>
				void drawTimeDivisions(juce::Graphics & g, juce::Rectangle<float> rect);

			template<typename ISA, typename Eval>
				void calculateFundamentalPeriod();

			template<typename ISA, typename Eval>
				void calculateTriggeringOffset();

			void resizeAudioStorage();

			template<typename ISA>
				void runPeakFilter();

			template<typename ISA, typename Eval>
				void analyseAndSetupState();

			template<typename ISA>
				void preprocessAudio(AFloat ** buffer, std::size_t numChannels, std::size_t & numSamples);

			template<typename ISA, class Analyzer>
				void executeSamplingWindows(AFloat ** buffer, std::size_t numChannels, std::size_t & numSamples);

			template<typename ISA>
				void audioProcessing(AFloat ** buffer, std::size_t numChannels, std::size_t numSamples);

			template<typename ISA>
				void paint2DGraphics(juce::Graphics & g);

			bool checkAndInformInvalidCombinations();

			void initPanelAndControls();

			/// <summary>
			/// Returns the combined gain value for the current frame.
			/// </summary>
			double getGain();

			void setLastMousePos(const juce::Point<float> position) noexcept;

			struct FilterStates
			{
				enum Entry
				{
					Slow = 0,
					Left = 0,
					Fast = 1,
					Right = 1
				};

				AudioStream::DataType envelope[2];

			} filters;

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
				bool isFrozen, antialias, diagnostics, dotSamples, customTrigger, overlayChannels, colourChannelsByFrequency, drawCursorTracker, isSuspended;
				float primitiveSize;

				double effectiveWindowSize;
				double windowTimeOffset;
				double beatDivision;
				double customTriggerFrequency;

				double autoGain, manualGain;

				double viewOffsets[4];
				std::int64_t transportPosition;
				juce::Colour colourBackground, colourGraph, colourPrimary, colourSecondary, colourTracker;

				EnvelopeModes envelopeMode;
				SubSampleInterpolation sampleInterpolation;
				OscilloscopeContent::TriggeringMode triggerMode;
				OscilloscopeContent::TimeMode timeMode;
				OscChannels channelMode;

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
			// unused.
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
				PreprocessingTriggerState preTriggerState{};
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

			struct ChannelData
			{
				static const std::size_t Bands = 3;
				typedef cpl::dsp::LinkwitzRileyNetwork<AFloat, Bands> Crossover;
				typedef cpl::GraphicsND::UPixel<cpl::GraphicsND::ComponentOrder::OpenGL> PixelType;
				typedef cpl::CLIFOStream<AFloat, 32> AudioBuffer;
				typedef cpl::CLIFOStream<PixelType, 32> ColourBuffer;

				struct Channel
				{
					AudioBuffer audioData;
					ColourBuffer colourData;
					Crossover::BandArray smoothFilters{};
					Crossover network;
					juce::Colour defaultKey;
				};

				Channel & defaultChannel()
				{
					return channels[0];
				}

				void resizeStorage(std::size_t samples, std::size_t capacity = -1)
				{
					if (capacity == -1)
						capacity = cpl::Math::nextPow2Inc(samples);

					for (auto & c : channels)
					{
						c.audioData.setStorageRequirements(samples, capacity);
						c.colourData.setStorageRequirements(samples, capacity);
					}

					for (auto & c : midSideColour)
					{
						c.setStorageRequirements(samples, capacity);
					}
				}

				void resizeChannels(std::size_t newChannels)
				{
					if (newChannels > channels.size())
					{
						auto original = channels.size();
						auto diff = newChannels - original;

						while (diff--)
						{
							channels.emplace_back();
						}

						if (original > 0)
							resizeStorage(channels.back().audioData.getSize(), channels.back().audioData.getCapacity());
					}
				}

				void tuneCrossOver(double lowCrossover, double highCrossover, double sampleRate)
				{
					networkCoeffs = Crossover::Coefficients::design({ static_cast<AFloat>(lowCrossover / sampleRate), static_cast<AFloat>(highCrossover / sampleRate) });
				}

				void tuneColourSmoothing(double milliseconds, double sampleRate)
				{
					smoothFilterPole = cpl::dsp::SmoothedParameterState<AFloat, 1>::design(milliseconds, sampleRate);
				}

				Crossover::Coefficients networkCoeffs;
				cpl::dsp::SmoothedParameterState<AFloat, 1>::PoleState smoothFilterPole;
				std::vector<Channel> channels{ 1 };
				ColourBuffer midSideColour[2];
				Crossover::BandArray midSideSmoothsFilters[2];
			};

			std::size_t medianPos;
			std::array<MedianData, MedianData::FilterSize> medianTriggerFilter;

			cpl::CMutex::Lockable bufferLock;
			ChannelData channelData;

			class DefaultKey;

			class SampleColourEvaluatorBase
			{
			public:

				typedef COscilloscope::ChannelData ChannelData;

				typedef ChannelData::AudioBuffer::ProxyView::const_iterator AudioIt;
				typedef ChannelData::ColourBuffer::ProxyView::const_iterator ColourIt;
				typedef ChannelData::AudioBuffer::ProxyView::value_type AudioT;
				typedef ChannelData::ColourBuffer::ProxyView::value_type ColourT;

			};

			template<OscChannels channelConfiguration, std::size_t ColourIndice>
				class SampleColourEvaluator;

			template<std::size_t ChannelIndex, std::size_t ColourIndice>
				class SimpleChannelEvaluator;

			template<std::size_t ChannelIndex, std::size_t ColourIndice, typename BinaryFunction>
				class MidSideEvaluatorBase;
		};

	};

#endif
