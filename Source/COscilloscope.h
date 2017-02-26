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

			COscilloscope(const std::string & nameId, AudioStream & data, ProcessorState * params);
			virtual ~COscilloscope();

			// Component overrides
			void onGraphicsRendering(juce::Graphics & g) override;
			void mouseWheelMove(const juce::MouseEvent& event, const juce::MouseWheelDetails& wheel) override;
			void mouseDoubleClick(const juce::MouseEvent& event) override;
			void mouseDrag(const juce::MouseEvent& event) override;
			void mouseUp(const juce::MouseEvent& event) override;
			void mouseDown(const juce::MouseEvent& event) override;
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
			double getGain();

		protected:

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

			template<typename V>
			void vectorGLRendering();

			// vector-accelerated drawing, rendering and processing
			template<typename V>
				void drawWavePlot(cpl::OpenGLRendering::COpenGLStack &);

			template<typename V>
				void drawWireFrame(juce::Graphics & g, juce::Rectangle<float> rect, float gain);

			template<typename V>
				void drawTimeDivisions(juce::Graphics & g, juce::Rectangle<float> rect, double horizontalGranularity);

			void calculateFundamentalPeriod();
			void calculateTriggeringOffset();
			void resizeAudioStorage();

			template<typename V>
				void runPeakFilter(const AudioStream::AudioBufferAccess &);

			template<typename V>
				void audioProcessing(typename cpl::simd::scalar_of<V>::type ** buffer, std::size_t numChannels, std::size_t numSamples);

			template<typename V>
				void paint2DGraphics(juce::Graphics & g);

			void initPanelAndControls();

			// guis and whatnot
			cpl::CBoxFilter<double, 60> avgFps;

			juce::MouseCursor displayCursor;

			// vars
			long long lastFrameTick, renderCycles;

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

			// contains non-atomic structures
			struct StateOptions
			{
				bool normalizeGain, isFrozen, antialias, diagnostics, dotSamples, customTrigger, overlayChannels, colourChannelsByFrequency;
				float primitiveSize;
				float envelopeCoeff;
				double effectiveWindowSize;
				double windowTimeOffset;
				double beatDivision;
				double customTriggerFrequency;

				double viewOffsets[4];
				std::int64_t transportPosition;
				juce::Colour colourBackground, colourGraph, colourPrimary, colourSecondary;
				cpl::ValueT envelopeGain;
				EnvelopeModes envelopeMode;
				SubSampleInterpolation sampleInterpolation;
				OscilloscopeContent::TriggeringMode triggerMode;
				OscilloscopeContent::TimeMode timeMode;
			} state;

			using VO = OscilloscopeContent::ViewOffsets;

			OscilloscopeContent * content;
			AudioStream & audioStream;
			//cpl::AudioBuffer audioStreamCopy;
			juce::Component * editor;
			// unused.
			std::unique_ptr<char> textbuf;
			unsigned long long processorSpeed; // clocks / sec
			juce::Point<float> lastMousePos;
			cpl::aligned_vector<std::complex<double>, 32> transformBuffer;
			cpl::aligned_vector<double, 16> temporaryBuffer;

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
			cpl::CLIFOStream<float, 16> lifoStream;

		};

	};

#endif
