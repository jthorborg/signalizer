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

		Interface for the vectorscope view

*************************************************************************************/

#ifndef SIGNALIZER_VECTORSCOPE_H
	#define SIGNALIZER_VECTORSCOPE_H

	#include "Signalizer.h"
	#include <cpl/Utility.h>
	#include <memory>
	#include <cpl/simd.h>
	#include "../Common/ConcurrentConfig.h"

	namespace cpl
	{
		namespace OpenGLRendering
		{
			class COpenGLStack;

		};
	};

	namespace Signalizer
	{

		class VectorScopeContent;

		class VectorScope final
			: public GraphicsWindow
			, private ParameterSet::RTListener
		{

		public:

			static const double higherAutoGainBounds;
			static const double lowerAutoGainBounds;

			VectorScope(
				std::shared_ptr<const SharedBehaviour>& globalBehaviour, 
				std::shared_ptr<const ConcurrentConfig>& config,
				std::shared_ptr<AudioStream::Output>& data, 
				std::shared_ptr<VectorScopeContent>& params
			);

			virtual ~VectorScope();

			// Component overrides
			void mouseWheelMove(const juce::MouseEvent& event, const juce::MouseWheelDetails& wheel) override;
			void mouseDoubleClick(const juce::MouseEvent& event) override;
			void mouseDrag(const juce::MouseEvent& event) override;
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

		protected:

			void paint2DGraphics(juce::Graphics & g, std::size_t numChannels);
			/// <summary>
			/// Handles all set flags in mtFlags.
			/// </summary>
			void handleFlagUpdates();

		private:

			struct FilterStates
			{
				enum Entry
				{
					Slow = 0,
					Left = 0,
					Fast = 1,
					Right = 1
				};

				cpl::relaxed_atomic<AudioStream::DataType> envelope[2];
				cpl::relaxed_atomic<AudioStream::DataType> balance[2][2];
				cpl::relaxed_atomic<AudioStream::DataType> phase[2];

			};

			struct Processor : public AudioStream::Listener
			{
				Signalizer::CriticalSection<std::vector<std::string>> channelNames;
				FilterStates filters{};
				cpl::relaxed_atomic<double> envelopeGain;
				cpl::relaxed_atomic<float> 
					stereoCoeff, 
					envelopeCoeff, 
					/// <summary>
					/// A constant factor slower than the stereoCoeff
					/// </summary>
					secondStereoFilterSpeed;

				cpl::relaxed_atomic<bool> isSuspended, normalizeGain;
				cpl::relaxed_atomic<EnvelopeModes> envelopeMode;

				/// <summary>
				/// Set this if the audio buffer window size was changed from somewhere else.
				/// </summary>
				cpl::ABoolFlag audioWindowWasResized;
				std::shared_ptr<const SharedBehaviour> globalBehaviour;

				void onStreamAudio(AudioStream::ListenerContext& source, AudioStream::DataType** buffer, std::size_t numChannels, std::size_t numSamples) override;
				void onStreamPropertiesChanged(AudioStream::ListenerContext& source, const AudioStream::AudioStreamInfo& before) override;

				template<typename ISA>
				void audioProcessing(AudioStream::DataType** buffer, std::size_t numChannels, std::size_t numSamples);

				Processor(std::shared_ptr<const SharedBehaviour>& behaviour)
					: globalBehaviour(behaviour)
					, secondStereoFilterSpeed(0.25f)
					, envelopeGain(1)
				{
				}
			};

			struct AudioDispatcher
			{
				template<typename ISA> static void dispatch(Processor& v, AFloat** buffer, std::size_t numChannels, std::size_t numSamples)
				{
					v.audioProcessing<ISA>(buffer, numChannels, numSamples);
				}
			};


			struct RenderingDispatcher
			{
				template<typename ISA> static void dispatch(VectorScope & v) { v.vectorGLRendering<ISA>(); }
			};


			void parameterChangedRT(cpl::Parameters::Handle localHandle, cpl::Parameters::Handle globalHandle, ParameterSet::BaseParameter * param) override;
			void deserialize(cpl::CSerializer::Builder & builder, cpl::Version version) override {};
			void serialize(cpl::CSerializer::Archiver & archive, cpl::Version version) override {};

			template<typename ISA>
				void vectorGLRendering();

			// vector-accelerated drawing, rendering and processing
			template<typename ISA>
				void drawPolarPlot(cpl::OpenGLRendering::COpenGLStack &, const AudioStream::AudioBufferAccess &, std::size_t);

			template<typename ISA>
				void drawRectPlot(cpl::OpenGLRendering::COpenGLStack &, const AudioStream::AudioBufferAccess &, std::size_t);

			template<typename ISA>
				void drawWireFrame(cpl::OpenGLRendering::COpenGLStack &);

			template<typename ISA>
				void drawGraphText(cpl::OpenGLRendering::COpenGLStack &, const AudioStream::AudioBufferAccess &);

			template<typename ISA>
				void drawStereoMeters(cpl::OpenGLRendering::COpenGLStack &, const AudioStream::AudioBufferAccess &);

			template<typename ISA>
				void runPeakFilter(const AudioStream::AudioBufferAccess &);

			void initPanelAndControls();

			struct Flags
			{
				cpl::ABoolFlag
					firstRun,
					/// <summary>
					/// Set this to resize the audio windows (like, when the audio window size (as in fft size) is changed.
					/// The argument to the resizing is the state.newWindowSize
					/// </summary>
					initiateWindowResize;
			} mtFlags;

			// contains non-atomic structures
			struct StateOptions
			{
				bool isPolar, isFrozen, fillPath, fadeHistory, antialias, diagnostics, drawLegend, scalePolar;
				float primitiveSize, rotation;
				juce::Colour colourBackground, colourWire, colourAxis, colourWaveform, colourMeter, colourWidget;
				cpl::ValueT userGain;
			} state;

			std::shared_ptr<VectorScopeContent> content;
			std::shared_ptr<const ConcurrentConfig> config;
			std::shared_ptr<AudioStream::Output> audioStream;
			std::shared_ptr<Processor> processor;
			std::shared_ptr<const SharedBehaviour> globalBehaviour;

			std::vector<std::unique_ptr<juce::OpenGLTexture>> textures;
		};

	};

#endif
