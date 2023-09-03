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

	file:Spectrum.h

		Interface for the spectrum view.

*************************************************************************************/

#ifndef SIGNALIZER_SPECTRUM_H
	#define SIGNALIZER_SPECTRUM_H

	#include "../Signalizer.h"
	#include <cpl/Common.h>
	#include <cpl/gui/CViews.h>
	#include <cpl/Utility.h>
	#include <cpl/dsp/CComplexResonator.h>
	#include <memory>
	#include <cpl/rendering/COpenGLImage.h>
	#include <cpl/special/AxisTools.h>
	#include <cpl/dsp/CPeakFilter.h>
	#include <queue>
	#include <cpl/lib/LockFreeQueue.h>
	#include <cpl/lib/BlockingLockFreeQueue.h>
	#include <vector>
	#include "SpectrumParameters.h"
	#include <cpl/dsp/SmoothedParameterState.h>
	#include "TransformConstant.h"
	#include "TransformPair.h"
	#include <cpl/lib/LockFreeDataQueue.h>

	#define USE_DUST_FFT

	namespace cpl
	{
		namespace OpenGLRendering
		{
			class COpenGLStack;

		};
	};


	namespace Signalizer
	{

		class Spectrum final
			: public GraphicsWindow
			, protected AudioStream::Listener
			, private ParameterSet::RTListener
		{
		public:
			typedef AudioStream::DataType fpoint;

			typedef TransformPair<float> TransformPair;
			typedef TransformPair::Constant Constant;
			typedef TransformPair::ProcessingType ProcessingType;

			struct DBRange
			{
				double low, high;
			};

			Spectrum(
				std::shared_ptr<const SharedBehaviour>& globalBehaviour,
				std::shared_ptr<const ConcurrentConfig>& config,
				std::shared_ptr<AudioStream::Output>& stream,
				std::shared_ptr<SpectrumContent>& params
			);
			virtual ~Spectrum();

			// Component overrides
			void mouseWheelMove(const juce::MouseEvent& event, const juce::MouseWheelDetails& wheel) override;
			void mouseDoubleClick(const juce::MouseEvent& event) override;
			void mouseDrag(const juce::MouseEvent& event) override;

			void resized() override;


			// OpenGLRender overrides
			void onOpenGLRendering() override;
			void initOpenGL() override;
			void closeOpenGL() override;

			// View overrides
			void suspend() override;
			void resume() override;
			void freeze() override;
			void unfreeze() override;
			void resetState() override;

			// interface


			void setDBs(double low, double high, bool updateControls = false);
			DBRange getDBs() const noexcept;

			/// <summary>
			/// Returns the working size of the audio buffer.
			/// It is guaranteed to be 0 > getWindowSize() <= audioStream.getAudioHistorySize()
			/// </summary>
			/// <returns></returns>
			std::size_t getWindowSize() const noexcept;
			void setWindowSize(std::size_t size);

		protected:
			typedef TransformPair::UComplex UComplex;

			struct RenderingDispatcher
			{
				template<typename ISA> static void dispatch(Spectrum & c) { c.vectorGLRendering<ISA>(); }
			};

            template<typename ISA>
                void vectorGLRendering();

			virtual void paint2DGraphics(juce::Graphics & g, const Constant& constant, /*const*/ TransformPair& primaryTransform);

			virtual void parameterChangedRT(cpl::Parameters::Handle localHandle, cpl::Parameters::Handle globalHandle, ParameterSet::BaseParameter * param) override;

			void calculateSpectrumColourRatios(Constant& constant);
		private:
			typedef cpl::GraphicsND::UPixel<cpl::GraphicsND::ComponentOrder::OpenGL> UPixel;
			typedef std::array<juce::Colour, SpectrumContent::LineGraphs::LineEnd> LineColours;
			typedef struct StreamState;
			typedef cpl::aligned_vector<UPixel, 16> FrameVector;
			typedef cpl::CLockFreeDataQueue<FrameVector> SFrameQueue;

			std::size_t getValidWindowSize(std::size_t in) const noexcept;

			/// <summary>
			/// Returns the number of needed channels required to process the current
			/// channel configuration.
			/// </summary>
			std::size_t getStateConfigurationChannels() const noexcept;

			/// <summary>
			/// Transforms state.audioBlobSizeMs into samples.
			/// </summary>
			std::size_t getBlobSamples() const noexcept;

			/// <summary>
			/// Returns the samplerate of the currently connected channels.
			/// </summary>
			/// <returns></returns>
			float getSampleRate() const noexcept;
			/// <summary>
			/// The number of filters used for algoritms not based on 2^n ffts
			/// This number is generally in relation to the number of pixels on the axis.
			/// </summary>
			int getNumFilters() const noexcept;
			/// <summary>
			/// The number of pixels/points in the frequency axis.
			/// </summary>
			int getAxisPoints() const noexcept;
			/// <summary>
			/// For a certain blob size and refresh rate, N frames has to be generated each display update.
			/// </summary>
			double getOptimalFramesPerUpdate() const noexcept;
			/// <summary>
			/// Returns an estimate of how many frames whom are ready to be rendered, through processNextSpectrumFrame()
			/// </summary>
			std::size_t getApproximateStoredFrames() const noexcept;
			/// <summary>
			/// Inits the UI.
			/// </summary>
			void initPanelAndControls();
			/// <summary>
			/// Handles audio re-allocations, tunings and such on the correct thread (you)
			/// Call this before any graphic rendering. If openGL is enabled, this MUST be called on the
			/// openGL thread.
			/// Every rendering-state change is handled in here to minimize duplicate heavy changes/resizes
			/// (certain operations imply others)
			/// </summary>
			void handleFlagUpdates(StreamState& sac);

			template<typename ISA>
				void renderColourSpectrum(const Constant& constant, TransformPair& transform, cpl::OpenGLRendering::COpenGLStack &);

			template<typename ISA>
				void renderTransformAsGraph(cpl::OpenGLRendering::COpenGLStack &, const TransformPair& transform, const LineColours& one, const LineColours& two);

			template<typename ISA>
				void renderLineGrid(cpl::OpenGLRendering::COpenGLStack&);

			void renderSpectrogramGrid(cpl::OpenGLRendering::COpenGLStack& ogs);

			void drawFrequencyTracking(juce::Graphics & g, const float fps, const Constant& constant, /*const*/ TransformPair& transform);

			/// <summary>
			/// Calculates the apparant worst-case scalloping loss given the current transform, size, view and window as a fraction.
			/// </summary>
			double getScallopingLossAtCoordinate(std::size_t coordinate, const Constant& constant);

			/// <summary>
			/// Some calculations rely on the view not changing so everything doesn't have to be recalculated constantly.
			/// This resets them, requesting a recalculation.
			/// </summary>
			void resetStaticViewAssumptions();

			void recalculateLegend(StreamState& cs);

			// TODO: technically, avoid any data races by making every member a std::atomic (or refactor everything that
			// changes state into altering flags instead (check out valueChanged()))
			struct StateOptions
			{
				bool isEditorOpen, isFrozen, isSuspended;
				bool antialias;
				bool isLinear;

				/// <summary>
				/// Interpolation method for discrete bins to continuous space
				/// </summary>
				SpectrumContent::BinInterpolation binPolation;

				/// <summary>
				/// Is the spectrum a horizontal device (line graph)
				/// or a vertical coloured spectrum?
				/// </summary>
				SpectrumContent::DisplayMode displayMode;

				/// <summary>
				/// The current selected algorithm that will digest the audio data into the display.
				/// </summary>
				SpectrumContent::TransformAlgorithm algo;
				/// <summary>
				/// How the incoming data is interpreted, channel-wise.
				/// </summary>
				SpectrumChannels configuration;

				SpectrumContent::ViewScaling viewScale;

				float primitiveSize;

				/// <summary>
				/// Describes the lower and higher limit of the dynamic range of the display.
				/// </summary>
				cpl::relaxed_atomic<double> dynRange[2];

				/// <summary>
				/// colourOne & two = colours for the main line graphs.
				/// graphColour = colour for the frequency & db grid.
				/// </summary>
				juce::Colour colourGrid, colourBackground, colourWidget;

				std::array<ColourRotation, SpectrumContent::LineGraphs::LineEnd> colourOne, colourTwo;

				/// <summary>
				/// The window size..
				/// </summary>
				std::size_t windowSize;

				/// <summary>
				/// The current view of the spectrum (might be zoomed, etc.).
				/// Use only on rendering thread
				/// </summary>
				cpl::Utility::Bounds<double> viewRect;

				/// <summary>
				/// Logarithmic displays has to start at a positive value.
				/// </summary>
				double minLogFreq = 10;

				cpl::weak_atomic<std::size_t> newWindowSize;
				cpl::weak_atomic<float> sampleRate;


				float alphaFloodFill;

				std::size_t axisPoints, numFilters, transformSize;

				SpectrumContent::LineGraphs frequencyTrackingGraph;

				ChangeVersion::Listener audioStreamChanged;
				double windowScale;

				bool drawLegend{ };
				LegendCache legend;

			} state;


			/// <summary>
			/// Set these flags and their status will be handled in the next handleFlagUpdates() call, which
			/// shall be called before any graphic rendering.
			/// </summary>
			struct Flags
			{
				cpl::ABoolFlag
					/// <summary>
					/// Set this to resize the audio windows (like, when the audio window size (as in fft size) is changed.
					/// The argument to the resizing is the state.newWindowSize
					/// </summary>
					initiateWindowResize,
					/// <summary>
					/// Set if anything in the audio stream changed
					/// </summary>
					/// <remarks>
					/// Can be non-atomic?
					/// </remarks>
					audioStreamChanged,
					/// <summary>
					/// Set this if the audio buffer window size was changed from somewhere else.
					/// </summary>
					audioWindowWasResized,
					/// <summary>
					///
					/// </summary>
					audioMemoryResize,
					workingMemoryResize,
					/// <summary>
					/// Set this when resizing the window
					/// </summary>
					resized,
					/// <summary>
					/// Set this when the view is changed, zoomed, whatever. Implicitly set by resized
					/// </summary>
					viewChanged,
					/// <summary>
					/// Set this to affect purely visual changes to the graph (like divisions and such)
					/// </summary>
					frequencyGraphChange,
					/// <summary>
					/// Recomputes the window
					/// </summary>
					windowKernelChange,
					/// <summary>
					/// see @resetState() override
					/// </summary>
					resetStateBuffers,
					/// <summary>
					/// Set when openGL context is (re)created
					/// </summary>
					openGLInitiation,
					/// <summary>
					/// Set when the dynamic range is changed.
					/// </summary>
					dynamicRangeChange,
					/// <summary>
					/// Set to update how the transforms are displayed (spectrum, graphs, etc?)
					/// </summary>
					displayModeChange,
					/// <summary>
					/// Set this to recalculate the slopes
					/// </summary>
					slopeMapChanged,
					mouseMove;
			} flags;

			struct StreamState
			{
				std::vector<TransformPair> pairs;
				Constant constant;
				ChangeVersion audioStreamChangeVersion;
				double streamLocalSampleRate;				
				std::vector<std::string> channelNames;
			};

			struct ProcessorShell : public AudioStream::Listener
			{
				std::shared_ptr<const SharedBehaviour> globalBehaviour;
				cpl::CLockFreeDataQueue<FrameVector> frameQueue;
				CriticalSection<StreamState> streamState;
				cpl::relaxed_atomic<bool> isSuspended;

				void onStreamAudio(AudioStream::ListenerContext& source, AudioStream::DataType** buffer, std::size_t numChannels, std::size_t numSamples) override;
				void onStreamPropertiesChanged(AudioStream::ListenerContext& source, const AudioStream::AudioStreamInfo& before) override;

				ProcessorShell(std::shared_ptr<const SharedBehaviour>& behaviour);
			};

			friend struct AudioDispatcher;

			std::shared_ptr<const SharedBehaviour> globalBehaviour;			
			std::shared_ptr<const ConcurrentConfig> config;
			cpl::OpenGLRendering::COpenGLImage oglImage;
			cpl::special::FrequencyAxis frequencyGraph, complexFrequencyGraph;
			cpl::special::DBMeterAxis dbGraph;

			// non-state variables
			std::shared_ptr<SpectrumContent> content;
			std::shared_ptr<ProcessorShell> processor;
			cpl::Utility::Bounds<double> oldViewRect;
			double scallopLoss;
			int lastPeak;

			/// <summary>
			/// Lenght/height after state updates.
			/// </summary>
			// TODO: Delete?
			int framePixelPosition;
			double oldWindowSize;
			double framesPerUpdate;
			std::vector<cpl::GraphicsND::UPixel<cpl::GraphicsND::ComponentOrder::OpenGL>> columnUpdate;

			/// <summary>
			/// The connected, incoming stream of data.
			/// </summary>
			std::shared_ptr<AudioStream::Output> audioStream;

			struct SmoothedPeakState
			{
				typedef cpl::dsp::SmoothedParameterState<double, 8> SmoothParam; 


				void clearNonNormals()
				{
					if (!std::isnormal(currentPeak))
						currentPeak = 0;
					if (!std::isnormal(frequencyState.getState()))
						frequencyState.reset();
					if (!std::isnormal(peakDBsState.getState()))
						peakDBsState.reset();
				}

				void process(double peakDBs, double frequency)
				{
					auto peak = std::pow(10, peakDBs / 20);

					currentPeak *= peakPole;

					if (peak > currentPeak)
					{
						currentPeak = 1.2 * (peak / peakPole);
						currentFrequency = frequency;
						currentPeakDBs = peakDBs;
					}

					
					frequencyState.process(filterPole, currentFrequency);
					peakDBsState.process(filterPole, currentPeakDBs);
				}

				void design(double ms, double sampleRate)
				{
					peakPole = SmoothParam::design(ms * 10, sampleRate);
					filterPole = SmoothParam::design(ms / 5, sampleRate);
				}

				double getPeakDBs()
				{
					return peakDBsState.getState();
				}

				double getFrequency()
				{
					return frequencyState.getState();
				}

			private:
				SmoothParam::PoleState peakPole {}, filterPole{};
				SmoothParam frequencyState, peakDBsState;
				double currentFrequency{}, currentPeakDBs{}, currentPeak{};

			} peakState;
		};

	};

#endif
