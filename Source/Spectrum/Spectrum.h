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

	#include "Signalizer.h"
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
		:
			public cpl::COpenGLView,
			protected AudioStream::Listener,
			private ParameterSet::RTListener
		{

		public:

			typedef UComplexFilter<AudioStream::DataType> UComplex;
			typedef AudioStream::DataType fpoint;
			typedef double fftType;

			class SFrameBuffer
			{
			public:

				typedef cpl::aligned_vector < UComplex, 32 > FrameVector;

				SFrameBuffer()
					: sampleBufferSize(), sampleCounter(), currentCounter(), frameQueue(10, 1000)
				{

				}
				std::size_t sampleBufferSize;
				std::size_t currentCounter;
				std::uint64_t sampleCounter;

				cpl::CLockFreeQueue<FrameVector *> frameQueue;
			};



			struct DBRange
			{
				double low, high;
			};

			Spectrum(const SharedBehaviour & globalBehaviour, const std::string & nameId, AudioStream & data, ProcessorState * state);
			virtual ~Spectrum();

			// Component overrides
			void mouseWheelMove(const juce::MouseEvent& event, const juce::MouseWheelDetails& wheel) override;
			void mouseDoubleClick(const juce::MouseEvent& event) override;
			void mouseDrag(const juce::MouseEvent& event) override;
			void mouseUp(const juce::MouseEvent& event) override;
			void mouseDown(const juce::MouseEvent& event) override;
			void mouseMove(const juce::MouseEvent& event) override;
			void mouseExit(const juce::MouseEvent & e) override;
			void mouseEnter(const juce::MouseEvent & e) override;

			void resized() override;


			// OpenGLRender overrides
			void onOpenGLRendering() override;
			void onGraphicsRendering(juce::Graphics & g) override;
			void initOpenGL() override;
			void closeOpenGL() override;

			// View overrides
			void suspend() override;
			void resume() override;
			void freeze() override;
			void unfreeze() override;
			void resetState() override;
			std::unique_ptr<juce::Component> createEditor();

			bool isEditorOpen() const;

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

			struct RenderingDispatcher
			{
				template<typename ISA> static void dispatch(Spectrum & c) { c.vectorGLRendering<ISA>(); }
			};

			struct AudioDispatcher
			{
				template<typename ISA> static void dispatch(Spectrum & c, AFloat ** buffer, std::size_t numChannels, std::size_t numSamples)
				{
					c.audioProcessing<ISA>(buffer, numChannels, numSamples);
				}
			};


            template<typename ISA>
                void vectorGLRendering();

			virtual void paint2DGraphics(juce::Graphics & g);

			virtual void parameterChangedRT(cpl::Parameters::Handle localHandle, cpl::Parameters::Handle globalHandle, ParameterSet::BaseParameter * param) override;

			bool onAsyncAudio(const AudioStream & source, AudioStream::DataType ** buffer, std::size_t numChannels, std::size_t numSamples) override;
			void onAsyncChangedProperties(const AudioStream & source, const AudioStream::AudioStreamInfo & before) override;


			/// <summary>
			/// Maps the current resonating system according to the current model (linear/logarithmic) and the current
			/// subsection of the complete spectrum such that a linear array of output data matches pixels 1:1, as well as
			/// formats the data into the filterResults array according to the channel mode (SpectrumChannels).
			///
			/// Call prepareTransform(), then doTransform(), then mapToLinearSpace().
			/// After the call to mapToLinearSpace, the results are written to getWorkingMemory().
			/// Returns the complex amount of filters processed.
			/// Needs exclusive access to audioResource.
			/// </summary>
			std::size_t mapToLinearSpace();

			/// <summary>
			/// Runs the transform (of any kind) results through potential post filters and other features, before displaying it.
			/// The transform will be rendered into filterResults after this.
			/// </summary>
			template<class InVector>
				void postProcessTransform(const InVector & transform, std::size_t size);

			/// <summary>
			/// Post processes the transform that will be interpreted according to what's selected.
			/// Needs exclusive access to audioResource.
			/// </summary>
			void postProcessStdTransform();

			/// <summary>
			/// Call this when something affects the view scaling, view size, mapping of frequencies, display modes etc.
			/// Set state accordingly first.
			/// </summary>
			void displayReordered();

			/// <summary>
			/// For some transform algorithms, it may be a no-op, but for others (like FFTs) that may need zero-padding
			/// or windowing, this is done here.
			///
			/// Call prepareTransform(), then doTransform(), then mapToLinearSpace()
			/// Needs exclusive access to audioResource.
			/// </summary>
			bool prepareTransform(const AudioStream::AudioBufferAccess & audio);

			/// <summary>
			/// For some transform algorithms, it may be a no-op, but for others (like FFTs) that may need zero-padding
			/// or windowing, this is done here.
			///
			/// This functions considers the additional arguments as more recent audio than the audio buffers (and as of such, considers numSamples less audio from
			/// the first argument).
			/// Needs exclusive access to audioResource.
			/// </summary>
			bool prepareTransform(const AudioStream::AudioBufferAccess & audio, fpoint ** preliminaryAudio, std::size_t numChannels, std::size_t numSamples);

			/// <summary>
			/// Again, some algorithms may not need this, but this ensures the transform is done after this call.
			///
			/// Call prepareTransform(), then doTransform(), then mapToLinearSpace()
			/// Needs exclusive access to audioResource.
			/// </summary>
			void doTransform();

			/// <summary>
			/// internally used for now.
			/// </summary>
			bool processNextSpectrumFrame();

			void calculateSpectrumColourRatios();
		private:

			void deserialize(cpl::CSerializer::Builder & builder, cpl::Version version) override;
			void serialize(cpl::CSerializer::Archiver & archive, cpl::Version version) override;

			std::size_t getValidWindowSize(std::size_t in) const noexcept;

			/// <summary>
			/// Returns the number of needed channels required to process the current
			/// channel configuration.
			/// </summary>
			std::size_t getStateConfigurationChannels() const noexcept;

			/// <summary>
			/// Transforms state.audioBlobSizeMs into samples.
			/// </summary>
			/// <returns></returns>
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
			/// <returns></returns>
			int getAxisPoints() const noexcept;
			/// <summary>
			/// For a certain blob size and refresh rate, N frames has to be generated each display update.
			/// </summary>
			/// <returns></returns>
			double getOptimalFramesPerUpdate() const noexcept;
			/// <summary>
			/// Returns an estimate of how many frames whom are ready to be rendered, through processNextSpectrumFrame()
			/// </summary>
			/// <returns></returns>
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
			void handleFlagUpdates();
			/// <summary>
			/// Computes the time-domain window kernel (see state.windowKernel) of getWindowSize() size.
			/// </summary>
			void computeWindowKernel();

			template<typename T>
				T * getAudioMemory()
				{
					return reinterpret_cast<T*>(audioMemory.data());
				}

			/// <summary>
			/// All inputs must be normalized. Scales the input to the display decibels, and runs it through peak filters.
			/// newVals = current vector of floats / doubles * 2 (complex), output from CSignalTransform::**dft() of size * 2
			/// for mode = left / merge / mid / side / right
			/// 	newVals is a complex vector of floats of size
			/// for mode = separate, mid&side
			/// 	newVals is a complex vector of floats of size * 2
			/// 	newVals[n * 2 + 0] = lreal
			/// 	newVals[n * 2 + 1] = limag
			/// 	newVals[n * 2 + size + 0] = rreal
			/// 	newVals[n * 2 + size + 1] = rimag
			/// for mode = phase
			/// 	newVals[n * 2 + 0] = mag
			/// 	newVals[n * 2 + 1] = phase cancellation(with 1 being totally cancelled)
			/// </summary>
			template<class V2>
				void mapAndTransformDFTFilters(SpectrumChannels type, const V2 & newVals, std::size_t size, double lowerFraction, double upperFraction, float clip);

			/// <summary>
			/// Returns the number of T elements available in the audio space buffer
			/// (as returned by getAudioMemory<T>())
			/// </summary>
			template<typename T>
				std::size_t getNumAudioElements() const noexcept;

			/// <summary>
			/// Gets the relay buffer for the channel.
			/// Note: you should call ensureRelayBufferSize before.
			/// </summary>
			fpoint * getRelayBufferChannel(std::size_t channel);
			/// <summary>
			/// Ensures the storage for the relay buffer.
			/// Not suited for real-time usage (allocates memory)
			/// </summary>
			void ensureRelayBufferSize(std::size_t channels, std::size_t numSamples);
			/// <summary>
			/// Returns the number of T elements available to be used as a FFT (power2 buffer size).
			/// Guaranteed to be less than getNumAudioElements.
			/// TODO: hoist into .inl file (other cases in this file as well)
			/// </summary>
			template<typename T>
				std::size_t getFFTSpace() const noexcept
				{
					auto size = audioMemory.size();
					return size ? (size - 1) / sizeof(T) : 0;
				}

			template<typename T>
				T * getWorkingMemory();

			template<typename T>
				std::size_t getNumWorkingElements() const noexcept;

			template<typename ISA>
				void renderColourSpectrum(cpl::OpenGLRendering::COpenGLStack &);

			template<typename ISA>
				void renderLineGraph(cpl::OpenGLRendering::COpenGLStack &);

			template<typename ISA>
				void audioProcessing(float ** buffer, std::size_t numChannels, std::size_t numSamples);

			template<typename ISA>
				void resonatingDispatch(float ** buffer, std::size_t numChannels, std::size_t numSamples);

			template<typename ISA>
				void addAudioFrame();

			/// <summary>
			/// Copies the state from the complex resonator into the output buffer.
			/// The output vector is assumed to accept index assigning of std::complex of fpoints.
			/// It is assumed the output vector can hold at least numChannels (of configuration) times
			/// numFilters.
			/// Output of channels are stored at numFilters offsets.
			///
			/// Returns the total number of complex samples copied into the output
			/// </summary>
			template<typename ISA, class Vector>
				std::size_t copyResonatorStateInto(Vector & output);

			void drawFrequencyTracking(juce::Graphics & g);

			/// <summary>
			/// Calculates the apparant worst-case scalloping loss given the current transform, size, view and window as a fraction.
			/// </summary>
			/// <param name="coordinate"></param>
			/// <returns></returns>
			double getScallopingLossAtCoordinate(std::size_t coordinate);

			/// <summary>
			/// Some calculations rely on the view not changing so everything doesn't have to be recalculated constantly.
			/// This resets them, requesting a recalculation.
			/// </summary>
			void resetStaticViewAssumptions();

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
				std::atomic<SpectrumContent::TransformAlgorithm> algo;
				/// <summary>
				/// How the incoming data is interpreted, channel-wise.
				/// </summary>
				SpectrumChannels configuration;

				SpectrumContent::ViewScaling viewScale;

				float primitiveSize;

				/// <summary>
				/// Describes the lower and higher limit of the dynamic range of the display.
				/// </summary>
				std::atomic<double> dynRange[2];

				/// <summary>
				/// colourOne & two = colours for the main line graphs.
				/// graphColour = colour for the frequency & db grid.
				/// </summary>
				juce::Colour colourGrid, colourBackground, colourTracker, colourOne[SpectrumContent::LineGraphs::LineEnd], colourTwo[SpectrumContent::LineGraphs::LineEnd];

				/// <summary>
				/// Colours for spectrum
				/// </summary>
				cpl::GraphicsND::UPixel<cpl::GraphicsND::ComponentOrder::OpenGL> colourSpecs[SpectrumContent::numSpectrumColours + 1];
				float normalizedSpecRatios[SpectrumContent::numSpectrumColours + 1];

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

				/// <summary>
				/// The window function applied to the input. Precomputed into windowKernel.
				/// </summary>
				std::atomic<cpl::dsp::WindowTypes> dspWindow;

				std::atomic<std::size_t> newWindowSize;
				std::atomic<float> sampleRate;
				/// <summary>
				/// internal testing flag
				/// </summary>
				bool iAuxMode;

				float alphaFloodFill;

				std::size_t axisPoints, numFilters;

				SpectrumContent::LineGraphs frequencyTrackingGraph;
			} state;


			/// <summary>
			/// Set these flags and their status will be handled in the next handleFlagUpdates() call, which
			/// shall be called before any graphic rendering.
			/// </summary>
			struct Flags
			{
				/// <summary>
				/// Dont set this! Set by the flag handler to assert state
				/// </summary>
				volatile bool internalFlagHandlerRunning;

				cpl::ABoolFlag
					/// <summary>
					/// Set this to resize the audio windows (like, when the audio window size (as in fft size) is changed.
					/// The argument to the resizing is the state.newWindowSize
					/// </summary>
					initiateWindowResize,
					/// <summary>
					/// Set if anything in the audio stream changed
					/// </summary>
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
					/// This flag will be true the first time handleFlagUpdates() is called
					/// </summary>
					firstChange,
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

			/// <summary>
			/// visual objects
			/// </summary>
			const SharedBehaviour & globalBehaviour;
			juce::MouseCursor displayCursor;
			cpl::OpenGLRendering::COpenGLImage oglImage;
			cpl::special::FrequencyAxis frequencyGraph, complexFrequencyGraph;
			cpl::special::DBMeterAxis dbGraph;
			cpl::CBoxFilter<double, 60> avgFps;

			// non-state variables
			SpectrumContent * content;
			unsigned long long processorSpeed; // clocks / sec
			double audioThreadUsage;
			double laggedFPS;
			juce::Point<float> lastMousePos;
			std::vector<std::unique_ptr<juce::OpenGLTexture>> textures;
			long long lastFrameTick, renderCycles;
			bool wasResized;
			cpl::Utility::Bounds<double> oldViewRect;
			std::atomic_bool hasMainThreadInitializedAudioStreamDependenant;
			std::atomic_bool isMouseInside;
			double scallopLoss;
			int lastPeak;
			/*struct NewChanges
			{
				std::atomic<DisplayMode> displayMode;
				std::atomic<std::size_t> windowSize;
				std::atomic<cpl::iCtrlPrec_t> divLimit;
				std::atomic<SpectrumChannels> configuration;
				std::atomic<cpl::iCtrlPrec_t> stretch;
				std::atomic<signed int> frequencyTrackingGraph;
				std::atomic<float> primitiveSize;
				std::atomic<float> alphaFloodFill;
				std::atomic<double> referenceTuning;
			} newc; */

			struct CurrentMouse
			{
				std::atomic<float>
					x, y;
			} cmouse;

			/// <summary>
			/// see cpl::dsp::windowScale
			/// </summary>
			fftType windowScale;

			/// <summary>
			/// Lenght/height after state updates.
			/// </summary>
			std::size_t relayWidth, relayHeight;
			int framePixelPosition;
			double oldWindowSize;
			int droppedAudioFrames;
			double framesPerUpdate;
			cpl::CPeakFilter<double> fpuFilter;
			std::vector<cpl::GraphicsND::UPixel<cpl::GraphicsND::ComponentOrder::OpenGL>> columnUpdate;

			struct LineGraphDesc
			{
				/// <summary>
				/// The peak filter coefficient, describing the decay rate of the filters.
				/// </summary>
				cpl::CPeakFilter<fpoint> filter;
				/// <summary>
				/// The'raw' formatted state output of the mapped transform algorithms.
				/// </summary>
				cpl::aligned_vector<UComplex, 32> states;
				/// <summary>
				/// The decay/peak-filtered and scaled outputs of the transforms,
				/// with each element corrosponding to a complex output pixel of getAxisPoints() size.
				/// Resized in displayReordered
				/// </summary>
				cpl::aligned_vector<UComplex, 32> results;

				void resize(std::size_t n)
				{
					states.resize(n); results.resize(n);
				}

				void zero() {
					std::memset(states.data(), 0, states.size() * sizeof(UComplex));
					std::memset(results.data(), 0, results.size() * sizeof(UComplex));
				}
			};
			// dsp objects
			std::array<LineGraphDesc, SpectrumContent::LineGraphs::LineEnd> lineGraphs;
			/// <summary>
			/// The complex resonator used for iir spectrums
			/// </summary>
			cpl::dsp::CComplexResonator<fpoint, 2> cresonator;
			/// <summary>
			/// An array, of numFilters size, with each element being the frequency for the filter of
			/// the corresponding logical display pixel unit.
			/// </summary>
			std::vector<fpoint> mappedFrequencies;
			/// <summary>
			/// The connected, incoming stream of data.
			/// </summary>
			AudioStream & audioStream;
			/// <summary>
			/// Temporary memory buffer for audio applications. Resized in setWindowSize (since the size is a function of the window size)
			/// </summary>
			cpl::aligned_vector<char, 32> audioMemory;
			/// <summary>
			/// used for 'real-time' audio processing, when the incoming buffer needs to be rearranged without changing it.
			/// </summary>
			struct RelayBuffer
			{
				cpl::aligned_vector<fpoint, 32> buffer;
				std::size_t samples, channels;
			} relay;

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
				SmoothParam::PoleState peakPole, filterPole;
				SmoothParam frequencyState, peakDBsState;
				double currentFrequency{}, currentPeakDBs{}, currentPeak{};

			} peakState;

			/// <summary>
			/// Temporary memory buffer for other applications.
			/// Resized in displayReordered
			/// </summary>
			cpl::aligned_vector<char, 32> workingMemory;
			/// <summary>
			/// The time-domain representation of the dsp-window applied to fourier transforms.
			/// </summary>
			cpl::aligned_vector<double, 32> windowKernel;

			cpl::aligned_vector<fpoint, 32> slopeMap;
			/// <summary>
			/// All audio processing not done in the audio thread (not real-time, async audio) must acquire this lock.
			/// Notice, you must always acquire this lock before accessing the audio buffers (should you intend to).
			/// </summary>
			cpl::CMutex::Lockable audioResource;

			SFrameBuffer sfbuf;
		};

	};

#endif
