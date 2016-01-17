
#ifndef _CSPECTRUM_H
	#define _CSPECTRUM_H

	#include <cpl/Common.h>
	#include <cpl/CAudioBuffer.h>
	#include <cpl/CViews.h>
	#include <cpl/GraphicComponents.h>
	#include <cpl/Utility.h>
	#include <cpl/gui/Controls.h>
	#include <cpl/dsp/CComplexResonator.h>
	#include <memory>
	#include <cpl/rendering/COpenGLImage.h>
	#include <cpl/CFrequencyGraph.h>
	#include <cpl/dsp/CPeakFilter.h>
	#include <cpl/CDBMeterGraph.h>
	#include <queue>
	#include <cpl/lib/LockFreeQueue.h>
	#include <cpl/lib/BlockingLockFreeQueue.h>
	#include <vector>
	#include <cpl/gui/CPresetWidget.h>
	#include "CommonSignalizer.h"
	#include <cpl/gui/CDSPWindowWidget.h>

	namespace cpl
	{
		namespace OpenGLEngine
		{
			class COpenGLStack;
			
		};
	};


	namespace Signalizer
	{



		
		class CSpectrum
		: 
			public cpl::COpenGLView, 
			protected cpl::CBaseControl::PassiveListener,
			protected cpl::CBaseControl::ValueFormatter,
			protected cpl::CBaseControl::Listener,
			protected AudioStream::Listener,
			protected juce::ComponentListener
		{

		public:

			static const std::size_t numSpectrumColours = 5;
			static const double minDBRange;
			static const double kMinDbs;
			static const double kMaxDbs;
			
			typedef UComplexFilter<AudioStream::DataType> UComplex;
			typedef AudioStream::DataType fpoint;
			typedef double fftType;

			class SFrameBuffer : public cpl::CMutex::Lockable
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

			enum class BinInterpolation
			{
				None,
				Linear,
				Lanczos
			};

			enum class DisplayMode
			{
				LineGraph,
				ColourSpectrum
			};

			enum class TransformAlgorithm
			{
				FFT, RSNT
			};

			enum class ViewScaling
			{
				Linear,
				Logarithmic
			};

			

			struct DBRange
			{
				double low, high;
			};

			CSpectrum(AudioStream & data);
			virtual ~CSpectrum();

			// Component overrides
			void mouseWheelMove(const MouseEvent& event, const MouseWheelDetails& wheel) override;
			void mouseDoubleClick(const MouseEvent& event) override;
			void mouseDrag(const MouseEvent& event) override;
			void mouseUp(const MouseEvent& event) override;
			void mouseDown(const MouseEvent& event) override;
			void resized() override;


			// OpenGLRender overrides
			void onOpenGLRendering() override;
			void onGraphicsRendering(Graphics & g) override;
			void initOpenGL() override;
			void closeOpenGL() override;

			// View overrides
			void suspend() override;
			void resume() override;
			void freeze() override;
			void unfreeze() override;
			void resetState() override;
			std::unique_ptr<juce::Component> createEditor() override;
			// CSerializer overrides
			void load(cpl::CSerializer::Builder & builder, long long int version) override;
			void save(cpl::CSerializer::Archiver & archive, long long int version) override;

			// cbasecontrol overrides
			void valueChanged(const cpl::CBaseControl *) override;
			bool valueChanged(cpl::CBaseControl *) override;
			bool stringToValue(const cpl::CBaseControl * ctrl, const std::string & buffer, cpl::iCtrlPrec_t & value) override;
			bool valueToString(const cpl::CBaseControl * ctrl, std::string & buffer, cpl::iCtrlPrec_t value) override;
			void onObjectDestruction(const cpl::CBaseControl::ObjectProxy & destroyedObject) override;
			bool isEditorOpen() const;

			// interface


			void setDBs(double low, double high, bool updateControls = false);
			void setDBs(DBRange &, bool updateControls = false);
			DBRange getDBs() const noexcept;

			/// <summary>
			/// Returns the working size of the audio buffer.
			/// It is guaranteed to be 0 > getWindowSize() <= audioStream.getAudioHistorySize()
			/// </summary>
			/// <returns></returns>
			std::size_t getWindowSize() const noexcept;
			void setWindowSize(std::size_t size);

		protected:

			void setTransformOptions();

			virtual void paint2DGraphics(juce::Graphics & g);

			bool onAsyncAudio(const AudioStream & source, float ** buffer, std::size_t numChannels, std::size_t numSamples) override;
			void onAsyncChangedProperties(const AudioStream & source, const AudioStream::AudioStreamInfo & before) override;
			void componentBeingDeleted(Component & 	component) override;

			/// <summary>
			/// Tunes the system accordingly to the view and current zoom.
			/// </summary>
			void mapFrequencies();
			/// <summary>
			/// Maps the current resonating system according to the current model (linear/logarithmic) and the current
			/// subsection of the complete spectrum such that a linear array of output data matches pixels 1:1, as well as 
			/// formats the data into the filterResults array according to the channel mode (ChannelConfiguration).
			/// 
			/// Call prepareTransform(), then doTransform(), then mapToLinearSpace(). 
			/// After the call to mapToLinearSpace, the results are written to getWorkingMemory().
			/// Returns the complex amount of filters processed.
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
			/// </summary>
			void prepareTransform(const AudioStream::AudioBufferAccess & audio);

			/// <summary>
			/// For some transform algorithms, it may be a no-op, but for others (like FFTs) that may need zero-padding
			/// or windowing, this is done here.
			/// 
			/// This functions considers the additional arguments as more recent audio than the audio buffers (and as of such, considers numSamples less audio from
			/// the first argument).
			/// </summary>
			void prepareTransform(const AudioStream::AudioBufferAccess & audio, fpoint ** preliminaryAudio, std::size_t numChannels, std::size_t numSamples);

			/// <summary>
			/// Again, some algorithms may not need this, but this ensures the transform is done after this call.
			/// 
			/// Call prepareTransform(), then doTransform(), then mapToLinearSpace()
			/// </summary>
			void doTransform();

			/// <summary>
			/// internally used for now.
			/// </summary>
			bool processNextSpectrumFrame();

			void calculateSpectrumColourRatios();
		private:

			void audioConsumerThread();

			std::size_t getValidWindowSize(std::size_t in) const noexcept;

			/// <summary>
			/// Returns the number of needed channels required to process the current
			/// channel configuration.
			/// </summary>
			std::size_t CSpectrum::getStateConfigurationChannels() const noexcept;

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
				T * getAudioMemory();

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
			/// Guaranteed to be < getNumAudioElements.
			/// </summary>
			template<typename T>
				std::size_t getFFTSpace() const noexcept;

			template<typename T>
				T * getWorkingMemory();

			template<typename T>
				std::size_t getNumWorkingElements() const noexcept;

			template<typename V>
				void renderColourSpectrum(cpl::OpenGLEngine::COpenGLStack &);

			template<typename V>
				void renderLineGraph(cpl::OpenGLEngine::COpenGLStack &);

			template<typename V>
				void audioProcessing(float ** buffer, std::size_t numChannels, std::size_t numSamples);

			template<typename V>
				void resonatingDispatch(float ** buffer, std::size_t numChannels, std::size_t numSamples);

			template<typename V>
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
			template<typename V, class Vector>
				std::size_t copyResonatorStateInto(Vector & output);

			// vars

			struct StateOptions
			{
				bool isEditorOpen, isFrozen, isSuspended;
				bool antialias;
				bool isLinear;
				
				/// <summary>
				/// Interpolation method for discrete bins to continuous space
				/// </summary>
				BinInterpolation binPolation;

				/// <summary>
				/// Is the spectrum a horizontal device (line graph)
				/// or a vertical coloured spectrum?
				/// </summary>
				DisplayMode displayMode;

				/// <summary>
				/// The current selected algorithm that will digest the audio data into the display.
				/// </summary>
				TransformAlgorithm algo;
				/// <summary>
				/// How the incoming data is interpreted, channel-wise.
				/// </summary>
				ChannelConfiguration configuration;

				ViewScaling viewScale;

				float primitiveSize;

				/// <summary>
				/// Describes the lower and higher limit of the dynamic range of the display.
				/// </summary>
				DBRange dynRange;

				/// <summary>
				/// colourOne & two = colours for the main line graphs.
				/// graphColour = colour for the frequency & db grid.
				/// </summary>
				juce::Colour colourOne, colourTwo, colourGrid, colourBackground;

				/// <summary>
				/// Colours for spectrum
				/// </summary>
				cpl::GraphicsND::UPixel<cpl::GraphicsND::ComponentOrder::OpenGL> colourSpecs[numSpectrumColours + 1];
				float normalizedSpecRatios[numSpectrumColours + 1];

				/// <summary>
				/// The window size..
				/// </summary>
				std::size_t windowSize;

				/// <summary>
				/// The current view of the spectrum (might be zoomed, etc.)
				/// </summary>
				cpl::Utility::Bounds<double> viewRect;

				/// <summary>
				/// Logarithmic displays has to start at a positive value.
				/// </summary>
				double minLogFreq;

				/// <summary>
				/// The window function applied to the input. Precomputed into windowKernel.
				/// </summary>
				cpl::dsp::WindowTypes dspWindow;

				/// <summary>
				/// internal testing flag
				/// </summary>
				bool iAuxMode;

				/// <summary>
				/// How often the audio stream is sampled for updating in the screen. Only for DisplayMode::ColourSpectrum modes.
				/// Effectively, this sets each horizontal pixel to contain this amount of miliseconds of information.
				/// </summary>
				double audioBlobSizeMs;

				/// <summary>
				/// For displaymode == ColourSpectrum, this value indicates how much devations in amount of frame updates are smoothed such that
				/// the display doesn't jitter too much.
				/// </summary>
				double bufferSmoothing;

				std::size_t axisPoints, numFilters;
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
					/// Set this to resize the audio windows (like, when the audio window size (as in fft size) is changed
					/// </summary>
					initiateWindowResize,
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
					firstChange;
			} flags;

			/// <summary>
			/// GUI elements
			/// </summary>
			juce::Component * editor;
			cpl::CComboBox kviewScaling, kalgorithm, kchannelConfiguration, kdisplayMode, kdspWindow, kbinInterpolation;
			cpl::CDSPWindowWidget kdspWin;
			cpl::CKnobSlider klowDbs, khighDbs, kdecayRate, kwindowSize, kpctForDivision, kblobSize, kframeUpdateSmoothing;
			cpl::CColourControl kline1Colour, kline2Colour, kgridColour, kbackgroundColour;
			cpl::CPresetWidget presetManager;
			cpl::CButton kdiagnostics, kfreeQ;
			/// <summary>
			/// Colour interpolators
			/// </summary>
			cpl::CColourControl kspecColours[5];
			cpl::CKnobSlider kspecRatios[5];
			/// <summary>
			/// visual objects
			/// </summary>

			juce::MouseCursor displayCursor;
			cpl::OpenGLEngine::COpenGLImage oglImage;
			cpl::CFrequencyGraph frequencyGraph;
			cpl::CDBMeterGraph dbGraph;
			cpl::CBoxFilter<double, 60> avgFps;
			

			// non-state variables
			unsigned long long processorSpeed; // clocks / sec
			double audioThreadUsage;
			juce::Point<float> lastMousePos;
			std::vector<std::unique_ptr<juce::OpenGLTexture>> textures;
			long long lastFrameTick, renderCycles;
			bool wasResized, isSuspended;
			cpl::Utility::Bounds<double> oldViewRect;
			DisplayMode newDisplayMode;
			std::atomic<std::size_t> newWindowSize;
			/// <summary>
			/// see cpl::dsp::windowScale
			/// </summary>
			fftType windowScale;

			/// <summary>
			/// Lenght/height after state updates.
			/// </summary>
			std::size_t relayWidth, relayHeight;
			int framePixelPosition;
			int droppedAudioFrames;
			double framesPerUpdate;
			cpl::CPeakFilter<double> fpuFilter;
			std::vector<cpl::GraphicsND::UPixel<cpl::GraphicsND::ComponentOrder::OpenGL>> columnUpdate;


			// dsp objects
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
			/// The peak filter coefficient, describing the decay rate of the filters.
			/// </summary>
			cpl::CPeakFilter<float> peakFilter;
			/// <summary>
			/// The'raw' formatted state output of the mapped transform algorithms.
			/// Resized in displayReordered
			/// </summary>
			cpl::aligned_vector<UComplexFilter<fpoint>, 32> filterStates;
			/// <summary>
			/// The decay/peak-filtered and scaled outputs of the transforms,
			/// with each element corrosponding to a complex output pixel of getAxisPoints() size.
			/// Resized in displayReordered
			/// </summary>
			cpl::aligned_vector<UComplexFilter<fpoint>, 32> filterResults;
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

			/// <summary>
			/// Temporary memory buffer for other applications.
			/// Resized in displayReordered
			/// </summary>
			cpl::aligned_vector<char, 32> workingMemory;
			/// <summary>
			/// The time-domain representation of the dsp-window applied to fourier transforms.
			/// </summary>
			cpl::aligned_vector<double, 32> windowKernel;

			/// <summary>
			/// All audio processing not done in the audio thread must acquire this lock. It is free to acquire and is 100% userspace.
			/// The only time it may be locked is during initialization/resets, in which case audio drop-outs
			/// is to be expected.
			/// </summary>
			//cpl::CMutex::Lockable audioStateLock;

			SFrameBuffer sfbuf;
		};
	
	};

#endif