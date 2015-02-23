#ifndef CSpectrum_H
#define CSpectrum_H

#include <cpl/Common.h>
#include <cpl/CLowPass.h>
#include <cpl/CAudioBuffer.h>
#include <cpl/GraphicComponents.h>
#include <cpl/ComponentContainers.h>
#include <cpl/Graphics.h>
#include <vector>
#include <memory>
#include <complex>
#include <cpl/MacroConstants.h>
#include <cpl/lib/AlignedAllocator.h>
#include <cpl/UserInterface.h>
#include <cpl/ffts.h>
#include <cpl/dsp/CSignalTransform.h>
#include <cpl/dsp/CPeakFilter.h>
#include <cpl/dsp/SDFTSystem.h>
#include <cpl/dsp/CComplexResonator.h>

namespace Signalizer
{
	using namespace juce;
	class Engine;
	
	typedef double ResonatorType;


	template<typename Scalar>
		union FTFilterResult
		{
			FTFilterResult() : real(0), imag(0) { }
			struct
			{
				Scalar real, imag;
			};
			struct
			{
				Scalar magnitude, phase;
			};
			struct
			{
				Scalar leftMagnitude, rightMagnitude;
			};
		};

	template<typename T>
		struct ResonatorSystem : public cpl::CMutex::Lockable
		{
			T * poles;
			T * c1;
			T * c2;
			T * x;
			T * y;
			T * real;
			T * imag;
			std::vector<T, cpl::AlignmentAllocator<T, 32u>> buffer;
		};

	template<typename T, std::size_t Vectors>
		struct SDFTSystem2 : public cpl::CMutex::Lockable
		{
			static const std::size_t numVectors = Vectors;
			T * realPoles[Vectors];
			T * imagPoles[Vectors];
			T * realState[Vectors];
			T * imagState[Vectors];
			unsigned * N;
			std::vector<T, cpl::AlignmentAllocator<T, 32u>> buffer;
		};

	class CSpectrum 
	: 
		public cpl::SubView, 
		private cpl::CCtrlListener, 
		private cpl::CButtonGroup<cpl::CRenderButton>::CMultiButtonCallback,
		private cpl::COpenGLView::COpenGLRenderBack,
		private cpl::MouseDelegate::MouseCallBack,
		protected cpl::CAudioListener
	{

	public:


		struct OptionButtons
		{
			enum Enum {
				drawCurve, threaded, parallel,
				end
			};
		};

		struct fftWindows
		{
			enum Enum {
				Hanning,
				Hamming,
				Rectangular,
				Blackmann,
				Triangular,
				Gaussian,
				end,
			};
		};

		struct Quality
		{
			enum Enum {
				Free,
				Bound,
				end
			};
		};

		struct Algorithm
		{
			enum Enum {
				FFT, RSNT, MQDFT,
				end
			};
		};

		struct ViewTypes
		{
			enum Enum {
				SpectralView, SpectroGram,
				end
			};
		};
		struct ScaleButtons
		{
			enum Enum {
				X, Y,
				end
			};
		};
		struct ChannelConfiguration
		{
			enum Enum {
				Left, Right, Merge, Phase, Separate,
				end
			};
		};

		struct DBRange
		{
			double low, high;
		};


		CSpectrum(AudioBuffer & data);
		virtual ~CSpectrum();

		void rearrange(int width, int height, bool callback = true);
		void paint(Graphics & g) override;

		// subview
		bool serialize(juce::MemoryBlock & data) override;
		bool restore(const void * data, std::size_t size) override;
		bool setFullScreenMode(bool toggle) override;
		void repaintMainContent() override;
		void setDBs(double low, double high, bool updateControls = false);
		void setDBs(DBRange &, bool updateControls = false);
		void setWindowSize(std::size_t smps);
		DBRange getDBs();

	protected:
		// openglrenderback
		virtual void renderOpenGL() override;
		virtual void renderNormal(juce::Graphics & g) override;
		// subview
		virtual void attachToOpenGL(juce::OpenGLContext & ctx) override
		{
			oglc = &ctx;
			renderWindow.attachToOpenGL(&ctx);

		}
		virtual bool audioCallback(float ** data, std::size_t numChannels, std::size_t numSamples);

		template<typename T>
			void resonate(float ** data, std::size_t numChannels, std::size_t numSamples);

		template<typename T>
			void wresonate(float ** data, std::size_t numChannels, std::size_t numSamples);

		// mouse events
		void cbMouseWheelMove(const MouseEvent& event, const MouseWheelDetails& wheel) override;
		void cbMouseDoubleClick(const MouseEvent& event) override;
		void cbMouseDrag(const MouseEvent& event) override;
		void cbMouseUp(const MouseEvent& event) override;
		void cbMouseDown(const MouseEvent& event) override;

		virtual void freeze() override;
		virtual void unfreeze() override;
		// cbasecontrol
		bool valueChanged(cpl::CBaseControl *) override;
		void resized() override;
		void rezoomed();

		double getGain();
		//using cpl::Graphics::RGBPixel;
		cpl::Graphics::RGBPixel getSTDColor();
		void buttonSelected(cpl::CBaseControl * c, int index) override;
		void buttonDeselected(cpl::CBaseControl * c, int index) override;


	

	private:
		
		void prepareTransform();
		void doTransform();
		void mapToLinearSpace();
		void mapResonatingSystem();
		void mapFrequencies();
		void aRFTMapping();
		void doFFT();
		void reallocBuffers();
		void initPanelAndControls();
		void computeWindowKernel();
		void setBandwidth();
		void viewChanged();
		int getWindowSize();
		inline float scale(float input, float min, float max);
		inline float invScale(float input, float min, float max);

		template<typename T>
		T * getAudioMemory()
		{
			return reinterpret_cast<T*>(memory.data());
		}

		template<typename T>
		std::size_t getNumAudioElements()
		{
			return memory.size() / sizeof(T);
		}


		template<typename T>
		T * getWorkingMemory()
		{
			return reinterpret_cast<T*>(workspace.data());
		}

		template<typename T>
		std::size_t getNumWorkingElements()
		{
			return workspace.size() / sizeof(T);
		}

		// guis and whatnot
		cpl::COpenGLView renderWindow;
		cpl::CBoxFilter<float, 60> avgFps;
		cpl::CControlPanel panel;
		cpl::CButtonGroup<cpl::CRenderButton> * kviewType;
		cpl::CButtonGroup<cpl::CRenderButton> * koptions;
		cpl::CButtonGroup<cpl::CRenderButton> * kalgorithm;
		cpl::CButtonGroup<cpl::CRenderButton> * kscaleButtons;
		cpl::CButtonGroup<cpl::CRenderButton> * kchannelConf;

		cpl::CKnobEx kaux, kaux2, kgain, kgraph, kdecay, klowDbs, khighDbs;
		cpl::CKnobEx kwindow;
		cpl::CColorKnob kstdColor;

		juce::MouseCursor displayCursor;
		
		// used for the realtime audio stream
		std::vector<char, cpl::AlignmentAllocator<char, 32u>> memory;
		// used for algorithms which is not in-place or temporary workspace memory
		std::vector<char, cpl::AlignmentAllocator<char, 32u>> workspace;
		juce::Image spectrumImg;
		std::unique_ptr<juce::Graphics> spectrumGraphics;
		std::vector<float> mappedFrequencies;
		cpl::MouseDelegate mouseListener;
		// vars
		long long lastFrameTick, renderCycles;
		bool firstResize, isFrozen;
		double lowestCoord, highestCoord;
		DBRange dbs;
		Algorithm::Enum algorithmType;
		Quality::Enum qSettings;
		ViewTypes::Enum viewType;

		float xoffset, yoffset, ydrag, xdrag;
		float expViewScale;
		int windowSize;
		double maxMag;
		unsigned long long processorSpeed; // clocks / sec
		unsigned numChannels;
		unsigned numFilters;
		unsigned imagePtr;
		bool debug;
		bool isXLog, isYLog;
		// data
		std::unique_ptr<char> textbuf;
		AudioBuffer & audioData;
		juce::Rectangle<int> displaySize;
		std::vector<FTFilterResult<float>, cpl::AlignmentAllocator<float, 32>> filterResults;
		std::vector<FTFilterResult<float>, cpl::AlignmentAllocator<float, 32>> filterStates;
		std::vector<double> windowKernel;
		cpl::dsp::CSignalTransform transformer;
		cpl::CPeakFilter<float> flt;
		ResonatorSystem<float> resonator;
		cpl::dsp::CComplexResonator<double, 1> sdft;
		std::vector<std::vector<ResonatorType, cpl::AlignmentAllocator<ResonatorType, 32>>> relayBuffer;

	};
	
	
	
	
	
};

#endif