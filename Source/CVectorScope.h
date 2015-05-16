
#ifndef _CVECTORSCOPE_H
	#define _CVECTORSCOPE_H

	#include <cpl/Common.h>
	#include <cpl/CAudioBuffer.h>
	#include <cpl/CViews.h>
	#include <cpl/GraphicComponents.h>
	#include <cpl/Utility.h>
	#include <cpl/gui/Controls.h>

	#include <memory>

	namespace cpl
	{
		namespace OpenGLEngine
		{
			class COpenGLStack;
			
		};
	};


	namespace Signalizer
	{
		
		template<typename T, std::size_t size>
			class LookupTable
			{
			public:
				typedef T Ty;
				static const std::size_t tableSize = size;
				
				inline Ty linearLookup(Ty dx) const noexcept
				{
					Ty scaled = dx * tableSize;
					std::size_t x1 = std::size_t(scaled);
					std::size_t x2 = x1 + 1;
					Ty fraction = scaled - x1;
					
					return table[x1] * (Ty(1) - fraction) + table[x2] * fraction;
				}
				
				Ty table[tableSize + 1];
			};
		
		template<typename T, std::size_t size>
			class QuarterCircleLut : public LookupTable<T, size>
			{
			public:
				QuarterCircleLut()
				{
					double increase = 1.0 / (size - 1);
					for (int i = 0; i < size; ++i)
					{
						
						// describe frist left upper part of circle
						this->table[i] = (T)std::sin(std::acos(1.0 - increase * i));
					}
					this->table[this->tableSize] = (T)1;
				}
			};

		enum class EnvelopeModes : int;
		
		class CVectorScope 
		: 
			public cpl::COpenGLView, 
			protected cpl::CBaseControl::PassiveListener,
			protected cpl::CBaseControl::ValueFormatter,
			protected cpl::CAudioListener,
			protected juce::ComponentListener,
			public cpl::Utility::CNoncopyable
		{

		public:

			static const double higherAutoGainBounds;
			static const double lowerAutoGainBounds;

			CVectorScope(cpl::AudioBuffer & data);
			virtual ~CVectorScope();

			// Component overrides
			void paint(Graphics & g) override;
			void mouseWheelMove(const MouseEvent& event, const MouseWheelDetails& wheel) override;
			void mouseDoubleClick(const MouseEvent& event) override;
			void mouseDrag(const MouseEvent& event) override;
			void mouseUp(const MouseEvent& event) override;
			void mouseDown(const MouseEvent& event) override;
			// OpenGLRender overrides
			void renderOpenGL() override;
			void initOpenGL() override;
			void closeOpenGL() override;
			// View overrides
			juce::Component * getWindow() override;
			void suspend() override;
			void resume() override;
			void freeze() override;
			void unfreeze() override;

			std::unique_ptr<juce::Component> createEditor() override;
			// CSerializer overrides
			void load(cpl::CSerializer::Builder & builder, long long int version) override;
			void save(cpl::CSerializer::Archiver & archive, long long int version) override;

			// cbasecontrol overrides
			void valueChanged(const cpl::CBaseControl *) override;
			bool stringToValue(const cpl::CBaseControl * ctrl, const std::string & buffer, cpl::iCtrlPrec_t & value) override;
			bool valueToString(const cpl::CBaseControl * ctrl, std::string & buffer, cpl::iCtrlPrec_t value) override;
			void onObjectDestruction(const cpl::CBaseControl::ObjectProxy & destroyedObject) override;
			bool isEditorOpen() const;
			double getGain();
		protected:
			bool audioCallback(cpl::CAudioSource & source, float ** buffer, std::size_t numChannels, std::size_t numSamples) override;
			void componentBeingDeleted(Component & 	component) override;

		private:

			// vector-accelerated drawing, rendering and processing
			template<typename V>
				void drawPolarPlot(cpl::OpenGLEngine::COpenGLStack &, const cpl::AudioBuffer &);

			template<typename V>
				void drawRectPlot(cpl::OpenGLEngine::COpenGLStack &, const cpl::AudioBuffer &);

			template<typename V>
				void drawWireFrame(cpl::OpenGLEngine::COpenGLStack &);

			template<typename V>
				void drawGraphText(cpl::OpenGLEngine::COpenGLStack &, const cpl::AudioBuffer &);

			template<typename V>
				void drawStereoMeters(cpl::OpenGLEngine::COpenGLStack &, const cpl::AudioBuffer &);

			template<typename V>
				void runPeakFilter(cpl::AudioBuffer & buffer, std::size_t samples);

			template<typename V>
				void audioProcessing(float ** buffer, std::size_t numChannels, std::size_t numSamples);

			void setGainAsFraction(double newFraction);
			double mapScaleToFraction(double dbs);
			void initPanelAndControls();

			// guis and whatnot
			cpl::CBoxFilter<double, 60> avgFps;

			cpl::CButton kantiAlias, kfadeOld, kdrawLines, kdiagnostics;
			cpl::CKnobSlider kwindow, krotation, kgain, kprimitiveSize, kenvelopeSmooth, kstereoSmooth;
			cpl::CColourControl kdrawingColour, kgraphColour, kbackgroundColour, kskeletonColour, kmeterColour;
			cpl::CTransformWidget ktransform;
			cpl::CComboBox kopMode, kenvelopeMode;
			juce::MouseCursor displayCursor;

			// vars
			long long lastFrameTick, renderCycles;

			cpl::iCtrlPrec_t envelopeGain;

			struct FilterStates
			{
				float envelope[2];
				float balance[2][2];
				float phase[2];
				
			} filters;

			struct StateOptions
			{
				bool isPolar, normalizeGain, isFrozen, fillPath, fadeHistory, antialias, isEditorOpen, diagnostics, doStereoMeasurements;
				float primitiveSize, rotation;
				float stereoCoeff;
				float envelopeCoeff;
				float secondStereoFilterSpeed;
				juce::Colour colourBackground, colourWire, colourGraph, colourDraw, colourMeter;
				EnvelopeModes envelopeMode;
			} state;
			
			// data
			cpl::AudioBuffer & audioStream;
			cpl::AudioBuffer audioStreamCopy;
			cpl::Utility::LazyPointer<QuarterCircleLut<GLfloat, 128>> circleData;
			juce::Component * editor;
			
			// unused.
			std::unique_ptr<char> textbuf;
			unsigned long long processorSpeed; // clocks / sec
			juce::Point<float> lastMousePos;
			std::vector<std::unique_ptr<juce::OpenGLTexture>> textures;

		};
	
	};

#endif