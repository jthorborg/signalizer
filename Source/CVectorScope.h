
#ifndef _CVECTORSCOPE_H
	#define _CVECTORSCOPE_H

	#include "CommonSignalizer.h"
	#include <cpl/GraphicComponents.h>
	#include <cpl/Utility.h>
	#include <cpl/gui/Controls.h>
	#include <cpl/gui/CPresetWidget.h>
	#include <memory>
	#include <cpl/simd.h>
	#include <cpl/lib/CFIFOEventSystem.h>

	namespace cpl
	{
		namespace OpenGLRendering
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
						
						// describe first left upper part of circle
						// maybe use the parabola like any sane person
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
			protected AudioStream::Listener,
			protected juce::ComponentListener
		{

		public:

			static const double higherAutoGainBounds;
			static const double lowerAutoGainBounds;

			CVectorScope(AudioStream & data);
			virtual ~CVectorScope();

			// Component overrides
			void onGraphicsRendering(Graphics & g) override;
			void mouseWheelMove(const MouseEvent& event, const MouseWheelDetails& wheel) override;
			void mouseDoubleClick(const MouseEvent& event) override;
			void mouseDrag(const MouseEvent& event) override;
			void mouseUp(const MouseEvent& event) override;
			void mouseDown(const MouseEvent& event) override;
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

			std::unique_ptr<juce::Component> createEditor() override;

			// cbasecontrol overrides
			void valueChanged(const cpl::CBaseControl *) override;
			bool stringToValue(const cpl::CBaseControl * ctrl, const std::string & buffer, cpl::iCtrlPrec_t & value) override;
			bool valueToString(const cpl::CBaseControl * ctrl, std::string & buffer, cpl::iCtrlPrec_t value) override;
			void onObjectDestruction(const cpl::CBaseControl::ObjectProxy & destroyedObject) override;
			bool isEditorOpen() const;
			double getGain();
		protected:
			bool onAsyncAudio(const AudioStream & source, AudioStream::DataType ** buffer, std::size_t numChannels, std::size_t numSamples) override;
			void onAsyncChangedProperties(const AudioStream & source, const AudioStream::AudioStreamInfo & before) override;
			void componentBeingDeleted(Component & 	component) override;

			virtual void paint2DGraphics(juce::Graphics & g);

			/// <summary>
			/// Handles all set flags in mtFlags.
			/// Make sure the audioStream is locked while doing this.
			/// </summary>
			virtual void handleFlagUpdates();
		private:

			void deserialize(cpl::CSerializer::Builder & builder, long long int version) override;
			void serialize(cpl::CSerializer::Archiver & archive, long long int version) override;

			template<typename V>
				void vectorGLRendering();

			// vector-accelerated drawing, rendering and processing
			template<typename V>
				void drawPolarPlot(cpl::OpenGLRendering::COpenGLStack &, const AudioStream::AudioBufferAccess &);

			template<typename V>
				void drawRectPlot(cpl::OpenGLRendering::COpenGLStack &, const AudioStream::AudioBufferAccess &);

			template<typename V>
				void drawWireFrame(cpl::OpenGLRendering::COpenGLStack &);

			template<typename V>
				void drawGraphText(cpl::OpenGLRendering::COpenGLStack &, const AudioStream::AudioBufferAccess &);

			template<typename V>
				void drawStereoMeters(cpl::OpenGLRendering::COpenGLStack &, const AudioStream::AudioBufferAccess &);

			template<typename V>
				void runPeakFilter(const AudioStream::AudioBufferAccess &);

			template<typename V>
				void audioProcessing(typename cpl::simd::scalar_of<V>::type ** buffer, std::size_t numChannels, std::size_t numSamples);

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
			cpl::CPresetWidget kpresets;
			juce::MouseCursor displayCursor;

			// vars
			long long lastFrameTick, renderCycles;

			cpl::iCtrlPrec_t envelopeGain;

			struct FilterStates
			{
				AudioStream::DataType envelope[2];
				AudioStream::DataType balance[2][2];
				AudioStream::DataType phase[2];
				
			} filters;

			struct Flags
			{
				cpl::ABoolFlag
					firstRun,
					/// <summary>
					/// Set this if the audio buffer window size was changed from somewhere else.
					/// </summary>
					audioWindowWasResized;
			} mtFlags;

			struct StateOptions
			{
				bool isPolar, normalizeGain, isFrozen, fillPath, fadeHistory, antialias, isEditorOpen, diagnostics, doStereoMeasurements;
				float primitiveSize, rotation;
				float stereoCoeff;
				float envelopeCoeff;
				/// <summary>
				/// A constant factor slower than the stereoCoeff
				/// </summary>
				float secondStereoFilterSpeed;
				juce::Colour colourBackground, colourWire, colourGraph, colourDraw, colourMeter;
				EnvelopeModes envelopeMode;
			} state;
			
			// data
			AudioStream & audioStream;
			//cpl::AudioBuffer audioStreamCopy;
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