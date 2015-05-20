
#ifndef _CSPECTRUM_H
	#define _CSPECTRUM_H

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
	
		
		class CSpectrum
		: 
			public cpl::COpenGLView, 
			protected cpl::CBaseControl::PassiveListener,
			protected cpl::CBaseControl::ValueFormatter,
			protected cpl::CAudioListener,
			protected juce::ComponentListener,
			public cpl::Utility::CNoncopyable
		{

		public:

			CSpectrum(cpl::AudioBuffer & data);
			virtual ~CSpectrum();

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

		protected:

			bool audioCallback(cpl::CAudioSource & source, float ** buffer, std::size_t numChannels, std::size_t numSamples) override;
			void componentBeingDeleted(Component & 	component) override;

		private:

			void initPanelAndControls();

			// guis and whatnot
			cpl::CBoxFilter<double, 60> avgFps;
			cpl::CButton kdiagnostics;
			juce::MouseCursor displayCursor;

			// vars
			long long lastFrameTick, renderCycles;

			struct StateOptions
			{
				bool isEditorOpen, isFrozen, isSuspended;
				bool antialias;

				float primitiveSize;
				juce::Colour colourBackground;
			} state;
			
			// data
			cpl::AudioBuffer & audioStream;
			cpl::AudioBuffer audioStreamCopy;
			juce::Component * editor;
			
			// unused.
			unsigned long long processorSpeed; // clocks / sec
			juce::Point<float> lastMousePos;
			std::vector<std::unique_ptr<juce::OpenGLTexture>> textures;

		};
	
	};

#endif