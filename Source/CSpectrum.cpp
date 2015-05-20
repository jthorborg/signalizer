
#include "CSpectrum.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include "SignalizerDesign.h"
#include <cpl/rendering/OpenGLRasterizers.h>
#include <cpl/simd.h>

namespace Signalizer
{

	std::unique_ptr<juce::Component> CSpectrum::createEditor()
	{
		auto content = new Signalizer::CContentPage();
	
		editor = content;
		editor->addComponentListener(this);
		state.isEditorOpen = editor ? true : false;
		return std::unique_ptr<juce::Component>(content);

	}

	CSpectrum::CSpectrum(cpl::AudioBuffer & data)
	:
		audioStream(data),
		processorSpeed(0), 
		audioStreamCopy(2),
		lastFrameTick(0),
		lastMousePos(),
		editor(nullptr),
		state()
	{
		setOpaque(true);
		processorSpeed = juce::SystemStats::getCpuSpeedInMegaherz();
		initPanelAndControls();
		listenToSource(data[0]);
	}

	void CSpectrum::componentBeingDeleted(Component & component)
	{
		if (&component == editor)
		{
			editor = nullptr;
			state.isEditorOpen = false;
		}
	}


	void CSpectrum::suspend()
	{


	}

	void CSpectrum::resume()
	{


	}

	CSpectrum::~CSpectrum()
	{
		notifyDestruction();
		if (editor)
			editor->removeComponentListener(this);
	}

	void CSpectrum::initPanelAndControls()
	{

		setMouseCursor(juce::MouseCursor::DraggingHandCursor);
	}

	void CSpectrum::save(cpl::CSerializer::Archiver & archive, long long int version)
	{

	}

	void CSpectrum::load(cpl::CSerializer::Builder & builder, long long int version)
	{


	}



	void CSpectrum::freeze()
	{
		state.isFrozen = true;
		std::vector<cpl::CMutex> locks;
		// lock streams firstly if we are synced.
		if (isSynced)
		{
			for (auto & buffer : audioStream)
			{
				locks.emplace_back(buffer);
			}
		}

		for (unsigned i = 0; i < audioStream.size(); ++i)
		{
			audioStream[i].clone(audioStreamCopy[i]);
		}
	}

	void CSpectrum::unfreeze()
	{
		state.isFrozen = false;
	}


	void CSpectrum::mouseWheelMove(const MouseEvent& event, const MouseWheelDetails& wheel)
	{

	}
	void CSpectrum::mouseDoubleClick(const MouseEvent& event)
	{

	}
	void CSpectrum::mouseDrag(const MouseEvent& event)
	{

	}
	void CSpectrum::mouseUp(const MouseEvent& event)
	{

	}
	void CSpectrum::mouseDown(const MouseEvent& event)
	{
		lastMousePos = event.position;
	}

	bool CSpectrum::valueToString(const cpl::CBaseControl * ctrl, std::string & buffer, cpl::iCtrlPrec_t value)
	{
		char buf[200];
		return false;
	}


	bool CSpectrum::stringToValue(const cpl::CBaseControl * ctrl, const std::string & buffer, cpl::iCtrlPrec_t & value)
	{

		return false;
	}
	void CSpectrum::onObjectDestruction(const cpl::CBaseControl::ObjectProxy & destroyedObject)
	{
		// hmmm.....
	}

	void CSpectrum::valueChanged(const cpl::CBaseControl * ctrl)
	{

	}

	template<typename V>
		void CSpectrum::audioProcessing(float ** buffer, std::size_t numChannels, std::size_t numSamples)
		{
			using namespace cpl::simd;

			if (numChannels != 2)
				return;




		}

	bool CSpectrum::audioCallback(cpl::CAudioSource & source, float ** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		return false;
	}

};
