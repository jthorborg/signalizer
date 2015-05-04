
#include "CVectorScope.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include "SignalizerDesign.h"
#include <cpl/rendering/OpenGLRasterizers.h>
#include <cpl/simd.h>

namespace Signalizer
{
	static std::vector<std::string> OperationalModeNames = {"Lissajous", "Polar"};
	static std::vector<std::string> EnvelopeModeNames = {"None", "RMS", "Peak Decay"};
	
	enum class EnvelopeModes : int
	{
		None,
		RMS,
		PeakDecay
	};
	
	static const double lowerAutoGainBounds = cpl::Math::dbToFraction(-120.0);
	static const double higherAutoGainBounds = cpl::Math::dbToFraction(120.0);

	std::unique_ptr<juce::Component> CVectorScope::createEditor()
	{
		auto content = new Signalizer::CContentPage();

		if (auto page = content->addPage("Settings", "icons/svg/wrench.svg"))
		{
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&ktransform, 0);
				page->addSection(section, "Transform");
			}
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&krotation, 0);
				section->addControl(&kwindow, 1);
				section->addControl(&kgain, 0);
				section->addControl(&kenvelopeSmooth, 1);
				section->addControl(&kenvelopeMode, 1);
				section->addControl(&kopMode, 0);
				page->addSection(section, "Utility");
				//
			}
		}

		if (auto page = content->addPage("Rendering", "icons/svg/painting.svg"))
		{
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kantiAlias, 0);
				section->addControl(&kfadeOld, 1);
				section->addControl(&kdrawLines, 2);
				section->addControl(&kdiagnostics, 3);
				page->addSection(section, "Options");
			}
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kdrawingColour, 0);
				section->addControl(&kgraphColour, 0);
				section->addControl(&kbackgroundColour, 0);
				section->addControl(&kskeletonColour, 0);
				section->addControl(&kprimitiveSize, 1);
				page->addSection(section, "Look");
			}
		}

		editor = content;
		editor->addComponentListener(this);

		return std::unique_ptr<juce::Component>(content);

	}

	CVectorScope::CVectorScope(cpl::AudioBuffer & data)
	:
		audioStream(data),
		kwindow("Window size", cpl::CKnobSlider::ControlType::ms),
		krotation("Wave Z-rotation"),
		kgain("Input gain"),
		kgraphColour("Graph colour"),
		kbackgroundColour("Background colour"),
		kdrawingColour("Drawing colour"),
		kskeletonColour("Skeleton colour"),
		kprimitiveSize("Primitive size"),
		kenvelopeSmooth("Env. smooth time", cpl::CKnobSlider::ControlType::ms),
		processorSpeed(0), 
		audioStreamCopy(2),
		lastFrameTick(0), 
		isFrozen(false),
		lastMousePos(),
		normalizeGain(false),
		envelopeSmooth(0.99998),
		envelopeFilters{},
		envelopeGain(1),
		editor(nullptr)
	{
		setOpaque(true);
		textbuf = std::unique_ptr<char>(new char[300]);
		processorSpeed = juce::SystemStats::getCpuSpeedInMegaherz();
		initPanelAndControls();
		listenToSource(data[0]);
	}

	void CVectorScope::componentBeingDeleted(Component & component)
	{
		if (&component == editor)
		{
			editor = nullptr;
		}
	}

	bool CVectorScope::isEditorOpen() const
	{
		return editor ? true : false;
	}

	void CVectorScope::suspend()
	{


	}

	void CVectorScope::resume()
	{


	}

	juce::Component * CVectorScope::getWindow()
	{
		return this;
	}

	CVectorScope::~CVectorScope()
	{
		notifyDestruction();
		if (editor)
			editor->removeComponentListener(this);
	}

	void CVectorScope::initPanelAndControls()
	{
		// listeners
		kwindow.bAddPassiveChangeListener(this);
		kwindow.bAddFormatter(this);
		kgain.bAddFormatter(this);
		kgain.bAddPassiveChangeListener(this);
		krotation.bAddFormatter(this);
		kprimitiveSize.bAddFormatter(this);
		kenvelopeMode.bAddPassiveChangeListener(this);
		kopMode.bAddPassiveChangeListener(this);
		kenvelopeSmooth.bAddPassiveChangeListener(this);
		// buttons n controls
		kantiAlias.bSetTitle("Antialias");
		kantiAlias.setToggleable(true);
		kfadeOld.bSetTitle("Fade older points");
		kfadeOld.setToggleable(true);
		kdrawLines.bSetTitle("Interconnect samples");
		kdrawLines.setToggleable(true);
		kdiagnostics.bSetTitle("Diagnostics");
		kdiagnostics.setToggleable(true);
		kenvelopeMode.bSetTitle("Auto-gain mode");
		kenvelopeMode.setValues(EnvelopeModeNames);
		
		kenvelopeMode.bSetValue(0);
		kopMode.bSetValue(0);
		
		// design
		kopMode.setValues(OperationalModeNames);
		kopMode.bSetTitle("Operational mode");
		
		
		// descriptions.
		kwindow.bSetDescription("The size of the displayed time window.");
		kgain.bSetDescription("How much the input (x,y) is scaled (or the input gain)" \
							 " - additional transform that only affects the waveform, and not the graph");
		krotation.bSetDescription("The amount of degrees to rotate the waveform around the origin (z-rotation)"\
			" - additional transform that only affects the waveform, and not the graph.");
		kantiAlias.bSetDescription("Antialiases rendering (if set - see global settings for amount). May slow down rendering.");
		kfadeOld.bSetDescription("If set, gradually older samples will be faded linearly.");
		kdrawLines.bSetDescription("If set, interconnect samples linearly.");
		kdrawingColour.bSetDescription("The main colour to paint with.");
		kgraphColour.bSetDescription("The colour of the graph.");
		kbackgroundColour.bSetDescription("The background colour of the view.");
		kdiagnostics.bSetDescription("Toggle diagnostic information in top-left corner.");
		kskeletonColour.bSetDescription("The colour of the box skeleton indicating the OpenGL camera clip box.");
		kprimitiveSize.bSetDescription("The size of the rendered primitives (eg. lines or points).");
		kenvelopeMode.bSetDescription("Monitors the audio stream and automatically scales the input gain such that it approaches unity intensity (envelope following).");
		kenvelopeSmooth.bSetDescription("Responsiveness - or the time it takes for the envelope follower to decay.");
		kopMode.bSetDescription("Changes the presentation of the data - Lissajous is the classic XY mode on oscilloscopes, while the polar mode is a wrapped circle of the former.");
		// design
		setMouseCursor(juce::MouseCursor::DraggingHandCursor);
	}

	void CVectorScope::save(cpl::CSerializer::Archiver & archive, long long int version)
	{
		archive << kwindow;
		archive << kgain;
		archive << krotation;
		archive << kantiAlias;
		archive << kfadeOld;
		archive << kdiagnostics;
		archive << kdrawLines;
		archive << kgraphColour;
		archive << kbackgroundColour;
		archive << kdrawingColour;
		archive << ktransform.getTransform3D();
		archive << kskeletonColour;
		archive << kprimitiveSize;
		archive << kenvelopeMode;
	}

	void CVectorScope::load(cpl::CSerializer::Builder & builder, long long int version)
	{
		builder >> kwindow;
		builder >> kgain;
		builder >> krotation;
		builder >> kantiAlias;
		builder >> kfadeOld;
		builder >> kdiagnostics;
		builder >> kdrawLines;
		builder >> kgraphColour;
		builder >> kbackgroundColour;
		builder >> kdrawingColour;
		builder >> ktransform.getTransform3D();
		builder >> kskeletonColour;
		builder >> kprimitiveSize;
		builder >> kenvelopeMode;
		
		ktransform.syncEditor();

	}



	void CVectorScope::freeze()
	{
		isFrozen = true;
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

	void CVectorScope::unfreeze()
	{
		isFrozen = false;
	}


	void CVectorScope::mouseWheelMove(const MouseEvent& event, const MouseWheelDetails& wheel)
	{
		auto amount = wheel.deltaY;
		auto normalizedSign = amount > 0 ? 1.0f : -1.0f;
		if (event.mods.isCtrlDown())
		{
			// increase gain
			kgain.bSetValue(kgain.bGetValue() + amount / 20);
		}
		else /* zoom graph */
		{
			
			auto mouseXCoord = (float(event.x) / getWidth()) * 2 - 1;
			auto mouseYCoord = (float(event.y) / getHeight()) * 2 - 1;

			auto & matrix = ktransform.getTransform3D();
			//matrix.scale.x += amount;
			//matrix.scale.y += amount;
			auto actualAmount = (1 + amount / 5) * matrix.scale.z;
			auto deltaIncrease = (actualAmount - matrix.scale.z) / matrix.scale.z;
			matrix.scale.z = actualAmount;

			matrix.scale.x *= 1 + deltaIncrease;
			matrix.scale.y *= 1 + deltaIncrease;
			//matrix.position.x *=  1 + normalizedSign * deltaIncrease;
			//matrix.position.y *= 1 + normalizedSign * deltaIncrease;
			ktransform.syncEditor();

		}
		
	}
	void CVectorScope::mouseDoubleClick(const MouseEvent& event)
	{
		
		if (event.mods.isLeftButtonDown())
		{
			// reset all zooming, offsets etc. when doubleclicking left
			kgain.bSetValue(0.5f);
			auto & matrix = ktransform.getTransform3D();
			matrix.position.y = matrix.position.x = 0;
			ktransform.syncEditor();
		}
	}
	void CVectorScope::mouseDrag(const MouseEvent& event)
	{
		auto & matrix = ktransform.getTransform3D();
		auto factor = float(getWidth()) / getHeight();
		auto deltaDifference = event.position - lastMousePos;
		if (event.mods.isCtrlDown())
		{
			matrix.rotation.y = std::fmod(deltaDifference.x * 0.3f + matrix.rotation.y, 360.f);
			matrix.rotation.x = std::fmod(deltaDifference.y * 0.3f + matrix.rotation.x, 360.f);
		}
		else
		{
			matrix.position.x += deltaDifference.x / 500.f;
			matrix.position.y += factor * -deltaDifference.y / 500.f;
		}
		ktransform.syncEditor();
		
		
		lastMousePos = event.position;
	}
	void CVectorScope::mouseUp(const MouseEvent& event)
	{

	}
	void CVectorScope::mouseDown(const MouseEvent& event)
	{
		lastMousePos = event.position;
	}

	bool CVectorScope::valueToString(const cpl::CBaseControl * ctrl, std::string & buffer, cpl::iCtrlPrec_t value)
	{
		char buf[200];
		if (ctrl == &kgain)
		{
			sprintf(buf, "%.2f dB", cpl::Math::fractionToDB(getGain()));
			buffer = buf;
			return true;
		}
		else if (ctrl == &kwindow)
		{
			auto bufLength = cpl::Math::round<int>( value * 1000);
			sprintf(buf, "%d ms", bufLength);
			buffer = buf;
			return true;
		}
		else if (ctrl == &krotation)
		{
			sprintf(buf, "%.2f degs", value * 360);
			buffer = buf;
			return true;
		}
		else if (ctrl == &kprimitiveSize)
		{
			sprintf(buf, "%.2f pts", value * 10);
			buffer = buf;
			return true;
		}
		return false;
	}

	void CVectorScope::setGainAsFraction(double newFraction)
	{
		auto newValue = cpl::Math::UnityScale::Inv::exp(newFraction, lowerAutoGainBounds, higherAutoGainBounds);
		kgain.bSetValue(cpl::Math::confineTo(newValue, 0, 1));
	}

	double CVectorScope::getGain()
	{
		if (normalizeGain)
			return envelopeGain;
		return cpl::Math::UnityScale::exp(kgain.bGetValue(), lowerAutoGainBounds, higherAutoGainBounds);
	}

	double CVectorScope::mapScaleToFraction(double valueInDBs)
	{
		return cpl::Math::UnityScale::Inv::exp(cpl::Math::dbToFraction(valueInDBs), lowerAutoGainBounds, higherAutoGainBounds);
	}

	bool CVectorScope::stringToValue(const cpl::CBaseControl * ctrl, const std::string & buffer, cpl::iCtrlPrec_t & value)
	{

		if (ctrl == &kgain)
		{
			char * endPtr(nullptr);
			cpl::iCtrlPrec_t newVal = strtod(buffer.c_str(), &endPtr);
			if (endPtr > buffer.data())
			{
				value = mapScaleToFraction(newVal);
				return true;
			}
		}
		else if (ctrl == &kwindow)
		{
			double newValue;
			char * endPtr = nullptr;
			newValue = strtod(buffer.c_str(), &endPtr);
			if (endPtr > 0)
			{
				value = cpl::Math::confineTo(newValue / 1000.0, 0.0, 1.0);
				return true;
			}
		}
		else if (ctrl == &krotation)
		{
			double newValue;
			char * endPtr = nullptr;
			newValue = strtod(buffer.c_str(), &endPtr);
			if (endPtr > 0)
			{
				value = cpl::Math::confineTo(fmod(newValue, 360) / 360.0, 0.0, 1.0);
				return true;
			}
		}
		else if (ctrl == &kprimitiveSize)
		{
			double newValue;
			char * endPtr = nullptr;
			newValue = strtod(buffer.c_str(), &endPtr);
			if (endPtr > 0)
			{
				value = cpl::Math::confineTo(newValue / 10.0, 0.0, 1.0);
				return true;
			}
		}
		return false;
	}
	void CVectorScope::onObjectDestruction(const cpl::CBaseControl::ObjectProxy & destroyedObject)
	{
		// hmmm.....
	}

	void CVectorScope::valueChanged(const cpl::CBaseControl * ctrl)
	{
		if (ctrl == &kwindow)
		{
			auto bufLength = ctrl->bGetValue() * 1000;
			for (auto & buffer : audioStream)
			{
				buffer.setLength(bufLength);
			}
			return;
		}
		else if (ctrl == &kenvelopeMode)
		{
			envelopeMode = kenvelopeMode.getZeroBasedSelIndex<EnvelopeModes>();
			normalizeGain = envelopeMode != EnvelopeModes::None;
			if (!normalizeGain)
				envelopeGain = kgain.bGetValue();
		}
		else if(ctrl == &kenvelopeSmooth)
		{
			envelopeSmooth = std::exp(-1.0 / (kenvelopeSmooth.bGetValue() * audioStream[0].sampleRate));
		}
		else if (ctrl == &kgain)
		{
			if (!normalizeGain)
				envelopeGain = kgain.bGetValue();
		}
	}
	template<typename V>
	void CVectorScope::processAudioBuffer(float ** buffers, std::size_t channels, std::size_t samples)
	{
		/*
		using namespace cpl;
		using namespace cpl::simd;

		auto const loopIncrement = elements_of<V>::value;

		auto const smoothPole = set1<V>(envelopeSmooth * loopIncrement);
		auto leftFilter = 0;
		suitable_container<V> startLeft, startRight;
		double smoothPole = envelopeSmooth;

		for (int i = 0; i < loopIncrement; ++i)
		{
			startLeft[i] = envelopeFilters[0] * smoothPole;
			startRight[i] = envelopeFilter[1] * smoothPole;
			smoothPole *= smoothPole;
		}


		if (isPeakDetection)
		{
			// doesn't matter much if a remainder is missing..
			for (std::size_t i = 0; i < samples; i += loopIncrement)
			{

			}
		}
		*/

	}


	bool CVectorScope::audioCallback(cpl::CAudioSource & source, float ** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		if (normalizeGain)
		{
			const bool isPeakDetection = true;

			if (numChannels != 2)
				return false;

			double currentEnvelope = 0;

			if (envelopeMode == EnvelopeModes::PeakDecay)
			{

				for (std::size_t i = 0; i < numSamples; ++i)
				{
					// square here.
					double input = buffer[0][i] * buffer[0][i];

					if (input > envelopeFilters[0])
						envelopeFilters[0] = input;
					else
						envelopeFilters[0] *= envelopeSmooth;
					input = buffer[1][i] * buffer[1][i];

					if (input > envelopeFilters[1])
						envelopeFilters[1] = input;
					else
						envelopeFilters[1] *= envelopeSmooth;

				}
				currentEnvelope = 1.0 / std::max(std::sqrt(envelopeFilters[0]), std::sqrt(envelopeFilters[1]));
			}
			else
			{
				for (std::size_t i = 0; i < numSamples; ++i)
				{
					// square here.
					double input = buffer[0][i] * buffer[0][i];
					envelopeFilters[0] = input + envelopeSmooth * (envelopeFilters[0] - input);
					input = buffer[1][i] * buffer[1][i];
					envelopeFilters[1] = input + envelopeSmooth * (envelopeFilters[1] - input);
				}
				
				currentEnvelope = 1.0 / (2 * std::max(std::sqrt(envelopeFilters[0]), std::sqrt(envelopeFilters[1])));

			}



			if (std::isnormal(currentEnvelope))
			{
				envelopeGain = cpl::Math::confineTo(currentEnvelope, lowerAutoGainBounds, higherAutoGainBounds);
			}
		}
		return false;
	}


};
