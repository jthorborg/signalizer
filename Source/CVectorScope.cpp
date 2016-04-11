
#include "CVectorScope.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include "SignalizerDesign.h"
#include <cpl/rendering/OpenGLRasterizers.h>
#include <cpl/simd.h>
#include <cpl/LexicalConversion.h>

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
	
	enum class OperationalModes : int
	{
		Lissajous,
		Polar
	};
	
	const double CVectorScope::lowerAutoGainBounds = cpl::Math::dbToFraction(-120.0);
	const double CVectorScope::higherAutoGainBounds = cpl::Math::dbToFraction(120.0);

	std::unique_ptr<juce::Component> CVectorScope::createEditor()
	{
		auto content = new Signalizer::CContentPage();

		if (auto page = content->addPage("Settings", "icons/svg/gear.svg"))
		{
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&ktransform, 0);
				page->addSection(section, "Transform");
			}
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kenvelopeMode, 0);
				section->addControl(&kenvelopeSmooth, 0);
				section->addControl(&kgain, 0);
				
				section->addControl(&kopMode, 1);
				section->addControl(&kstereoSmooth, 1);


				section->addControl(&krotation, 0);
				section->addControl(&kwindow, 1);


				page->addSection(section, "Utility");
				//
			}
		}

		if (auto page = content->addPage("Rendering", "icons/svg/brush.svg"))
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
				section->addControl(&kmeterColour, 1);
				section->addControl(&kprimitiveSize, 1);
				page->addSection(section, "Look");
			}
		}
		if (auto page = content->addPage("Utility", "icons/svg/wrench.svg"))
		{
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kpresets, 0);
				page->addSection(section, "Presets");
			}
		}
		editor = content;
		editor->addComponentListener(this);
		state.isEditorOpen = editor ? true : false;
		return std::unique_ptr<juce::Component>(content);

	}

	CVectorScope::CVectorScope(AudioStream & data)
	:
		COpenGLView("Vectorscope view"),
		audioStream(data),
		kwindow("Window size", cpl::CKnobSlider::ControlType::ms),
		krotation("Wave Z-rotation"),
		kgain("Input gain"),
		kgraphColour("Graph colour"),
		kbackgroundColour("Background colour"),
		kdrawingColour("Drawing colour"),
		kskeletonColour("Skeleton colour"),
		kprimitiveSize("Primitive size"),
		kmeterColour("Meter colour"),
		kenvelopeSmooth("Env. window", cpl::CKnobSlider::ControlType::ms),
		kstereoSmooth("Stereo window", cpl::CKnobSlider::ControlType::ms),
		processorSpeed(0), 
		//audioStreamCopy(2),
		lastFrameTick(0),
		lastMousePos(),
		envelopeGain(1),
		editor(nullptr),
		state(),
		filters(),
		kpresets(this, "vectorscope"),
		oldWindowSize(-1)
	{
		mtFlags.firstRun = true;
		state.secondStereoFilterSpeed = 0.25f;
		state.doStereoMeasurements = true;
		setOpaque(true);
		textbuf = std::unique_ptr<char>(new char[300]);
		processorSpeed = juce::SystemStats::getCpuSpeedInMegaherz();
		initPanelAndControls();
		listenToSource(audioStream);
	}

	void CVectorScope::componentBeingDeleted(Component & component)
	{
		if (&component == editor)
		{
			editor = nullptr;
			state.isEditorOpen = false;
		}
	}


	void CVectorScope::suspend()
	{
		oldWindowSize = kwindow.bGetValue();

	}

	void CVectorScope::resume()
	{
		if(oldWindowSize != -1)
		{
			kwindow.bSetValue(oldWindowSize);
			kwindow.bForceEvent(); // in case they (still) are the same
		}


	}

	juce::Component * CVectorScope::getWindow()
	{
		return this;
	}

	CVectorScope::~CVectorScope()
	{
		// TODO: detach listener
		notifyDestruction();
		if (editor)
			editor->removeComponentListener(this);
	}

	void CVectorScope::initPanelAndControls()
	{
		// listeners
		kwindow.bAddPassiveChangeListener(this);
		kgain.bAddPassiveChangeListener(this);
		kenvelopeMode.bAddPassiveChangeListener(this);
		kopMode.bAddPassiveChangeListener(this);
		kenvelopeSmooth.bAddPassiveChangeListener(this);
		krotation.bAddPassiveChangeListener(this);
		kantiAlias.bAddPassiveChangeListener(this);
		kfadeOld.bAddPassiveChangeListener(this);
		kdrawLines.bAddPassiveChangeListener(this);
		kdrawingColour.bAddPassiveChangeListener(this);
		kgraphColour.bAddPassiveChangeListener(this);
		kskeletonColour.bAddPassiveChangeListener(this);
		kbackgroundColour.bAddPassiveChangeListener(this);
		kprimitiveSize.bAddPassiveChangeListener(this);
		kdiagnostics.bAddPassiveChangeListener(this);
		kopMode.bAddPassiveChangeListener(this);
		kstereoSmooth.bAddPassiveChangeListener(this);
		kmeterColour.bAddPassiveChangeListener(this);
		// formatters
		kwindow.bAddFormatter(this);
		kgain.bAddFormatter(this);
		kprimitiveSize.bAddFormatter(this);
		krotation.bAddFormatter(this);


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
		kmeterColour.bSetDescription("The colour of the stereo meters (balance and phase)");
		kprimitiveSize.bSetDescription("The size of the rendered primitives (eg. lines or points).");
		kenvelopeMode.bSetDescription("Monitors the audio stream and automatically scales the input gain such that it approaches unity intensity (envelope following).");
		kenvelopeSmooth.bSetDescription("Responsiveness (RMS window size) - or the time it takes for the envelope follower to decay.");
		kopMode.bSetDescription("Changes the presentation of the data - Lissajous is the classic XY mode on oscilloscopes, while the polar mode is a wrapped circle of the former.");
		kstereoSmooth.bSetDescription("Responsiveness (RMS window size) - or the time it takes for the stereo meters to follow.");

		// design
		setMouseCursor(juce::MouseCursor::DraggingHandCursor);
	}

	void CVectorScope::serialize(cpl::CSerializer::Archiver & archive, cpl::Version version)
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
		archive << kenvelopeSmooth;
		archive << kopMode;
		archive << kstereoSmooth;
		archive << kmeterColour;
	}

	void CVectorScope::deserialize(cpl::CSerializer::Builder & builder, cpl::Version version)
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
		builder >> kenvelopeSmooth;
		builder >> kopMode;
		builder >> kstereoSmooth;
		builder >> kmeterColour;

		ktransform.syncEditor();

	}

	void CVectorScope::handleFlagUpdates()
	{
		bool firstRun = false;

		if (mtFlags.firstRun.cas())
		{
			firstRun = true;
		}

		if (mtFlags.audioWindowWasResized.cas())
		{
			if (!firstRun)
			{
				auto capacity = audioStream.getAudioHistoryCapacity();
				cpl::GUIUtils::MainEvent(*this, 
					[=] 
					{
						if (capacity == 0)
							kwindow.bSetInternal(0);
						else
							kwindow.bSetInternal(double(audioStream.getAudioHistorySize()) / capacity);
						kwindow.bRedraw();
					}
				);
			}
		}

	}

	void CVectorScope::freeze()
	{
		state.isFrozen = true;
	}

	void CVectorScope::unfreeze()
	{
		state.isFrozen = false;
	}


	void CVectorScope::mouseWheelMove(const MouseEvent& event, const MouseWheelDetails& wheel)
	{
		auto amount = wheel.deltaY;
		if (event.mods.isCtrlDown())
		{
			// increase gain
			kgain.bSetValue(kgain.bGetValue() + amount / 20);
		}
		else /* zoom graph */
		{

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
			auto dbVal = cpl::Math::UnityScale::exp(kgain.bGetValue(), lowerAutoGainBounds, higherAutoGainBounds);
			sprintf(buf, "%.2f dB", cpl::Math::fractionToDB(dbVal));
			buffer = buf;
			return true;
		}
		else if (ctrl == &kwindow)
		{
			auto sampleRate = audioStream.getInfo().sampleRate.load(std::memory_order_acquire);
			auto bufLength = cpl::Math::round<int>(value * (1000 * audioStream.getAudioHistoryCapacity() / sampleRate));
			if (sampleRate > 0)
				sprintf(buf, "%d ms", bufLength);
			else
				sprintf(buf, "sample rate error");
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
		return envelopeGain;
	}

	bool CVectorScope::handleMessage(cpl::CMessageSystem::CoalescedMessage & m)
	{
		return false;
	}

	double CVectorScope::mapScaleToFraction(double valueInDBs)
	{
		return cpl::Math::UnityScale::Inv::exp(cpl::Math::dbToFraction(valueInDBs), lowerAutoGainBounds, higherAutoGainBounds);
	}

	bool CVectorScope::stringToValue(const cpl::CBaseControl * ctrl, const std::string & buffer, cpl::iCtrlPrec_t & value)
	{

		if (ctrl == &kgain)
		{
			double newValue;
			if (cpl::lexicalConversion(buffer, newValue))
			{
				value = mapScaleToFraction(newValue);
				return true;
			}
		}
		else if (ctrl == &kwindow)
		{
			double newValue;

			if (buffer.find_first_of("ms") != std::string::npos)
			{
				if (cpl::lexicalConversion(buffer, newValue))
				{
					value = cpl::Math::confineTo(newValue / (1000 * audioStream.getAudioHistoryCapacity() / audioStream.getInfo().sampleRate), 0.0, 1.0);
					return true;
				}
			}
			else
			{
				if (cpl::lexicalConversion(buffer, newValue))
				{
					value = cpl::Math::confineTo(newValue / audioStream.getAudioHistoryCapacity(), 0.0, 1.0);
					return true;
				}
			}

		}
		else if (ctrl == &krotation)
		{
			double newValue;
			if (cpl::lexicalConversion(buffer, newValue))
			{
				value = cpl::Math::confineTo(fmod(newValue, 360) / 360.0, 0.0, 1.0);
				return true;
			}
		}
		else if (ctrl == &kprimitiveSize)
		{
			double newValue;
			if (cpl::lexicalConversion(buffer, newValue))
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
			struct RetryResizer
			{
				RetryResizer(CVectorScope * h) : handle(h) {};
				CVectorScope * handle;

				void operator()()
				{
					auto currentCapacity = handle->audioStream.getAudioHistoryCapacity();
					if (currentCapacity > 0)
					{
						auto bufLength = cpl::Math::round<std::size_t>(handle->kwindow.bGetValue() * handle->audioStream.getAudioHistoryCapacity());
						handle->audioStream.setAudioHistorySize(bufLength);
					}
					else
					{
						cpl::GUIUtils::FutureMainEvent(200, RetryResizer(handle), handle);
					}

				}
			};
			RetryResizer(this)();
			return;
		}
		else if (ctrl == &kenvelopeMode)
		{
			state.envelopeMode = kenvelopeMode.getZeroBasedSelIndex<EnvelopeModes>();
			state.normalizeGain = state.envelopeMode != EnvelopeModes::None;
			// TODO: Consider if needed
			/*if (!state.normalizeGain)
				envelopeGain = kgain.bGetValue();*/
		}
		else if(ctrl == &kenvelopeSmooth)
		{
			state.envelopeCoeff = std::exp(-1.0 / (kenvelopeSmooth.bGetValue() * audioStream.getInfo().sampleRate));
		}
		else if (ctrl == &kstereoSmooth)
		{
			state.stereoCoeff = std::exp(-1.0 / (kstereoSmooth.bGetValue() * audioStream.getInfo().sampleRate));
		}
		else if (ctrl == &kgain)
		{
			envelopeGain = cpl::Math::UnityScale::exp(kgain.bGetValue(), lowerAutoGainBounds, higherAutoGainBounds);
		}
		else if(ctrl == &kopMode)
		{
			state.isPolar = kopMode.getZeroBasedSelIndex<OperationalModes>() == OperationalModes::Polar;
		}
		else if(ctrl == &kantiAlias)
		{
			state.antialias = ctrl->bGetBoolState();
		}
		else if(ctrl == &kfadeOld)
		{
			state.fadeHistory = ctrl->bGetBoolState();
		}
		else if(ctrl == &kdrawLines)
		{
			state.fillPath = ctrl->bGetBoolState();
		}
		else if(ctrl == &kdiagnostics)
		{
			state.diagnostics = ctrl->bGetBoolState();
		}
		else if(ctrl == &krotation)
		{
			state.rotation = (float)ctrl->bGetValue();
		}
		else if(ctrl == &kprimitiveSize)
		{
			state.primitiveSize = (float)ctrl->bGetValue();
		}
		else if(ctrl == &kdrawingColour)
		{
			state.colourDraw = kdrawingColour.getControlColourAsColour();
		}
		else if(ctrl == &kgraphColour)
		{
			state.colourGraph = kgraphColour.getControlColourAsColour();
		}
		else if(ctrl == &kskeletonColour)
		{
			state.colourWire = kskeletonColour.getControlColourAsColour();
		}
		else if(ctrl == &kbackgroundColour)
		{
			state.colourBackground = kbackgroundColour.getControlColourAsColour();
		}
		else if (ctrl == &kmeterColour)
		{
			state.colourMeter = kmeterColour.getControlColourAsColour();
		}
	}

	template<typename V>
		void CVectorScope::audioProcessing(typename cpl::simd::scalar_of<V>::type ** buffer, std::size_t numChannels, std::size_t numSamples)
		{
			using namespace cpl::simd;
			typedef typename scalar_of<V>::type T;
			if (numChannels != 2)
				return;

			// if all options are turned off, just return.
			if ((!state.normalizeGain && state.envelopeMode != EnvelopeModes::RMS) && !state.doStereoMeasurements)
				return;

			T filterEnv[2] = { filters.envelope[0], filters.envelope[1] };
			T stereoPoles[2] = { state.stereoCoeff, std::pow(state.stereoCoeff, state.secondStereoFilterSpeed) };

			const T cosineRotation = (T)std::cos(M_PI * 135 / 180);
			const T sineRotation = (T)std::sin(M_PI * 135 / 180);
			const V vMatrixReal = set1<V>(cosineRotation);
			const V vMatrixImag = set1<V>(sineRotation);

			const V vDummyAngle = set1<V>((T)(M_PI * 0.25));
			const V vZero = zero<V>();

			auto const loopIncrement = elements_of<V>::value;

			// ensure a perfect multiple and no buffer overrun
			numSamples -= numSamples & (loopIncrement - 1);

			suitable_container<V> outputPhases;

			for (std::size_t n = 0; n < numSamples; n += loopIncrement)
			{
				// some of the heavy math done in vector mode..
				const V vLeft = loadu<V>(buffer[0] + n);
				const V vRight = loadu<V>(buffer[1] + n); 

				// rotate 235 degrees...
				const V vX = vLeft * vMatrixReal - vRight * vMatrixImag;
				const V vY = vRight * vMatrixImag + vLeft * vMatrixReal;

				// check if both axes are zero
				const V vZeroAxes = vand((V)(vX == vZero), (V)(vY == vZero));

				// compute the phase angle and replace the zero vector elements with the dummy angle (to avoid nans, +/-infs are defined)
				const V vRadians = atan(vY / vX);
				const V vPhaseAngle = vselect(vDummyAngle, vRadians, vZeroAxes);

				outputPhases = vPhaseAngle;
				// the phase angle is discontinuous around PI, so we take the cosine
				// to avoid the discontinuity and give a small smoothing
				// for some arcane fucking reason, the phase response of this function IN THIS CONTEXT is wrong
				// testing reveals no problems, it just doesn't work right here. Uncomment if you find a fix.
#pragma message cwarn("Errornous cosine output here.")
				//outputPhases = cos(outputPhases * set1<V>(2));

				// scalar code segment for recursive IIR filters..
				for (std::size_t i = 0; i < loopIncrement; ++i)
				{
					auto const z = n + i;
					// collect squared inputs (really just a cheap abs)
					const auto lSquared = buffer[0][z] * buffer[0][z];
					const auto rSquared = buffer[1][z] * buffer[1][z];
					
					// average envelope 
					filterEnv[0] = lSquared + state.envelopeCoeff * (filterEnv[0] - lSquared);
					filterEnv[1] = rSquared + state.envelopeCoeff * (filterEnv[1] - rSquared);

					// balance average source
					filters.balance[0][0] = lSquared + stereoPoles[0] * (filters.balance[0][0] - lSquared);
					filters.balance[0][1] = rSquared + stereoPoles[0] * (filters.balance[0][1] - rSquared);
					filters.balance[1][0] = lSquared + stereoPoles[1] * (filters.balance[1][0] - lSquared);
					filters.balance[1][1] = rSquared + stereoPoles[1] * (filters.balance[1][1] - rSquared);

					// phase averaging
					// see larger comment above.
					outputPhases[i] = cos(outputPhases[i] * consts<T>::two);
					filters.phase[0] = outputPhases[i] + stereoPoles[0] * (filters.phase[0] - outputPhases[i]);
					filters.phase[1] = outputPhases[i] + stereoPoles[1] * (filters.phase[1] - outputPhases[i]);
				}


			}

			double currentEnvelope = 1.0 / (2 * std::max(std::sqrt(filterEnv[0]), std::sqrt(filterEnv[1])));

			// store calculated envelope
			if (state.envelopeMode == EnvelopeModes::RMS && state.normalizeGain)
			{
				// only update filters if this mode is on.
				filters.envelope[0] = filterEnv[0]; 
				filters.envelope[1] = filterEnv[1];

				if (std::isnormal(currentEnvelope))
				{
					envelopeGain = cpl::Math::confineTo(currentEnvelope, lowerAutoGainBounds, higherAutoGainBounds);
					cpl::GUIUtils::MainEvent(*this, [=] { if(state.isEditorOpen) setGainAsFraction(envelopeGain); });
				}
			}

		}

	bool CVectorScope::onAsyncAudio(const AudioStream & source, AudioStream::DataType ** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		switch (cpl::simd::max_vector_capacity<float>())
		{
		case 8:
		#ifdef CPL_COMPILER_SUPPORTS_AVX
				audioProcessing<cpl::Types::v8sf>(buffer, numChannels, numSamples);
				break;
		#endif
		case 4:
			audioProcessing<cpl::Types::v4sf>(buffer, numChannels, numSamples);
			break;
		default:
			audioProcessing<float>(buffer, numChannels, numSamples);
			break;
		}
		return false;
	}

	void CVectorScope::onAsyncChangedProperties(const AudioStream & source, const AudioStream::AudioStreamInfo & before)
	{
		mtFlags.audioWindowWasResized = true;
	}

};
