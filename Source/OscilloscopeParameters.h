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
 
	file:COscilloscope.h

		Interface for the oscilloscope parameters
 
*************************************************************************************/

#ifndef SIGNALIZER_COSCILLOSCOPEPARAMETERS_H
	#define SIGNALIZER_COSCILLOSCOPEPARAMETERS_H

	#include "CommonSignalizer.h"
	#include "SignalizerDesign.h"

	namespace Signalizer
	{

		class OscilloscopeContent final
			: public cpl::Parameters::UserContent
			, public ProcessorState
			, public ParameterSet::UIListener
		{
		public:

			static std::size_t const LookaheadSize = 8192;
			static std::size_t const InterpolationKernelSize = 11;

			enum class TriggeringMode
			{
				None,
				Spectral,
				//ZeroCrossing,
				end
			};

			enum class TimeMode
			{
				Time, Cycles, Beats
			};

			template<typename ParameterView>
			class WindowSizeTransformatter : public AudioHistoryTransformatter<ParameterView>
			{
			public:
				WindowSizeTransformatter(AudioStream & audioStream, std::size_t auxLookahead, Mode mode = Milliseconds)
					: AudioHistoryTransformatter(audioStream, mode), lookahead(auxLookahead), timeMode(TimeMode::Time)
				{

				}

				void setTimeModeFromUI(TimeMode newMode)
				{
					timeMode = newMode;
					param->updateFromUINormalized(param->getValueNormalized());
				}

			private:

				virtual void onAsyncChangedProperties(const Stream & changedSource, const typename Stream::AudioStreamInfo & before) override
				{
					// TODO: what if oldCapacity == 0?
					const auto oldFraction = param->getValueNormalized();
					auto oldCapacity = lastCapacity.load(std::memory_order_relaxed);
					auto beforeCapacity = before.audioHistoryCapacity.load(std::memory_order_acquire);
					if (oldCapacity == 0)
						oldCapacity = beforeCapacity;

					const auto newCapacity = changedSource.getInfo().audioHistoryCapacity.load(std::memory_order_relaxed);

					if (newCapacity > 0)
						lastCapacity.store(newCapacity, std::memory_order_relaxed);

					if (oldCapacity == 0 || newCapacity == 0)
					{
						param->updateFromProcessorNormalized(oldFraction, cpl::Parameters::UpdateFlags::All & ~cpl::Parameters::UpdateFlags::RealTimeSubSystem);
					}
					else
					{
						const auto sampleSizeBefore = oldCapacity * oldFraction;
						const auto newFraction = sampleSizeBefore / newCapacity;
						if (oldFraction != newFraction || beforeCapacity == 0)
							param->updateFromProcessorNormalized(newFraction, cpl::Parameters::UpdateFlags::All & ~cpl::Parameters::UpdateFlags::RealTimeSubSystem);
					}
				}

				virtual bool format(const ValueType & val, std::string & buf) override
				{
					char buffer[100];
					switch (timeMode)
					{
						case TimeMode::Cycles:
						{
							sprintf_s(buffer, u8"%.2f (%.2f r)", val, cpl::simd::consts<ValueType>::tau * val);
							buf = buffer;
							return true;
						}
						case TimeMode::Beats:
						{
							sprintf_s(buffer, "1/%.0f", val);
							buf = buffer;
							return true;
						}
						case TimeMode::Time: return AudioHistoryTransformatter<ParameterView>::format(val, buf);
					}
				}

				virtual bool interpret(const std::string & buf, ValueType & val) override
				{
					ValueType collectedValue;

					if (cpl::lexicalConversion(buf, collectedValue))
					{
						bool notSamples = true;
						if (timeMode != TimeMode::Time)
						{
							if (timeMode == TimeMode::Cycles)
							{
								if (buf.find("r") != std::string::npos)
								{
									collectedValue /= cpl::simd::consts<ValueType>::tau;
								}
							}
							else if (timeMode == TimeMode::Beats)
							{
								if (buf.find("bars") != std::string::npos)
								{
									collectedValue /= 4;
								}
							}
							val = collectedValue;
							return true;
						}
						else if (buf.find("s") != std::string::npos && (notSamples = buf.find("smps") == std::string::npos))
						{
							if (buf.find("ms") != std::string::npos)
							{
								collectedValue /= 1000;
							}
							collectedValue *= stream.getInfo().sampleRate.load(std::memory_order_relaxed);
						}
						else
						{
							// assume value is in miliseconds
							if (m == Milliseconds && notSamples)
							{
								collectedValue /= 1000;
								collectedValue *= stream.getInfo().sampleRate.load(std::memory_order_relaxed);
							}
						}

						val = collectedValue;
						return true;

					}

					return false;

				}

				virtual ValueType transform(ValueType val) const noexcept override
				{
					switch (timeMode)
					{
						case TimeMode::Cycles:
						{
							return cpl::Math::UnityScale::exp<ValueType>(val, 1, 32);
						}
						case TimeMode::Beats:
						{
							return cpl::Math::nextPow2Inc(cpl::Math::round<std::size_t>(cpl::Math::UnityScale::exp<ValueType>(1 - val, 1, 128)));
						}
						case TimeMode::Time:
						{
							const auto minExponential = 100;
							const auto capacity = stream.getAudioHistoryCapacity();

							const auto top = capacity;
							const auto expSamples = cpl::Math::UnityScale::exp<ValueType>(val, minExponential, top);
							const auto rescaled = cpl::Math::UnityScale::linear<ValueType>(cpl::Math::UnityScale::Inv::linear<ValueType>(expSamples, minExponential, top), 2, top);
							return rescaled;
						}
					}

				}


				virtual ValueType normalize(ValueType val) const noexcept override
				{
					switch (timeMode)
					{
						case TimeMode::Cycles:
						{
							return cpl::Math::UnityScale::Inv::exp<ValueType>(val, 1, 128);
						}
						case TimeMode::Beats:
						{
							return cpl::Math::UnityScale::Inv::exp<ValueType>(1 - val, 1, 32);
						}
						case TimeMode::Time:
						{
							const auto minExponential = 100;
							const auto capacity = stream.getAudioHistoryCapacity();
							const auto top = capacity;
							const auto linear = cpl::Math::UnityScale::Inv::linear<ValueType>(val, 1, top);
							const auto expSamples = cpl::Math::UnityScale::linear<ValueType>(linear, minExponential, top);

							const auto normalized = cpl::Math::UnityScale::Inv::exp<ValueType>(expSamples, minExponential, top);
							return normalized;
						}
					}
				}

				std::atomic<Scaling> scale;
				std::size_t lookahead;
				TimeMode timeMode;
			};

			class OscilloscopeController 
				: public CContentPage
			{
			public:

				OscilloscopeController(OscilloscopeContent & parentValue)
					: parent(parentValue)
					, kantiAlias(&parentValue.antialias)
					, kdiagnostics(&parentValue.diagnostics)
					, kwindow(&parentValue.windowSize)
					, kgain(&parentValue.inputGain)
					, kprimitiveSize(&parentValue.primitiveSize)
					, kenvelopeSmooth(&parentValue.envelopeWindow)
					, kdrawingColour(&parentValue.drawingColour)
					, kgraphColour(&parentValue.graphColour)
					, kbackgroundColour(&parentValue.backgroundColour)
					, kskeletonColour(&parentValue.skeletonColour)
					, ktransform(&parentValue.transform)
					, kenvelopeMode(&parentValue.autoGain.param)
					, kpresets(&valueSerializer, "oscilloscope")
					, ksubSampleInterpolationMode(&parentValue.subSampleInterpolation.param)
					, kpctForDivision(&parentValue.pctForDivision)
					, kchannelConfiguration(&parentValue.channelConfiguration.param)
					, ktriggerPhaseOffset(&parentValue.triggerPhaseOffset)
					, ktriggerMode(&parentValue.triggerMode.param)
					, ktimeMode(&parentValue.timeMode.param)
					, kdotSamples(&parentValue.dotSamples)
					, kcolourSmoothingTime(&parentValue.colourSmoothing)

					, editorSerializer(
						*this, 
						[](auto & oc, auto & se, auto version) { oc.serializeEditorSettings(se, version); },
						[](auto & oc, auto & se, auto version) { oc.deserializeEditorSettings(se, version); }
					)
					, valueSerializer(
						*this,
						[](auto & oc, auto & se, auto version) { oc.serializeAll(se, version); },
						[](auto & oc, auto & se, auto version) { oc.deserializeAll(se, version); }
					)
				{
					initControls();
					initUI();
				}

				cpl::SafeSerializableObject & getEditorSO() override { return editorSerializer; }

				~OscilloscopeController()
				{
					notifyDestruction();
				}

				void initControls()
				{
					kwindow.bSetTitle("Window size");
					kgain.bSetTitle("Input gain");
					kgraphColour.bSetTitle("Graph colour");
					kbackgroundColour.bSetTitle("Backg. colour");
					kdrawingColour.bSetTitle("Drawing colour");
					kskeletonColour.bSetTitle("Skeleton colour");
					kprimitiveSize.bSetTitle("Primitive size");
					kenvelopeSmooth.bSetTitle("Env. window");
					ksubSampleInterpolationMode.bSetTitle("Sample interpolation");
					kpctForDivision.bSetTitle("Grid div. space");
					kchannelConfiguration.bSetTitle("Channel conf.");
					ktriggerMode.bSetTitle("Trigger mode");
					ktriggerPhaseOffset.bSetTitle("Trigger phase");
					ktimeMode.bSetTitle("Time mode");

					// buttons n controls
					kantiAlias.setSingleText("Antialias");
					kantiAlias.setToggleable(true);
					kdiagnostics.setSingleText("Diagnostics");
					kdiagnostics.setToggleable(true);
					kenvelopeMode.bSetTitle("Auto-gain mode");
					kdotSamples.bSetTitle("Dot samples");
					kdotSamples.setToggleable(true);

					// descriptions.
					kwindow.bSetDescription("The size of the displayed time window.");
					kgain.bSetDescription("How much the input (x,y) is scaled (or the input gain)" \
						" - additional transform that only affects the waveform, and not the graph");
					kantiAlias.bSetDescription("Antialiases rendering (if set - see global settings for amount). May slow down rendering.");
					kdrawingColour.bSetDescription("The main colour to paint with.");
					kgraphColour.bSetDescription("The colour of the graph.");
					kbackgroundColour.bSetDescription("The background colour of the view.");
					kdiagnostics.bSetDescription("Toggle diagnostic information in top-left corner.");
					kskeletonColour.bSetDescription("The colour of the box skeleton indicating the OpenGL camera clip box.");
					kprimitiveSize.bSetDescription("The size of the rendered primitives (eg. lines or points).");
					kenvelopeMode.bSetDescription("Monitors the audio stream and automatically scales the input gain such that it approaches unity intensity (envelope following).");
					kenvelopeSmooth.bSetDescription("Responsiveness (RMS window size) - or the time it takes for the envelope follower to decay.");
					ksubSampleInterpolationMode.bSetDescription("Controls how point samples are interpolated to wave forms");
					kpctForDivision.bSetDescription("The minimum amount of free space that triggers a recursed frequency grid division; smaller values draw more frequency divisions.");
					kchannelConfiguration.bSetDescription("Select how the audio channels are interpreted.");
					ktriggerMode.bSetDescription("Select a mode for triggering waveforms - i.e. syncing them to the grid");
					ktriggerPhaseOffset.bSetDescription("A custom +/- full-circle offset for the phase on triggering");
					ktimeMode.bSetDescription("Specifies the working units of the time display");
					kdotSamples.bSetDescription("Marks sample positions when drawing subsampled interpolated lines");
				}

				void initUI()
				{
					if (auto page = addPage("Settings", "icons/svg/gear.svg"))
					{
						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&ktransform, 0);
							page->addSection(section, "Transform");
						}
						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&kenvelopeMode, 0);
							section->addControl(&kenvelopeSmooth, 1);
							
							section->addControl(&kchannelConfiguration, 0);

							section->addControl(&kgain, 1);

							section->addControl(&kwindow, 0);
							section->addControl(&ktimeMode, 1);
							section->addControl(&kpctForDivision, 0);
							section->addControl(&kcolourSmoothingTime, 1);

							page->addSection(section, "Utility");
						}
						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&ktriggerMode, 0);
							section->addControl(&ktriggerPhaseOffset, 1);

							page->addSection(section, "Triggering");
						}
					}

					if (auto page = addPage("Rendering", "icons/svg/brush.svg"))
					{
						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&kantiAlias, 0);
							section->addControl(&kdiagnostics, 1);
							section->addControl(&kdotSamples, 2);
							page->addSection(section, "Options");
						}
						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{							
							section->addControl(&kprimitiveSize, 0);
							section->addControl(&ksubSampleInterpolationMode, 1);

							section->addControl(&kdrawingColour, 0);
							section->addControl(&kgraphColour, 1);
							section->addControl(&kbackgroundColour, 0);
							section->addControl(&kskeletonColour, 1);

							page->addSection(section, "Look");
						}
					}

					if (auto page = addPage("Utility", "icons/svg/wrench.svg"))
					{
						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&kpresets, 0);
							page->addSection(section, "Presets");
						}
					}
				}

			private:

				void serializeEditorSettings(cpl::CSerializer::Archiver & archive, cpl::Version version)
				{
					archive << kwindow;
					archive << kgain;
					archive << kantiAlias;
					archive << kdiagnostics;
					archive << kgraphColour;
					archive << kbackgroundColour;
					archive << kdrawingColour;
					archive << ktransform;
					archive << kskeletonColour;
					archive << kprimitiveSize;
					archive << kenvelopeMode;
					archive << kenvelopeSmooth;
					archive << ksubSampleInterpolationMode;
					archive << kpctForDivision;
					archive << kchannelConfiguration;
					archive << kpctForDivision;
					archive << ktriggerPhaseOffset;
					archive << ktriggerMode;
					archive << ktimeMode;
					archive << kdotSamples;
					archive << kcolourSmoothingTime;
				}

				void deserializeEditorSettings(cpl::CSerializer::Archiver & builder, cpl::Version version)
				{
					// in general, controls should never restore values. However, older versions
					// of Signalizer does exactly this, so to keep backwards-compatibility, we 
					// can obtain the preset values through this.
					cpl::Serialization::ScopedModifier m(cpl::CSerializer::Modifiers::RestoreValue, version < cpl::Version(0, 2, 8));
					builder << m;

					builder >> kwindow;
					builder >> kgain;
					builder >> kantiAlias;
					builder >> kdiagnostics;
					builder >> kgraphColour;
					builder >> kbackgroundColour;
					builder >> kdrawingColour;
					builder >> ktransform;
					builder >> kskeletonColour;
					builder >> kprimitiveSize;
					builder >> kenvelopeMode;
					builder >> kenvelopeSmooth;
					builder >> ksubSampleInterpolationMode;
					builder >> kpctForDivision;
					builder >> kchannelConfiguration;
					builder >> kpctForDivision;
					builder >> ktriggerPhaseOffset;
					builder >> ktriggerMode;
					builder >> ktimeMode;
					builder >> kdotSamples;
					builder >> kcolourSmoothingTime;
				}


				// entrypoints for completely storing values and settings in independant blobs (the preset widget)
				void serializeAll(cpl::CSerializer::Archiver & archive, cpl::Version version)
				{
					if (version < cpl::Version(0, 2, 8))
					{
						// presets from < 0.2.8 only store editor settings with values
						serializeEditorSettings(archive, version);
					}
					else
					{
						// store parameter and editor settings separately
						serializeEditorSettings(archive.getContent("Editor"), version);
						archive.getContent("Parameters") << parent;
					}

				}

				// entrypoints for completely storing values and settings in independant blobs (the preset widget)
				void deserializeAll(cpl::CSerializer::Builder & builder, cpl::Version version)
				{
					if (version < cpl::Version(0, 2, 8))
					{
						// presets from < 0.2.8 only store editor settings with values
						deserializeEditorSettings(builder, version);
					}
					else
					{
						// store parameter and editor settings separately
						deserializeEditorSettings(builder.getContent("Editor"), version);
						builder.getContent("Parameters") >> parent;
					}
				}

				cpl::CButton kantiAlias, kdiagnostics, kdotSamples;
				cpl::CValueKnobSlider kwindow, kgain, kprimitiveSize, kenvelopeSmooth, kpctForDivision, ktriggerPhaseOffset, kcolourSmoothingTime;
				cpl::CColourControl kdrawingColour, kgraphColour, kbackgroundColour, kskeletonColour;
				cpl::CTransformWidget ktransform;
				cpl::CValueComboBox kenvelopeMode, ksubSampleInterpolationMode, kchannelConfiguration, ktriggerMode, ktimeMode;
				cpl::CPresetWidget kpresets;

				OscilloscopeContent & parent;

				SSOSurrogate<OscilloscopeController>
					editorSerializer,
					valueSerializer;
			};

			OscilloscopeContent(std::size_t offset, bool shouldCreateShortNames, SystemView system)
				: systemView(system)
				, parameterSet("Oscilloscope", "OS.", system.getProcessor(), static_cast<int>(offset))
				, audioHistoryTransformatter(system.getAudioStream(), LookaheadSize, audioHistoryTransformatter.Milliseconds)

				, dbRange(cpl::Math::dbToFraction(-120.0), cpl::Math::dbToFraction(120.0))
				, windowRange(0, 1000)
				, degreeRange(0, 360)
				, ptsRange(0.01, 10)
				, phaseRange(-180, 180)
				, reverseUnitRange(1, 0)
				, colourSmoothRange(1, 1000)

				, msFormatter("ms")
				, degreeFormatter("degs")
				, ptsFormatter("pts")

				, autoGain("AutoGain")
				, envelopeWindow("EnvWindow", windowRange, msFormatter)
				, inputGain("InputGain", dbRange, dbFormatter)
				, windowSize("WindowSize", audioHistoryTransformatter, audioHistoryTransformatter)
				, antialias("AntiAlias", boolRange, boolFormatter)
				, diagnostics("Diagnostics", boolRange, boolFormatter)
				, primitiveSize("PixelSize", ptsRange, ptsFormatter)
				, subSampleInterpolation("SampleIntp")
				, pctForDivision("PctDiv", unityRange, pctFormatter)
				, channelConfiguration("ChConf")
				, triggerPhaseOffset("TrgPhase", phaseRange, degreeFormatter)
				, triggerMode("TrgMode")
				, timeMode("TimeMode")
				, dotSamples("DotSmps", boolRange, boolFormatter)
				, colourSmoothing("ColSmth", colourSmoothRange, msFormatter)

				, colourBehavior()
				, drawingColour(colourBehavior, "Draw.")
				, graphColour(colourBehavior, "Graph.")
				, backgroundColour(colourBehavior, "BackG.")
				, skeletonColour(colourBehavior, "Skelt.")

				, tsfBehaviour()
				, transform(tsfBehaviour)

			{
				viewOffsets.emplace_back("ViewLeft", unityRange, basicFormatter);
				viewOffsets.emplace_back("ViewTop", unityRange, basicFormatter);
				viewOffsets.emplace_back("ViewRight", reverseUnitRange, basicFormatter);
				viewOffsets.emplace_back("ViewBottom", reverseUnitRange, basicFormatter);

				autoGain.fmt.setValues({ "None", "RMS", "Peak decay" });
				subSampleInterpolation.fmt.setValues({ "None", "Rectangular", "Linear", "Lanczos 5" });
				channelConfiguration.fmt.setValues({ "Left", "Right", "Mid/Merge", "Side", "Separate", "Mid+Side"});
				triggerMode.fmt.setValues({ "None", "Spectral" /*, "Zero-crossings" */});
				timeMode.fmt.setValues({ "Time", "Cycles", "Beats" });
				// order matters
				auto singleParameters = { 
					&autoGain.param, 
					&envelopeWindow,
					&inputGain, 
					&windowSize, 
					&antialias,
					&diagnostics, 
					&primitiveSize,
					&subSampleInterpolation.param,
					&channelConfiguration.param,
					&pctForDivision,
					&triggerPhaseOffset,
					&triggerMode.param,
					&timeMode.param,
					&dotSamples,
					&colourSmoothing
				};

				for (auto sparam : singleParameters)
				{
					parameterSet.registerSingleParameter(sparam->generateUpdateRegistrator());
				}

				for (auto & v : viewOffsets)
				{
					parameterSet.registerSingleParameter(v.generateUpdateRegistrator());
				}

				parameterSet.registerParameterBundle(&drawingColour, drawingColour.getBundleName());
				parameterSet.registerParameterBundle(&graphColour, graphColour.getBundleName());
				parameterSet.registerParameterBundle(&backgroundColour, backgroundColour.getBundleName());
				parameterSet.registerParameterBundle(&skeletonColour, skeletonColour.getBundleName());
				parameterSet.registerParameterBundle(&transform, "3D.");

				parameterSet.seal();

				postParameterInitialization();

				timeMode.param.getParameterView().addListener(this);
			}

			~OscilloscopeContent()
			{
				timeMode.param.getParameterView().removeListener(this);
			}

			void parameterChangedUI(cpl::Parameters::Handle localHandle, cpl::Parameters::Handle globalHandle, ParameterSet::ParameterView * parameter)
			{
				if (parameter == &timeMode.param.getParameterView())
				{
					audioHistoryTransformatter.setTimeModeFromUI(timeMode.param.getAsTEnum<TimeMode>());
				}
			}

			virtual std::unique_ptr<StateEditor> createEditor() override
			{
				return std::make_unique<OscilloscopeController>(*this);
			}

			virtual ParameterSet & getParameterSet() override
			{
				return parameterSet;
			}

			virtual void serialize(cpl::CSerializer::Archiver & archive, cpl::Version v) override
			{
				archive << windowSize;
				archive << inputGain;
				archive << antialias;
				archive << diagnostics;
				archive << graphColour;
				archive << backgroundColour;
				archive << drawingColour;
				archive << transform;
				archive << skeletonColour;
				archive << primitiveSize;
				archive << autoGain.param;
				archive << envelopeWindow;
				archive << subSampleInterpolation.param;
				archive << channelConfiguration.param;
				archive << pctForDivision;
				archive << triggerPhaseOffset;
				archive << triggerMode.param;
				archive << timeMode.param;
				for (auto && v : viewOffsets)
				{
					archive << v;
				}
				archive << dotSamples;
				archive << colourSmoothing;
			}

			virtual void deserialize(cpl::CSerializer::Builder & builder, cpl::Version v) override
			{
				builder >> windowSize;
				builder >> inputGain;
				builder >> antialias;
				builder >> diagnostics;
				builder >> graphColour;
				builder >> backgroundColour;
				builder >> drawingColour;
				builder >> transform;
				builder >> skeletonColour;
				builder >> primitiveSize;
				builder >> autoGain.param;
				builder >> envelopeWindow;
				builder >> subSampleInterpolation.param;
				builder >> channelConfiguration.param;
				builder >> pctForDivision;
				builder >> triggerPhaseOffset;
				builder >> triggerMode.param;
				builder >> timeMode.param;
				for (auto & v : viewOffsets)
				{
					builder >> v;
				}
				builder >> dotSamples;
				builder >> colourSmoothing;
			}

			WindowSizeTransformatter<ParameterSet::ParameterView> audioHistoryTransformatter;
			SystemView systemView;
			ParameterSet parameterSet;

			cpl::UnitFormatter<double>
				msFormatter,
				degreeFormatter,
				ptsFormatter;

			cpl::PercentageFormatter<double>
				pctFormatter;

			cpl::DBFormatter<double> dbFormatter;
			cpl::BooleanFormatter<double> boolFormatter;

			cpl::BasicFormatter<double> basicFormatter;
			cpl::BooleanRange<double> boolRange;

			cpl::ExponentialRange<double> dbRange, colourSmoothRange;

			cpl::LinearRange<double>
				ptsRange,
				windowRange,
				degreeRange,
				phaseRange,
				reverseUnitRange;

			cpl::UnityRange<double> unityRange;
			
			enum ViewOffsets
			{
				Left, Top, Right, Bottom, end
			};

			cpl::ParameterValue<ParameterSet::ParameterView>
				envelopeWindow,
				inputGain,
				windowSize,
				antialias,
				diagnostics,
				primitiveSize,
				pctForDivision,
				triggerPhaseOffset,
				dotSamples,
				colourSmoothing;

			std::vector<cpl::ParameterValue<ParameterSet::ParameterView>> viewOffsets;

			ChoiceParameter
				autoGain,
				channelConfiguration,
				subSampleInterpolation,
				triggerMode,
				timeMode;


			cpl::ParameterColourValue<ParameterSet::ParameterView>::SharedBehaviour colourBehavior;

			cpl::ParameterColourValue<ParameterSet::ParameterView>
				drawingColour,
				graphColour,
				backgroundColour,
				skeletonColour;

			cpl::ParameterTransformValue<ParameterSet::ParameterView>::SharedBehaviour<ParameterSet::ParameterView::ValueType> tsfBehaviour;

			cpl::ParameterTransformValue<ParameterSet::ParameterView> transform;

		private:

			void postParameterInitialization()
			{
				audioHistoryTransformatter.initialize(windowSize.getParameterView());
			}
		};
	
	};

#endif