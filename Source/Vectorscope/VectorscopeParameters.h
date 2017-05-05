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
 
	file:VectorScope.h

		Interface for the vectorscope view parameters
 
*************************************************************************************/

#ifndef SIGNALIZER_CVECTORSCOPEPARAMETERS_H
	#define SIGNALIZER_CVECTORSCOPEPARAMETERS_H

	#include "Signalizer.h"

	namespace Signalizer
	{

		class VectorScopeContent final
			: public cpl::Parameters::UserContent
			, public ProcessorState
		{
		public:

			class VectorScopeController 
				: public CContentPage
			{
			public:

				VectorScopeController(VectorScopeContent & parentValue)
					: parent(parentValue)
					, kantiAlias(&parentValue.antialias)
					, kfadeOld(&parentValue.fadeOlderPoints)
					, kdrawLines(&parentValue.interconnectSamples)
					, kdiagnostics(&parentValue.diagnostics)
					, kwindow(&parentValue.windowSize)
					, krotation(&parentValue.waveZRotation)
					, kgain(&parentValue.inputGain)
					, kprimitiveSize(&parentValue.primitiveSize)
					, kenvelopeSmooth(&parentValue.envelopeWindow)
					, kstereoSmooth(&parentValue.stereoWindow)
					, kdrawingColour(&parentValue.drawingColour)
					, kgraphColour(&parentValue.graphColour)
					, kbackgroundColour(&parentValue.backgroundColour)
					, kskeletonColour(&parentValue.skeletonColour)
					, kmeterColour(&parentValue.meterColour)
					, ktransform(&parentValue.transform)
					, kopMode(&parentValue.operationalMode.param)
					, kenvelopeMode(&parentValue.autoGain.param)
					, kpresets(&valueSerializer, "vectorscope")
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

				~VectorScopeController()
				{
					notifyDestruction();
				}

				void initControls()
				{
					kwindow.bSetTitle("Window size");
					krotation.bSetTitle("Wave Z-rotation");
					kgain.bSetTitle("Input gain");
					kgraphColour.bSetTitle("Graph colour");
					kbackgroundColour.bSetTitle("Backg. colour");
					kdrawingColour.bSetTitle("Drawing colour");
					kskeletonColour.bSetTitle("Skeleton colour");
					kprimitiveSize.bSetTitle("Primitive size");
					kmeterColour.bSetTitle("Meter colour");
					kenvelopeSmooth.bSetTitle("Env. window");
					kstereoSmooth.bSetTitle("Stereo window");

					// buttons n controls
					kantiAlias.setSingleText("Antialias");
					kantiAlias.setToggleable(true);
					kfadeOld.setSingleText("Fade older points");
					kfadeOld.setToggleable(true);
					kdrawLines.setSingleText("Interconnect samples");
					kdrawLines.setToggleable(true);
					kdiagnostics.setSingleText("Diagnostics");
					kdiagnostics.setToggleable(true);
					kenvelopeMode.bSetTitle("Auto-gain mode");

					// design
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

					if (auto page = addPage("Rendering", "icons/svg/brush.svg"))
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
					archive << krotation;
					archive << kantiAlias;
					archive << kfadeOld;
					archive << kdiagnostics;
					archive << kdrawLines;
					archive << kgraphColour;
					archive << kbackgroundColour;
					archive << kdrawingColour;
					archive << ktransform;
					archive << kskeletonColour;
					archive << kprimitiveSize;
					archive << kenvelopeMode;
					archive << kenvelopeSmooth;
					archive << kopMode;
					archive << kstereoSmooth;
					archive << kmeterColour;
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
					builder >> krotation;
					builder >> kantiAlias;
					builder >> kfadeOld;
					builder >> kdiagnostics;
					builder >> kdrawLines;
					builder >> kgraphColour;
					builder >> kbackgroundColour;
					builder >> kdrawingColour;
					builder >> ktransform;
					builder >> kskeletonColour;
					builder >> kprimitiveSize;
					builder >> kenvelopeMode;
					builder >> kenvelopeSmooth;
					builder >> kopMode;
					builder >> kstereoSmooth;
					builder >> kmeterColour;
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

				cpl::CButton kantiAlias, kfadeOld, kdrawLines, kdiagnostics;
				cpl::CValueKnobSlider kwindow, krotation, kgain, kprimitiveSize, kenvelopeSmooth, kstereoSmooth;
				cpl::CColourControl kdrawingColour, kgraphColour, kbackgroundColour, kskeletonColour, kmeterColour;
				cpl::CTransformWidget ktransform;
				cpl::CValueComboBox kopMode, kenvelopeMode;
				cpl::CPresetWidget kpresets;

				VectorScopeContent & parent;

				SSOSurrogate<VectorScopeController>
					editorSerializer,
					valueSerializer;
			};

			VectorScopeContent(std::size_t offset, bool shouldCreateShortNames, SystemView system)
				: systemView(system)
				, parameterSet("Vectorscope", "VS.", system.getProcessor(), static_cast<int>(offset))
				, audioHistoryTransformatter(system.getAudioStream(), audioHistoryTransformatter.Milliseconds)

				, dbRange(cpl::Math::dbToFraction(-120.0), cpl::Math::dbToFraction(120.0))
				, windowRange(0, 1000)
				, degreeRange(0, 360)
				, ptsRange(0.01, 10)

				, msFormatter("ms")
				, degreeFormatter("degs")
				, ptsFormatter("pts")

				, autoGain("AutoGain")
				, operationalMode("OpMode")
				, envelopeWindow("EnvWindow", windowRange, msFormatter)
				, stereoWindow("StereoWindow", windowRange, msFormatter)
				, inputGain("InputGain", dbRange, dbFormatter)
				, windowSize("WindowSize", audioHistoryTransformatter, audioHistoryTransformatter)
				, waveZRotation("WaveZ", degreeRange, degreeFormatter)
				, antialias("AntiAlias", boolRange, boolFormatter)
				, fadeOlderPoints("FadeOld", boolRange, boolFormatter)
				, interconnectSamples("Interconnect", boolRange, boolFormatter)
				, diagnostics("Diagnostics", boolRange, boolFormatter)
				, primitiveSize("PixelSize", ptsRange, ptsFormatter)


				, colourBehaviour()
				, drawingColour(colourBehaviour, "Draw.")
				, graphColour(colourBehaviour, "Graph.")
				, backgroundColour(colourBehaviour, "BackG.")
				, skeletonColour(colourBehaviour, "Skelt.")
				, meterColour(colourBehaviour, "Meter.")

				, tsfBehaviour()
				, transform(tsfBehaviour)
			{
				operationalMode.fmt.setValues({ "Lissajous", "Polar" });
				autoGain.fmt.setValues({ "None", "RMS", "Peak decay" });

				auto singleParameters = { 
					&autoGain.param, &operationalMode.param, &envelopeWindow, &stereoWindow,
					&inputGain, &windowSize, &waveZRotation, &antialias,
					&fadeOlderPoints, &interconnectSamples, &diagnostics, &primitiveSize,
				};

				for (auto sparam : singleParameters)
				{
					parameterSet.registerSingleParameter(sparam->generateUpdateRegistrator());
				}

				parameterSet.registerParameterBundle(&drawingColour, drawingColour.getBundleName());
				parameterSet.registerParameterBundle(&graphColour, graphColour.getBundleName());
				parameterSet.registerParameterBundle(&backgroundColour, backgroundColour.getBundleName());
				parameterSet.registerParameterBundle(&skeletonColour, skeletonColour.getBundleName());
				parameterSet.registerParameterBundle(&meterColour, meterColour.getBundleName());
				parameterSet.registerParameterBundle(&transform, "3D.");

				parameterSet.seal();

				postParameterInitialization();
			}

			virtual std::unique_ptr<StateEditor> createEditor() override
			{
				return std::make_unique<VectorScopeController>(*this);
			}

			virtual ParameterSet & getParameterSet() override
			{
				return parameterSet;
			}

			virtual void serialize(cpl::CSerializer::Archiver & archive, cpl::Version v) override
			{
				archive << windowSize;
				archive << inputGain;
				archive << waveZRotation;
				archive << antialias;
				archive << fadeOlderPoints;
				archive << diagnostics;
				archive << interconnectSamples;
				archive << graphColour;
				archive << backgroundColour;
				archive << drawingColour;
				archive << transform;
				archive << skeletonColour;
				archive << primitiveSize;
				archive << autoGain.param;
				archive << envelopeWindow;
				archive << operationalMode.param;
				archive << stereoWindow;
				archive << meterColour;
			}

			virtual void deserialize(cpl::CSerializer::Builder & builder, cpl::Version v) override
			{
				builder >> windowSize;
				builder >> inputGain;
				builder >> waveZRotation;
				builder >> antialias;
				builder >> fadeOlderPoints;
				builder >> diagnostics;
				builder >> interconnectSamples;
				builder >> graphColour;
				builder >> backgroundColour;
				builder >> drawingColour;
				builder >> transform;
				builder >> skeletonColour;
				builder >> primitiveSize;
				builder >> autoGain.param;
				builder >> envelopeWindow;
				builder >> operationalMode.param;
				builder >> stereoWindow;
				builder >> meterColour;
			}

			AudioHistoryTransformatter<ParameterSet::ParameterView> audioHistoryTransformatter;
			SystemView systemView;
			ParameterSet parameterSet;

			cpl::UnitFormatter<double>
				msFormatter,
				degreeFormatter,
				ptsFormatter;

			cpl::DBFormatter<double> dbFormatter;
			cpl::BooleanFormatter<double> boolFormatter;

			cpl::BooleanRange<double> boolRange;

			cpl::ExponentialRange<double> dbRange;

			cpl::LinearRange<double>
				ptsRange,
				windowRange,
				degreeRange;

			cpl::UnityRange<double> unityRange;


			cpl::ParameterValue<ParameterSet::ParameterView>
				envelopeWindow,
				stereoWindow,
				inputGain,
				windowSize,
				waveZRotation,
				antialias,
				fadeOlderPoints,
				interconnectSamples,
				diagnostics,
				primitiveSize;

			ChoiceParameter
				autoGain,
				operationalMode;

			cpl::ParameterColourValue<ParameterSet::ParameterView>::SharedBehaviour colourBehaviour;

			cpl::ParameterColourValue<ParameterSet::ParameterView>
				drawingColour,
				graphColour,
				backgroundColour,
				skeletonColour,
				meterColour;

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