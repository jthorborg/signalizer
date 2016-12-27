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
		{
		public:

			class OscilloscopeController 
				: public CContentPage
				// the preset widget controls the complete serialization logic,
				// for outside serialization we implementent specific access instead
				// to only serialize GUI settings.
				, private cpl::SafeSerializableObject
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
					, kenvelopeMode(&parentValue.autoGain)
					, kpresets(this, "oscilloscope")
					, ksubSampleInterpolationMode(&parentValue.subSampleInterpolation)
				{
					initControls();
					initUI();
				}

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


					// buttons n controls
					kantiAlias.setSingleText("Antialias");
					kantiAlias.setToggleable(true);
					kdiagnostics.setSingleText("Diagnostics");
					kdiagnostics.setToggleable(true);
					kenvelopeMode.bSetTitle("Auto-gain mode");

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
							section->addControl(&kdiagnostics, 1);
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

				void serializeEditorSettings(cpl::CSerializer::Archiver & archive, cpl::Version version) override
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
				}

				void deserializeEditorSettings(cpl::CSerializer::Archiver & builder, cpl::Version version) override
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
				}

			private:

				// entrypoints for completely storing values and settings in independant blobs (the preset widget)
				void serialize(cpl::CSerializer::Archiver & archive, cpl::Version version) override
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
				void deserialize(cpl::CSerializer::Builder & builder, cpl::Version version) override
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

				cpl::CButton kantiAlias, kdiagnostics;
				cpl::CValueKnobSlider kwindow, kgain, kprimitiveSize, kenvelopeSmooth;
				cpl::CColourControl kdrawingColour, kgraphColour, kbackgroundColour, kskeletonColour;
				cpl::CTransformWidget ktransform;
				cpl::CValueComboBox kenvelopeMode, ksubSampleInterpolationMode;
				cpl::CPresetWidget kpresets;

				OscilloscopeContent & parent;
			};

			OscilloscopeContent(std::size_t offset, bool shouldCreateShortNames, SystemView system)
				: systemView(system)
				, parameterSet("Oscilloscope", "OS.", system.getProcessor(), static_cast<int>(offset))
				, audioHistoryTransformatter(system.getAudioStream(), audioHistoryTransformatter.Miliseconds)

				, dbRange(cpl::Math::dbToFraction(-120.0), cpl::Math::dbToFraction(120.0))
				, windowRange(0, 1000)
				, degreeRange(0, 360)
				, ptsRange(0.01, 10)

				, msFormatter("ms")
				, degreeFormatter("degs")
				, ptsFormatter("pts")
				, gainModeFormatter(gainModeTransformer)
				, opModeFormatter(opModeTransformer)
				, subSampleFormatter(subSampleTransformer)

				, autoGain("AutoGain", gainModeTransformer, gainModeFormatter)
				, envelopeWindow("EnvWindow", windowRange, msFormatter)
				, inputGain("InputGain", dbRange, dbFormatter)
				, windowSize("WindowSize", audioHistoryTransformatter, audioHistoryTransformatter)
				, antialias("AntiAlias", boolRange, boolFormatter)
				, diagnostics("Diagnostics", boolRange, boolFormatter)
				, primitiveSize("PixelSize", ptsRange, ptsFormatter)
				, subSampleInterpolation("SampleIntp", subSampleTransformer, subSampleFormatter)

				, colourBehavior()
				, drawingColour(colourBehavior, "Draw.")
				, graphColour(colourBehavior, "Graph.")
				, backgroundColour(colourBehavior, "BackG.")
				, skeletonColour(colourBehavior, "Skelt.")

				, tsfBehaviour()
				, transform(tsfBehaviour)
			{
				opModeFormatter.setValues({ "Lissajous", "Polar" });
				gainModeFormatter.setValues({ "None", "RMS", "Peak decay" });
				subSampleFormatter.setValues({ "None", "Rectangular", "Linear", "Lanczos 5" });


				auto singleParameters = { 
					&autoGain, &envelopeWindow,
					&inputGain, &windowSize, &subSampleInterpolation, &antialias,
					&diagnostics, &primitiveSize,
				};

				for (auto sparam : singleParameters)
				{
					parameterSet.registerSingleParameter(sparam->generateUpdateRegistrator());
				}

				parameterSet.registerParameterBundle(&drawingColour, drawingColour.getBundleName());
				parameterSet.registerParameterBundle(&graphColour, graphColour.getBundleName());
				parameterSet.registerParameterBundle(&backgroundColour, backgroundColour.getBundleName());
				parameterSet.registerParameterBundle(&skeletonColour, skeletonColour.getBundleName());
				parameterSet.registerParameterBundle(&transform, "3D.");

				parameterSet.seal();

				postParameterInitialization();
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
				archive << autoGain;
				archive << envelopeWindow;
				archive << subSampleInterpolation;
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
				builder >> autoGain;
				builder >> envelopeWindow;
				builder >> subSampleInterpolation;
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
			cpl::ChoiceFormatter<double>
				gainModeFormatter,
				opModeFormatter,
				subSampleFormatter;

			cpl::ChoiceTransformer<double>
				gainModeTransformer,
				opModeTransformer,
				subSampleTransformer;

			cpl::BooleanRange<double> boolRange;

			cpl::ExponentialRange<double> dbRange;

			cpl::LinearRange<double>
				ptsRange,
				windowRange,
				degreeRange;

			cpl::UnityRange<double> unityRange;


			cpl::ParameterValue<ParameterSet::ParameterView>
				autoGain,
				envelopeWindow,
				inputGain,
				windowSize,
				antialias,
				diagnostics,
				primitiveSize,
				subSampleInterpolation;

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