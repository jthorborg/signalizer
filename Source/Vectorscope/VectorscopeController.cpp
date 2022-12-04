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
 
	file:VectorscopeController.cpp
		
		Implementation of the editor for the Vectorscope
 
*************************************************************************************/

#include "Signalizer.h"
#include "../Common/SignalizerDesign.h"
#include "VectorscopeParameters.h"
#include "Vectorscope.h"

namespace Signalizer
{
	class VectorScopeController 
		: public CContentPage
	{
	public:

		VectorScopeController(VectorScopeContent& parentValue, std::shared_ptr<VectorScopeContent>&& shared)
			: parent(std::move(shared))
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
			, kwaveformColour(&parentValue.waveformColour)
			, kaxisColour(&parentValue.axisColour)
			, kbackgroundColour(&parentValue.backgroundColour)
			, kwireframeColour(&parentValue.wireframeColour)
			, kmeterColour(&parentValue.meterColour)
			, ktransform(&parentValue.transform)
			, kopMode(&parentValue.operationalMode.param)
			, kenvelopeMode(&parentValue.autoGain.param)
			, kshowLegend(&parentValue.showLegend)
			, kpresets(&valueSerializer, "vectorscope")
			, kwidgetColour(&parentValue.widgetColour)
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
			kaxisColour.bSetTitle("Axis colour");
			kbackgroundColour.bSetTitle("Backg. colour");
			kwaveformColour.bSetTitle("Audio colour");
			kwireframeColour.bSetTitle("Wireframe colour");
			kprimitiveSize.bSetTitle("Primitive size");
			kmeterColour.bSetTitle("Meter colour");
			kenvelopeSmooth.bSetTitle("Env. window");
			kstereoSmooth.bSetTitle("Stereo window");
			kwidgetColour.bSetTitle("Widget colour");

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
			kshowLegend.bSetTitle("Show legend");
			kshowLegend.setToggleable(true);

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
			kwaveformColour.bSetDescription("The main colour of the waveform audio.");
			kaxisColour.bSetDescription("The colour of the axis lines of the graph.");
			kbackgroundColour.bSetDescription("The background colour of the view.");
			kdiagnostics.bSetDescription("Toggle diagnostic information in top-left corner.");
			kwireframeColour.bSetDescription("The colour of the wireframe attached to the graph.");
			kmeterColour.bSetDescription("The colour of the stereo meters (balance and phase)");
			kprimitiveSize.bSetDescription("The size of the rendered primitives (eg. lines or points).");
			kenvelopeMode.bSetDescription("Monitors the audio stream and automatically scales the input gain such that it approaches unity intensity (envelope following).");
			kenvelopeSmooth.bSetDescription("Responsiveness (RMS window size) - or the time it takes for the envelope follower to decay.");
			kopMode.bSetDescription("Changes the presentation of the data - Lissajous is the classic XY mode on oscilloscopes, while the polar mode is a wrapped circle of the former.");
			kstereoSmooth.bSetDescription("Responsiveness (RMS window size) - or the time it takes for the stereo meters to follow.");
			kshowLegend.bSetDescription("Display a legend of the channels and assigned colours");
			kwidgetColour.bSetDescription("Colour of widgets on the screen (like legends)");
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
				}
			}

			if (auto page = addPage("Rendering", "icons/svg/brush.svg"))
			{
				if (auto section = new Signalizer::CContentPage::MatrixSection())
				{
					section->addControl(&kantiAlias, 0);
					section->addControl(&kfadeOld, 1);
					section->addControl(&kdrawLines, 2);
					section->addControl(&kshowLegend, 3);
					page->addSection(section, "Options");
				}
				if (auto section = new Signalizer::CContentPage::MatrixSection())
				{
					section->addControl(&kwaveformColour, 0);
					section->addControl(&kaxisColour, 0);
					section->addControl(&kbackgroundColour, 0);
					section->addControl(&kwireframeColour, 0);
					section->addControl(&kmeterColour, 1);
					section->addControl(&kwidgetColour, 1);
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

				if (auto section = new Signalizer::CContentPage::MatrixSection())
				{
					section->addControl(&kdiagnostics, 0);
					page->addSection(section, "Options");
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
			archive << kaxisColour;
			archive << kbackgroundColour;
			archive << kwaveformColour;
			archive << ktransform;
			archive << kwireframeColour;
			archive << kprimitiveSize;
			archive << kenvelopeMode;
			archive << kenvelopeSmooth;
			archive << kopMode;
			archive << kstereoSmooth;
			archive << kmeterColour;
			archive << kshowLegend;
			archive << kwidgetColour;
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
			builder >> kaxisColour;
			builder >> kbackgroundColour;
			builder >> kwaveformColour;
			builder >> ktransform;
			builder >> kwireframeColour;
			builder >> kprimitiveSize;
			builder >> kenvelopeMode;
			builder >> kenvelopeSmooth;
			builder >> kopMode;
			builder >> kstereoSmooth;
			builder >> kmeterColour;
			
			if (version > cpl::Version(0, 3, 6))
			{
				builder >> kshowLegend;
				builder >> kwidgetColour;
			}
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
				archive.getContent("Parameters") << *parent;
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
				builder.getContent("Parameters") >> *parent;
			}
		}

		cpl::CButton kantiAlias, kfadeOld, kdrawLines, kdiagnostics, kshowLegend;
		cpl::CValueKnobSlider kwindow, krotation, kgain, kprimitiveSize, kenvelopeSmooth, kstereoSmooth;
		cpl::CColourControl kwaveformColour, kaxisColour, kbackgroundColour, kwireframeColour, kmeterColour, kwidgetColour;
		cpl::CTransformWidget ktransform;
		cpl::CValueComboBox kopMode, kenvelopeMode;
		cpl::CPresetWidget kpresets;

		std::shared_ptr<VectorScopeContent> parent;

		SSOSurrogate<VectorScopeController>
			editorSerializer,
			valueSerializer;
	};
	

	std::unique_ptr<StateEditor> VectorScopeContent::createEditor()
	{
		return std::make_unique<VectorScopeController>(*this, shared_from_this());
	}

	std::unique_ptr<cpl::CSubView> VectorScopeContent::createView(
		std::shared_ptr<const SharedBehaviour>& globalBehaviour,
		std::shared_ptr<const ConcurrentConfig>& config,
		std::shared_ptr<AudioStream::Output>& stream
	)
	{
		return std::make_unique<VectorScope>(globalBehaviour, config, stream, shared_from_this());
	}
};
