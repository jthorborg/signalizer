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

	file:OscilloscopeController.cpp

		Editor for the oscilloscope content

*************************************************************************************/

#include "Signalizer.h"
#include "../Common/SignalizerDesign.h"
#include "OscilloscopeParameters.h"

namespace Signalizer
{
	class OscilloscopeController
		: public CContentPage
		, public ParameterSet::UIListener
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
			, kprimaryColour(&parentValue.getColour(0))
			, ksecondaryColour(&parentValue.getColour(1))
			, kgraphColour(&parentValue.graphColour)
			, kbackgroundColour(&parentValue.backgroundColour)
			, klowColour(&parentValue.lowColour)
			, kmidColour(&parentValue.midColour)
			, khighColour(&parentValue.highColour)
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
			, ktriggerOnCustomFrequency(&parentValue.triggerOnCustomFrequency)
			, kcustomFrequency(&parentValue.customTriggerFrequency)
			, koverlayChannels(&parentValue.overlayChannels)
			, kchannelColouring(&parentValue.channelColouring.param)
			, kcolourSmoothingTime(&parentValue.colourSmoothing)
			, kcursorTracker(&parentValue.cursorTracker)
			, ktrackerColour(&parentValue.trackerColour)
			, kfreqColourBlend(&parentValue.frequencyColouringBlend)
			, ktriggerHysteresis(&parentValue.triggerHysteresis)
			, ktriggerThreshold(&parentValue.triggerThreshold)
			, kshowLegend(&parentValue.showLegend)
			, ktriggerChannel(&parentValue.triggeringChannel)

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

			parent.getParameterSet().addUIListener(parent.timeMode.param.getParameterView().getHandle(), this);
			parent.getParameterSet().addUIListener(parent.triggerMode.param.getParameterView().getHandle(), this);

			enforceTriggeringModeCompability();
		}


		void parameterChangedUI(cpl::Parameters::Handle localHandle, cpl::Parameters::Handle globalHandle, ParameterSet::ParameterView * parameter) override
		{
			if(parameter == &parent.timeMode.param.getParameterView())
				enforceTriggeringModeCompability();
		}

		void enforceTriggeringModeCompability()
		{
			if (parent.timeMode.param.getAsTEnum<OscilloscopeContent::TimeMode>() == OscilloscopeContent::TimeMode::Cycles)
			{
				for (int i = 0; i < parent.triggerMode.tsf.getQuantization(); ++i)
				{
					ktriggerMode.setEnabledStateFor(i, i == cpl::enum_cast<int>(OscilloscopeContent::TimeMode::Cycles));
				}

				ktriggerMode.getValueReference().setTransformedValue(cpl::enum_cast<double>(OscilloscopeContent::TimeMode::Cycles));
			}
			else
			{
				for (int i = 0; i < parent.triggerMode.tsf.getQuantization(); ++i)
				{
					ktriggerMode.setEnabledStateFor(i, true);
				}
			}

		}

		cpl::SafeSerializableObject & getEditorSO() override { return editorSerializer; }

		~OscilloscopeController()
		{
			parent.getParameterSet().removeUIListener(parent.timeMode.param.getParameterView().getHandle(), this);
			parent.getParameterSet().removeUIListener(parent.triggerMode.param.getParameterView().getHandle(), this);
			notifyDestruction();
		}

		void initControls()
		{
			kwindow.bSetTitle("Window size");
			kgain.bSetTitle("Input gain");
			kgraphColour.bSetTitle("Graph colour");
			kbackgroundColour.bSetTitle("Backg. colour");
			kprimaryColour.bSetTitle("Primary colour");
			ksecondaryColour.bSetTitle("Secondary colour");
			klowColour.bSetTitle("Low band colour");
			kmidColour.bSetTitle("Mid band colour");
			khighColour.bSetTitle("High band colour");
			kprimitiveSize.bSetTitle("Primitive size");
			kenvelopeSmooth.bSetTitle("Env. window");
			ksubSampleInterpolationMode.bSetTitle("Sample interpolation");
			kpctForDivision.bSetTitle("Grid div. space");
			kchannelConfiguration.bSetTitle("Channel conf.");
			ktriggerMode.bSetTitle("Trigger mode");
			ktriggerPhaseOffset.bSetTitle("Trigger phase");
			ktimeMode.bSetTitle("Time mode");
			kcustomFrequency.bSetTitle("Custom trigger");
			kchannelColouring.bSetTitle("Channel colouring");
			kcolourSmoothingTime.bSetTitle("Colour smoothing");
			ktrackerColour.bSetTitle("Tracker colour");
			kfreqColourBlend.bSetTitle("Colour blend");
			ktriggerHysteresis.bSetTitle("Hysteresis");
			ktriggerThreshold.bSetTitle("Trigger thrshld");
			ktriggerChannel.bSetTitle("Trigger channel");
			// buttons n controls

			kantiAlias.setSingleText("Antialias");
			kantiAlias.setToggleable(true);
			kdiagnostics.setSingleText("Diagnostics");
			kdiagnostics.setToggleable(true);
			kenvelopeMode.bSetTitle("Auto-gain mode");
			kdotSamples.bSetTitle("Dot samples");
			kdotSamples.setToggleable(true);
			ktriggerOnCustomFrequency.bSetTitle("Trigger on custom");
			ktriggerOnCustomFrequency.setToggleable(true);
			koverlayChannels.bSetTitle("Overlay channels");
			koverlayChannels.setToggleable(true);
			kcursorTracker.bSetTitle("Cursor tracker");
			kcursorTracker.setToggleable(true);
			kshowLegend.bSetTitle("Show legend");
			kshowLegend.setToggleable(true);

			// descriptions.
			kwindow.bSetDescription("The size of the displayed time window.");
			kgain.bSetDescription("How much the input (x,y) is scaled (or the input gain)" \
				" - additional transform that only affects the waveform, and not the graph");
			kantiAlias.bSetDescription("Antialiases rendering (if set - see global settings for amount). May slow down rendering.");
			kprimaryColour.bSetDescription("Colour for the first channel");
			ksecondaryColour.bSetDescription("Colour for the second channel");
			kgraphColour.bSetDescription("The colour of the graph.");
			kbackgroundColour.bSetDescription("The background colour of the view.");
			kdiagnostics.bSetDescription("Toggle diagnostic information in top-left corner.");
			kprimitiveSize.bSetDescription("The size of the rendered primitives (eg. lines or points).");
			kenvelopeMode.bSetDescription("Monitors the audio stream and automatically scales the input gain such that it approaches unity intensity (envelope following).");
			kenvelopeSmooth.bSetDescription("Responsiveness (RMS window size) - or the time it takes for the envelope follower to decay.");
			ksubSampleInterpolationMode.bSetDescription("Controls how point samples are interpolated to wave forms");
			kpctForDivision.bSetDescription("The minimum amount of free space that triggers a recursed frequency grid division; smaller values draw more frequency divisions.");
			kchannelConfiguration.bSetDescription("Select how the audio channels are interpreted.");
			ktriggerMode.bSetDescription("Select a mode for triggering waveforms - i.e. syncing to frequency content, time or transition information");
			ktriggerPhaseOffset.bSetDescription("A custom +/- full-circle offset for the phase on triggering");
			ktimeMode.bSetDescription("Specifies the working units of the time display");
			kdotSamples.bSetDescription("Marks sample positions when drawing subsampled interpolated lines");
			ktriggerOnCustomFrequency.bSetDescription("If toggled, the oscilloscope will trigger on the specified custom frequency instead of autodetecting a fundamental");
			kcustomFrequency.bSetDescription("Specifies a custom frequency to trigger on - input units can be notes (like a#2), hz, radians (rads), beats, samples (smps) or period (ms)");
			klowColour.bSetDescription("Colour for the low frequency band");
			kmidColour.bSetDescription("Colour for the mid frequency band");
			khighColour.bSetDescription("Colour for the high frequency band");
			kchannelColouring.bSetDescription("Method for colouring of channels, static equals just the drawing colour while spectral paints with separate colours for each frequency band");
			kfreqColourBlend.bSetDescription("Blending between the static and spectral colouring of channels");
			koverlayChannels.bSetDescription("Toggle to paint multiple channels on top of each other, otherwise they are painted in separate views");
			kcolourSmoothingTime.bSetDescription("Smooths the colour variation over the period of time");
			kcursorTracker.bSetDescription("Enable to create a tracker at the cursor displaying (x,y) values");
			ktrackerColour.bSetDescription("Colour of the cursor tracker");
			ktriggerHysteresis.bSetDescription("The hysteresis of the triggering function defines an opaque measure of how resistant the trigger is to change");
			ktriggerThreshold.bSetDescription("The triggering function will not consider any candidates below the threshold");
			kshowLegend.bSetDescription("Display a legend of the signals/colours");
			ktriggerChannel.bSetDescription("Adjust the channel used for triggering when in separate channel mode");
		}

		void initUI()
		{
			if (auto page = addPage("Settings", "icons/svg/gear.svg"))
			{
				/* if (auto section = new Signalizer::CContentPage::MatrixSection())
				{
					section->addControl(&ktransform, 0);
					page->addSection(section, "Transform");
				} */

				if (auto section = new Signalizer::CContentPage::MatrixSection())
				{
					section->addControl(&koverlayChannels, 0);
					section->addControl(&kcursorTracker, 1);
					section->addControl(&kshowLegend, 2);
					page->addSection(section, "Options");
				}

				if (auto section = new Signalizer::CContentPage::MatrixSection())
				{
					section->addControl(&kgain, 0);
					section->addControl(&kchannelConfiguration, 1);

					section->addControl(&kenvelopeSmooth, 0);
					section->addControl(&kenvelopeMode, 1);

					section->addControl(&kpctForDivision, 0);

					page->addSection(section, "Utility");
				}
				if (auto section = new Signalizer::CContentPage::MatrixSection())
				{
					section->addControl(&kwindow, 0);
					section->addControl(&ktimeMode, 1);
					section->addControl(&ktriggerMode, 0);
					section->addControl(&ktriggerPhaseOffset, 1);

					section->addControl(&ktriggerThreshold, 0);
					section->addControl(&ktriggerHysteresis, 1);

					section->addControl(&kcustomFrequency, 0);
					section->addControl(&ktriggerOnCustomFrequency, 1);

					section->addControl(&ktriggerChannel, 0);

					page->addSection(section, "Spatial");
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

					section->addControl(&kgraphColour, 0);
					section->addControl(&kbackgroundColour, 1);
					section->addControl(&ktrackerColour, 0);

					page->addSection(section, "Look");
				}

				if (auto section = new Signalizer::CContentPage::MatrixSection())
				{
					section->addControl(&kcolourSmoothingTime, 0);
					section->addControl(&kchannelColouring, 1);

					section->addControl(&kprimaryColour, 0);
					section->addControl(&ksecondaryColour, 1);

					section->addControl(&kfreqColourBlend, 0);
					section->addControl(&klowColour, 1);
					section->addControl(&kmidColour, 0);
					section->addControl(&khighColour, 1);

					page->addSection(section, "Spectral colouring");
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
			archive << kprimaryColour;
			archive << ktransform;
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
			archive << ktriggerOnCustomFrequency;
			archive << kcustomFrequency;
			archive << koverlayChannels;
			archive << kchannelColouring;
			archive << klowColour << kmidColour << khighColour;
			archive << ksecondaryColour;
			archive << kcolourSmoothingTime;
			archive << kcursorTracker;
			archive << ktrackerColour;
			archive << kfreqColourBlend;
			archive << ktriggerHysteresis;
			archive << ktriggerThreshold;
			archive << kshowLegend;
			archive << ktriggerChannel;
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
			builder >> kprimaryColour;
			builder >> ktransform;
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
			builder >> ktriggerOnCustomFrequency;
			builder >> kcustomFrequency;
			builder >> koverlayChannels;
			builder >> kchannelColouring;
			builder >> klowColour >> kmidColour >> khighColour;
			builder >> ksecondaryColour;
			builder >> kcolourSmoothingTime;

			if (version > cpl::Version(0, 3, 1))
			{
				builder >> kcursorTracker;
				builder >> ktrackerColour;
				builder >> kfreqColourBlend;
			}

			if (version > cpl::Version(0, 3, 2))
			{
				builder >> ktriggerHysteresis;
				builder >> ktriggerThreshold;
			}

			if (version > cpl::Version(0, 3, 3))
			{
				builder >> kshowLegend;
			}

			if (version > cpl::Version(0, 3, 4))
			{
				builder >> ktriggerChannel;
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

		cpl::CButton kantiAlias, kdiagnostics, kdotSamples, ktriggerOnCustomFrequency, koverlayChannels, kcursorTracker, kshowLegend;
		cpl::CValueInputControl kcustomFrequency;
		cpl::CValueKnobSlider
			kwindow, kgain, kprimitiveSize, kenvelopeSmooth, kpctForDivision, ktriggerPhaseOffset, kcolourSmoothingTime, kfreqColourBlend,
			ktriggerHysteresis, ktriggerThreshold, ktriggerChannel;
		cpl::CColourControl kprimaryColour, ksecondaryColour, kgraphColour, kbackgroundColour, klowColour, kmidColour, khighColour, ktrackerColour;
		cpl::CTransformWidget ktransform;
		cpl::CValueComboBox kenvelopeMode, ksubSampleInterpolationMode, kchannelConfiguration, ktriggerMode, ktimeMode, kchannelColouring;
		cpl::CPresetWidget kpresets;

		OscilloscopeContent & parent;

		SSOSurrogate<OscilloscopeController>
			editorSerializer,
			valueSerializer;
	};

	std::unique_ptr<StateEditor> OscilloscopeContent::createEditor()
	{
		return std::make_unique<OscilloscopeController>(*this);
	}

};
