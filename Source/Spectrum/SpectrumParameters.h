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

	file:SpectrumParameters.h

		Interface for the vectorscope view

*************************************************************************************/

#ifndef SIGNALIZER_SPECTRUMPARAMETERS_H
	#define SIGNALIZER_SPECTRUMPARAMETERS_H

	#include "Signalizer.h"

	namespace Signalizer
	{

		class SpectrumContent final
			: public cpl::Parameters::UserContent
			, public ProcessorState
		{
		public:

			enum LineGraphs
			{
				None = -2, Transform = -1, LineMain = 0, LineSecond, LineEnd
			};

			enum class BinInterpolation
			{
				None,
				Linear,
				Lanczos
			};

			enum class DisplayMode
			{
				LineGraph,
				ColourSpectrum
			};

			enum class TransformAlgorithm
			{
				FFT, RSNT
			};

			enum class ViewScaling
			{
				Linear,
				Logarithmic
			};

			static const std::size_t numSpectrumColours = 5;
			static constexpr double kMinDbs = -24 * 16;
			// the maximum level of dbs to display
			static constexpr double kMaxDbs = 24 * 4;

			class SpectrumController final
				: public CContentPage
				, private cpl::ValueEntityBase::Listener
			{
			public:

				SpectrumController(SpectrumContent & parentValue)
					: parent(parentValue)
					, kdspWin(&parentValue.dspWin)
					, kslope(&parentValue.slope)
					
					, kviewScaling(&parentValue.viewScaling.param)
					, kalgorithm(&parentValue.algorithm.param)
					, kchannelConfiguration(&parentValue.channelConfiguration.param)
					, kdisplayMode(&parentValue.displayMode.param)
					, kbinInterpolation(&parentValue.binInterpolation.param)
					, kfrequencyTracker(&parentValue.frequencyTracker.param)
					, ktrackerColour(&parentValue.trackerColour)
					, ktrackerSmoothing(&parentValue.trackerSmoothing)

					, klowDbs(&parentValue.lowDbs)
					, khighDbs(&parentValue.highDbs)
					, kwindowSize(&parentValue.windowSize)
					, kpctForDivision(&parentValue.pctForDivision)
					, kblobSize(&parentValue.blobSize)
					, kframeUpdateSmoothing(&parentValue.frameUpdateSmoothing)
					, kspectrumStretching(&parentValue.spectrumStretching)
					, kprimitiveSize(&parentValue.primitiveSize)
					, kfloodFillAlpha(&parentValue.floodFillAlpha)
					, kreferenceTuning(&parentValue.referenceTuning)
					, kdiagnostics(&parentValue.diagnostics)
					, kfreeQ(&parentValue.freeQ)

					// TODO: Figure out automatic way to initialize N array in constructor
					, kgridColour(&parentValue.gridColour)
					, kbackgroundColour(&parentValue.backgroundColour)

					/*, kspecColours {
						{ &parentValue.specColours[0] },{ &parentValue.specColours[1] },{ &parentValue.specColours[2] },
						{ &parentValue.specColours[3] },{ &parentValue.specColours[4] }
					}

					, kspecRatios{
						{ &parentValue.specRatios[0] }, { &parentValue.specRatios[1] }, { &parentValue.specRatios[2] }, { &parentValue.specRatios[3] }, { &parentValue.specRatios[4] }
					}

					, klines{
						{ { &parentValue.lines[0].decay }, { &parentValue.lines[0].colourOne }, { &parentValue.lines[0].colourTwo } },
						{ { &parentValue.lines[1].decay }, { &parentValue.lines[1].colourOne }, { &parentValue.lines[1].colourTwo } },
					} */

					, presetManager(&valueSerializer, "spectrum")

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
					for(std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
					{
						klines.emplace_back(std::make_unique<LineControl>(&parentValue.lines[i].decay, &parentValue.lines[i].colourOne, &parentValue.lines[i].colourTwo));
					}

					for(std::size_t i = 0; i < numSpectrumColours; ++i)
					{
						kspecColours.emplace_back(std::make_unique<cpl::CColourControl>(&parentValue.specColours[i]));
						kspecRatios.emplace_back(std::make_unique<cpl::CValueKnobSlider>(&parentValue.specRatios[i]));
					}


					initControls();
					initUI();
					setTransformOptions();
				}

				cpl::SafeSerializableObject & getEditorSO() override { return editorSerializer; }

				~SpectrumController()
				{
					parent.algorithm.param.removeListener(this);
					parent.displayMode.param.removeListener(this);
					notifyDestruction();
				}

				virtual void valueEntityChanged(ValueEntityListener * sender, cpl::ValueEntityBase * value) override
				{
					if (value == &parent.displayMode.param || value == &parent.algorithm.param)
					{
						setTransformOptions();
					}
				}

				void setTransformOptions()
				{
					auto algo = cpl::enum_cast<TransformAlgorithm>(parent.algorithm.param.getTransformedValue());
					auto dispMode = cpl::enum_cast<DisplayMode>(parent.displayMode.param.getTransformedValue());

					if (algo == TransformAlgorithm::FFT)
					{
						kdspWin.setWindowOptions(cpl::CDSPWindowWidget::ChoiceOptions::All);
					}
					else if (algo == TransformAlgorithm::RSNT)
					{
						kdspWin.setWindowOptions(cpl::CDSPWindowWidget::ChoiceOptions::FiniteDFTWindows);
					}

					if (dispMode == DisplayMode::ColourSpectrum)
					{
						// disable all multichannel configurations
						for (std::size_t i = 0; i < (size_t)SpectrumChannels::End; i++)
						{
							if (i >(size_t)SpectrumChannels::OffsetForMono)
							{
								kchannelConfiguration.setEnabledStateFor(i, false);
							}
						}
					}
					else
					{
						// enable them.
						for (std::size_t i = 0; i < (size_t)SpectrumChannels::End; i++)
						{
							kchannelConfiguration.setEnabledStateFor(i, true);
						}
					}
				}

				void initControls()
				{
					parent.algorithm.param.addListener(this);
					parent.displayMode.param.addListener(this);

					for (std::size_t i = 0; i < numSpectrumColours; ++i)
					{
						kspecColours[i]->bSetTitle("Spectrum " + std::to_string(i + 1));
						kspecRatios[i]->bSetTitle("Gradient ratio " + std::to_string(i + 1));
						kspecRatios[i]->bSetDescription("How large a part of the gradient this colour occupies (the 4 values are normalized).");
						kspecColours[i]->bSetDescription("The four colours (together with the background colour) represents the linear interpolating colour function for the intensity found in the graph, for the colour spectrum, composing a gradient.");
					}

					// ------ titles -----------
					kviewScaling.bSetTitle("Graph scale");
					kalgorithm.bSetTitle("Transform algorithm");
					kchannelConfiguration.bSetTitle("Channel conf.");
					kdisplayMode.bSetTitle("Display mode");
					kfrequencyTracker.bSetTitle("Frequency tracking");
					kframeUpdateSmoothing.bSetTitle("Upd. smoothing");
					kbinInterpolation.bSetTitle("Bin interpolation");
					klowDbs.bSetTitle("Lower limit");
					khighDbs.bSetTitle("Upper limit");
					kwindowSize.bSetTitle("Window size");
					kfreeQ.setSingleText("Unbound Q");
					kdiagnostics.setSingleText("Diagnostics");
					kdiagnostics.setToggleable(true);
					kfreeQ.setToggleable(true);
					kspectrumStretching.bSetTitle("Spectrum stretch");
					kprimitiveSize.bSetTitle("Primitive size");
					kfloodFillAlpha.bSetTitle("Flood fill %");
					kgridColour.bSetTitle("Grid colour");
					kbackgroundColour.bSetTitle("Backg. colour");
					kreferenceTuning.bSetTitle("A4 ref. tuning");
					kpctForDivision.bSetTitle("Grid div. space");
					kblobSize.bSetTitle("Update speed");
					ktrackerSmoothing.bSetTitle("Tracker smooth");
					ktrackerColour.bSetTitle("Tracker colour");

					// ------ descriptions -----
					kviewScaling.bSetDescription("Set the scale of the frequency-axis of the coordinate system.");
					kalgorithm.bSetDescription("Select the algorithm used for transforming the incoming audio data.");
					kchannelConfiguration.bSetDescription("Select how the audio channels are interpreted.");
					kdisplayMode.bSetDescription("Select how the information is displayed; line graphs are updated each frame while the colour spectrum maintains the previous history.");
					kbinInterpolation.bSetDescription("Choice of interpolation for transform algorithms that produce a discrete set of values instead of an continuous function.");
					kdiagnostics.bSetDescription("Toggle diagnostic information in top-left corner.");
					klowDbs.bSetDescription("The lower limit of the displayed dynamic range.");
					khighDbs.bSetDescription("The upper limit of the displayed dynamic range");
					kwindowSize.bSetDescription("The window size of the audio data, affects time/frequency resolution.");
					kgridColour.bSetDescription("The colour of the dB/frequency grid.");
					kbackgroundColour.bSetDescription("The colour of the background.");
					kpctForDivision.bSetDescription("The minimum amount of free space that triggers a recursed frequency grid division; smaller values draw more frequency divisions.");
					kblobSize.bSetDescription("Controls how much audio data a horizontal unit represents; effectively controls the update rate of the colour spectrum.");
					kframeUpdateSmoothing.bSetDescription("Reduces jitter in spectrum updates at the (possible) expense of higher graphical latency.");
					kfreeQ.bSetDescription("Frees the quality factor from being bounded by the window size for transforms that support it. "
						"Although it (possibly) makes response time slower, it also makes the time/frequency resolution exact, and is a choice for analyzing static material.");
					kspectrumStretching.bSetDescription("Stretches the spectrum horizontally, emulating a faster update rate (useful for transforms which are not continuous).");
					kfrequencyTracker.bSetDescription("Specifies which pair of graphs that is evaluated for nearby peak estimations.");
					kprimitiveSize.bSetDescription("The size of the rendered primitives (eg. lines or points).");
					kfloodFillAlpha.bSetDescription("For line graphs, add a flood fill of the same colour for each line with the following alpha %");
					kreferenceTuning.bSetDescription("Reference tuning for A4; used when converting to/from musical notes and frequencies");
					ktrackerSmoothing.bSetDescription("Eliminates small fluctuations and holds analysis values in the frequency tracker for a longer time");
					ktrackerColour.bSetDescription("Colour of the frequency tracker");

					klines[LineGraphs::LineMain]->colourOne.bSetDescription("The colour of the first channel of the main graph.");
					klines[LineGraphs::LineMain]->colourTwo.bSetDescription("The colour of the second channel of the main graph.");

					klines[LineGraphs::LineMain]->colourOne.bSetTitle("Graph 1 colour");
					klines[LineGraphs::LineMain]->colourTwo.bSetTitle("Graph 2 colour");

					klines[LineGraphs::LineMain]->decay.bSetTitle("Main decay");
					klines[LineGraphs::LineMain]->decay.bSetDescription("Decay rate of the main graph channels; allows the graph to decay more slowly, but still reacting to peaks.");

					for (std::size_t i = LineGraphs::LineMain + 1; i < LineGraphs::LineEnd; ++i)
					{
						auto graphNumber = std::to_string(i);
						klines[i]->colourOne.bSetTitle("Aux 1 colour");
						klines[i]->colourTwo.bSetTitle("Aux 2 colour");
						klines[i]->decay.bSetTitle("Aux " + graphNumber + " decay");
						klines[i]->decay.bSetDescription("Decay rate of auxillary graph " + graphNumber + " channels; allows the graph to decay more slowly, but still reacting to peaks.");
						klines[i]->colourOne.bSetDescription("The colour of the first channel of auxillary graph " + graphNumber + ".");
						klines[i]->colourTwo.bSetDescription("The colour of the second channel of auxillary graph " + graphNumber + ".");
					}

				}

				void initUI()
				{
					if (auto page = addPage("Settings", "icons/svg/gear.svg"))
					{
						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&kviewScaling, 0);
							section->addControl(&kchannelConfiguration, 0);
							section->addControl(&kdisplayMode, 1);
							section->addControl(&kfrequencyTracker, 1);
							page->addSection(section);
						}
						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&klowDbs, 1);
							section->addControl(&khighDbs, 0);
							section->addControl(&kblobSize, 0);
							section->addControl(&kwindowSize, 1);
							section->addControl(&kpctForDivision, 0);
							section->addControl(&kspectrumStretching, 1);
							page->addSection(section);
						}
						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
							{
								section->addControl(&klines[i]->decay, i & 1);
							}
							page->addSection(section);
						}
					}
					if (auto page = addPage("Algorithm", "icons/svg/formulae.svg"))
					{
						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&kalgorithm, 0);
							section->addControl(&kbinInterpolation, 1);
							page->addSection(section);
						}
						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&kdspWin, 0);
							page->addSection(section);
						}

						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&kslope, 0);
							page->addSection(section);
						}

						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&kfreeQ, 0);
							page->addSection(section);
						}

					}

					if (auto page = addPage("Rendering", "icons/svg/brush.svg"))
					{
						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&kgridColour, 0);
							section->addControl(&kbackgroundColour, 1);
							section->addControl(&ktrackerColour, 0);
							page->addSection(section);
						}
						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
							{
								section->addControl(&klines[i]->colourOne, 0);
								section->addControl(&klines[i]->colourTwo, 1);
							}
							page->addSection(section);
						}

						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							for (std::size_t i = 0; i < numSpectrumColours; ++i)
							{
								section->addControl(kspecColours[i].get(), 0);
								section->addControl(kspecRatios[i].get(), 1);
							}
							page->addSection(section);
						}


					}
					if (auto page = addPage("Utility", "icons/svg/wrench.svg"))
					{
						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&presetManager, 0);
							page->addSection(section);
						}

						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&kframeUpdateSmoothing, 0);
							section->addControl(&kprimitiveSize, 1);
							section->addControl(&kfloodFillAlpha, 0);
							section->addControl(&kreferenceTuning, 1);
							section->addControl(&ktrackerSmoothing, 0);
							section->addControl(&kdiagnostics, 1);
							page->addSection(section);
						}
					}
				}

			private:

				void serializeEditorSettings(cpl::CSerializer::Archiver & archive, cpl::Version version)
				{
					archive << kviewScaling;
					archive << kalgorithm;
					archive << kchannelConfiguration;
					archive << kdisplayMode;
					archive << khighDbs;
					archive << klowDbs;
					archive << kwindowSize;
					archive << kpctForDivision;

					for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
					{
						archive << klines[i]->colourOne;
						archive << klines[i]->colourTwo;
						archive << klines[i]->decay;
					}

					archive << kgridColour;
					archive << kblobSize;
					archive << kbackgroundColour;
					archive << kframeUpdateSmoothing;

					for (std::size_t i = 0; i < numSpectrumColours; ++i)
					{
						archive << *kspecColours[i];
						archive << *kspecRatios[i];
					}

					archive << kbinInterpolation;
					archive << cpl::Serialization::Reserve(sizeof(double) * 2);
					archive << kdspWin;
					archive << kfreeQ;
					archive << kspectrumStretching;
					archive << kfrequencyTracker;
					archive << kprimitiveSize;
					archive << kfloodFillAlpha;
					archive << kslope;
					archive << kreferenceTuning;
					archive << ktrackerSmoothing;
					archive << ktrackerColour;
				}

				void deserializeEditorSettings(cpl::CSerializer::Archiver & builder, cpl::Version version)
				{
					// in general, controls should never restore values. However, older versions
					// of Signalizer does exactly this, so to keep backwards-compatibility, we
					// can obtain the preset values through this.
					cpl::Serialization::ScopedModifier m(cpl::CSerializer::Modifiers::RestoreValue, version < cpl::Version(0, 2, 8));
					builder << m;

					builder >> kviewScaling;
					builder >> kalgorithm;
					builder >> kchannelConfiguration;
					builder >> kdisplayMode;
					// set high first, so low isn't capped
					builder >> khighDbs;
					builder >> klowDbs;
					builder >> kwindowSize;
					builder >> kpctForDivision;

					for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
					{
						builder >> klines[i]->colourOne;
						builder >> klines[i]->colourTwo;
						builder >> klines[i]->decay;
					}

					builder >> kgridColour;
					builder >> kblobSize;
					builder >> kbackgroundColour;
					builder >> kframeUpdateSmoothing;

					for (std::size_t i = 0; i < numSpectrumColours; ++i)
					{
						builder >> *kspecColours[i];
						builder >> *kspecRatios[i];
					}

					builder >> kbinInterpolation;
					// view rect
					double x1, x2;
					builder >> x1 >> x2;
					if (builder.getModifier(cpl::CSerializer::Modifiers::RestoreValue))
					{
						// TODO: Fill in view rect
					}

					builder >> kdspWin;
					builder >> kfreeQ;
					builder >> kspectrumStretching;
					builder >> kfrequencyTracker;
					builder >> kprimitiveSize;
					builder >> kfloodFillAlpha;

					if (version > cpl::Version::fromParts(0, 2, 6))
					{
						builder >> kslope;
						builder >> kreferenceTuning;
					}

					if (version >= cpl::Version(0, 3, 1))
					{
						builder >> ktrackerSmoothing >> ktrackerColour;
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

				cpl::CValueComboBox
					kviewScaling,
					kalgorithm,
					kchannelConfiguration,
					kdisplayMode,
					kbinInterpolation,
					kfrequencyTracker;

				cpl::CDSPWindowWidget kdspWin;
				cpl::CPowerSlopeWidget kslope;
				cpl::CValueKnobSlider
					klowDbs,
					khighDbs,
					kwindowSize,
					kpctForDivision,
					kblobSize,
					kframeUpdateSmoothing,
					kspectrumStretching,
					kprimitiveSize,
					kfloodFillAlpha,
					kreferenceTuning,
					ktrackerSmoothing;

				cpl::CColourControl kgridColour, kbackgroundColour, ktrackerColour;

				struct LineControl
				{
					LineControl(cpl::ValueEntityBase * decayValue, cpl::ColourValue * one, cpl::ColourValue * two) : decay(decayValue), colourOne(one), colourTwo(two) {}
					cpl::CValueKnobSlider decay;
					cpl::CColourControl colourOne, colourTwo;
				};

				// TODO: turn into array once aggregrate initialization of member arrays doesn't require present copy or move constructors
				std::vector<std::unique_ptr<LineControl>> klines;
				std::vector<std::unique_ptr<cpl::CColourControl>> kspecColours;
				std::vector<std::unique_ptr<cpl::CValueKnobSlider>> kspecRatios;

				cpl::CPresetWidget presetManager;
				cpl::CButton kdiagnostics, kfreeQ;

				SpectrumContent & parent;

				SSOSurrogate<SpectrumController>
					editorSerializer,
					valueSerializer;
			};

			SpectrumContent(std::size_t offset, bool shouldCreateShortNames, SystemView system)
				: systemView(system)
				, parameterSet("Spectrum", "SC.", system.getProcessor(), static_cast<int>(offset))
				, audioHistoryTransformatter(system.getAudioStream(), audioHistoryTransformatter.Samples)

				, dynamicRange(kMinDbs, kMaxDbs)
				, literalDBFormatter("dB")
				, reverseUnitRange(1, 0)
				, unitRange(0, 1)
				, blobRange(0.5, 1000)
				, spectrumStretchRange(1, 20)
				, primitiveRange(0.01, 10)
				, referenceRange(220, 880)
				, smoothCappedRange(0, 0.996)
				, trackerSmoothRange(0, 1000)

				, dspWin(windowBehavior, "DSPWin")
				, slope(slopeBehavior, "Slope")

				, viewScaling("VScale")
				, algorithm("Algo")
				, channelConfiguration("ChConf")
				, displayMode("DispMode")
				, binInterpolation("BinInt")
				, frequencyTracker("FTracker")

				, lowDbs("LowDBs", dynamicRange, literalDBFormatter)
				, highDbs("HighDBs", dynamicRange, literalDBFormatter)
				, windowSize("WindowSize", audioHistoryTransformatter, audioHistoryTransformatter)
				, pctForDivision("PctDiv", unitRange, basicFormatter)
				, blobSize("BlobSize", blobRange, msFormatter)
				, frameUpdateSmoothing("FrmSmooth", smoothCappedRange, basicFormatter)
				, spectrumStretching("SpcStretch", spectrumStretchRange, basicFormatter)
				, primitiveSize("PrtSize", primitiveRange, basicFormatter)
				, floodFillAlpha("AlphaFill", unitRange, pctFormatter)
				, referenceTuning("RefTuning", referenceRange, basicFormatter)
				, viewLeft("ViewLeft", unitRange, basicFormatter)
				, viewRight("ViewRight", reverseUnitRange, basicFormatter)
				, diagnostics("Diagnostics", boolRange, boolFormatter)
				, freeQ("FreeQ", boolRange, boolFormatter)
				, trackerSmoothing("TrckSmth", trackerSmoothRange, msFormatter)

				, colourBehaviour()

				// TODO: Figure out automatic way to initialize N array in constructor
				, gridColour(colourBehaviour, "Grid.")
				, backgroundColour(colourBehaviour, "Bck.")
				, trackerColour(colourBehaviour, "Trck.")
				, specColours { { colourBehaviour , "Grdnt1."}, { colourBehaviour , "Grdnt2." }, { colourBehaviour , "Grdnt3." }, { colourBehaviour , "Grdnt4." }, { colourBehaviour , "Grdnt5." } }

				, specRatios{
					{ "GradRatio1", unitRange, basicFormatter },
					{ "GradRatio2", unitRange, basicFormatter },
					{ "GradRatio3", unitRange, basicFormatter },
					{ "GradRatio4", unitRange, basicFormatter },
					{ "GradRatio5", unitRange, basicFormatter }
				}

				, lines {
					{ { "Grph1.Decay", unitRange, dbSecFormatter }, { colourBehaviour , "Grph1.1."}, { colourBehaviour , "Grph1.2." } },
					{ { "Grph2.Decay", unitRange, dbSecFormatter }, { colourBehaviour , "Grph2.1." }, { colourBehaviour , "Grph2.1." } }
				}

			{
				dbSecFormatter.setUnit("dB/s");

				viewScaling.fmt.setValues({ "Linear", "Logarithmic" });
				algorithm.fmt.setValues({ "FFT", "Resonator" });
				channelConfiguration.fmt.setValues({ "Left", "Right", "Mid/Merge", "Side", "Phase", "Separate", "Mid+Side", "Complex" });
				displayMode.fmt.setValues({ "Line graph", "Colour spectrum" });
				binInterpolation.fmt.setValues({ "None", "Linear", "Lanczos" });

				std::vector<std::string> frequencyTrackingOptions;

				frequencyTrackingOptions.push_back("None");
				frequencyTrackingOptions.push_back("Transform");
				for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
				{
					if (i == LineGraphs::LineMain)
						frequencyTrackingOptions.push_back("Main graph");
					else
						frequencyTrackingOptions.push_back("Aux graph " + std::to_string(i));
				}

				frequencyTracker.fmt.setValues(frequencyTrackingOptions);

				auto singleParameters = {
					&lowDbs, &highDbs, &windowSize, &pctForDivision, &blobSize, &frameUpdateSmoothing, &spectrumStretching,
					&primitiveSize, &floodFillAlpha, &referenceTuning, &viewLeft, &viewRight, &diagnostics, &freeQ, &trackerSmoothing
				};

				for (auto sparam : singleParameters)
				{
					parameterSet.registerSingleParameter(sparam->generateUpdateRegistrator());
				}

				for (auto sparam : { &viewScaling, &algorithm, &channelConfiguration, &displayMode, &binInterpolation, &frequencyTracker })
				{
					parameterSet.registerSingleParameter(sparam->param.generateUpdateRegistrator());
				}

				for (auto & sparam : specRatios)
				{
					parameterSet.registerSingleParameter(sparam.generateUpdateRegistrator());
				}

				auto regBundle = [&](cpl::Parameters::BundleUpdate<ParameterSet::ParameterView> & bundle, const std::string & n) { parameterSet.registerParameterBundle(&bundle, n); };

				regBundle(dspWin, "DWin.");
				regBundle(slope, "Slope.");
				regBundle(gridColour, gridColour.getBundleName());
				regBundle(backgroundColour, backgroundColour.getBundleName());
				regBundle(trackerColour, trackerColour.getBundleName());

				for (std::size_t i = 0; i < std::extent<decltype(specColours)>::value; ++i)
					regBundle(specColours[i], specColours[i].getBundleName());

				for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
				{
					parameterSet.registerSingleParameter(lines[i].decay.generateUpdateRegistrator());
					regBundle(lines[i].colourOne, lines[i].colourOne.getBundleName());
					regBundle(lines[i].colourTwo, lines[i].colourTwo.getBundleName());
				}

				parameterSet.seal();

				postParameterInitialization();

			}

			virtual std::unique_ptr<StateEditor> createEditor() override
			{
				return std::make_unique<SpectrumController>(*this);
			}

			virtual ParameterSet & getParameterSet() override
			{
				return parameterSet;
			}

			virtual void serialize(cpl::CSerializer::Archiver & archive, cpl::Version v) override
			{
				archive << viewScaling.param;
				archive << algorithm.param;
				archive << channelConfiguration.param;
				archive << displayMode.param;
				archive << highDbs;
				archive << lowDbs;
				archive << windowSize;
				archive << pctForDivision;

				for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
				{
					archive << lines[i].colourOne;
					archive << lines[i].colourTwo;
					archive << lines[i].decay;
				}

				archive << gridColour;
				archive << blobSize;
				archive << backgroundColour;
				archive << frameUpdateSmoothing;

				for (std::size_t i = 0; i < numSpectrumColours; ++i)
				{
					archive << specColours[i];
					archive << specRatios[i];
				}

				archive << binInterpolation.param;
				archive << viewLeft;
				archive << viewRight;
				archive << dspWin;
				archive << freeQ;
				archive << spectrumStretching;
				archive << frequencyTracker.param;
				archive << primitiveSize;
				archive << floodFillAlpha;
				archive << slope;
				archive << referenceTuning;
				archive << audioHistoryTransformatter;

				archive << trackerSmoothing << trackerColour;
			}

			virtual void deserialize(cpl::CSerializer::Builder & builder, cpl::Version v) override
			{
				builder >> viewScaling.param;
				builder >> algorithm.param;
				builder >> channelConfiguration.param;
				builder >> displayMode.param;
				// set high first, so low isn't capped
				builder >> highDbs;
				builder >> lowDbs;
				builder >> windowSize;
				builder >> pctForDivision;

				for (std::size_t i = 0; i < LineGraphs::LineEnd; ++i)
				{
					builder >> lines[i].colourOne;
					builder >> lines[i].colourTwo;
					builder >> lines[i].decay;
				}

				builder >> gridColour;
				builder >> blobSize;
				builder >> backgroundColour;
				builder >> frameUpdateSmoothing;

				for (std::size_t i = 0; i < numSpectrumColours; ++i)
				{
					builder >> specColours[i];
					builder >> specRatios[i];
				}

				builder >> binInterpolation.param;
				builder >> viewLeft;
				builder >> viewRight;
				builder >> dspWin;
				builder >> freeQ;
				builder >> spectrumStretching;
				builder >> frequencyTracker.param;
				builder >> primitiveSize;
				builder >> floodFillAlpha;
				builder >> slope;
				builder >> referenceTuning;
				builder >> audioHistoryTransformatter;

				if (v >= cpl::Version(0, 3, 1))
				{
					builder >> trackerSmoothing >> trackerColour;
				}
			}

			SystemView systemView;
			ParameterSet parameterSet;
			Signalizer::AudioHistoryTransformatter<ParameterSet::ParameterView> audioHistoryTransformatter;

			typedef cpl::ParameterValue<ParameterSet::ParameterView> Parameter;



			cpl::ParameterWindowDesignValue<ParameterSet::ParameterView> dspWin;
			cpl::ParameterWindowDesignValue<ParameterSet::ParameterView>::SharedBehaviour windowBehavior;

			cpl::ParameterPowerSlopeValue<ParameterSet::ParameterView> slope;
			cpl::ParameterPowerSlopeValue<ParameterSet::ParameterView>::SharedBehaviour slopeBehavior;

			ChoiceParameter
				viewScaling,
				algorithm,
				channelConfiguration,
				displayMode,
				binInterpolation,
				frequencyTracker;

			Parameter
				lowDbs,
				highDbs,
				windowSize,
				pctForDivision,
				/// <summary>
				/// How often the audio stream is sampled for updating in the screen. Only for DisplayMode::ColourSpectrum modes.
				/// Effectively, this sets each horizontal pixel to contain this amount of miliseconds of information.
				/// </summary>
				blobSize,
				/// <summary>
				/// For displaymode == ColourSpectrum, this value indicates how much devations in amount of frame updates are smoothed such that
				/// the display doesn't jitter too much.
				/// </summary>
				frameUpdateSmoothing,
				spectrumStretching,
				primitiveSize,
				floodFillAlpha,
				referenceTuning,
				viewLeft,
				viewRight,
				freeQ,
				diagnostics,
				specRatios[numSpectrumColours],
				trackerSmoothing;

			cpl::ParameterColourValue<ParameterSet::ParameterView>
				gridColour,
				backgroundColour,
				specColours[numSpectrumColours],
				trackerColour;

			struct LineControl
			{
				Parameter decay;
				cpl::ParameterColourValue<ParameterSet::ParameterView> colourOne, colourTwo;
			} lines[LineGraphs::LineEnd];

			cpl::BooleanRange<double> boolRange;
			cpl::BooleanFormatter<double> boolFormatter;

			cpl::DBFormatter<SFloat> dbSecFormatter;
			cpl::UnitFormatter<SFloat> literalDBFormatter;
			cpl::UnitFormatter<SFloat> msFormatter;
			cpl::PercentageFormatter<SFloat> pctFormatter;
			cpl::BasicFormatter<SFloat> basicFormatter;

			cpl::LinearRange<SFloat>
				dynamicRange,
				reverseUnitRange,
				unitRange,
				smoothCappedRange,
				spectrumStretchRange,
				primitiveRange,
				referenceRange,
				trackerSmoothRange;

			cpl::ExponentialRange<SFloat>
				blobRange;

			cpl::ParameterColourValue<ParameterSet::ParameterView>::SharedBehaviour colourBehaviour;


		private:

			void postParameterInitialization()
			{
				audioHistoryTransformatter.initialize(windowSize.getParameterView());
			}
		};

	};

#endif
