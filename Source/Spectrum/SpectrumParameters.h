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

			virtual std::unique_ptr<StateEditor> createEditor() override;

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
