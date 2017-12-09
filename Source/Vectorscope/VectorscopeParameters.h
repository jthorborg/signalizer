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

			virtual std::unique_ptr<StateEditor> createEditor() override;

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