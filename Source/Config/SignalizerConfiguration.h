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

	file:SignalizerConfiguration.h

		A (possibly) unversioned file allowing to configure Signalizer (components)
		to your project. Example: Here you can change whether Signalizer internally
		uses floats or doubles, the threaded parameter model and such.

		Signalizer requires this file is directly includeable, but it can be anywhere.

*************************************************************************************/


#ifndef SIGNALIZER_CONFIGURATION_SIGNALIZER_H
	#define SIGNALIZER_CONFIGURATION_SIGNALIZER_H

	#include <cpl/CAudioStream.h>
	#include <cpl/infrastructure/parameters/ParameterSystem.h>
	#include <cpl/infrastructure/values/Values.h>

	namespace Signalizer
	{
		/// <summary>
		/// Floating-point type used for parameters etc. in Signalizer
		/// </summary>
		typedef double SFloat;
		/// <summary>
		/// Floating point type used for audio
		/// </summary>
		typedef float AFloat;
		/// <summary>
		/// Floating point type used for parameters of the host system
		/// </summary>
		typedef float PFloat;

		typedef cpl::FormattedParameter<SFloat, cpl::ThreadedParameter<SFloat>> Parameter;
		typedef cpl::ParameterGroup<SFloat, PFloat, Parameter> ParameterSet;


		// TODO: Figure out why sizes around 256 causes buffer overruns
		typedef cpl::CAudioStream<AFloat, 64> AudioStream;


	};
#endif
