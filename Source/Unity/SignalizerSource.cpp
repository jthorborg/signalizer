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

	file:LumpedSignalizer.cpp

		Direct include for a unity-build-embed version of Signalizer. Use the headers
		as you normally would.

*************************************************************************************/

#include "SignalizerConfiguration.h"

#include "../Oscilloscope/Oscilloscope.cpp"
#include "../Oscilloscope/OscilloscopeController.cpp"
#include "../Oscilloscope/OscilloscopeRendering.cpp"

#include "../Spectrum/Spectrum.cpp"
#include "../Spectrum/SpectrumController.cpp"
#include "../Spectrum/SpectrumDSP.cpp"
#include "../Spectrum/SpectrumRendering.cpp"

#include "../Vectorscope/Vectorscope.cpp"
#include "../Vectorscope/VectorscopeController.cpp"
#include "../Vectorscope/VectorscopeRendering.cpp"