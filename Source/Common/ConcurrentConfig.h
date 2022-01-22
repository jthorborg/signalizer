/*************************************************************************************

	Signalizer - cross-platform audio visualization plugin - v. 0.x.y

	Copyright (C) 2021 Janus Lynggaard Thorborg (www.jthorborg.com)

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

	file:ConcurrentConfig.h

		Read-only atomic access to properties linked to this processor

*************************************************************************************/


#ifndef SIGNALIZER_CONCURRENTCONFIG_H
	#define SIGNALIZER_CONCURRENTCONFIG_H

	#include <cpl/Common.h>
	#include <cpl/lib/weak_atomic.h>

	namespace Signalizer
	{
		struct ConcurrentConfig
		{
			cpl::relaxed_atomic<double> sampleRate;
			cpl::relaxed_atomic<std::size_t> historySize, historyCapacity;
			cpl::relaxed_atomic<double> bpm;
			cpl::relaxed_atomic<unsigned> numChannels;
		};
	}

#endif
