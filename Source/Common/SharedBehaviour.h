/*************************************************************************************
 
	Signalizer - cross-platform audio visualization plugin - v. 0.x.y
 
	Copyright (C) 2017 Janus Lynggaard Thorborg (www.jthorborg.com)
 
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
 
	file:SharedBehaviour.h

		Shared/global options for views.
 
*************************************************************************************/

#ifndef SIGNALIZER_SHAREDBEHAVIOUR_H
#define SIGNALIZER_SHAREDBEHAVIOUR_H

#include <cpl/lib/weak_atomic.h>
#include "cpl/profiling/Profiling.h"

namespace Signalizer
{
	class SharedBehaviour
	{
	public:

		SharedBehaviour(
			std::shared_ptr<cpl::Profiling::Lane> realtimeLane, 
			std::shared_ptr<cpl::Profiling::Lane> asyncDspLane,
			std::shared_ptr<cpl::Profiling::Lane> renderingLane
		)
			: realtime(realtimeLane)
			, async(asyncDspLane)
			, rendering(renderingLane)
		{

		}

		// TODO: Don't expose these in a global. Move overlay rendering to MainEditor instead.
		cpl::Profiling::Lane& getRealtimeLane() const { return *realtime; }
		cpl::Profiling::Lane& getAsyncDSPLane() const { return *async; };
		std::shared_ptr<cpl::Profiling::Lane> getRenderingLane() const { return rendering; };

		cpl::relaxed_atomic<bool>
			hideWidgetsOnMouseExit = false,
			stopProcessingOnSuspend = false,
			showLegend = false;

	private:
		mutable std::shared_ptr<cpl::Profiling::Lane>
			realtime, async, rendering;
	};
};


#endif