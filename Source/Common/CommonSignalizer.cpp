/*************************************************************************************

	Signalizer - cross-platform audio visualization plugin - v. 0.x.y

	Copyright (C) 2023 Janus Lynggaard Thorborg (www.jthorborg.com)

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

	file:CommonSignalizer.cpp

		Implementation of CommonSignalizer.h

*************************************************************************************/

#include "CommonSignalizer.h"

namespace Signalizer
{

	CriticalSection<Assumptions>& getFailedAssumptions()
	{
		static CriticalSection<Assumptions> failedAssumptions;
		return failedAssumptions;
	}

	void revealExceptionLog()
	{
		juce::File f(cpl::GetExceptionLogFilePath());

		if (f.existsAsFile())
			f.revealToUser();
	}

	bool triggerNonTerminalAssumption(const char* assumption, const char* file, const int line, const char* function)
	{
		std::size_t seed = 0;

		// boost::hash_combine
		auto hash = [](std::size_t& seed, auto thing)
		{
			auto hasher = std::hash<decltype(thing)>();
			seed ^= hasher(thing) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		};

		hash(seed, assumption); hash(seed, file); hash(seed, line); hash(seed, function);

		auto assumptions = getFailedAssumptions().lock();

		// only report once
		if (assumptions->set.count(seed))
			return false;

		assumptions->set.insert(seed);

		std::string message =
			::cpl::programInfo.name + " (" + ::cpl::programInfo.version.toString() + "): in " + file + ":" + ::std::to_string(line) + " -> " + function + ":\n" +
			"Runtime assumption \"" + assumption + "\" failed.";

		CPL_BREAKIFDEBUGGED();
		CPL_DEBUGOUT((message + "\n").c_str());
		cpl::LogException(message);

		assumptions->current.emplace_back(std::move(message));

		return false;
	}
}
