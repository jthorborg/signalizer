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

	file:ProcessorPersistentState.h

		Parameters specific to operation of the plugin that are not exposed to the host.

*************************************************************************************/

#ifndef SIGNALIZER_PROCESSORPERSISTENTSTATE_H
	#define SIGNALIZER_PROCESSORPERSISTENTSTATE_H

	#include "Signalizer.h"
	#include "../Common/CommonSignalizer.h"

	namespace Signalizer
	{
		class ProcessorPersistentState : public cpl::CSerializer::Serializable
		{
		public:

			cpl::ValueEntityBase* getGraphSerializationValue() { return &graphSerializationChoiceValue; }

			ProcessorPersistentState()
				: graphSerializationChoiceValue(&graphTransformer, &graphFormatter)
				, graphFormatter(graphTransformer)
			{
				graphFormatter.setValues({ "Full", "Ignore in session", "Ignore always" });
			}

			bool shouldSerializeGraph()
			{
				return getGraphChoice() == GraphSerializationChoice::Full;
			}

			const char* getGraphSerializationDescription()
			{
				return
					"Advanced options for altering how sidechaining information is restored when saved as a preset from within a host.\n\n"
					"Full is normal operation.\n\n"
					"Ignore in session means any presets generated until this plugin is reloaded will be restored without sidechaining information "
					"- this is useful for generating presets that should not restore connections to other plugins initially, but afterwards operate normally.\n\n"
					"Ignore always means Signalizer should never attempt to save or restore sidechain information, and always assume a default layout of inputs to outputs. "
					"You can use this as a workaround or if you never use the sidechaining feature anyway and experience issues.";
			}

		private:

			enum class GraphSerializationChoice : std::int8_t
			{
				/// <summary>
				/// We always serialize the full graph when the host asks us to
				/// </summary>
				Full,
				/// <summary>
				/// We do not serialize the host graph, but we
				/// do not propagate that state to the next time we (or a clone) is loaded.
				/// This is useful for creating host presets that should function normally afterwards
				/// </summary>
				IgnoreSession,
				/// <summary>
				/// Same as <see cref="IgnoreSession"/>, except we perpetually stay ignorant of the host graph.
				/// This is useful as a workaround or if you want to create presets that work with the sidechaining feature.
				/// </summary>
				IgnoreAlways
			};

			GraphSerializationChoice getGraphChoice()
			{
				return graphSerializationChoiceValue.getAsTEnum<GraphSerializationChoice>();
			}

			void setGraphChoice(GraphSerializationChoice choice)
			{
				graphSerializationChoiceValue.setAsTEnum(choice);
			}

			void serialize(cpl::CSerializer::Archiver& ar, cpl::Version version) override 
			{ 
				auto current = getGraphChoice();

				if (current == GraphSerializationChoice::IgnoreSession)
					current = GraphSerializationChoice::Full;

				ar << current;
			}

			void deserialize(cpl::CSerializer::Builder& builder, cpl::Version version) override
			{
				GraphSerializationChoice choice;
				builder >> choice;

				CPL_RUNTIME_ASSERTION(choice != GraphSerializationChoice::IgnoreSession);

				setGraphChoice(choice);
			}

			cpl::ChoiceTransformer<cpl::ValueT> graphTransformer;
			cpl::ChoiceFormatter<cpl::ValueT> graphFormatter;

			cpl::SelfcontainedValue<cpl::ChoiceTransformer<cpl::ValueT>, cpl::ChoiceFormatter<cpl::ValueT>> graphSerializationChoiceValue;
		};
	};
#endif
