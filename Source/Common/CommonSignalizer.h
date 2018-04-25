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

	file:CommonSignalizer.h

		Interface for common code and types used in Signalizer.

*************************************************************************************/


#ifndef SIGNALIZER_COMMON_SIGNALIZER_H
	#define SIGNALIZER_COMMON_SIGNALIZER_H

	#include <cpl/Common.h>
	#include <complex>
	#include "SignalizerConfiguration.h"

	namespace Signalizer
	{
		enum class EnvelopeModes : int
		{
			None,
			RMS,
			PeakDecay
		};

		enum class SubSampleInterpolation : int
		{
			None,
			Rectangular,
			Linear,
			Lanczos
		};

		class StateEditor;
		class ProcessorState;
		class SystemView;

		class ProcessorState
			: public cpl::SafeSerializableObject
		{
		public:

			virtual std::unique_ptr<StateEditor> createEditor() = 0;
			virtual ParameterSet & getParameterSet() = 0;

			virtual ~ProcessorState() {}

		};

		class SystemView
		{
		public:

			SystemView(AudioStream & audioStream, ParameterSet::AutomatedProcessor & automatedProcessor)
				: stream(audioStream), processor(automatedProcessor)
			{

			}

			ParameterSet::AutomatedProcessor & getProcessor() noexcept { return processor; }
			AudioStream & getAudioStream() noexcept { return stream; }

		private:
			AudioStream & stream;
			ParameterSet::AutomatedProcessor & processor;
		};

		struct ChoiceParameter
		{
			cpl::ParameterValue<ParameterSet::ParameterView> param;
			cpl::ChoiceFormatter<SFloat> fmt;
			cpl::ChoiceTransformer<SFloat> tsf;

			ChoiceParameter(const std::string & name)
				: param(name, tsf, fmt)
				, fmt(tsf)
			{
			}
		};

		template<typename ParameterView>
		class AudioHistoryTransformatter
			: public ParameterView::ParameterType::Transformer
			, public ParameterView::ParameterType::Formatter
			, public cpl::CSerializer::Serializable
			, protected AudioStream::Listener
		{
		public:

			typedef typename ParameterView::ParameterType::ValueType ValueType;

			enum Mode
			{
				Samples,
				Milliseconds
			};

			enum Scaling
			{
				Linear,
				Exponential
			};

			AudioHistoryTransformatter(AudioStream & audioStream, Mode mode = Milliseconds)
				: param(nullptr), stream(audioStream), m(mode), lastCapacity(0)
			{

			}

			void setModeFromUI(Mode newMode)
			{
				if (m != newMode)
				{
					m = newMode;
					// updates displays, even though the internal fraction stays the same.
					this->param->updateFromUINormalized(param->getValueNormalized());
				}
			}

			void initialize(ParameterView & view)
			{
				param = &view;
				listenToSource(stream, true);
			}

			void serialize(cpl::CSerializer::Archiver & ar, cpl::Version v) override
			{
				ar << lastCapacity.load(std::memory_order_relaxed);
			}

			void deserialize(cpl::CSerializer::Builder & ar, cpl::Version v) override
			{
				std::uint64_t val;
				ar >> val;
				lastCapacity.store(val, std::memory_order_relaxed);
			}

			~AudioHistoryTransformatter()
			{
				detachFromSource();
			}

		protected:

			virtual void onAsyncChangedProperties(const Stream & changedSource, const typename Stream::AudioStreamInfo & before) override
			{
				// TODO: what if oldCapacity == 0?
				const auto oldFraction = param->getValueNormalized();
				auto oldCapacity = lastCapacity.load(std::memory_order_relaxed);
				auto beforeCapacity = before.audioHistoryCapacity.load(std::memory_order_acquire);
				if (oldCapacity == 0)
					oldCapacity = beforeCapacity;

				const auto newCapacity = changedSource.getInfo().audioHistoryCapacity.load(std::memory_order_relaxed);

				if (newCapacity > 0)
					lastCapacity.store(newCapacity, std::memory_order_relaxed);

				if(oldCapacity == 0 || newCapacity == 0)
				{
					param->updateFromProcessorNormalized(oldFraction, cpl::Parameters::UpdateFlags::All & ~cpl::Parameters::UpdateFlags::RealTimeSubSystem);
				}
				else
				{
					const auto sampleSizeBefore = oldCapacity * oldFraction;
					const auto newFraction = sampleSizeBefore / newCapacity;
					if (oldFraction != newFraction || beforeCapacity == 0)
						param->updateFromProcessorNormalized(newFraction, cpl::Parameters::UpdateFlags::All & ~cpl::Parameters::UpdateFlags::RealTimeSubSystem);
				}
			}

			virtual bool format(const ValueType & val, std::string & buf) override
			{
				char buffer[100];

				if (m == Milliseconds)
				{
					sprintf_s(buffer, "%.2f ms", 1000 * val / stream.getInfo().sampleRate.load(std::memory_order_relaxed));
				}
				else
				{
					sprintf_s(buffer, "%.0f smps", val);
				}

				buf = buffer;
				return true;
			}

			virtual bool interpret(const cpl::string_ref buf, ValueType & val) override
			{
				ValueType collectedValue;


				if (cpl::lexicalConversion(buf, collectedValue))
				{
					bool notSamples = true;

					if (buf.find("s") != std::string::npos && (notSamples = buf.find("smps") == std::string::npos))
					{
						if (buf.find("ms") != std::string::npos)
						{
							collectedValue /= 1000;
						}
						collectedValue *= stream.getInfo().sampleRate.load(std::memory_order_relaxed);
					}
					else
					{
						// assume value is in miliseconds
						if (m == Milliseconds && notSamples)
						{
							collectedValue /= 1000;
							collectedValue *= stream.getInfo().sampleRate.load(std::memory_order_relaxed);
						}
					}

					val = collectedValue;
					return true;

				}

				return false;

			}

			virtual ValueType transform(ValueType val) const noexcept override
			{
				auto samples = std::round(val * stream.getAudioHistoryCapacity());

				/* if (m == Milliseconds)
				{
					samples /= stream.getInfo().sampleRate.load(std::memory_order_relaxed);
					samples *= 1000;
				} */


				return samples;
			}


			virtual ValueType normalize(ValueType val) const noexcept override
			{
				/* if (m == Milliseconds)
				{
					val /= stream.getInfo().sampleRate.load(std::memory_order_relaxed);
					val *= 1000;
				} */

				return val / stream.getAudioHistoryCapacity();
			}

			ParameterView * param;
			AudioStream & stream;
			Mode m;
			/// <summary>
			/// Represents the last capacity change we were informated about,
			/// and thus currently represents the magnitude scale of the fraction.
			/// </summary>
			std::atomic<std::uint64_t> lastCapacity;
		};

		typedef std::unique_ptr<ProcessorState>(*ParameterCreater)(std::size_t offset, bool createShortNames, SystemView system);

		enum class OscChannels
		{
			/// <summary>
			/// Only the left channel will be analyzed and displayed.
			/// </summary>
			Left,
			/// <summary>
			/// Only the right channel will be analyzed and displayed.
			/// </summary>
			Right,
			/// <summary>
			/// Left and right will be merged (added) together and processed
			/// in mono mode, equivalent to mid in M/S processing
			/// </summary>
			Merge,
			Mid = Merge,
			/// <summary>
			/// The difference between the two channels will be processed.
			/// Equal to left - right, and is equivalent to side in M/S processing
			/// </summary>
			Side,
			/// <summary>
			/// If the configuration is over this value,
			/// the processing requires more than one channel.
			/// </summary>
			OffsetForMono = Side,
			/// <summary>
			/// Both channels are displayed.
			/// </summary>
			Separate,
			/// <summary>
			/// First channel is mid, second channel is side.
			/// </summary>
			MidSide,
			End
		};

		enum class SpectrumChannels
		{
			/// <summary>
			/// Only the left channel will be analyzed and displayed.
			/// </summary>
			Left,
			/// <summary>
			/// Only the right channel will be analyzed and displayed.
			/// </summary>
			Right,
			/// <summary>
			/// Left and right will be merged (added) together and processed
			/// in mono mode, equivalent to mid in M/S processing
			/// </summary>
			Merge,
			Mid = Merge,
			/// <summary>
			/// The difference between the two channels will be processed.
			/// Equal to left - right, and is equivalent to side in M/S processing
			/// </summary>
			Side,
			/// <summary>
			/// If the configuration is over this value,
			/// the processing requires more than one channel.
			/// </summary>
			OffsetForMono = Side,
			/// <summary>
			/// Both channels will be processed seperately, and the average of the magnitude will be displayed,
			/// together with a graph of the scaled phase cancellation.
			/// </summary>
			Phase,
			/// <summary>
			/// Both channels are displayed.
			/// </summary>
			Separate,
			/// <summary>
			/// First channel is mid, second channel is side.
			/// </summary>
			MidSide,
			/// <summary>
			/// Channels 1 and 2 are interpreted as a complex sequence of real and imaginary numbers
			/// </summary>
			Complex,
			End
		};

		template<typename Scalar>
		union UComplexFilter
		{
			UComplexFilter() : real(0), imag(0) { }

			UComplexFilter(const std::complex<Scalar> & c)
				: real(c.real()), imag(c.imag())
			{
			}

			UComplexFilter & operator = (const std::complex<Scalar> & c) noexcept
			{
				real = c.real();
				imag = c.imag();
				return *this;
			}

			struct
			{
				Scalar real, imag;
			};
			struct
			{
				Scalar magnitude, phase;
			};
			struct
			{
				Scalar leftMagnitude, rightMagnitude;
			};

			UComplexFilter operator * (Scalar left) const noexcept
			{
				UComplexFilter ret;
				ret.real = real * left;
				ret.imag = imag * left;
				return ret;
			}

			UComplexFilter operator + (const UComplexFilter & left) const noexcept
			{
				UComplexFilter ret;
				ret.real = left.real + real;
				ret.imag = left.imag + imag;
				return ret;
			}

			operator std::complex<Scalar>() const noexcept
			{
				return std::complex<Scalar>(real, imag);
			}
		};

		template<typename T>
		struct UniqueHandle
		{
		public:

			UniqueHandle(std::unique_ptr<T> && ref)
			{
				handle.reset(ref.release());
				handle.get_deleter().doDelete = true;
			}

			UniqueHandle<T> & operator = (std::unique_ptr<T> && ref)
			{
				handle.reset(ref.release());
				handle.get_deleter().doDelete = true;
				ref.get_deleter().doDelete = false;
				return *this;
			}

			UniqueHandle<T> weakCopy()
			{
				UniqueHandle<T> ret;
				ret.handle.reset(handle.get());
				ret.handle.get_deleter().doDelete = false;
				return std::move(ret);
			}

			T * operator -> () const noexcept { return *get(); }
			T * operator -> () noexcept { return *get(); }

			T * get() noexcept { return handle.get(); }
			T * get() const noexcept { return handle.get(); }

			/// <summary>
			/// Transforms ownership of the contained object. Throws if it isn't owned.
			/// </summary>
			T * acquire()
			{
				if (!handle.get_deleter().doDelete)
					CPL_RUNTIME_EXCEPTION("UniqueHandle asked to release ownership of something it doesn't own");
				handle.get_deleter().doDelete = false;
				return handle.release();
			}

			/// <summary>
			/// Deletes the object if owned, and clears afterwards.
			/// </summary>
			void forget() noexcept
			{
				handle = nullptr;
			}

			/// <summary>
			/// May leak memory. Removes reference to any object, but doesn't delete it.
			/// </summary>
			void clear() noexcept
			{
				handle.release();
			}

			static UniqueHandle<T> null() { return {}; }

		private:

			UniqueHandle()
			{
				handle.get_deleter().doDelete = false;
			}

			std::unique_ptr<T, cpl::Utility::MaybeDelete<T>> handle;
		};



		/// <summary>
		/// Provides a optionally lazily loaded instance of some object, where you can (de)serialize the state independently of the instance
		/// </summary>
		template<typename T>
		class DecoupledStateObject
		{
			static_assert(std::is_base_of<cpl::DestructionNotifier, T>::value, "State object must be able to notify of destruction");

		public:

			typedef std::function<void(T &, cpl::CSerializer::Archiver &, cpl::Version)> FSerializer;
			typedef std::function<void(T &, cpl::CSerializer::Builder &, cpl::Version)> FDeserializer;
			typedef std::function<std::unique_ptr<T>()> FGenerator;

			DecoupledStateObject(FGenerator generatorFunc, FSerializer serializerFunc, FDeserializer deserializerFunc)
				: generator(generatorFunc)
				, serializer(serializerFunc)
				, deserializer(deserializerFunc)
				, objectDeathHook(*this)
				, cachedObject(UniqueHandle<T>::null())
			{

			}

			FGenerator replaceGenerator(FGenerator generatorFunc)
			{
				auto old = generator;
				generator = generatorFunc;
				return old;
			}

			FSerializer replaceSerializer(FSerializer serializerFunc)
			{
				auto old = serializer;
				serializer = serializerFunc;
				return old;
			}

			FDeserializer replaceDeserializer(FDeserializer deserializerFunc)
			{

				auto old = deserializer;
				deserializer = deserializerFunc;
				return old;
			}

			UniqueHandle<T> getUnique()
			{
				UniqueHandle<T> ret = hasCached() ? std::move(cachedObject) : create();
				cachedObject = ret.weakCopy();
				return std::move(ret);
			}

			UniqueHandle<T> getCached() { if (!hasCached()) cachedObject = create(); return cachedObject.weakCopy(); }
			bool hasCached() const noexcept { return cachedObject.get() != nullptr; }

			void setState(cpl::CSerializer::Builder & builder, cpl::Version v)
			{
				if (hasCached())
				{
					builder.setMasterVersion(v);
					deserializeState(*getCached().get(), &builder);
				}
				else
				{
					state = builder;
					state.setMasterVersion(v);
				}
			}

			/// <summary>
			/// If serialized state is refreshed, the version will be updated
			/// </summary>
			/// <returns></returns>
			const cpl::CSerializer::Builder & getState()
			{
				if (hasCached())
				{
					serializeState(*getCached().get());
				}
				return state;
			}

		private:

			void onObjectDestruction(cpl::DestructionNotifier * obj)
			{
				CPL_RUNTIME_ASSERTION(obj);
				CPL_RUNTIME_ASSERTION(obj == cachedObject.get());

				serializeState(*getCached().get());

				cachedObject.clear();
			}

			void serializeState(T & obj)
			{
				state.clear();
				state.setMasterVersion(cpl::programInfo.version);
				serializer(obj, state, cpl::programInfo.version);
			}

			void deserializeState(T & obj, cpl::CSerializer::Builder * optionalExternalState = nullptr)
			{
				auto & usableState = optionalExternalState ? *optionalExternalState : state;
				deserializer(obj, usableState, usableState.getLocalVersion());
			}

			UniqueHandle<T> create()
			{
				auto uref = generator();
				if (!state.isEmpty())
					deserializeState(*uref);
				objectDeathHook.listenToObject(uref.get());
				return std::move(uref);
			}

			class DestructionDelegate
				: public cpl::DestructionNotifier::EventListener
			{
			public:
				DestructionDelegate(DecoupledStateObject<T> & parentToNotify) : parent(parentToNotify) { }

				void listenToObject(cpl::DestructionNotifier * notifier)
				{
					notif = notifier;
					notifier->addEventListener(this);
				}

				virtual void onServerDestruction(cpl::DestructionNotifier * v) override
				{
					parent.onObjectDestruction(v);
					hasDied = true;
				}

				~DestructionDelegate()
				{
					if (!hasDied && notif)
					{
						notif->removeEventListener(this);
					}
				}

			private:
				bool hasDied = false;
				DecoupledStateObject<T> & parent;
				cpl::DestructionNotifier * notif = nullptr;
			};

			FGenerator generator;
			FSerializer serializer;
			FDeserializer deserializer;


			cpl::CSerializer state;
			UniqueHandle<T> cachedObject;
			DestructionDelegate objectDeathHook;
		};

		template<typename Parent>
			class SSOSurrogate : public cpl::SafeSerializableObject
			{
			public:
				typedef std::function<void(Parent &, cpl::CSerializer::Archiver &, cpl::Version)> FSerializer;
				typedef std::function<void(Parent &, cpl::CSerializer::Builder &, cpl::Version)> FDeserializer;

				SSOSurrogate(Parent & parent, FSerializer serializer, FDeserializer deserializer)
					: parent(parent), serializer(serializer), deserializer(deserializer) {}
			private:
				void serialize(cpl::CSerializer::Archiver & ar, cpl::Version version) override
				{
					serializer(parent, ar, version);
				}

				void deserialize(cpl::CSerializer::Builder & b, cpl::Version version) override
				{
					deserializer(parent, b, version);
				}

				Parent & parent;

				FSerializer serializer;
				FDeserializer deserializer;
			};

		struct ParameterMap
		{
			void insert(std::pair<std::string, std::unique_ptr<ProcessorState>> entry)
			{
				parameterSets.push_back(&entry.second->getParameterSet());
				map.emplace_back(std::move(entry));
			}

			ParameterSet::ParameterView * findParameter(cpl::Parameters::Handle index)
			{
				for (std::size_t i = 0; i < map.size(); ++i)
					if (auto param = parameterSets[i]->findParameter(index))
						return param;

				CPL_RUNTIME_EXCEPTION("Parameter index from host is out of bounds");
			}

			ParameterSet * getSet(const std::string & s) noexcept
			{
				for (std::size_t i = 0; i < map.size(); ++i)
				{
					if (map[i].first == s)
						return parameterSets[i];
				}

				return nullptr;
			}

			ParameterSet * getSet(std::size_t i) noexcept
			{
				return parameterSets[i];
			}

			ProcessorState * getState(std::size_t i) noexcept
			{
				return map[i].second.get();
			}


			ProcessorState * getState(const std::string & s) noexcept
			{
				for (std::size_t i = 0; i < map.size(); ++i)
				{
					if (map[i].first == s)
						return map[i].second.get();
				}

				return nullptr;
			}

			std::size_t numSetsAndState() const noexcept
			{
				return parameterSets.size();
			}

			std::size_t numParams() const noexcept 
			{
				std::size_t r(0);
				for (auto & i : parameterSets)
					r += i->size();
				return r;
			}

		private:

			std::vector<ParameterSet *> parameterSets;
			std::vector<std::pair<std::string, std::unique_ptr<ProcessorState>>> map;
		};

		template<typename NameVector, typename ColourVector>
		inline void PaintLegend(juce::Graphics& g, juce::Colour front, juce::Colour back, juce::Point<float> position, const NameVector& names, const ColourVector& colours, std::size_t count)
		{
			juce::GlyphArrangement arrangement;
			constexpr auto offset = 5;
			constexpr auto strokeSize = 50;
			position.y += offset * 3;
			position.x += offset;
			auto startingY = position.y;
			auto lineHeight = g.getCurrentFont().getHeight();

			for (std::size_t i = 0; i < count; ++i)
			{
				arrangement.addLineOfText(g.getCurrentFont(), names[i], position.x, position.y);
				position.y += offset + g.getCurrentFont().getHeight();
			}

			auto bounds = arrangement.getBoundingBox(0, -1, true).reduced(-offset).withTrimmedRight(-strokeSize);

			g.setColour(back);
			g.fillRoundedRectangle(bounds, offset);
			g.setColour(front);
			g.drawRoundedRectangle(bounds, offset, 1);

			arrangement.draw(g);


			for (std::size_t i = 0; i < count; ++i)
			{
				auto y = startingY + i * (offset + lineHeight) - lineHeight * 0.33f;
				g.setColour(colours[i]);
				g.drawLine(bounds.getRight() - strokeSize, y, bounds.getRight() - offset, y, 2.0f);
			}
		}

	};
	namespace std
	{
		template<typename T>
			inline T abs(const Signalizer::UComplexFilter<T> & f)
			{
				return sqrt(f.real * f.real + f.imag * f.imag);
			}

	}
#endif
