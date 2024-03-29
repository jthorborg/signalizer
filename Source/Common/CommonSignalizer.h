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
	#include <cpl/lib/weak_atomic.h>
	#include "ConcurrentConfig.h"
	#include <mutex>
	#include <cpl/AudioStream.h>
	#include <cpl/gui/CViews.h>
	#include <chrono>
	#include <variant>
	#include <cpl/Exceptions.h>
	#include <functional>

	namespace cpl
	{
		class CSubView;
	}

	namespace Signalizer
	{
		typedef std::make_signed<std::size_t>::type ssize_t;

		template<class T>
		using vector_set = std::set<T>;

		typedef std::int32_t PinInt;

		struct DirectedPortPair
		{
			PinInt Source;
			PinInt Destination;

			inline friend bool operator < (DirectedPortPair a, DirectedPortPair b)
			{
				return std::tie(a.Destination, a.Source) < std::tie(b.Destination, b.Source);
			}
		};

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
		class SystemView;
		class SharedBehaviour;
		struct ConcurrentConfig;
		
		class ProcessorState
			: public cpl::SafeSerializableObject
		{
		public:

			virtual std::unique_ptr<StateEditor> createEditor() = 0;
			virtual std::unique_ptr<cpl::CSubView> createView(
				std::shared_ptr<const SharedBehaviour> globalBehaviour,
				std::shared_ptr<const ConcurrentConfig> config,
				std::shared_ptr<AudioStream::Output>& stream
			) = 0;

			virtual ParameterSet& getParameterSet() = 0;
			virtual const char* getName() = 0;
			virtual void onPresentationStreamCreated(std::shared_ptr<AudioStream::Output>& output) {}

			virtual ~ProcessorState() {}
		};

		template<typename Derived>
		class ProcessorStreamState
			: public ProcessorState
			, public AudioStream::Listener
		{
		public:

			void onPresentationStreamCreated(std::shared_ptr<AudioStream::Output>& output) override
			{ 
				auto base = static_cast<Derived*>(this)->shared_from_this();

				if (auto old = previousOutput.lock())
					old->removeListener(base);

				output->addListener(base);
			}

		protected:
			std::weak_ptr<AudioStream::Output> previousOutput;
		};

		class GraphicsWindow : public cpl::COpenGLView
		{
		protected:

			GraphicsWindow(std::string name) 
				: COpenGLView(std::move(name)) 
			{
			
			}

			struct CurrentMouse
			{
				void setFromPoint(juce::Point<float> mousePosition)
				{
					x = mousePosition.x;
					y = mousePosition.y;
				}

				juce::Point<float> getPoint()
				{
					return juce::Point<float>(x, y);
				}

				cpl::relaxed_atomic<float>
					x, y;
			} originMouse, currentMouse;

			cpl::relaxed_atomic<bool> isMouseInside;
			juce::MouseCursor displayCursor;

			// guis and whatnot
			cpl::CBoxFilter<float, 64> avgDelta;
			cpl::CBoxFilter<float, 64> avgFrame;


			// Mouse overrides
			void mouseMove(const juce::MouseEvent& event) override
			{
				// TODO: implement beginChangeGesture()
				currentMouse.setFromPoint(event.position);
			}

			void mouseDrag(const juce::MouseEvent& event) override
			{
				// TODO: implement beginChangeGesture()
				currentMouse.setFromPoint(event.position);
			}

			void mouseDown(const juce::MouseEvent& event) override
			{
				// TODO: implement endChangeGesture()
				originMouse.setFromPoint(event.position);
			}

			void mouseExit(const juce::MouseEvent& e)  override
			{
				isMouseInside = false;
			}

			void mouseEnter(const juce::MouseEvent& e) override
			{
				isMouseInside = true;
			}

			// Graphics overrides
			void onGraphicsRendering(juce::Graphics& g) override
			{
				// do software rendering
				if (!isOpenGL())
				{
					g.fillAll(juce::Colours::black);
					g.setColour(juce::Colours::white);
					g.drawText("Enable OpenGL in settings to use the " + CView::getName(), getLocalBounds(), juce::Justification::centred);
				}
			}

			void postFrame()
			{
				avgDelta.setNext(openGlEndToEndTime());
				avgFrame.setNext(openGLFrameTime());
			}

			double getAverageFrameTime()
			{
				return avgFrame.getAverage();
			}

			double getAverageDelta()
			{
				return avgDelta.getAverage();
			}

			void computeAverageStats(float& fps, float& usagePercent)
			{
				const auto frame = avgFrame.getAverage();
				const auto delta = avgDelta.getAverage();

				usagePercent = 100 * (frame / delta);
				fps = 1.0f / delta;
			}
		};

		class SystemView
		{
		public:

			SystemView(std::shared_ptr<const ConcurrentConfig> config, ParameterSet::AutomatedProcessor& automatedProcessor)
				: config(config), processor(automatedProcessor)
			{

			}

			ParameterSet::AutomatedProcessor& getProcessor() noexcept { return processor; }
			std::shared_ptr<const ConcurrentConfig>& getConcurrentConfig() noexcept { return config; }

		private:

			std::shared_ptr<const ConcurrentConfig> config;			
			ParameterSet::AutomatedProcessor& processor;
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

			AudioHistoryTransformatter(std::shared_ptr<const ConcurrentConfig>& config, Mode mode = Milliseconds)
				: param(nullptr)
				, m(mode)
				, lastCapacity(config->historyCapacity)
				, sampleRate(config->sampleRate)
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
				// TODO: This can leave scope before us.
				param = &view;
			}

			void serialize(cpl::CSerializer::Archiver & ar, cpl::Version v) override
			{
				std::uint64_t val = lastCapacity;
				ar << val;
			}

			void deserialize(cpl::CSerializer::Builder & ar, cpl::Version v) override
			{
				std::uint64_t val;
				ar >> val;
				lastCapacity = val;
			}

			void onStreamPropertiesChanged(const AudioStream::ListenerContext& changedSource, const AudioStream::AudioStreamInfo& before)
			{
				// TODO: what if oldCapacity == 0?
				const auto oldFraction = param->getValueNormalized();
				auto oldCapacity = lastCapacity;
				auto beforeCapacity = before.audioHistoryCapacity;
				if (oldCapacity == 0)
					oldCapacity = beforeCapacity;

				const auto& info = changedSource.getInfo();
				const auto newCapacity = info.audioHistoryCapacity;
				sampleRate = info.sampleRate;

				if (newCapacity > 0)
					lastCapacity = newCapacity;

				if (oldCapacity == 0 || newCapacity == 0)
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
			
		protected:

			virtual bool format(const ValueType & val, std::string & buf) override
			{
				char buffer[100];

				if (m == Milliseconds)
				{
					cpl::sprintfs(buffer, "%.2f ms", 1000 * val / sampleRate);
				}
				else
				{
					cpl::sprintfs(buffer, "%.0f smps", val);
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
						collectedValue *= sampleRate;
					}
					else
					{
						// assume value is in miliseconds
						if (m == Milliseconds && notSamples)
						{
							collectedValue /= 1000;
							collectedValue *= sampleRate;
						}
					}

					val = collectedValue;
					return true;

				}

				return false;
			}

			virtual ValueType transform(ValueType val) const noexcept override
			{
				if (lastCapacity == 0)
				{
					// Don't really want to transform something if we don't know the capacity.
					CPL_BREAKIFDEBUGGED();
				}

				auto samples = std::round(val * lastCapacity);

				/* if (m == Milliseconds)
				{
					samples /= stream.getInfo().sampleRate.load(std::memory_order_relaxed);
					samples *= 1000;
				} */


				return samples;
			}


			virtual ValueType normalize(ValueType val) const noexcept override
			{
				if (lastCapacity == 0)
				{
					// Don't really want to normalize something if we don't know the capacity.
					CPL_BREAKIFDEBUGGED();
				}

				/* if (m == Milliseconds)
				{
					val /= stream.getInfo().sampleRate.load(std::memory_order_relaxed);
					val *= 1000;
				} */

				return val / lastCapacity;
			}

			ParameterView* param;
			Mode m;
			/// <summary>
			/// Represents the last capacity change we were informated about,
			/// and thus currently represents the magnitude scale of the fraction.
			/// </summary>
			cpl::relaxed_atomic<std::uint64_t> lastCapacity;
			cpl::relaxed_atomic<double> sampleRate;
		};

		typedef std::shared_ptr<ProcessorState>(*ContentCreater)(std::size_t offset, SystemView& system);

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
			void insert(std::pair<std::string, std::shared_ptr<ProcessorState>> entry)
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

			// TODO: Move these into private access?
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

			std::shared_ptr<ProcessorState>& getState(std::size_t i) noexcept
			{
				return map[i].second;
			}

			std::shared_ptr<ProcessorState>& getState(const std::string & s) noexcept
			{
				for (std::size_t i = 0; i < map.size(); ++i)
				{
					if (map[i].first == s)
						return map[i].second;
				}

				CPL_UNREACHABLE();
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
			std::vector<std::pair<std::string, std::shared_ptr<ProcessorState>>> map;
		};

		struct ColourRotation
		{
			ColourRotation() : base(), size(), stereo() {}

			ColourRotation(juce::Colour base, std::size_t size, bool stereo)
				: base(base), size(size), stereo(stereo)
			{

			}

			juce::Colour operator [](std::size_t index) const noexcept
			{
				if (stereo)
					index &= ~0x1ull;

				return base.withRotatedHue(index / size);
			}

			juce::Colour getBase() const noexcept
			{
				return base;
			}

			friend bool operator != (const ColourRotation& a, const ColourRotation& b)
			{
				return a.base != b.base || a.size != b.size || a.stereo != b.stereo;
			}


		private:
			juce::Colour base;
			float size;
			bool stereo;
		};

		/// <summary>
		/// Not atomic.
		/// </summary>
		struct ChangeVersion
		{
			struct Listener
			{
			public:
				bool consumeChanges(const ChangeVersion& v)
				{
					CPL_RUNTIME_ASSERTION(version <= v.version);

					if (v.version > version)
					{
						version = v.version;
						return true;
					}

					return false;
				}

			private:
				int version {};
			};

			void bump() 
			{
				version++;
			}

		private:
			int version{};
		};

		struct LegendCache
		{
			static constexpr int offset = 5;
			static constexpr int strokeSize = 50;

			typedef std::variant<juce::Colour, std::pair<juce::Colour, juce::Colour>, juce::ColourGradient> ColourItem;

			void setFont(juce::Font fontToUse)
			{
				font = fontToUse;
			}

			void reset(juce::Point<float> positionToUse)
			{
				position = positionToUse;

				position.y += offset * 3;
				position.x += offset;
				startingY = position.y;

				arrangement.clear();
				colours.clear();
			}

			void addLine(const juce::String& text, ColourItem&& item)
			{
				arrangement.addLineOfText(font, text, position.x, position.y);
				colours.emplace_back(std::move(item));
				position.y += offset + font.getHeight();
			}

			void addLine(const juce::String& text, juce::Colour colour)
			{
				addLine(text, ColourItem(colour));
			}

			void addLine(const juce::String& text, juce::Colour a, juce::Colour b)
			{
				addLine(text, ColourItem(std::make_pair(a, b)));
			}

			void addLine(const juce::String& text, juce::ColourGradient&& gradient)
			{
				addLine(text, ColourItem(std::move(gradient)));
			}

			void paint(juce::Graphics& g, juce::Colour front, juce::Colour back)
			{
				auto bounds = arrangement.getBoundingBox(0, -1, true).reduced(-offset).withTrimmedRight(-strokeSize);
				auto lineHeight = font.getHeight();

				g.setColour(back);
				g.fillRoundedRectangle(bounds, offset);
				g.setColour(front);
				g.drawRoundedRectangle(bounds, offset, 1);

				arrangement.draw(g);

				for (std::size_t i = 0; i < colours.size(); ++i)
				{
					auto y = startingY + i * (offset + lineHeight) - lineHeight * 0.33f;

					if (auto arg = std::get_if<juce::Colour>(&colours[i]))
					{
						g.setColour(*arg);
						g.drawLine(bounds.getRight() - strokeSize, y, bounds.getRight() - offset, y, 3.0f);
					}
					else if (auto arg = std::get_if<std::pair<juce::Colour, juce::Colour>>(&colours[i]))
					{
						g.setColour(arg->first);
						g.drawLine(bounds.getRight() - strokeSize, y, bounds.getRight() - offset - strokeSize / 2, y, 3.0f);
						g.setColour(arg->second);
						g.drawLine(bounds.getRight() - strokeSize / 2, y, bounds.getRight() - offset, y, 3.0f);
					}
					else if (auto arg = std::get_if<juce::ColourGradient>(&colours[i]))
					{
						arg->point1 = { bounds.getRight() - strokeSize, y };
						arg->point2 = { bounds.getRight() - offset, y };
						g.setGradientFill(*arg);
						g.drawLine(bounds.getRight() - strokeSize, y, bounds.getRight() - offset, y, 3.0f);
					}
				}
			}

		private:

			juce::GlyphArrangement arrangement;
			std::vector<ColourItem> colours;
			juce::Font font;
			juce::Point<float> position;
			float startingY;
		};

		template<typename T, class Mutex = std::mutex>
		class CriticalSection
		{
		public:

			class Access
			{
				friend class CriticalSection<T, Mutex>;

			public:
				
				T* operator -> () noexcept
				{
					return &data;
				}

				T& operator * () noexcept
				{
					return data;
				}

				Access(const Access& other) = delete;
				Access& operator =(const Access& other) = delete;

				Access(Access&& other) = default;
				Access& operator =(Access&& other) = default;

			private:

				Access(CriticalSection<T, Mutex>& data)
					: data(data.data)
					, lock(data.mutex)
				{

				}

				T& data;
				std::lock_guard<Mutex> lock;
			};

			Access lock()
			{
				return { *this };
			}

			template <typename... Args>
			CriticalSection(Args&& ...args)
				: data{ std::forward<Args>(args)... }
			{
			}

		private:
			T data;
			Mutex mutex;
		};

		class FloatColour : public std::array<float, 3>
		{
		public:
			FloatColour(const juce::Colour& colour)
			{
				this->operator[](0) = colour.getFloatRed();
				this->operator[](1) = colour.getFloatGreen();
				this->operator[](2) = colour.getFloatBlue();
			}

			FloatColour()
			{
				this->fill(0);
			}

			juce::Colour toJuceColour() const noexcept
			{
				return juce::Colour::fromFloatRGBA(
					this->operator[](0),
					this->operator[](1),
					this->operator[](2),
					1
				);
			}
		};

		template<typename T>
		inline bool assignAndChanged(T& source, const T& operand)
		{
			if (!(source != operand))
				return false;

			source = operand;
			return true;
		}

#define NONTERMINAL_ASSUMPTION(condition) ((condition) || triggerNonTerminalAssumption(#condition, __FILE__, __LINE__, __func__))

		struct Assumptions
		{
			typedef std::vector<std::string> Messages;
			Messages current;
			std::set<std::size_t> set;

			Messages takeAllMessages()
			{
				Messages messages;
				std::swap(current, messages);

				return messages;
			}
		};

		CriticalSection<Assumptions>& getFailedAssumptions();
		bool triggerNonTerminalAssumption(const char* assumption, const char* file, const int line, const char* function);
		void revealExceptionLog();
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
