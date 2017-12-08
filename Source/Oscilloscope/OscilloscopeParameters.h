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

	file:Oscilloscope.h

		Interface for the oscilloscope parameters

*************************************************************************************/

#ifndef SIGNALIZER_COSCILLOSCOPEPARAMETERS_H
	#define SIGNALIZER_COSCILLOSCOPEPARAMETERS_H

	#include "Signalizer.h"

	namespace Signalizer
	{

		class OscilloscopeContent final
			: public cpl::Parameters::UserContent
			, public ProcessorState
			, public ParameterSet::UIListener
		{
		public:

			static constexpr std::size_t LookaheadSize = 8192;
			static constexpr std::size_t InterpolationKernelSize = 10;

			enum class TriggeringMode
			{
				None,
				Spectral,
				Window,
				EnvelopeHold,
				ZeroCrossing,
				end
			};

			enum class TimeMode
			{
				Time, Cycles, Beats
			};

			enum class ColourMode
			{
				Static, SpectralEnergy
			};

			template<typename ParameterView>
			class WindowSizeTransformatter : public AudioHistoryTransformatter<ParameterView>
			{
			public:

			    typedef typename AudioHistoryTransformatter<ParameterView>::Mode Mode;
			    typedef typename AudioHistoryTransformatter<ParameterView>::Stream Stream;
			    typedef typename AudioHistoryTransformatter<ParameterView>::ValueType ValueType;
			    typedef typename AudioHistoryTransformatter<ParameterView>::Scaling Scaling;

				WindowSizeTransformatter(AudioStream & audioStream, std::size_t auxLookahead, Mode mode = Mode::Milliseconds)
					: AudioHistoryTransformatter<ParameterView>(audioStream, mode)
					, lookahead(auxLookahead)
					, timeMode(TimeMode::Time)
				{

				}

				void setTimeModeFromUI(TimeMode newMode)
				{
					timeMode = newMode;
					this->param->updateFromUINormalized(this->param->getValueNormalized());
				}

			private:

				virtual void onAsyncChangedProperties(const Stream & changedSource, const typename Stream::AudioStreamInfo & before) override
				{
					// TODO: what if oldCapacity == 0?
					const auto oldFraction = this->param->getValueNormalized();
					auto oldCapacity = this->lastCapacity.load(std::memory_order_relaxed);
					auto beforeCapacity = before.audioHistoryCapacity.load(std::memory_order_acquire);
					if (oldCapacity == 0)
						oldCapacity = beforeCapacity;

					const auto newCapacity = changedSource.getInfo().audioHistoryCapacity.load(std::memory_order_relaxed);

					if (newCapacity > 0)
						this->lastCapacity.store(newCapacity, std::memory_order_relaxed);

					if (oldCapacity == 0 || newCapacity == 0)
					{
						this->param->updateFromProcessorNormalized(oldFraction, cpl::Parameters::UpdateFlags::All & ~cpl::Parameters::UpdateFlags::RealTimeSubSystem);
					}
					else
					{
						const auto sampleSizeBefore = oldCapacity * oldFraction;
						const auto newFraction = sampleSizeBefore / newCapacity;
						if (oldFraction != newFraction || beforeCapacity == 0)
							this->param->updateFromProcessorNormalized(newFraction, cpl::Parameters::UpdateFlags::All & ~cpl::Parameters::UpdateFlags::RealTimeSubSystem);
					}
				}

				virtual bool format(const ValueType & val, std::string & buf) override
				{
					char buffer[100];
					switch (timeMode)
					{
						case TimeMode::Cycles:
						{
							sprintf_s(buffer, u8"%.2f (%.2f r)", val, cpl::simd::consts<ValueType>::tau * val);
							buf = buffer;
							return true;
						}
						case TimeMode::Beats:
						{
							sprintf_s(buffer, "1/%.0f", val);
							buf = buffer;
							return true;
						}
						default: case TimeMode::Time: return AudioHistoryTransformatter<ParameterView>::format(val, buf);
					}
				}

				virtual bool interpret(const cpl::string_ref buf, ValueType & val) override
				{
					ValueType collectedValue;

					if (cpl::lexicalConversion(buf, collectedValue))
					{
						bool notSamples = true;
						if (timeMode != TimeMode::Time)
						{
							if (timeMode == TimeMode::Cycles)
							{
								if (buf.find("r") != std::string::npos)
								{
									collectedValue /= cpl::simd::consts<ValueType>::tau;
								}
							}
							else if (timeMode == TimeMode::Beats)
							{
								std::size_t numeratorEnd = 0;
								if ((numeratorEnd = buf.find('/')) != std::string::npos)
								{

									ValueType numerator, denominator;

									if (cpl::lexicalConversion(std::string(buf.begin() + numeratorEnd + 1, buf.end()), denominator)
										&& cpl::lexicalConversion(std::string(buf.begin(), buf.begin() + numeratorEnd), numerator))
									{
										collectedValue = numerator / denominator;
									}
									else
									{
										return false;
									}

								}

								if (buf.find("bars") != std::string::npos)
								{
									collectedValue /= 4;
								}
								// stored as reciprocal.
								collectedValue = 1 / collectedValue;
							}
							val = collectedValue;
							return true;
						}
						else if (buf.find("s") != std::string::npos && (notSamples = buf.find("smps") == std::string::npos))
						{
							if (buf.find("ms") != std::string::npos)
							{
								collectedValue /= 1000;
							}
							collectedValue *= this->stream.getInfo().sampleRate.load(std::memory_order_relaxed);
						}
						else
						{
							// assume value is in miliseconds
							if (this->m == Mode::Milliseconds && notSamples)
							{
								collectedValue /= 1000;
								collectedValue *= this->stream.getInfo().sampleRate.load(std::memory_order_relaxed);
							}
						}

						val = collectedValue;
						return true;

					}

					return false;

				}

				virtual ValueType transform(ValueType val) const noexcept override
				{
					switch (timeMode)
					{
						case TimeMode::Cycles:
						{
							return cpl::Math::UnityScale::exp<ValueType>(val, 1, 32);
						}
						case TimeMode::Beats:
						{
							return cpl::Math::nextPow2Inc(cpl::Math::round<std::size_t>(cpl::Math::UnityScale::exp<ValueType>(1 - val, 1, 32)));
						}
						default: case TimeMode::Time:
						{
							const auto minExponential = 100;
							const auto capacity = this->stream.getAudioHistoryCapacity();

							const auto top = capacity;
							const auto expSamples = cpl::Math::UnityScale::exp<ValueType>(val, minExponential, top);
							const auto rescaled = cpl::Math::UnityScale::linear<ValueType>(cpl::Math::UnityScale::Inv::linear<ValueType>(expSamples, minExponential, top), 2, top);
							return rescaled;
						}
					}

				}


				virtual ValueType normalize(ValueType val) const noexcept override
				{
					switch (timeMode)
					{
						case TimeMode::Cycles:
						{
							return cpl::Math::UnityScale::Inv::exp<ValueType>(val, 1, 32);
						}
						case TimeMode::Beats:
						{
							return 1 - cpl::Math::UnityScale::Inv::exp<ValueType>(val, 1, 32);
						}
						default: case TimeMode::Time:
						{
							const auto minExponential = 100;
							const auto capacity = this->stream.getAudioHistoryCapacity();
							const auto top = capacity;
							const auto linear = cpl::Math::UnityScale::Inv::linear<ValueType>(val, 2, top);
							const auto expSamples = cpl::Math::UnityScale::linear<ValueType>(linear, minExponential, top);

							const auto normalized = cpl::Math::UnityScale::Inv::exp<ValueType>(expSamples, minExponential, top);
							return normalized;
						}
					}
				}

				std::atomic<Scaling> scale;
				std::size_t lookahead;
				TimeMode timeMode;
			};

			template<typename ParameterView>
			class LinearHzFormatter : public ParameterView::ParameterType::Formatter
			{
			public:

				typedef typename ParameterView::ParameterType::ValueType ValueType;

				LinearHzFormatter(AudioStream & as)
					: stream(as)
				{
					setTuningFromA4();
				}

				void setTuningFromA4(double hz = 440)
				{
					a4InHz = hz;
				}

				virtual bool format(const ValueType & val, std::string & buf) override
				{
					char buffer[100];
					sprintf_s(buffer, "%.5f Hz", val);
					buf = buffer;
					return true;
				}

				virtual bool interpret(const cpl::string_ref buf, ValueType & val) override
				{
					ValueType contained;

					// try to parse it as a note firstly:
					int octave;
					char tone;
					char hash;
					static const std::pair<char, int> notes[] = { { 'a', 0 },{ 'b', 2 },{ 'c', -9 },{ 'd', -7 },{ 'e', -5 },{ 'f', -4 },{ 'g', -2 } };

					if (std::sscanf(buf.c_str(), "%c%d", &tone, &octave) == 2 && !isdigit(tone))
					{
						tone = tolower(tone);
						if (tone >= 'a' && tone <= 'g')
						{
							auto offset = notes[tone - 'a'].second;
							auto note = octave * 12 + offset;
							val = a4InHz * std::pow(2, (note - 48) / 12.0);
							return true;
						}
						else
						{
							return false;
						}
					}
					else if (std::sscanf(buf.c_str(), "%c%c%d", &tone, &hash, &octave) == 3 && !isdigit(tone))
					{
						tone = tolower(tone);
						if (tone >= 'a' && tone <= 'g')
						{
							auto offset = notes[tone - 'a'].second;
							if (hash == '#') offset++;
							else if (tolower(hash) == 'b') offset--;
							auto note = octave * 12 + offset;
							val = a4InHz * std::pow(2, (note - 48) / 12.0);
							return true;
						}
						else
						{
							return false;
						}
					}
					else if (cpl::lexicalConversion(buf, contained))
					{

						if (buf.find("smps") != std::string::npos)
						{
							contained = stream.getAudioHistorySamplerate() / contained;
						}
						else if (buf.find("ms") != std::string::npos)
						{
							contained = 1.0 / (contained / 1000);
						}
						else if (buf.find("r") != std::string::npos)
						{
							contained = (contained / (2 * cpl::simd::consts<ValueType>::pi)) * stream.getAudioHistorySamplerate();
						}
						else if (buf.find("b") != std::string::npos)
						{
							// TODO: Tecnically illegal to acquire the async playhead here -
							contained = (contained * stream.getASyncPlayhead().getBPM()) / 60;
						}

						val = contained;

						return true;
					}

					return false;
				}


			private:
				double a4InHz;
				const AudioStream & stream;
			};

			class OscilloscopeController
				: public CContentPage
				, public ParameterSet::UIListener
			{
			public:

				OscilloscopeController(OscilloscopeContent & parentValue)
					: parent(parentValue)
					, kantiAlias(&parentValue.antialias)
					, kdiagnostics(&parentValue.diagnostics)
					, kwindow(&parentValue.windowSize)
					, kgain(&parentValue.inputGain)
					, kprimitiveSize(&parentValue.primitiveSize)
					, kenvelopeSmooth(&parentValue.envelopeWindow)
					, kprimaryColour(&parentValue.primaryColour)
					, ksecondaryColour(&parentValue.secondaryColour)
					, kgraphColour(&parentValue.graphColour)
					, kbackgroundColour(&parentValue.backgroundColour)
					, klowColour(&parentValue.lowColour)
					, kmidColour(&parentValue.midColour)
					, khighColour(&parentValue.highColour)
					, ktransform(&parentValue.transform)
					, kenvelopeMode(&parentValue.autoGain.param)
					, kpresets(&valueSerializer, "oscilloscope")
					, ksubSampleInterpolationMode(&parentValue.subSampleInterpolation.param)
					, kpctForDivision(&parentValue.pctForDivision)
					, kchannelConfiguration(&parentValue.channelConfiguration.param)
					, ktriggerPhaseOffset(&parentValue.triggerPhaseOffset)
					, ktriggerMode(&parentValue.triggerMode.param)
					, ktimeMode(&parentValue.timeMode.param)
					, kdotSamples(&parentValue.dotSamples)
					, ktriggerOnCustomFrequency(&parentValue.triggerOnCustomFrequency)
					, kcustomFrequency(&parentValue.customTriggerFrequency)
					, koverlayChannels(&parentValue.overlayChannels)
					, kchannelColouring(&parentValue.channelColouring.param)
					, kcolourSmoothingTime(&parentValue.colourSmoothing)
					, kcursorTracker(&parentValue.cursorTracker)
					, ktrackerColour(&parentValue.trackerColour)
					, kfreqColourBlend(&parentValue.frequencyColouringBlend)
					, ktriggerHysteresis(&parentValue.triggerHysteresis)
					, ktriggerThreshold(&parentValue.triggerThreshold)

					, editorSerializer(
						*this,
						[](auto & oc, auto & se, auto version) { oc.serializeEditorSettings(se, version); },
						[](auto & oc, auto & se, auto version) { oc.deserializeEditorSettings(se, version); }
					)
					, valueSerializer(
						*this,
						[](auto & oc, auto & se, auto version) { oc.serializeAll(se, version); },
						[](auto & oc, auto & se, auto version) { oc.deserializeAll(se, version); }
					)
				{
					initControls();
					initUI();

					parent.getParameterSet().addUIListener(parent.timeMode.param.getParameterView().getHandle(), this);
					parent.getParameterSet().addUIListener(parent.triggerMode.param.getParameterView().getHandle(), this);

					enforceTriggeringModeCompability();
				}


				void parameterChangedUI(cpl::Parameters::Handle localHandle, cpl::Parameters::Handle globalHandle, ParameterSet::ParameterView * parameter) override
				{
					if(parameter == &parent.timeMode.param.getParameterView())
						enforceTriggeringModeCompability();
				}

				void enforceTriggeringModeCompability()
				{
					if (parent.timeMode.param.getAsTEnum<TimeMode>() == TimeMode::Cycles)
					{
						for (int i = 0; i < parent.triggerMode.tsf.getQuantization(); ++i)
						{
							ktriggerMode.setEnabledStateFor(i, i == cpl::enum_cast<int>(TimeMode::Cycles));
						}

						ktriggerMode.getValueReference().setTransformedValue(cpl::enum_cast<double>(TimeMode::Cycles));
					}
					else
					{
						for (int i = 0; i < parent.triggerMode.tsf.getQuantization(); ++i)
						{
							ktriggerMode.setEnabledStateFor(i, true);
						}
					}

				}

				cpl::SafeSerializableObject & getEditorSO() override { return editorSerializer; }

				~OscilloscopeController()
				{
					parent.getParameterSet().removeUIListener(parent.timeMode.param.getParameterView().getHandle(), this);
					parent.getParameterSet().removeUIListener(parent.triggerMode.param.getParameterView().getHandle(), this);
					notifyDestruction();
				}

				void initControls()
				{
					kwindow.bSetTitle("Window size");
					kgain.bSetTitle("Input gain");
					kgraphColour.bSetTitle("Graph colour");
					kbackgroundColour.bSetTitle("Backg. colour");
					kprimaryColour.bSetTitle("Primary colour");
					ksecondaryColour.bSetTitle("Secondary colour");
					klowColour.bSetTitle("Low band colour");
					kmidColour.bSetTitle("Mid band colour");
					khighColour.bSetTitle("High band colour");
					kprimitiveSize.bSetTitle("Primitive size");
					kenvelopeSmooth.bSetTitle("Env. window");
					ksubSampleInterpolationMode.bSetTitle("Sample interpolation");
					kpctForDivision.bSetTitle("Grid div. space");
					kchannelConfiguration.bSetTitle("Channel conf.");
					ktriggerMode.bSetTitle("Trigger mode");
					ktriggerPhaseOffset.bSetTitle("Trigger phase");
					ktimeMode.bSetTitle("Time mode");
					kcustomFrequency.bSetTitle("Custom trigger");
					kchannelColouring.bSetTitle("Channel colouring");
					kcolourSmoothingTime.bSetTitle("Colour smoothing");
					ktrackerColour.bSetTitle("Tracker colour");
					kfreqColourBlend.bSetTitle("Colour blend");
					ktriggerHysteresis.bSetTitle("Hysteresis");
					ktriggerThreshold.bSetTitle("Trigger thrshld");
					// buttons n controls

					kantiAlias.setSingleText("Antialias");
					kantiAlias.setToggleable(true);
					kdiagnostics.setSingleText("Diagnostics");
					kdiagnostics.setToggleable(true);
					kenvelopeMode.bSetTitle("Auto-gain mode");
					kdotSamples.bSetTitle("Dot samples");
					kdotSamples.setToggleable(true);
					ktriggerOnCustomFrequency.bSetTitle("Trigger on custom");
					ktriggerOnCustomFrequency.setToggleable(true);
					koverlayChannels.bSetTitle("Overlay channels");
					koverlayChannels.setToggleable(true);
					kcursorTracker.bSetTitle("Cursor tracker");
					kcursorTracker.setToggleable(true);


					// descriptions.
					kwindow.bSetDescription("The size of the displayed time window.");
					kgain.bSetDescription("How much the input (x,y) is scaled (or the input gain)" \
						" - additional transform that only affects the waveform, and not the graph");
					kantiAlias.bSetDescription("Antialiases rendering (if set - see global settings for amount). May slow down rendering.");
					kprimaryColour.bSetDescription("Colour for the first channel");
					ksecondaryColour.bSetDescription("Colour for the second channel");
					kgraphColour.bSetDescription("The colour of the graph.");
					kbackgroundColour.bSetDescription("The background colour of the view.");
					kdiagnostics.bSetDescription("Toggle diagnostic information in top-left corner.");
					kprimitiveSize.bSetDescription("The size of the rendered primitives (eg. lines or points).");
					kenvelopeMode.bSetDescription("Monitors the audio stream and automatically scales the input gain such that it approaches unity intensity (envelope following).");
					kenvelopeSmooth.bSetDescription("Responsiveness (RMS window size) - or the time it takes for the envelope follower to decay.");
					ksubSampleInterpolationMode.bSetDescription("Controls how point samples are interpolated to wave forms");
					kpctForDivision.bSetDescription("The minimum amount of free space that triggers a recursed frequency grid division; smaller values draw more frequency divisions.");
					kchannelConfiguration.bSetDescription("Select how the audio channels are interpreted.");
					ktriggerMode.bSetDescription("Select a mode for triggering waveforms - i.e. syncing to frequency content, time or transition information");
					ktriggerPhaseOffset.bSetDescription("A custom +/- full-circle offset for the phase on triggering");
					ktimeMode.bSetDescription("Specifies the working units of the time display");
					kdotSamples.bSetDescription("Marks sample positions when drawing subsampled interpolated lines");
					ktriggerOnCustomFrequency.bSetDescription("If toggled, the oscilloscope will trigger on the specified custom frequency instead of autodetecting a fundamental");
					kcustomFrequency.bSetDescription("Specifies a custom frequency to trigger on - input units can be notes (like a#2), hz, radians (rads), beats, samples (smps) or period (ms)");
					klowColour.bSetDescription("Colour for the low frequency band");
					kmidColour.bSetDescription("Colour for the mid frequency band");
					khighColour.bSetDescription("Colour for the high frequency band");
					kchannelColouring.bSetDescription("Method for colouring of channels, static equals just the drawing colour while spectral paints with separate colours for each frequency band");
					kfreqColourBlend.bSetDescription("Blending between the static and spectral colouring of channels");
					koverlayChannels.bSetDescription("Toggle to paint multiple channels on top of each other, otherwise they are painted in separate views");
					kcolourSmoothingTime.bSetDescription("Smooths the colour variation over the period of time");
					kcursorTracker.bSetDescription("Enable to create a tracker at the cursor displaying (x,y) values");
					ktrackerColour.bSetDescription("Colour of the cursor tracker");
					ktriggerHysteresis.bSetDescription("The hysteresis of the triggering function defines an opaque measure of how resistant the trigger is to change");
					ktriggerThreshold.bSetDescription("The triggering function will not consider any candidates below the threshold");
				}

				void initUI()
				{
					if (auto page = addPage("Settings", "icons/svg/gear.svg"))
					{
						/* if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&ktransform, 0);
							page->addSection(section, "Transform");
						} */

						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&koverlayChannels, 0);
							section->addControl(&kcursorTracker, 1);
							page->addSection(section, "Options");
						}

						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&kgain, 0);
							section->addControl(&kchannelConfiguration, 1);

							section->addControl(&kenvelopeSmooth, 0);
							section->addControl(&kenvelopeMode, 1);

							section->addControl(&kpctForDivision, 0);

							page->addSection(section, "Utility");
						}
						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&kwindow, 0);
							section->addControl(&ktimeMode, 1);
							section->addControl(&ktriggerMode, 0);
							section->addControl(&ktriggerPhaseOffset, 1);

							section->addControl(&ktriggerThreshold, 0);
							section->addControl(&ktriggerHysteresis, 1);

							section->addControl(&kcustomFrequency, 0);
							section->addControl(&ktriggerOnCustomFrequency, 1);

							page->addSection(section, "Spatial");
						}
					}

					if (auto page = addPage("Rendering", "icons/svg/brush.svg"))
					{
						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&kantiAlias, 0);
							section->addControl(&kdiagnostics, 1);
							section->addControl(&kdotSamples, 2);
							page->addSection(section, "Options");
						}
						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&kprimitiveSize, 0);
							section->addControl(&ksubSampleInterpolationMode, 1);

							section->addControl(&kgraphColour, 0);
							section->addControl(&kbackgroundColour, 1);
							section->addControl(&ktrackerColour, 0);

							page->addSection(section, "Look");
						}

						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&kcolourSmoothingTime, 0);
							section->addControl(&kchannelColouring, 1);

							section->addControl(&kprimaryColour, 0);
							section->addControl(&ksecondaryColour, 1);

							section->addControl(&kfreqColourBlend, 0);
							section->addControl(&klowColour, 1);
							section->addControl(&kmidColour, 0);
							section->addControl(&khighColour, 1);

							page->addSection(section, "Spectral colouring");
						}
					}

					if (auto page = addPage("Utility", "icons/svg/wrench.svg"))
					{
						if (auto section = new Signalizer::CContentPage::MatrixSection())
						{
							section->addControl(&kpresets, 0);
							page->addSection(section, "Presets");
						}
					}
				}

			private:

				void serializeEditorSettings(cpl::CSerializer::Archiver & archive, cpl::Version version)
				{
					archive << kwindow;
					archive << kgain;
					archive << kantiAlias;
					archive << kdiagnostics;
					archive << kgraphColour;
					archive << kbackgroundColour;
					archive << kprimaryColour;
					archive << ktransform;
					archive << kprimitiveSize;
					archive << kenvelopeMode;
					archive << kenvelopeSmooth;
					archive << ksubSampleInterpolationMode;
					archive << kpctForDivision;
					archive << kchannelConfiguration;
					archive << kpctForDivision;
					archive << ktriggerPhaseOffset;
					archive << ktriggerMode;
					archive << ktimeMode;
					archive << kdotSamples;
					archive << ktriggerOnCustomFrequency;
					archive << kcustomFrequency;
					archive << koverlayChannels;
					archive << kchannelColouring;
					archive << klowColour << kmidColour << khighColour;
					archive << ksecondaryColour;
					archive << kcolourSmoothingTime;
					archive << kcursorTracker;
					archive << ktrackerColour;
					archive << kfreqColourBlend;
					archive << ktriggerHysteresis;
					archive << ktriggerThreshold;
				}

				void deserializeEditorSettings(cpl::CSerializer::Archiver & builder, cpl::Version version)
				{
					// in general, controls should never restore values. However, older versions
					// of Signalizer does exactly this, so to keep backwards-compatibility, we
					// can obtain the preset values through this.
					cpl::Serialization::ScopedModifier m(cpl::CSerializer::Modifiers::RestoreValue, version < cpl::Version(0, 2, 8));
					builder << m;

					builder >> kwindow;
					builder >> kgain;
					builder >> kantiAlias;
					builder >> kdiagnostics;
					builder >> kgraphColour;
					builder >> kbackgroundColour;
					builder >> kprimaryColour;
					builder >> ktransform;
					builder >> kprimitiveSize;
					builder >> kenvelopeMode;
					builder >> kenvelopeSmooth;
					builder >> ksubSampleInterpolationMode;
					builder >> kpctForDivision;
					builder >> kchannelConfiguration;
					builder >> kpctForDivision;
					builder >> ktriggerPhaseOffset;
					builder >> ktriggerMode;
					builder >> ktimeMode;
					builder >> kdotSamples;
					builder >> ktriggerOnCustomFrequency;
					builder >> kcustomFrequency;
					builder >> koverlayChannels;
					builder >> kchannelColouring;
					builder >> klowColour >> kmidColour >> khighColour;
					builder >> ksecondaryColour;
					builder >> kcolourSmoothingTime;

					if (version > cpl::Version(0, 3, 1))
					{
						builder >> kcursorTracker;
						builder >> ktrackerColour;
						builder >> kfreqColourBlend;
					}

					if (version > cpl::Version(0, 3, 2))
					{
						builder >> ktriggerHysteresis;
						builder >> ktriggerThreshold;
					}

				}


				// entrypoints for completely storing values and settings in independant blobs (the preset widget)
				void serializeAll(cpl::CSerializer::Archiver & archive, cpl::Version version)
				{
					if (version < cpl::Version(0, 2, 8))
					{
						// presets from < 0.2.8 only store editor settings with values
						serializeEditorSettings(archive, version);
					}
					else
					{
						// store parameter and editor settings separately
						serializeEditorSettings(archive.getContent("Editor"), version);
						archive.getContent("Parameters") << parent;
					}

				}

				// entrypoints for completely storing values and settings in independant blobs (the preset widget)
				void deserializeAll(cpl::CSerializer::Builder & builder, cpl::Version version)
				{
					if (version < cpl::Version(0, 2, 8))
					{
						// presets from < 0.2.8 only store editor settings with values
						deserializeEditorSettings(builder, version);
					}
					else
					{
						// store parameter and editor settings separately
						deserializeEditorSettings(builder.getContent("Editor"), version);
						builder.getContent("Parameters") >> parent;
					}
				}

				cpl::CButton kantiAlias, kdiagnostics, kdotSamples, ktriggerOnCustomFrequency, koverlayChannels, kcursorTracker;
				cpl::CValueInputControl kcustomFrequency;
				cpl::CValueKnobSlider
					kwindow, kgain, kprimitiveSize, kenvelopeSmooth, kpctForDivision, ktriggerPhaseOffset, kcolourSmoothingTime, kfreqColourBlend,
					ktriggerHysteresis, ktriggerThreshold;
				cpl::CColourControl kprimaryColour, ksecondaryColour, kgraphColour, kbackgroundColour, klowColour, kmidColour, khighColour, ktrackerColour;
				cpl::CTransformWidget ktransform;
				cpl::CValueComboBox kenvelopeMode, ksubSampleInterpolationMode, kchannelConfiguration, ktriggerMode, ktimeMode, kchannelColouring;
				cpl::CPresetWidget kpresets;

				OscilloscopeContent & parent;

				SSOSurrogate<OscilloscopeController>
					editorSerializer,
					valueSerializer;
			};

			OscilloscopeContent(std::size_t parameterOffset, bool shouldCreateShortNames, SystemView system)
				: systemView(system)
				, parameterSet("Oscilloscope", "OS.", system.getProcessor(), static_cast<int>(parameterOffset))
				, audioHistoryTransformatter(system.getAudioStream(), LookaheadSize, audioHistoryTransformatter.Milliseconds)

				, dbRange(cpl::Math::dbToFraction(-120.0), cpl::Math::dbToFraction(120.0))
				, windowRange(0, 1000)
				, degreeRange(0, 360)
				, ptsRange(0.01, 10)
				, phaseRange(-180, 180)
				, reverseUnitRange(1, 0)
				, customTriggerRange(5, 48000)
				, colourSmoothRange(0.001, 1000)
				, triggerThresholdRange(0, 4)
				, msFormatter("ms")
				, degreeFormatter("degs")
				, ptsFormatter("pts")
				, customTriggerFormatter(system.getAudioStream())

				, autoGain("AutoGain")
				, envelopeWindow("EnvWindow", windowRange, msFormatter)
				, inputGain("InputGain", dbRange, dbFormatter)
				, windowSize("WindowSize", audioHistoryTransformatter, audioHistoryTransformatter)
				, antialias("AntiAlias", boolRange, boolFormatter)
				, diagnostics("Diagnostics", boolRange, boolFormatter)
				, primitiveSize("PixelSize", ptsRange, ptsFormatter)
				, subSampleInterpolation("SampleIntp")
				, pctForDivision("PctDiv", unityRange, pctFormatter)
				, channelConfiguration("ChConf")
				, triggerPhaseOffset("TrgPhase", phaseRange, degreeFormatter)
				, triggerMode("TrgMode")
				, timeMode("TimeMode")
				, dotSamples("DotSmps", boolRange, boolFormatter)
				, triggerOnCustomFrequency("CustomTrg", boolRange, boolFormatter)
				, customTriggerFrequency("TrgFreq", customTriggerRange, customTriggerFormatter)
				, overlayChannels("Overlay", boolRange, boolFormatter)
				, channelColouring("Colouring")
				, colourSmoothing("ColSmooth", colourSmoothRange, msFormatter)
				, cursorTracker("CursorTrck", unityRange, boolFormatter)
				, frequencyColouringBlend("FColBlend", unityRange, pctFormatter)
				, triggerHysteresis("TrgHstrs", unityRange, pctFormatter)
				, triggerThreshold("TrgThrhold", triggerThresholdRange, dbFormatter)

				, colourBehaviour()
				, primaryColour(colourBehaviour, "Prim.")
				, secondaryColour(colourBehaviour, "Sec.")
				, graphColour(colourBehaviour, "Graph.")
				, backgroundColour(colourBehaviour, "BackG.")
				, lowColour(colourBehaviour, "Low.")
				, midColour(colourBehaviour, "Mid.")
				, highColour(colourBehaviour, "High.")
				, trackerColour(colourBehaviour, "Trckr.")
				, tsfBehaviour()
				, transform(tsfBehaviour)

			{
				viewOffsets.emplace_back("ViewLeft", unityRange, basicFormatter);
				viewOffsets.emplace_back("ViewTop", unityRange, basicFormatter);
				viewOffsets.emplace_back("ViewRight", reverseUnitRange, basicFormatter);
				viewOffsets.emplace_back("ViewBottom", reverseUnitRange, basicFormatter);

				autoGain.fmt.setValues({ "None", "RMS", "Peak decay" });
				subSampleInterpolation.fmt.setValues({ "None", "Rectangular", "Linear", "Lanczos" });
				channelConfiguration.fmt.setValues({ "Left", "Right", "Mid", "Side", "Separate", "Mid+Side"});
				triggerMode.fmt.setValues({ "None", "Spectral", "Window", "Envelope" , "Zero-crossing"});
				timeMode.fmt.setValues({ "Time", "Cycles", "Beats" });
				channelColouring.fmt.setValues({ "Static", "Spectral energy" });

				// order matters
				auto singleParameters = {
					&autoGain.param,
					&envelopeWindow,
					&inputGain,
					&windowSize,
					&antialias,
					&diagnostics,
					&primitiveSize,
					&subSampleInterpolation.param,
					&channelConfiguration.param,
					&pctForDivision,
					&triggerPhaseOffset,
					&triggerMode.param,
					&timeMode.param,
					&dotSamples,
					&triggerOnCustomFrequency,
					&customTriggerFrequency,
					&overlayChannels,
					&channelColouring.param,
					&colourSmoothing,
					&cursorTracker,
					&frequencyColouringBlend,
					&triggerHysteresis,
					&triggerThreshold
				};

				for (auto sparam : singleParameters)
				{
					parameterSet.registerSingleParameter(sparam->generateUpdateRegistrator());
				}

				for (auto & v : viewOffsets)
				{
					parameterSet.registerSingleParameter(v.generateUpdateRegistrator());
				}

				for (auto cparam : { &primaryColour, &secondaryColour, &graphColour, &backgroundColour, &lowColour, &midColour, &highColour, &trackerColour })
				{
					parameterSet.registerParameterBundle(cparam, cparam->getBundleName());
				}

				parameterSet.registerParameterBundle(&transform, "3D.");

				parameterSet.seal();
				audioHistoryTransformatter.initialize(windowSize.getParameterView());
				timeMode.param.getParameterView().addListener(this);
			}

			~OscilloscopeContent()
			{
				timeMode.param.getParameterView().removeListener(this);
			}

			void parameterChangedUI(cpl::Parameters::Handle localHandle, cpl::Parameters::Handle globalHandle, ParameterSet::ParameterView * parameter)
			{
				if (parameter == &timeMode.param.getParameterView())
				{
					audioHistoryTransformatter.setTimeModeFromUI(timeMode.param.getAsTEnum<TimeMode>());
				}
			}

			virtual std::unique_ptr<StateEditor> createEditor() override
			{
				return std::make_unique<OscilloscopeController>(*this);
			}

			virtual ParameterSet & getParameterSet() override
			{
				return parameterSet;
			}

			virtual void serialize(cpl::CSerializer::Archiver & archive, cpl::Version version) override
			{
				archive << windowSize;
				archive << inputGain;
				archive << antialias;
				archive << diagnostics;
				archive << graphColour;
				archive << backgroundColour;
				archive << primaryColour;
				archive << transform;
				archive << primitiveSize;
				archive << autoGain.param;
				archive << envelopeWindow;
				archive << subSampleInterpolation.param;
				archive << channelConfiguration.param;
				archive << pctForDivision;
				archive << triggerPhaseOffset;
				archive << triggerMode.param;
				archive << timeMode.param;
				for (auto && v : viewOffsets)
				{
					archive << v;
				}
				archive << dotSamples;
				archive << triggerOnCustomFrequency;
				archive << customTriggerFrequency;
				archive << overlayChannels;
				archive << channelColouring.param;
				archive << lowColour << midColour << highColour;
				archive << secondaryColour;
				archive << colourSmoothing;
				archive << cursorTracker;
				archive << trackerColour;
				archive << frequencyColouringBlend;
				archive << triggerHysteresis;
				archive << triggerThreshold;
			}

			virtual void deserialize(cpl::CSerializer::Builder & builder, cpl::Version version) override
			{
				builder >> windowSize;
				builder >> inputGain;
				builder >> antialias;
				builder >> diagnostics;
				builder >> graphColour;
				builder >> backgroundColour;
				builder >> primaryColour;
				builder >> transform;
				builder >> primitiveSize;
				builder >> autoGain.param;
				builder >> envelopeWindow;
				builder >> subSampleInterpolation.param;
				builder >> channelConfiguration.param;
				builder >> pctForDivision;
				builder >> triggerPhaseOffset;
				builder >> triggerMode.param;
				builder >> timeMode.param;

				for (auto & v : viewOffsets)
				{
					builder >> v;
				}

				builder >> dotSamples;
				builder >> triggerOnCustomFrequency;
				builder >> customTriggerFrequency;
				builder >> overlayChannels;
				builder >> channelColouring.param;
				builder >> lowColour >> midColour >> highColour;
				builder >> secondaryColour;
				builder >> colourSmoothing;

				if (version >= cpl::Version(0, 3, 1))
				{
					builder >> cursorTracker;
					builder >> trackerColour;
					builder >> frequencyColouringBlend;
				}

				if (version >= cpl::Version(0, 3, 2))
				{
					builder >> triggerHysteresis;
					builder >> triggerThreshold;
				}
			}

			WindowSizeTransformatter<ParameterSet::ParameterView> audioHistoryTransformatter;
			LinearHzFormatter<ParameterSet::ParameterView> customTriggerFormatter;
			SystemView systemView;
			ParameterSet parameterSet;

			cpl::UnitFormatter<double>
				msFormatter,
				degreeFormatter,
				ptsFormatter;

			cpl::PercentageFormatter<double>
				pctFormatter;

			cpl::DBFormatter<double> dbFormatter;
			cpl::BooleanFormatter<double> boolFormatter;

			cpl::BasicFormatter<double> basicFormatter;
			cpl::BooleanRange<double> boolRange;

			cpl::ExponentialRange<double> dbRange, colourSmoothRange;

			cpl::LinearRange<double>
				ptsRange,
				windowRange,
				degreeRange,
				phaseRange,
				reverseUnitRange,
				customTriggerRange,
				triggerThresholdRange;

			cpl::UnityRange<double> unityRange;

			enum ViewOffsets
			{
				Left, Top, Right, Bottom, end
			};

			cpl::ParameterValue<ParameterSet::ParameterView>
				envelopeWindow,
				inputGain,
				windowSize,
				antialias,
				diagnostics,
				primitiveSize,
				pctForDivision,
				triggerPhaseOffset,
				dotSamples,
				triggerOnCustomFrequency,
				customTriggerFrequency,
				overlayChannels,
				colourSmoothing,
				cursorTracker,
				frequencyColouringBlend,
				triggerHysteresis,
				triggerThreshold;

			std::vector<cpl::ParameterValue<ParameterSet::ParameterView>> viewOffsets;

			ChoiceParameter
				autoGain,
				channelConfiguration,
				subSampleInterpolation,
				triggerMode,
				timeMode,
				channelColouring;

			cpl::ParameterColourValue<ParameterSet::ParameterView>::SharedBehaviour colourBehaviour;

			cpl::ParameterColourValue<ParameterSet::ParameterView>
				primaryColour,
				graphColour,
				secondaryColour,
				backgroundColour,
				lowColour, midColour, highColour,
				trackerColour;

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
