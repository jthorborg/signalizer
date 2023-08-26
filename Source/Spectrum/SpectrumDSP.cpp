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

	file:SpectrumDSP.cpp

		Implementation of the dsp code for the spectrum.

*************************************************************************************/

#include "Spectrum.h"
#include <cpl/ffts.h>
#include <cpl/system/SysStats.h>
#include <cpl/lib/LockFreeDataQueue.h>
#include <cpl/stdext.h>
#include <cpl/JobSystem.h>
#include <execution>
#include "TransformDSP.inl"

namespace Signalizer
{
	Spectrum::ProcessorShell::ProcessorShell(std::shared_ptr<const SharedBehaviour>& behaviour)
		: globalBehaviour(behaviour), frameQueue(10)
	{
	}

	std::size_t Spectrum::getBlobSamples() const noexcept
	{
		return std::max<std::size_t>(10, static_cast<std::size_t>(content->blobSize.getTransformedValue() * 0.001 * getSampleRate()));
	}

	void Spectrum::resetState()
	{
		flags.resetStateBuffers = true;
	}

	bool Spectrum::processNextSpectrumFrame(const Constant& constant, TransformPair& transform)
	{
		SFrameQueue::ElementAccess access;
		if (processor->frameQueue.popElement(access))
		{
			FrameVector& curFrame(*access.getData());

			std::size_t numFilters = getNumFilters();

			// the size will be zero for a couple of frames, if there's some messing around with window sizes
			// or we get audio running before anything is actually initiated.
			if (curFrame.size() != 0)
			{
				if (curFrame.size() == numFilters)
				{
					transform.postProcessTransform(constant, cpl::as_uarray(curFrame).reinterpret<UComplex::Scalar>());
				}
				else
				{
					// linearly interpolate bins. if we win the cpu-lottery one day, change this to sinc.
					std::vector<UComplex> tempSpace(numFilters);

					// interpolation factor.
					UComplex::Scalar wspToNext = (curFrame.size() - 1) / UComplex::Scalar(std::max<std::size_t>(1, numFilters));

					for (std::size_t n = 0; n < numFilters; ++n)
					{
						auto y2 = n * wspToNext;
						auto x = static_cast<std::size_t>(y2);
						auto yFrac = y2 - x;
						tempSpace[n] = curFrame[x] * (UComplex::Scalar(1) - yFrac) + curFrame[x + 1] * yFrac;
					}
					
					transform.postProcessTransform(constant, cpl::as_uarray(tempSpace).reinterpret<UComplex::Scalar>());
				}
			}

			return true;
		}
		return false;
	}

	struct AudioDispatcher
	{
		template<typename ISA> static void dispatch(Spectrum::ProcessorShell& shell, AudioStream::ListenerContext& source, AudioStream::DataType** buffer, std::size_t numChannels, std::size_t numSamples)
		{
			constexpr std::size_t concurrency = 4;

			CPL_RUNTIME_ASSERTION(numChannels % 2 == 0);

			auto access = shell.streamState.lock();

			if (numChannels < 2 || access->pairs.empty())
				return;

			CPL_RUNTIME_ASSERTION((numChannels / 2) == access->pairs.size());
			CPL_RUNTIME_ASSERTION(source.getNumDeferredSamples() == 0);

			// conditional lock.
			std::optional<AudioStream::AudioBufferAccess> aba;
			if (access->constant.algo == SpectrumContent::TransformAlgorithm::FFT)
			{
				aba.emplace(source.getAudioBufferViews());
			}

			cpl::jobs::parallel_for(
				numChannels / 2,
				[&](auto i)
				{
					std::optional<Spectrum::TransformPair::AudioPair> views;
					if (access->constant.algo == SpectrumContent::TransformAlgorithm::FFT)
					{
						views = Spectrum::TransformPair::AudioPair{ aba->getView(i * 2), aba->getView(i * 2 + 1) };
					}
					
					access->pairs[i].audioEntryPoint<ISA>(
						access->constant, 
						views,
						{ buffer[i * 2], buffer[i * 2 + 1] },
						numSamples
					);
				}
			);

			if (access->constant.displayMode == SpectrumContent::DisplayMode::ColourSpectrum)
			{
				auto& state = access->pairs[0];

				for (std::size_t i = 0; i < state.sfbuf.size(); ++i)
				{
					Spectrum::SFrameQueue::ElementAccess access;

					// shouldn't be possible.
					if (!shell.frameQueue.acquireFreeElement<true, false>(access))
						break;

					auto& buffer = *access.getData();
					buffer = std::move(state.sfbuf[i]);
				}

				for (auto& pair : access->pairs)
					pair.sfbuf.clear();
			}
		}
	};

	void Spectrum::ProcessorShell::onStreamAudio(AudioStream::ListenerContext& source, AudioStream::DataType** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		if (isSuspended && globalBehaviour->stopProcessingOnSuspend)
			return;

		cpl::simd::dynamic_isa_dispatch<ProcessingType, AudioDispatcher>(*this, source, buffer, numChannels, numSamples);
	}

	void Spectrum::ProcessorShell::onStreamPropertiesChanged(AudioStream::ListenerContext& source, const AudioStream::AudioStreamInfo& before)
	{
		auto access = streamState.lock();

		access->audioStreamChangeVersion.bump();
		access->streamLocalSampleRate = source.getInfo().sampleRate;

		access->pairs.resize(source.getInfo().channels / 2);
	}

	float Spectrum::getSampleRate() const noexcept
	{
		return state.sampleRate;
	}

	int Spectrum::getAxisPoints() const noexcept
	{
		return static_cast<int>(state.axisPoints);
	}

	int Spectrum::getNumFilters() const noexcept
	{
		return getAxisPoints();
	}

	double Spectrum::getOptimalFramesPerUpdate() const noexcept
	{
#pragma message cwarn("collect this somewhere.")
		const double monitorRefreshRate = 60.0;
		auto res = double(isOpenGL() ? (monitorRefreshRate / getSwapInterval()) : refreshRate) / getBlobSamples();
		assert(std::isnormal(res));
		return res;
	}

	std::size_t Spectrum::getApproximateStoredFrames() const noexcept
	{
#pragma message cwarn("fix this to include channels, other processing methods.. etc.")
		return processor->frameQueue.enqueuededElements();
	}

	double Spectrum::getScallopingLossAtCoordinate(std::size_t coordinate, const Constant& constant)
	{
		auto ret = 0.6366; // default absolute worst case (equivalent to sinc(0.5), ie. rectangular windows

		if (state.axisPoints < 3)
			return ret;

		if (state.displayMode == SpectrumContent::DisplayMode::LineGraph)
		{
#pragma message cwarn("Fix this to use direct values")
			auto& value = content->dspWin;
			auto type = value.getWindowType();
			auto alpha = value.getAlpha();
			auto beta = value.getBeta();
			auto symmetry = value.getWindowShape();

			bool canOvershoot = false;
			auto sampleRate = getSampleRate();
			double normalizedBandwidth = 0, fractionateScallopLoss = normalizedBandwidth;


			auto safeIndex = cpl::Math::confineTo<std::size_t>(coordinate, 0, state.axisPoints - 2);
			if (state.algo == SpectrumContent::TransformAlgorithm::RSNT)
			{
				if (state.viewScale == SpectrumContent::ViewScaling::Linear)
				{

					normalizedBandwidth = state.windowSize * std::abs((double)constant.mapFrequency(safeIndex + 1) - constant.mapFrequency(safeIndex)) / sampleRate;

					normalizedBandwidth = std::min(0.5, normalizedBandwidth);
				}

				fractionateScallopLoss = cpl::dsp::windowScallopLoss(type, 4, normalizedBandwidth, cpl::dsp::Windows::Shape::Periodic, alpha, beta);
				// resonators have, per definition, at least 3 dB bandwidth, so number is equal to 10^(-3/20)
				fractionateScallopLoss = std::min(fractionateScallopLoss, 0.70794578438413791080221494218931);
			}
			else if (state.algo == SpectrumContent::TransformAlgorithm::FFT)
			{
				normalizedBandwidth = 0.5;

				if (state.binPolation == SpectrumContent::BinInterpolation::Lanczos && state.frequencyTrackingGraph != SpectrumContent::LineGraphs::Transform)
				{
					normalizedBandwidth = state.windowSize * std::abs((double)constant.mapFrequency(safeIndex + 1) - constant.mapFrequency(safeIndex)) / sampleRate;
					if (normalizedBandwidth < 0.5)
					{
						canOvershoot = true;
					}
					else
					{
						normalizedBandwidth = 0.5;
					}

				}

				fractionateScallopLoss = cpl::dsp::windowScallopLoss(type, 4, normalizedBandwidth, symmetry, alpha, beta);
			}

			ret = fractionateScallopLoss;
		}
		return ret;
	}
};
