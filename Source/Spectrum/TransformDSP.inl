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

	file:TransformDSP.inl

		Implementation for template functions of TransformPair.h

*************************************************************************************/

#ifndef SIGNALIZER_TRANSFORM_DSP_INL
#define SIGNALIZER_TRANSFORM_DSP_INL

#include "TransformPair.h"
#include "TransformConstant.h"

namespace Signalizer
{
	template<typename T>
	inline bool TransformPair<T>::prepareTransform(const Constant& constant, const AudioPair& views)
	{
		{
			// we need the buffers to be same size, and at least equal or greater in size of ours (cant fill in information).
			// this is a very rare condition that can be solved by locking the audio access during the flags update and this
			// call, however to avoid unnecessary locks we skip a frame instead once in a while.
			if (views[0].size() != views[1].size() || views[0].size() < constant.windowSize)
				return false;

			// can't underflow
			std::size_t offset = views[0].size() - constant.windowSize;

			switch (constant.algo)
			{
			case SpectrumContent::TransformAlgorithm::FFT:
			{
				auto buffer = getAudioMemory<std::complex<T>>(constant.transformSize);
				std::size_t channel = 1;
				std::size_t i = 0;

				switch (constant.configuration)
				{
				case SpectrumChannels::Left:
					channel = 0;
				case SpectrumChannels::Right:
				{
					for (std::size_t indice = 0; indice < AudioStream::bufferIndices; ++indice)
					{
						std::size_t range = views[channel].getItRange(indice);
						auto it = views[channel].getItIndex(indice);

						if (range > offset)
						{
							range -= offset;
							it += offset;

							while (range--)
							{
								buffer[i] = *it++ * constant.windowKernel[i];
								i++;
							}

							offset = 0;
						}
						else
						{
							offset -= range;
						}

					}
					break;
				}
				case SpectrumChannels::Merge:
				{
					for (std::size_t indice = 0; indice < AudioStream::bufferIndices; ++indice)
					{
						std::size_t range = views[0].getItRange(indice);
						auto left = views[0].getItIndex(indice);
						auto right = views[1].getItIndex(indice);

						if (range > offset)
						{
							range -= offset;
							left += offset;
							right += offset;

							while (range--)
							{
								buffer[i] = (*left++ + *right++) * constant.windowKernel[i] * 0.5f;
								i++;
							}
							offset = 0;
						}
						else
						{
							offset -= range;
						}
					}
					break;
				}
				case SpectrumChannels::Side:
				{
					for (std::size_t indice = 0; indice < AudioStream::bufferIndices; ++indice)
					{
						std::size_t range = views[0].getItRange(indice);
						auto left = views[0].getItIndex(indice);
						auto right = views[1].getItIndex(indice);

						if (range > offset)
						{
							range -= offset;
							left += offset;
							right += offset;

							while (range--)
							{
								buffer[i] = (*left++ - *right++) * constant.windowKernel[i] * 0.5f;
								i++;
							}

							offset = 0;
						}
						else
						{
							offset -= range;
						}
					}
					break;
				}
				case SpectrumChannels::MidSide:
				{
					for (std::size_t indice = 0; indice < AudioStream::bufferIndices; ++indice)
					{
						std::size_t range = views[0].getItRange(indice);
						auto left = views[0].getItIndex(indice);
						auto right = views[1].getItIndex(indice);

						if (range > offset)
						{
							range -= offset;
							left += offset;
							right += offset;

							while (range--)
							{
								buffer[i] = std::complex<T>
									(
										(*left + *right) * constant.windowKernel[i] * 0.5f,
										(*left - *right) * constant.windowKernel[i] * 0.5f
										);
								left++;
								right++;
								i++;
							}

							offset = 0;
						}
						else
						{
							offset -= range;
						}
					}
					break;
				}
				case SpectrumChannels::Phase:
				case SpectrumChannels::Separate:
				case SpectrumChannels::Complex:
				{
					for (std::size_t indice = 0; indice < AudioStream::bufferIndices; ++indice)
					{
						std::size_t range = views[0].getItRange(indice);
						auto left = views[0].getItIndex(indice);
						auto right = views[1].getItIndex(indice);

						if (range > offset)
						{
							range -= offset;
							left += offset;
							right += offset;

							while (range--)
							{
								buffer[i] = std::complex<T>
								{
									*left++ * constant.windowKernel[i],
									*right++ * constant.windowKernel[i]
								};
								i++;
							}
							offset = 0;
						}
						else
						{
							offset -= range;
						}
					}
					break;
				}
				}
				//zero-pad until buffer is filled

				for (size_t pad = i; pad < constant.transformSize; ++pad)
				{
					buffer[pad] = (T)0;
				}

				break;
			}
			break;
			}
		}
		return true;
	}

	template<typename T>
	inline bool TransformPair<T>::prepareTransform(const Constant& constant, const AudioPair& views, std::array<AFloat*, 2> preliminaryAudio, std::size_t numSamples)
	{
		{

			// we need the buffers to be same size, and at least equal or greater in size of ours (cant fill in information).
			// this is a very rare condition that can be solved by locking the audio access during the flags update and this
			// call, however to avoid unnecessary locks we skip a frame instead once in a while.
			if (views[0].size() != views[1].size() || views[0].size() < constant.windowSize)
				return false;

			// extra discarded samples in case the incoming buffer is larger than our own
			std::size_t extraDiscardedSamples = views[0].size() - constant.windowSize;

			switch (constant.algo)
			{
			case SpectrumContent::TransformAlgorithm::FFT:
			{
				auto buffer = getAudioMemory<std::complex<T>>(constant.transformSize);
				std::size_t channel = 1;
				std::size_t i = 0;
				std::size_t stop = std::min(numSamples, constant.windowSize);

				std::size_t offset = stop + extraDiscardedSamples;
				auto sizeToStopAt = constant.windowSize - offset;

				switch (constant.configuration)
				{
				case SpectrumChannels::Left:
					channel = 0;
				case SpectrumChannels::Right:
				{
					// get start from buffers - first indice is a special case.
					// notice the buffers are reversed in time, so we pull old data first with an offset,
					// and then fill in the new
					for (std::size_t indice = 0; indice < AudioStream::bufferIndices; ++indice)
					{
						std::size_t range = views[channel].getItRange(indice);
						auto it = views[channel].getItIndex(indice);

						if (range > offset)
						{

							// advance our current progress (from preliminary audio)
							// notice this only has effect for the first buffer coming in here.
							range -= offset;
							it += offset;

							while (range-- && i < sizeToStopAt)
							{
								buffer[i] = *it++ * constant.windowKernel[i];
								i++;
							}

							offset = 0;

						}
						else
						{
							offset -= range;
						}
					}


					// process preliminary
					for (std::size_t k = 0; k < stop; ++i, k++)
					{
						buffer[i] = preliminaryAudio[channel][k] * constant.windowKernel[i];
					}


					break;
				}
				case SpectrumChannels::Merge:
				{
					for (std::size_t indice = 0; indice < AudioStream::bufferIndices; ++indice)
					{
						std::size_t range = views[channel].getItRange(indice);
						auto it = views[channel].getItIndex(indice);
						auto left = views[0].getItIndex(indice);
						auto right = views[1].getItIndex(indice);
						if (range > offset)
						{

							range -= offset;
							it += offset;

							while (range-- && i < sizeToStopAt)
							{
								buffer[i] = (*left++ + *right++) * constant.windowKernel[i] * 0.5f;
								i++;
							}

							offset = 0;

						}
						else
						{
							offset -= range;
						}
					}

					for (std::size_t k = 0; k < stop; ++i, k++)
					{
						buffer[i] = (preliminaryAudio[0][k] + preliminaryAudio[1][k]) * constant.windowKernel[i] * (T)0.5;
					}

					break;
				}
				case SpectrumChannels::Side:
				{
					for (std::size_t indice = 0; indice < AudioStream::bufferIndices; ++indice)
					{
						std::size_t range = views[channel].getItRange(indice);
						auto it = views[channel].getItIndex(indice);
						auto left = views[0].getItIndex(indice);
						auto right = views[1].getItIndex(indice);
						if (range > offset)
						{

							range -= offset;
							it += offset;

							while (range-- && i < sizeToStopAt)
							{
								buffer[i] = (*left++ - *right++) * constant.windowKernel[i] * (T)0.5;
								i++;
							}

							offset = 0;

						}
						else
						{
							offset -= range;
						}
					}

					for (std::size_t k = 0; k < stop; ++i, k++)
					{
						buffer[i] = (preliminaryAudio[0][k] - preliminaryAudio[1][k]) * constant.windowKernel[i] * (T)0.5;
					}

					break;
				}
				case SpectrumChannels::MidSide:
				{
					for (std::size_t indice = 0; indice < AudioStream::bufferIndices; ++indice)
					{
						std::size_t range = views[channel].getItRange(indice);
						auto it = views[channel].getItIndex(indice);
						auto left = views[0].getItIndex(indice);
						auto right = views[1].getItIndex(indice);
						if (range > offset)
						{

							range -= offset;
							it += offset;

							while (range-- && i < sizeToStopAt)
							{
								buffer[i] = std::complex<T>
									(
										(*left + *right) * constant.windowKernel[i] * (T)0.5,
										(*left - *right) * constant.windowKernel[i] * (T)0.5
										);
								left++;
								right++;
								i++;
							}

							offset = 0;

						}
						else
						{
							offset -= range;
						}
					}

					for (std::size_t k = 0; k < stop; ++i, k++)
					{
						buffer[i] = std::complex<T>
							(
								(preliminaryAudio[0][k] + preliminaryAudio[1][k]) * constant.windowKernel[i] * (T)0.5,
								(preliminaryAudio[0][k] - preliminaryAudio[1][k]) * constant.windowKernel[i] * (T)0.5
								);
					}

					break;
				}
				case SpectrumChannels::Phase:
				case SpectrumChannels::Separate:
				case SpectrumChannels::Complex:
				{

					for (std::size_t indice = 0; indice < AudioStream::bufferIndices; ++indice)
					{
						std::size_t range = views[channel].getItRange(indice);
						auto it = views[channel].getItIndex(indice);
						auto left = views[0].getItIndex(indice);
						auto right = views[1].getItIndex(indice);
						if (range > offset)
						{

							range -= offset;
							it += offset;

							while (range-- && i < sizeToStopAt)
							{
								buffer[i] = std::complex<T>
								{
									*left++ * constant.windowKernel[i],
									*right++ * constant.windowKernel[i]
								};
								i++;
							}

							offset = 0;

						}
						else
						{
							offset -= range;
						}
					}

					for (std::size_t k = 0; k < stop; ++i, k++)
					{
						buffer[i] = std::complex<T>
							(
								preliminaryAudio[0][k] * constant.windowKernel[i],
								preliminaryAudio[1][k] * constant.windowKernel[i]
								);
					}

					break;
				}
				}
				//zero-pad until buffer is filled
				for (size_t pad = i; pad < constant.transformSize; ++pad)
				{
					buffer[pad] = 0;
				}

				break;
			}
			break;
			}
		}
		return true;
	}

	template<typename T>
	inline void TransformPair<T>::doTransform(const Constant& constant)
	{
		switch (constant.algo)
		{
		case SpectrumContent::TransformAlgorithm::FFT:
		{
			if (constant.transformSize > 0)
				constant.fft.forward(
					getAudioMemory<std::complex<T>>(constant.transformSize),
					getAudioMemory<std::complex<T>>(constant.transformSize),
					getWork<std::complex<T>>(constant.transformSize)
				);
			break;
		}
		}
	}

	template<typename T>
	inline void TransformPair<T>::mapToLinearSpace(const Constant& constant)
	{
		using namespace cpl;

		switch (constant.algo)
		{
		case SpectrumContent::TransformAlgorithm::FFT:
		{
			const auto lanczosFilterSize = 5;
			cpl::ssize_t bin = 0, oldBin = 0, maxLBin, maxRBin = 0;
			Types::fsint_t N = static_cast<Types::fsint_t>(constant.transformSize);

			// we rely on mapping indexes, so we need N > 2 at least.
			if (N < 3 || constant.sampleRate < 1)
				return;

			std::size_t numBins = N >> 1;
			auto const topFrequency = constant.sampleRate / 2;
			auto const freqToBin = static_cast<T>(numBins / topFrequency);

			std::complex<T> leftMax, rightMax;

			T maxLMag, maxRMag, newLMag, newRMag;

			// complex transform results, N + 1 size
			auto csf = getAudioMemory<std::complex<T>>(constant.transformSize + 1);
			// buffer for single results, numPoints * 2 size
			auto wsp = getWork<T>(constant.axisPoints * 2);
			// buffer for complex results, numPoints size
			auto csp = getWork<std::complex<T>>(constant.axisPoints * 2);

			// this will make scaling correct regardless of amount of zero-padding
			// notice the 0.5: fft's of size 32 will output 16 for exact frequency bin matches,
			// so we halve the reciprocal scaling factor to normalize the size.
			auto const invSize = static_cast<T>(constant.windowKernelScale / (constant.windowSize * 0.5));

			switch (constant.configuration)
			{
			case SpectrumChannels::Left:
			case SpectrumChannels::Right:
			case SpectrumChannels::Merge:
			case SpectrumChannels::Side:
			{
				oldBin = constant.mappedFrequencies[0] * freqToBin;

				// the DC (0) and nyquist bin are NOT 'halved' due to the symmetric nature of the fft,
				// so halve these:
				csf[0] *= 0.5;
				csf[N >> 1] *= 0.5;

				// TODO: Vectorize
				for (std::size_t i = 0; i < numBins; ++i)
				{
					csf[i] = std::abs(csf[i]);
				}

				double fftBandwidth = 1.0 / numBins;
				//double pxlBandwidth = 1.0 / numPoints;
				cpl::Types::fint_t x = 0;
				switch (constant.binPolation)
				{
				case SpectrumContent::BinInterpolation::None:
					for (x = 0; x < constant.axisPoints - 1; ++x)
					{
						// bandwidth of the filter for this 'line', or point
						double bwForLine = (constant.mappedFrequencies[x + 1] - constant.mappedFrequencies[x]) / topFrequency;
						// as long as the bandwidth is smaller than our fft resolution, we interpolate the points
						// otherwise, break out and sample the max values of the bins inside the bandwidth
						if (bwForLine > fftBandwidth)
							break;
						// +0.5 to centerly space bins.
						auto index = Math::confineTo((std::size_t)(constant.mappedFrequencies[x] * freqToBin + 0.5), 0, numBins);
						csp[x] = invSize * csf[index];
					}
					break;
				case SpectrumContent::BinInterpolation::Linear:
					for (x = 0; x < constant.axisPoints - 1; ++x)
					{
						double bwForLine = (constant.mappedFrequencies[x + 1] - constant.mappedFrequencies[x]) / topFrequency;
						if (bwForLine > fftBandwidth)
							break;

						csp[x] = invSize * dsp::linearFilter<std::complex<T>>(csf, constant.mappedFrequencies[x] * freqToBin);
					}
					break;
				case SpectrumContent::BinInterpolation::Lanczos:
					for (x = 0; x < constant.axisPoints - 1; ++x)
					{

						double bwForLine = (constant.mappedFrequencies[x + 1] - constant.mappedFrequencies[x]) / topFrequency;
						if (bwForLine > fftBandwidth)
							break;

						csp[x] = invSize * dsp::lanczosFilter<std::complex<T>, true>(csf, constant.mappedFrequencies[x] * freqToBin, lanczosFilterSize);
					}
					break;
				default:
					break;
				}

				oldBin = constant.mappedFrequencies[x] * freqToBin;

				for (; x < constant.axisPoints; ++x)
				{
					maxLMag = maxRMag = newLMag = newRMag = 0;

					bin = static_cast<std::size_t>(constant.mappedFrequencies[x] * freqToBin);
#ifdef DEBUG
					if ((std::size_t)bin >= constant.transformSize)
						CPL_RUNTIME_EXCEPTION("Corrupt frequency mapping!");
#endif
					maxRBin = maxLBin = bin;

					auto diff = bin - oldBin;
					auto counter = diff ? 1 : 0;
					// here we loop over all the bins that is mapped for a single coordinate
					do
					{
						auto offset = oldBin + counter;
						newLMag = Math::square(csf[offset]);
						// select highest number in this chunk for display. Not exactly correct, though.
						if (newLMag > maxLMag)
						{
							maxLBin = oldBin + counter;
							maxLMag = newLMag;
						}
						counter++;
						diff--;
					} while (diff > 0);

					csp[x] = invSize * csf[maxLBin];

					oldBin = bin;
				}

				break;
			}
			case SpectrumChannels::Phase:
			{
				// two-for-one pass, first channel is 0... N/2 -1, second is N/2 .. N -1
				dsp::separateTransformsIPL(csf);

				// fix up DC and nyquist bins (see previous function documentation)
				csf[N] = csf[0].imag() * 0.5;
				csf[0] = csf[0].real() * 0.5;
				csf[N >> 1] *= 0.5;
				csf[(N >> 1) - 1] *= 0.5;

				// TODO: rewrite this.
				// firstly, we have to do a phase cancellation pass.
				// this is because magnitude interpolation is wrong for
				// complex vectors, it needs to be done on magnitude.
				// however, phase calculation needs to be done on vectors.
				// so we first do a phase cancellation pass, and afterwards
				// absolute the lower bottom of the transform, that needs to be interpolated.

				// The index of the transform, where the bandwidth is higher than mapped pixels (so no more interpolation is needed)
				// TODO: This can be calculated from view mapping scale and N pixels.
				std::size_t bandWidthBreakingPoint = constant.axisPoints;

				double fftBandwidth = 1.0 / numBins;
				//double pxlBandwidth = 1.0 / numPoints;
				cpl::Types::fint_t x = 0;
				switch (constant.binPolation)
				{
				case SpectrumContent::BinInterpolation::Linear:
				{
					// phase pass
					for (x = 0; x < constant.axisPoints - 1; ++x)
					{
						double bwForLine = (constant.mappedFrequencies[x + 1] - constant.mappedFrequencies[x]) / topFrequency;
						if (bwForLine > fftBandwidth)
						{
							bandWidthBreakingPoint = x;
							break;
						}

						auto iLeft = dsp::linearFilter<std::complex<T>>(csf, constant.mappedFrequencies[x] * freqToBin);
						auto iRight = dsp::linearFilter<std::complex<T>>(csf, N - (constant.mappedFrequencies[x] * freqToBin));

						auto cancellation = invSize * std::sqrt(Math::square(iLeft + iRight));
						auto mid = invSize * (std::abs(iLeft) + std::abs(iRight));

						wsp[x * 2 + 1] = T(1) - (mid > 0 ? (cancellation / mid) : 0);
					}

					// normalize vectors we need to magnitude-interpolate.
					/* auto stop = std::size_t(bandWidthBreakingPoint * pointsToBin);
					for (std::size_t i = 1; i < stop; ++i)
					{
						csf[i] = std::abs(csf[i]);
						csf[N - i] = std::abs(csf[N - i]);
					} */

					std::size_t normalizedPosition = 0U;

					// magnitude pass.
					// TODO: we aren't calculating the last point, because that means we have to normalize one position further
					// to ensure correct magnitude interpretation. That will however nullify the phase response in the next pass,
					// since normalization eliminates phase information.

					const std::size_t linearFilterSize = 1;

					for (x = 0; x < bandWidthBreakingPoint; ++x)
					{
						// fractional bin position
						auto binPosition = constant.mappedFrequencies[x] * freqToBin;

						// normalize vectors we need to magnitude-interpolate.
						while ((binPosition + linearFilterSize) > normalizedPosition && x < bandWidthBreakingPoint - (linearFilterSize + 1))
						{
							csf[normalizedPosition] = std::abs(csf[normalizedPosition]);
							csf[N - normalizedPosition] = std::abs(csf[N - normalizedPosition]);
							normalizedPosition++;
						}

						auto iLeft = dsp::linearFilter<std::complex<T>>(csf, binPosition);
						auto iRight = dsp::linearFilter<std::complex<T>>(csf, N - binPosition);
						// TODO: abs not really needed, because of normalization.
						wsp[x * 2] = invSize * (std::abs(iLeft) + std::abs(iRight));

					}
					break;
				}
				case SpectrumContent::BinInterpolation::Lanczos:
				{
					for (x = 0; x < constant.axisPoints - 1; ++x)
					{

						double bwForLine = (constant.mappedFrequencies[x + 1] - constant.mappedFrequencies[x]) / topFrequency;
						if (bwForLine > fftBandwidth)
						{
							bandWidthBreakingPoint = x;
							break;
						}

						auto iLeft = dsp::lanczosFilter<std::complex<T>, true>(csf, constant.mappedFrequencies[x] * freqToBin, lanczosFilterSize);
						auto iRight = dsp::lanczosFilter<std::complex<T>, true>(csf, N - (constant.mappedFrequencies[x] * freqToBin), lanczosFilterSize);

						auto cancellation = invSize * std::sqrt(Math::square(iLeft + iRight));
						auto mid = invSize * (std::abs(iLeft) + std::abs(iRight));

						wsp[x * 2 + 1] = T(1) - (mid > 0 ? (cancellation / mid) : 0);
					}

					std::size_t normalizedPosition = 0U;

					// magnitude pass.
					// TODO: we aren't calculating the last point, because that means we have to normalize one position further
					// to ensure correct magnitude interpretation. That will however nullify the phase response in the next pass,
					// since normalization eliminates phase information.


					for (x = 0; x < bandWidthBreakingPoint; ++x)
					{
						// fractional bin position
						auto binPosition = constant.mappedFrequencies[x] * freqToBin;

						// normalize vectors we need to magnitude-interpolate.
						while ((binPosition + lanczosFilterSize) > normalizedPosition && x < bandWidthBreakingPoint - (lanczosFilterSize + 1))
						{
							csf[normalizedPosition] = std::abs(csf[normalizedPosition]);
							csf[N - normalizedPosition] = std::abs(csf[N - normalizedPosition]);
							normalizedPosition++;
						}

						auto iLeft = dsp::lanczosFilter<std::complex<T>, true>(csf, constant.mappedFrequencies[x] * freqToBin, lanczosFilterSize);
						auto iRight = dsp::lanczosFilter<std::complex<T>, true>(csf, N - (constant.mappedFrequencies[x] * freqToBin), lanczosFilterSize);

						wsp[x * 2] = invSize * (std::abs(iLeft) + std::abs(iRight));
					}

					break;
				}
				default:
					for (x = 0; x < constant.axisPoints - 1; ++x)
					{
						// bandwidth of the filter for this 'line', or point
						double bwForLine = (constant.mappedFrequencies[x + 1] - constant.mappedFrequencies[x]) / topFrequency;
						// as long as the bandwidth is smaller than our fft resolution, we interpolate the points
						// otherwise, break out and sample the max values of the bins inside the bandwidth
						if (bwForLine > fftBandwidth)
							break;

						// +0.5 to centerly space bins.
						auto index = Math::confineTo((std::size_t)(constant.mappedFrequencies[x] * freqToBin + 0.5), 0, numBins - 1);
						auto iLeft = csf[index];
						auto iRight = csf[N - index];

						auto cancellation = invSize * std::abs(iLeft + iRight);
						auto mid = invSize * (std::abs(iLeft) + std::abs(iRight));

						wsp[x * 2] = mid;

						wsp[x * 2 + 1] = T(1) - (mid > 0 ? (cancellation / mid) : 0);
					}
					break;
				}


				// the process after interpolation is much simpler, as we dont have to account
				// for wrongly interpolation of phase-mangled vectors.
				if (x < constant.axisPoints)
					oldBin = constant.mappedFrequencies[x] * freqToBin;



				for (; x < constant.axisPoints; ++x)
				{
					std::size_t maxBin = 0;
					T maxValue = 0, newMag = 0;
					bin = static_cast<std::size_t>(constant.mappedFrequencies[x] * freqToBin);
#ifdef DEBUG
					if ((std::size_t)bin >= constant.transformSize)
						CPL_RUNTIME_EXCEPTION("Corrupt frequency mapping!");
#endif

					auto diff = bin - oldBin;
					auto counter = diff ? 1 : 0;
					// here we loop over all the bins that is mapped for a single coordinate
					do
					{
						auto offset = oldBin + counter;
						newMag = std::max(Math::square(csf[offset]), Math::square(csf[N - offset]));
						// select highest number in this chunk for display. Not exactly correct, though.
						if (newMag > maxValue)
						{
							maxValue = newMag;
							maxBin = oldBin + counter;
						}
						counter++;
						diff--;
					} while (diff > 0);

					leftMax = csf[maxBin];
					rightMax = csf[N - maxBin];

					auto interference = invSize * std::abs(leftMax + rightMax);
					auto mid = invSize * (std::abs(leftMax) + std::abs(rightMax));
					auto cancellation = interference / mid;
					wsp[x * 2] = mid;
					wsp[x * 2 + 1] = T(1) - (mid > 0 ? cancellation : 0);

					oldBin = bin;
				}

				break;
			}
			case SpectrumChannels::Separate:
			case SpectrumChannels::MidSide:
			{
				// two-for-one pass, first channel is 0... N/2 -1, second is N/2 .. N -1
				dsp::separateTransformsIPL(csf.slice(0, N));

				// fix up DC and nyquist bins (see previous function documentation)
				csf[N] = csf[0].imag() * 0.5;
				csf[0] = csf[0].real() * 0.5;
				csf[N >> 1] *= 0.5;
				csf[(N >> 1) - 1] *= 0.5;

				for (decltype(N) i = 1; i < N; ++i)
				{
					csf[i] = std::abs(csf[i]);
				}

				// The index of the transform, where the bandwidth is higher than mapped pixels (so no more interpolation is needed)
				// TODO: This can be calculated from view mapping scale and N pixels.
				std::size_t bandWidthBreakingPoint = constant.axisPoints - 1;

				double fftBandwidth = 1.0 / numBins;
				//double pxlBandwidth = 1.0 / numPoints;
				cpl::Types::fint_t x;
				switch (constant.binPolation)
				{
				case SpectrumContent::BinInterpolation::Linear:

					// phase pass
					for (x = 0; x < constant.axisPoints - 1; ++x)
					{
						double bwForLine = (constant.mappedFrequencies[x + 1] - constant.mappedFrequencies[x]) / topFrequency;
						if (bwForLine > fftBandwidth)
						{
							bandWidthBreakingPoint = x;
							break;
						}

						auto iLeft = dsp::linearFilter<std::complex<T>>(csf, constant.mappedFrequencies[x] * freqToBin);
						auto iRight = dsp::linearFilter<std::complex<T>>(csf, N - (constant.mappedFrequencies[x] * freqToBin));

						csp[x] = invSize * iLeft;
						csp[constant.axisPoints + x] = invSize * iRight;
					}

					break;
				case SpectrumContent::BinInterpolation::Lanczos:
					for (x = 0; x < constant.axisPoints - 1; ++x)
					{

						double bwForLine = (constant.mappedFrequencies[x + 1] - constant.mappedFrequencies[x]) / topFrequency;
						if (bwForLine > fftBandwidth)
						{
							bandWidthBreakingPoint = x;
							break;
						}

						auto iLeft = dsp::lanczosFilter<std::complex<T>, true>(csf, constant.mappedFrequencies[x] * freqToBin, 5);
						auto iRight = dsp::lanczosFilter<std::complex<T>, true>(csf, N - (constant.mappedFrequencies[x] * freqToBin), lanczosFilterSize);

						csp[x] = invSize * iLeft;
						csp[constant.axisPoints + x] = invSize * iRight;
					}

					break;
				default:
					for (x = 0; x < constant.axisPoints - 1; ++x)
					{
						// bandwidth of the filter for this 'line', or point
						double bwForLine = (constant.mappedFrequencies[x + 1] - constant.mappedFrequencies[x]) / topFrequency;
						// as long as the bandwidth is smaller than our fft resolution, we interpolate the points
						// otherwise, break out and sample the max values of the bins inside the bandwidth
						if (bwForLine > fftBandwidth)
							break;

						// +0.5 to centerly space bins.
						auto index = Math::confineTo((std::size_t)(constant.mappedFrequencies[x] * freqToBin + 0.5), 0, numBins);

						csp[x] = invSize * csf[index];
						csp[constant.axisPoints + x] = invSize * csf[N - index];
					}
					break;
				}


				// the process after interpolation is much simpler, as we dont have to account
				// for wrongly interpolation of phase-mangled vectors.

				oldBin = constant.mappedFrequencies[x] * freqToBin;

				for (; x < constant.axisPoints; ++x)
				{
					maxLMag = maxRMag = newLMag = newRMag = 0;

					bin = static_cast<std::size_t>(constant.mappedFrequencies[x] * freqToBin);
#ifdef DEBUG
					if ((std::size_t)bin >= constant.transformSize)
						CPL_RUNTIME_EXCEPTION("Corrupt frequency mapping!");
#endif
					maxRBin = maxLBin = bin;

					auto diff = bin - oldBin;
					auto counter = diff ? 1 : 0;
					// here we loop over all the bins that is mapped for a single coordinate
					do
					{
						auto offset = oldBin + counter;
						//offset <<= 1;
						newLMag = Math::square(csf[offset]);
						newRMag = Math::square(csf[N - offset]);
						// select highest number in this chunk for display. Not exactly correct, though.
						if (newLMag > maxLMag)
						{
							maxLBin = oldBin + counter;
							maxLMag = newLMag;
						}

						if (newRMag > maxRMag)
						{
							maxRBin = N - (oldBin + counter);
							maxRMag = newRMag;
						}

						counter++;
						diff--;
					} while (diff > 0);

					csp[x] = invSize * csf[maxLBin];
					csp[constant.axisPoints + x] = invSize * csf[maxRBin];
					oldBin = bin;
				}
			}
			break;
			case SpectrumChannels::Complex:
			{
				// two-for-one pass, first channel is 0... N/2 -1, second is N/2 .. N -1

				// fix up DC and nyquist bins (see previous function documentation)
				//csf[N] = csf[0].imag() * 0.5;
				csf[0] *= (T)0.5;

				double fftBandwidth = 1.0 / (numBins * 2);
				//double pxlBandwidth = 1.0 / numPoints;
				cpl::Types::fint_t x = 0;

				for (decltype(N) i = 1; i < N; ++i)
				{
					csf[i] = std::abs(csf[i]);
				}

				while (x < constant.axisPoints)
				{

					switch (constant.binPolation)
					{
					case SpectrumContent::BinInterpolation::Linear:
					{
						for (; x < constant.axisPoints; ++x)
						{
							if (x != constant.axisPoints - 1)
							{
								double bwForLine = (constant.mappedFrequencies[x + 1] - constant.mappedFrequencies[x]) / topFrequency;
								if (bwForLine > fftBandwidth)
									break;
							}
							csp[x] = invSize * dsp::linearFilter<std::complex<T>>(csf, constant.mappedFrequencies[x] * freqToBin);
						}
						break;
					}
					case SpectrumContent::BinInterpolation::Lanczos:
					{
						for (; x < constant.axisPoints; ++x)
						{
							if (x != constant.axisPoints - 1)
							{
								double bwForLine = (constant.mappedFrequencies[x + 1] - constant.mappedFrequencies[x]) / topFrequency;
								if (bwForLine > fftBandwidth)
									break;
							}

							csp[x] = invSize * dsp::lanczosFilter<std::complex<T>, true>(csf, constant.mappedFrequencies[x] * freqToBin, lanczosFilterSize);
						}

						break;
					}
					default:
						for (; x < constant.axisPoints; ++x)
						{
							if (x != constant.axisPoints - 1)
							{
								double bwForLine = (constant.mappedFrequencies[x + 1] - constant.mappedFrequencies[x]) / topFrequency;
								if (bwForLine > fftBandwidth)
									break;
							}
							// +0.5 to centerly space bins.
							auto index = Math::confineTo((std::size_t)(constant.mappedFrequencies[x] * freqToBin + 0.5), 0, N);
							csp[x] = invSize * csf[index];
						}
						break;
					}
					if (x != constant.axisPoints)
						oldBin = constant.mappedFrequencies[x] * freqToBin;

					for (; x < constant.axisPoints; ++x)
					{
						maxLMag = maxRMag = newLMag = newRMag = 0;

						bin = static_cast<std::size_t>(constant.mappedFrequencies[x] * freqToBin);
#ifdef DEBUG
						if ((std::size_t)bin >= constant.transformSize)
							CPL_RUNTIME_EXCEPTION("Corrupt frequency mapping!");
#endif
						maxRBin = maxLBin = bin;
						if (x != constant.axisPoints - 1)
						{
							double bwForLine = (constant.mappedFrequencies[x + 1] - constant.mappedFrequencies[x]) / topFrequency;
							if (bwForLine < fftBandwidth)
								break;
						}

						auto diff = bin - oldBin;
						auto counter = diff ? 1 : 0;
						// here we loop over all the bins that is mapped for a single coordinate
						do
						{
							auto offset = oldBin + counter;
							//offset <<= 1;
							newLMag = Math::square(csf[offset]);
							// select highest number in this chunk for display. Not exactly correct, though.
							if (newLMag > maxLMag)
							{
								maxLBin = oldBin + counter;
								maxLMag = newLMag;
							}

							counter++;
							diff--;
						} while (diff > 0);

						csp[x] = invSize * csf[maxLBin];
						oldBin = bin;
					}
				}
			}

			break;
			}
			break;
		}
		case SpectrumContent::TransformAlgorithm::RSNT:
		{
			auto configurationChannels = constant.getStateConfigurationChannels();
			auto wsp = getWork<std::complex<float>>(constant.resonator.getNumFilters() * configurationChannels);
			std::size_t filtersPerChannel = copyResonatorStateInto<T>(constant, constant.dspWindow, wsp, configurationChannels) / configurationChannels;

			switch (constant.configuration)
			{
			case SpectrumChannels::Phase:
			{
				for (std::size_t x = 0; x < filtersPerChannel; ++x)
				{

					auto iLeft = wsp[x];
					auto iRight = wsp[x + filtersPerChannel];

					auto cancellation = std::sqrt(Math::square(iLeft + iRight));
					auto mid = std::abs(iLeft) + std::abs(iRight);


					wsp[x] = std::complex<float>(mid, AFloat(1) - (mid > 0 ? (cancellation / mid) : 0));

				}

				break;
			}
			// rest of cases do not need any handling
			}

		}
		break;
		}
	}

	template<typename T>
	template<typename ISA>
	inline void TransformPair<T>::addAudioFrame(const Constant& constant)
	{
		mapToLinearSpace(constant);
		postProcessStdTransform(constant);

		constexpr auto i = SpectrumContent::LineGraphs::LineMain;

		auto frame = ComplexSpectrumFrame(lineGraphs[i].results.begin(), lineGraphs[i].results.begin() + constant.axisPoints /* channels ? */);
		sfbuf.emplace_back(std::move(frame));
	}

	template<typename T>
	template<typename V, class Vector>
	inline std::size_t TransformPair<T>::copyResonatorStateInto(const Constant& constant, cpl::dsp::WindowTypes windowType, Vector& output, std::size_t outChannels)
	{
		auto numResFilters = constant.resonator.getNumFilters();
		// casts from std::complex<T> * to T * which is well-defined.
		// TODO: std::ranges
		auto buf = reinterpret_cast<AFloat*>(cpl::data(output));
		cresonator.getWholeWindowedState<V>(constant.resonator, windowType, buf, outChannels, numResFilters);

		return numResFilters << (outChannels - 1);
	}

	template<typename T>
	template<typename ISA>
	inline void TransformPair<T>::audioEntryPoint(const Constant& constant, const std::optional<AudioPair>& views, std::array<AFloat*, 2> buffer, std::size_t numSamples, bool historyMayBeDeferred)
	{
		if (constant.displayMode == SpectrumContent::DisplayMode::ColourSpectrum)
		{
			std::int64_t n = numSamples;
			std::size_t offset = 0;

			while (n > 0)
			{
				std::int64_t numRemainingSamples = processedSamplesSinceLastFrame > constant.sampleBufferSize ? 0 : constant.sampleBufferSize - processedSamplesSinceLastFrame;
				const auto availableSamples = numRemainingSamples + std::min(std::int64_t(0), n - numRemainingSamples);

				// do some resonation
				if (constant.algo == SpectrumContent::TransformAlgorithm::RSNT)
				{
					resonatingDispatch<ISA>(constant, { buffer[0] + offset, buffer[1] + offset }, availableSamples);
				}

				processedSamplesSinceLastFrame += availableSamples;

				if (processedSamplesSinceLastFrame >= constant.sampleBufferSize)
				{
					bool transformReady = true;
					if (constant.algo == SpectrumContent::TransformAlgorithm::FFT)
					{
						// the abstract timeline consists of the old data in the audio stream, with the following audio presented in this function.
						// thus, the more we include of the buffer ('offbuf') the newer the data segment gets.
						if ((transformReady = prepareTransform(constant, *views, buffer, !historyMayBeDeferred ? availableSamples + offset : 0)))
							doTransform(constant);
					}

					if (transformReady)
						addAudioFrame<ISA>(constant);

					processedSamplesSinceLastFrame = 0;
				}

				offset += availableSamples;
				n -= availableSamples;
			}
		}
		else if (constant.algo == SpectrumContent::TransformAlgorithm::RSNT)
		{
			resonatingDispatch<ISA>(constant, buffer, numSamples);
		}

	}

	template<typename T>
	template<typename ISA>
	inline void TransformPair<T>::resonatingDispatch(const Constant& constant, std::array<AFloat*, 2> buffer, std::size_t numSamples)
	{
		typedef AFloat TReso;

		auto one = [&](auto fn)
		{
			auto work = this->getWork<T>(numSamples);

			for (std::size_t i = 0; i < numSamples; ++i)
			{
				work[i] = static_cast<T>(fn(buffer[0][i], buffer[1][i]));
			}

			return std::array{ work };
		};

		auto two = [&](auto fn)
		{
			auto work = this->getWork<T>(numSamples * 2);

			for (std::size_t i = 0; i < numSamples; ++i)
			{
				auto result = fn(buffer[0][i], buffer[1][i]);
				work[i] = static_cast<T>(result.first);
				work[i + numSamples] = static_cast<T>(result.second);
			}

			return std::array{ work.slice(0, numSamples), work.slice(numSamples, numSamples) };
		};

		// TODO: asserts?
		if (numSamples < 1)
			return;

		switch (constant.configuration)
		{
		case SpectrumChannels::Right:
		{
			auto work = one([](AFloat left, AFloat right) { return right; });
			cresonator.resonateReal<typename ISA::V>(constant.resonator, work, 1, numSamples);			
			break;
		}
		case SpectrumChannels::Left:
		{
			auto work = one([](AFloat left, AFloat right) { return left; });
			cresonator.resonateReal<typename ISA::V>(constant.resonator, work, 1, numSamples);
			break;
		}
		case SpectrumChannels::Mid:
		{
			auto work = one([](AFloat left, AFloat right) { return left + right; });
			cresonator.resonateReal<typename ISA::V>(constant.resonator, work, 1, numSamples);
			break;
		}
		case SpectrumChannels::Side:
		{
			auto work = one([](AFloat left, AFloat right) { return left - right; });
			cresonator.resonateReal<typename ISA::V>(constant.resonator, work, 1, numSamples);
			break;
		}
		case SpectrumChannels::MidSide:
		{
			auto work = two([](AFloat left, AFloat right) { return std::pair{ left - right, left + right }; });
			cresonator.resonateReal<typename ISA::V>(constant.resonator, work, 2, numSamples);
			break;
		}
		case SpectrumChannels::Phase:
		case SpectrumChannels::Separate:
		{
			auto work = two([](AFloat left, AFloat right) { return std::pair{ left, right }; });
			cresonator.resonateReal<typename ISA::V>(constant.resonator, work, 2, numSamples);
			break;
		}
		case SpectrumChannels::Complex:
		{
			auto work = two([](AFloat left, AFloat right) { return std::pair{ left, right }; });
			cresonator.resonateReal<typename ISA::V>(constant.resonator, work, 2, numSamples);
			break;
		}
		}
	}

	template<typename T>
	template<class V2>
	void TransformPair<T>::mapAndTransformDFTFilters(const Constant& constant, const V2& newVals, std::size_t size)
	{
		CPL_RUNTIME_ASSERTION(size == constant.axisPoints);

		for (std::size_t k = 0; k < lineGraphs.size(); ++k)
		{
			lineGraphs[k].resize(size);
		}

		double lowerFraction = cpl::Math::dbToFraction<double>(constant.lowDBs);
		double upperFraction = cpl::Math::dbToFraction<double>(constant.highDBs);

		auto deltaYRecip = static_cast<T>(1.0 / std::log(upperFraction / lowerFraction));
		auto minFracRecip = static_cast<T>(1.0 / lowerFraction);

		T lowerClip = static_cast<T>(constant.clipDB);

		switch (constant.configuration)
		{
		case SpectrumChannels::Left:
		case SpectrumChannels::Merge:
		case SpectrumChannels::Right:
		case SpectrumChannels::Side:
		case SpectrumChannels::Complex:
		{

			for (cpl::Types::fint_t i = 0; i < size; ++i)
			{

				auto newReal = newVals[i * 2];
				auto newImag = newVals[i * 2 + 1];
				// mag = abs(cmplx)
				auto magnitude = sqrt(newReal * newReal + newImag * newImag);


				for (std::size_t k = 0; k < lineGraphs.size(); ++k)
				{
					lineGraphs[k].states[i].magnitude *= constant.filter[k].pole;

					if (magnitude > lineGraphs[k].states[i].magnitude)
					{
						lineGraphs[k].states[i].magnitude = static_cast<T>(magnitude);
					}

					auto deltaX = constant.slopeMap[i] * lineGraphs[k].states[i].magnitude * minFracRecip;
					// deltaX mostly zero here - add simd check
					auto result = deltaX > 0 ? std::log(deltaX) * deltaYRecip : lowerClip;
					lineGraphs[k].results[i].magnitude = static_cast<T>(result);
					lineGraphs[k].results[i].phase = 0;
				}

			}
			break;
		}
		case SpectrumChannels::Separate:
		case SpectrumChannels::MidSide:
		{

			for (cpl::Types::fint_t i = 0; i < size; ++i)
			{

				auto lreal = newVals[i * 2];
				auto rreal = newVals[i * 2 + size * 2];
				auto limag = newVals[i * 2 + 1];
				auto rimag = newVals[i * 2 + size * 2 + 1];
				// mag = abs(cmplx)
				auto lmag = sqrt(lreal * lreal + limag * limag);
				auto rmag = sqrt(rreal * rreal + rimag * rimag);

				for (std::size_t k = 0; k < lineGraphs.size(); ++k)
				{
					lineGraphs[k].states[i].leftMagnitude *= constant.filter[k].pole;
					lineGraphs[k].states[i].rightMagnitude *= constant.filter[k].pole;

					if (lmag > lineGraphs[k].states[i].leftMagnitude)
					{
						lineGraphs[k].states[i].leftMagnitude = static_cast<T>(lmag);
					}
					if (rmag > lineGraphs[k].states[i].rightMagnitude)
					{
						lineGraphs[k].states[i].rightMagnitude = static_cast<T>(rmag);
					}
					// log10(y / _min) / log10(_max / _min);
					auto deltaLX = constant.slopeMap[i] * lineGraphs[k].states[i].leftMagnitude * minFracRecip;
					auto deltaRX = constant.slopeMap[i] * lineGraphs[k].states[i].rightMagnitude * minFracRecip;
					// deltaX mostly zero here - add simd check
					auto lResult = deltaLX > 0 ? std::log(deltaLX) * deltaYRecip : lowerClip;
					auto rResult = deltaRX > 0 ? std::log(deltaRX) * deltaYRecip : lowerClip;
					lineGraphs[k].results[i].leftMagnitude = static_cast<T>(lResult);
					lineGraphs[k].results[i].rightMagnitude = static_cast<T>(rResult);
				}
			}
			break;
		}
		case SpectrumChannels::Phase:
		{
			T phaseFilters[SpectrumContent::LineGraphs::LineEnd];

			for (std::size_t k = 0; k < lineGraphs.size(); ++k)
				phaseFilters[k] = std::pow<T>(constant.filter[k].pole, 0.3);

			for (cpl::Types::fint_t i = 0; i < size; ++i)
			{

				auto mag = newVals[i * 2];
				auto phase = newVals[i * 2 + 1];
				// mag = abs(cmplx)

				mag *= consts::half;

				for (std::size_t k = 0; k < lineGraphs.size(); ++k)
				{
					lineGraphs[k].states[i].magnitude *= constant.filter[k].pole;

					if (mag > lineGraphs[k].states[i].magnitude)
					{
						lineGraphs[k].states[i].magnitude = static_cast<T>(mag);
					}
					phase *= mag;

					lineGraphs[k].states[i].phase = phase + phaseFilters[k] * (lineGraphs[k].states[i].phase - phase);


					// log10(y / _min) / log10(_max / _min);
					// deltaX mostly zero here - add simd check
					auto deltaX = constant.slopeMap[i] * lineGraphs[k].states[i].magnitude * minFracRecip;
					auto deltaY = constant.slopeMap[i] * lineGraphs[k].states[i].phase * minFracRecip;
					// deltaX mostly zero here - add simd check
					auto result = deltaX > 0 ? std::log(deltaX) * deltaYRecip : lowerClip;
					lineGraphs[k].results[i].magnitude = static_cast<T>(result);
					lineGraphs[k].results[i].phase = static_cast<T>(deltaY > 0 ? std::log(deltaY) * deltaYRecip : lowerClip);
				}
			}
			break;
		}
		};
	}

	template<typename T>
	void TransformPair<T>::postProcessStdTransform(const Constant& constant)
	{
		mapAndTransformDFTFilters(constant, getTransformResult(constant), constant.axisPoints);
	}

	template<typename T>
	cpl::uarray<const T> TransformPair<T>::getTransformResult(const Constant& constant)
	{
		return getWork<T>(constant.axisPoints * constant.getStateConfigurationChannels() * 2);
	}
}

#endif
