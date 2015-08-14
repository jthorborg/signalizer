#include "CSpectrum.h"
#include <cpl/ffts.h>
#include <cpl/SysStats.h>
#include <cpl/lib/LockFreeDataQueue.h>

//#define _ZERO_PHASE

namespace Signalizer
{

	template<typename T>
		T * CSpectrum::getAudioMemory()
		{

			return reinterpret_cast<T*>(audioMemory.data());
		}

	template<typename T>
		std::size_t CSpectrum::getNumAudioElements() const noexcept
		{
			return audioMemory.size() / sizeof(T);
		}


	template<typename T>
		T * CSpectrum::getWorkingMemory()
		{
			return reinterpret_cast<T*>(workingMemory.data());
		}

	template<typename T>
		std::size_t CSpectrum::getNumWorkingElements() const noexcept
		{
			return workingMemory.size() / sizeof(T);
		}

	int CSpectrum::getBlobSamples() const noexcept
	{
		return state.audioBlobSizeMs * 0.001 * getSampleRate();
	}

	void CSpectrum::resetState()
	{
		flags.resetStateBuffers = true;
	}

	void CSpectrum::audioConsumerThread()
	{
		// signal that we exist.
		audioThreadIsInitiated.store(true);

		AudioFrame recv;
		int pops(20);
		std::vector<fpoint> audioInput[2];

		typedef decltype(cpl::Misc::ClockCounter()) ctime_t;

		// when it returns false, its time to quit this thread.
		while (audioFifo.popElementBlocking(recv))
		{
			ctime_t then = cpl::Misc::ClockCounter();

			// each time we get into here, it's very likely there's abunch of messages waiting.
			auto numExtraEntries = audioFifo.enqueuededElements();

			// always resize queue before emptying
			if (pops++ > 10)
			{
				audioFifo.grow(0, true, 0.3f, 2);
				pops = 0;
			}

			for (auto & ai : audioInput)
				ai.resize(AudioFrame::capacity * (1 + numExtraEntries)); // worst case possible

			audioInput[0].insert(audioInput[0].begin(), recv.begin(), recv.end());
			std::size_t numSamples = recv.size;

			for (size_t i = 0; i < numExtraEntries; i++)
			{
				if (!audioFifo.popElementBlocking(recv))
					return;
				audioInput[0].insert(audioInput[0].begin() + numSamples, recv.begin(), recv.end());
				numSamples += recv.size;
			}

			// process

			if (!flags.firstChange)
			{
				float * buffers[2] = { audioInput[0].data(), audioInput[1].data() };

				auto & cpuinfo = cpl::SysStats::CProcessorInfo::instance();

				_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

				if (cpuinfo.test(cpuinfo.AVX))
				{
					// size of magic here to support floats and doubles through changing the typename
					audioProcessing<typename cpl::simd::vector_of<fpoint, 8 / (sizeof(fpoint) / 4)>::type>(buffers, 2, numSamples);
				}
				else if (cpuinfo.test(cpuinfo.SSE2))
				{
					audioProcessing<typename cpl::simd::vector_of<fpoint, 4 / (sizeof(fpoint) / 4)>::type>(buffers, 2, numSamples);
				}
				else
				{
					audioProcessing<typename cpl::simd::vector_of<fpoint, 1 / (sizeof(fpoint) / 4)>::type>(buffers, 2, numSamples);
				}
			}


			// clear buffers.


			ctime_t elapsed = cpl::Misc::ClockCounter() - then;


			double fractionOfSecond = numSamples / getSampleRate();
			double cpuUsedOfSecond = (elapsed * 1e-4) / processorSpeed;

			double cpuUsedCurrently = cpuUsedOfSecond / fractionOfSecond;

			// lowpass it a bit.
			//this->audioThreadUsage = cpuUsedCurrently + 0.9 * (this->audioThreadUsage - cpuUsedCurrently);
			this->audioThreadUsage = cpuUsedCurrently;
		}

	}

	void CSpectrum::computeWindowKernel()
	{
		size_t sampleSize = getWindowSize();
		std::size_t i = 0;


		switch (state.dspWindow)
		{
		default:
			case cpl::dsp::WindowTypes::Hann:
			{
				for (; i < sampleSize; ++i)
				{
					windowKernel[i] = 1.0 - std::cos(TAU * i / (sampleSize));
				}
				break;
			}
			case cpl::dsp::WindowTypes::Rectangular:
			case cpl::dsp::WindowTypes::Hamming:
			case cpl::dsp::WindowTypes::FlatTop:
			case cpl::dsp::WindowTypes::Blackman:
			case cpl::dsp::WindowTypes::Triangular:
			case cpl::dsp::WindowTypes::Parzen:
			case cpl::dsp::WindowTypes::Nuttall:
			case cpl::dsp::WindowTypes::BlackmanNutall:
			case cpl::dsp::WindowTypes::BlackmanHarris:
			case cpl::dsp::WindowTypes::Gaussian:
			case cpl::dsp::WindowTypes::Slepian:
			case cpl::dsp::WindowTypes::DolphChebyshev:
			case cpl::dsp::WindowTypes::Kaiser:
			case cpl::dsp::WindowTypes::Ultraspherical:
			case cpl::dsp::WindowTypes::End:
				for (; i < sampleSize; ++i)
				{
					windowKernel[i] = 1.0;
				}
			break;
		}


		size_t fullSize = getNumAudioElements<std::complex<double>>();
		// zero-padding
		for (; i < fullSize; ++i)
		{
			windowKernel[i] = 0.0;
		}
		// uncomment if you want to use window for convolution
		//cpl::dsp::CSignalTransform::sfft(reinterpret_cast<double*>(windowKernel.data()), fullSize);
	}



	void CSpectrum::prepareTransform()
	{
		// we can 'freeze' the window pretty effectively, by not allowing new input
		// and still be able to zoom around. 
		if (state.isFrozen)
			return;
		const bool locking = shouldSynchronize();
		auto size = getWindowSize(); // the size of the transform, containing samples
									 // the quantized (to next power of 2) samples of this transform
									 // that is, the size + additional zero-padding
		auto fullSize = getNumAudioElements<std::complex<double>>();

		auto const channelConfiguration = kchannelConfiguration.getZeroBasedSelIndex<ChannelConfiguration>();

		auto numChannels = 1 + (channelConfiguration > ChannelConfiguration::Merge);

		{
			// acquire locks on audio buffegrs
			std::vector<cpl::CMutex> locks;
			if (locking)
				for (auto & buf : audioStream)
					locks.emplace_back(buf);

			// note here, that since our audiobuffer is a circular buffer,
			// index[0] is actually the last input, and index[size - 1] is the first
			// while the fft doesn't care whether the input is reversed or not, 
			// the varying Q-transforms do care, since they vary the amount of samples
			// considered.

			switch (state.algo)
			{
			case TransformAlgorithm::FFT:
			{
				auto buffer = getAudioMemory<std::complex<double>>();
				std::size_t channel = 1;
				std::size_t halfSize = fullSize >> 1;
				std::size_t halfPadding = (fullSize - size) >> 1;
				std::size_t halfInputSize = size >> 1;
				std::size_t i = 0, offset = halfSize + halfPadding;


				/*
					Most people zero-pad the end of the input buffer to the FFT.
					This is okay if you don't do any peak/phase interpolation. In reality,
					you should shift the input half a period + number of padding in the buffer,
					and apply the zeroes _exactly_ in the middle. This preserves the phases, which are
					critical for interpolation methods. See more here :
					https://ccrma.stanford.edu/~jos/parshl/Filling_FFT_Input_Buffer.html
				*/

				switch (channelConfiguration)
				{
				case ChannelConfiguration::Left:
					channel = 0;
				case ChannelConfiguration::Right:
				{
#ifdef _ZERO_PHASE
					for (; i < halfInputSize; ++i)
					{
						buffer[i + offset] = audioStream[channel].singleCheckAccess(i) * windowKernel[i];
					}

					offset = halfInputSize;
					for (; i < size; ++i)
					{
						buffer[i - offset] = audioStream[channel].singleCheckAccess(i) * windowKernel[i];
					}
					break;
#else
					for (; i < size; ++i)
					{
						buffer[i] = audioStream[channel].singleCheckAccess(i) * windowKernel[i];
					}
#endif
				}
				case ChannelConfiguration::Phase:
				case ChannelConfiguration::Separate:
				{
					for (; i < halfInputSize; ++i)
					{
						buffer[i + offset] = std::complex<double>
						(
							audioStream[0].singleCheckAccess(i)  * windowKernel[i],
							audioStream[1].singleCheckAccess(i)  * windowKernel[i]
						);
					}

					offset = halfInputSize;
					for (; i < size; ++i)
					{
						buffer[i - offset] = std::complex<double>
						(
							audioStream[0].singleCheckAccess(i) * windowKernel[i],
							audioStream[1].singleCheckAccess(i) * windowKernel[i]
						);
					}
					break;
				}
				case ChannelConfiguration::Merge:
					for (; i < halfInputSize; ++i)
					{
						buffer[i + offset] = (audioStream[0].singleCheckAccess(i) + audioStream[1].singleCheckAccess(i)) * windowKernel[i] * 0.5;
					}

					offset = halfInputSize;
					for (; i < size; ++i)
					{
						buffer[i - offset] = (audioStream[0].singleCheckAccess(i) + audioStream[1].singleCheckAccess(i)) * windowKernel[i] * 0.5;;
					}
				}

				/*
				zero-pad until buffer is filled
				*/
#ifdef _ZERO_PHASE
				offset = halfSize + halfPadding;
				for (size_t pad = halfSize - halfPadding; pad < offset; ++pad)
				{
					buffer[pad] = 0;
				}
#else
				for (size_t pad = i; pad < fullSize; ++pad)
				{
					buffer[pad] = 0;
				}
#endif

				break;
			}
			/*// this case is different from CDFT, since the input musnt be windowed
			case Algorithm::MQDFT:
			{
				auto buffer = getAudioMemory<float>();
				size = getWindowSize();
				std::size_t channel = 1;
				winOsc.reset(size, 1, M_PI / 2); // offset phase by 90 degrees to create a cosine instead.
												 // this is used for reversing the index.
				auto const N = size - 1;
				switch (channelConfiguration)
				{
				case ChannelConfiguration::Left:
					channel = 0;
				case ChannelConfiguration::Right:
					for (unsigned i = 0; i < size; ++i)
					{
						buffer[N - i] = audioData[channel].singleCheckAccess(i);
					}
					break;
				case ChannelConfiguration::Phase:
				case ChannelConfiguration::Separate:
					for (unsigned n = 0; n < numChannels; ++n)
					{
						for (unsigned i = 0; i < size; ++i)
						{
							buffer[(N - i + n * size)] = audioData[n].singleCheckAccess(i);
						}
					}
					break;
				case ChannelConfiguration::Merge:
					for (unsigned i = 0; i < size; ++i)
					{
						buffer[N - i] = (audioData[0].singleCheckAccess(i) + audioData[1].singleCheckAccess(i)) * 0.5;
					}
					break;
				}
				break;
			}*/
			}
		}
	}

	void CSpectrum::doTransform()
	{
		auto const channelConfiguration = kchannelConfiguration.getZeroBasedSelIndex<ChannelConfiguration>();
		switch (state.algo)
		{
		case TransformAlgorithm::FFT:
		{
			auto const numSamples = getNumAudioElements<std::complex<double>>();
			switch (channelConfiguration)
			{
			case ChannelConfiguration::Left:
			case ChannelConfiguration::Right:
			case ChannelConfiguration::Merge:
				cpl::signaldust::DustFFT_fwdDa(getAudioMemory<double>(), numSamples);
				break;
			case ChannelConfiguration::Phase:
			case ChannelConfiguration::Separate:
				cpl::signaldust::DustFFT_fwdDa(getAudioMemory<double>(), numSamples);

			}

			break;
		}
		/*case Algorithm::MQDFT:
		{
			auto size = getWindowSize();
			switch (channelConfiguration)
			{
			case ChannelConfiguration::Left:
			case ChannelConfiguration::Right:
			case ChannelConfiguration::Merge:
				transformer.mqdft<1>(getAudioMemory<float>(), size);
				break;
			case ChannelConfiguration::Phase:
			case ChannelConfiguration::Separate:
				transformer.mqdft<2>(getAudioMemory<float>(), size);

			}
			break;
		}
		// Resonators dont need preparation as they are continious and real time
		case Algorithm::RSNT:
			break;*/
		}
	}

	template<typename Ty>
		struct DualComplex
		{
			typedef Ty type;

			std::complex<type> val[2];

		};
	template<typename Ty>
		inline DualComplex<Ty> getZFromNFFT(Ty * tsf, std::size_t idx, std::size_t N)
		{
			idx <<= 1;
			N <<= 1;
			Ty x1 = tsf[idx];
			Ty x2 = tsf[N - idx];
			Ty y1 = tsf[idx + 1];
			Ty y2 = tsf[N - idx + 1];

			DualComplex<Ty> ret;
			//ret.val[0] = std::complex<Ty>((x1 + x2) * 0.5, (y1 + y2) * 0.5);
			//ret.val[1] = std::complex<Ty>((y1 - y2) * 0.5, -(x1 - x2) * 0.5);

			ret.val[0] = std::complex<Ty>((x1 + x2) * 0.5, (y1 - y2) * 0.5);
			ret.val[1] = std::complex<Ty>((y1 + y2) * 0.5, -(x1 - x2) * 0.5);

			return ret;
		}


	/*

		all inputs must be scaled.

		oldVals = state. First time must be zero. vector of floats of size. changed during call
		newVals = current vector of floats * 2 (complex), output from CSignalTransform::**dft() of size * 2
		output = vector of single floats, logarithmically mapped to 0 - 1 range, of size

		oldVals will be changed during the call.

			for mode = left
				newVals is a complex vector of floats of size
			for mode = merge, phase, dual
				newVals is a complex vector of floats of size * 2, ChannelConfiguration::Separated channels
				newVals[n * 2 + 0] = lreal
				newVals[n * 2 + 1] = limag
				newVals[n * 2 + size + 0] = rreal
				newVals[n * 2 + size + 1] = rimag

	*/
	template<CSpectrum::ChannelConfiguration type, class Scalar, class V2>
		void mapAndTransformDFTFilters(UComplexFilter<Scalar> * __RESTRICT__ oldVals, 
			const V2 & newVals, UComplexFilter<Scalar> * __RESTRICT__ output, std::size_t size,
			float lowDbs, float highDbs, cpl::CPeakFilter<Scalar> filter)
		{

			double lowerFraction = cpl::Math::dbToFraction<double>(lowDbs);
			double upperFraction = cpl::Math::dbToFraction<double>(highDbs);
			auto deltaYRecip = static_cast<Scalar>(1.0 / log(upperFraction / lowerFraction));
			auto minFracRecip = static_cast<Scalar>(1.0 / lowerFraction);
			auto halfRecip = Scalar(0.5);

			switch (type)
			{
			case CSpectrum::ChannelConfiguration::Left:
			case CSpectrum::ChannelConfiguration::Merge:
			case CSpectrum::ChannelConfiguration::Right:
			{

				for (cpl::Types::fint_t i = 0; i < size; ++i)
				{

					auto newReal = newVals[i * 2];
					auto newImag = newVals[i * 2 + 1];
					// mag = abs(cmplx)
					auto magnitude = sqrt(newReal * newReal + newImag * newImag);

					oldVals[i].magnitude *= filter.pole;

					if (magnitude > oldVals[i].magnitude)
					{
						oldVals[i].magnitude = magnitude;
					} 
					// log10(y / _min) / log10(_max / _min);

					auto deltaX = oldVals[i].magnitude * minFracRecip;
					// deltaX mostly zero here - add simd check
					auto result = deltaX ? log(deltaX) * deltaYRecip : 0;
					output[i].magnitude = result;
					output[i].phase = 0;
				}
				break;
			}
			case CSpectrum::ChannelConfiguration::Separate:
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

					oldVals[i].leftMagnitude *= filter.pole;
					oldVals[i].rightMagnitude *= filter.pole;

					if (lmag > oldVals[i].leftMagnitude)
					{
						oldVals[i].leftMagnitude = lmag;
					} 
					if (rmag > oldVals[i].rightMagnitude)
					{
						oldVals[i].rightMagnitude = rmag;
					}
					// log10(y / _min) / log10(_max / _min);
					auto deltaLX = oldVals[i].leftMagnitude * minFracRecip;
					auto deltaRX = oldVals[i].rightMagnitude * minFracRecip;
					// deltaX mostly zero here - add simd check
					auto lResult = deltaLX ? log(deltaLX) * deltaYRecip : 0;
					auto rResult = deltaRX ? log(deltaRX) * deltaYRecip : 0;
					output[i].leftMagnitude = lResult;
					output[i].rightMagnitude = rResult;
				}
				break;
			}
			case CSpectrum::ChannelConfiguration::Phase:
			{
				// make filter for phases slower since it's not in logdomain
				auto const phaseFilter = std::pow(filter.pole, 0.1);
				for (cpl::Types::fint_t i = 0; i < size; ++i)
				{

					auto lreal = newVals[i * 2];
					auto rreal = newVals[i * 2 + size * 2];
					auto limag = newVals[i * 2 + 1];
					auto rimag = newVals[i * 2 + size * 2 + 1];
					// mag = abs(cmplx)

					auto preal = lreal + rreal;
					auto pimag = limag + rimag;

					/*auto preal = rreal + limag;
					auto pimag = lreal + rimag;*/
					auto phase = sqrt(preal * preal + pimag * pimag);
					auto lmag = sqrt(lreal * lreal + limag * limag);
					auto rmag = sqrt(rreal * rreal + rimag * rimag);
					auto mag = lmag + rmag;
					// here we must avoid division by zero - either inject some noise, or compare

					if (mag != (decltype(mag))0)
						phase = 1 - phase / mag;
					else
						phase = 1;
					mag *= halfRecip;
					oldVals[i].magnitude *= filter.pole;
					if (mag > oldVals[i].magnitude)
					{
						oldVals[i].magnitude = mag;
					} 

					oldVals[i].phase = phase + phaseFilter * (oldVals[i].phase - phase);

					output[i].phase = oldVals[i].phase;
					// log10(y / _min) / log10(_max / _min);
					// deltaX mostly zero here - add simd check
					auto deltaX = oldVals[i].magnitude * minFracRecip;
					// deltaX mostly zero here - add simd check
					auto result = deltaX ? log(deltaX) * deltaYRecip : 0;
					output[i].magnitude = result;
				}
				break;
			}
			};


		}


	void CSpectrum::mapToLinearSpace()
	{
		auto const numPoints = getAxisPoints();
		auto const & dbRange = getDBs();
		std::size_t numFilters = getNumFilters();

		switch (state.algo)
		{
		case TransformAlgorithm::FFT:
		{
			int bin = 0, oldBin = 0, maxLBin, maxRBin = 0;
			std::size_t numBins = getNumAudioElements<std::complex<double>>() / 2;

			std::size_t N = numBins << 1;
			auto const topFrequency = getSampleRate() / 2;
			auto const freqToBin = double(numBins) / topFrequency;



			typedef double ftype;

			std::complex<ftype> leftMax, rightMax;

			ftype maxLMag, maxRMag, newLMag, newRMag;
			ftype * tsf = getAudioMemory<ftype>();
			std::complex<ftype> * csf = getAudioMemory<std::complex<ftype>>();
			ftype * wsp = getWorkingMemory<ftype>();


			// this will make scaling correct regardless of amount of zero-padding
			// notice the 0.5: fft's of size 32 will output 16 for exact frequency bin matches,
			// so we halve the reciprocal scaling factor to normalize the size.
			auto const invSize = 1.0 / (getWindowSize() * 0.5);
			
			// the DC (0) and nyquist bin are NOT 'halved' due to the symmetric nature of the fft,
			// so halve these:
			csf[0] *= 0.5;
			csf[N >> 2] *= 0.5;

			switch (state.configuration)
			{
			case ChannelConfiguration::Left:
			case ChannelConfiguration::Right:
			case ChannelConfiguration::Merge:
			{
				oldBin = mappedFrequencies[0] * freqToBin;

				

				
				// following is needed if 'bad' zero-padding is present
				// bad zeropadding (not in the middle) screws with the phase
				// so we replace the complex number with the absolute, so the interpolation
				// works. pretty bad though.

#pragma message cwarn("SIMD")
				for (int i = 0; i < numBins; ++i)
				{
					csf[i] = std::abs(csf[i]);
				}
				

				double fftBandwidth = 1.0 / numBins;
				//double pxlBandwidth = 1.0 / numPoints;
				cpl::Types::fint_t x;
				switch (state.binPolation)
				{
				case BinInterpolation::None:
					for (x = 0; x < numPoints - 1; ++x)
					{
						// bandwidth of the filter for this 'line', or point
						double bwForLine = (mappedFrequencies[x + 1] - mappedFrequencies[x]) / topFrequency;
						// as long as the bandwidth is smaller than our fft resolution, we interpolate the points
						// otherwise, break out and sample the max values of the bins inside the bandwidth
						if (bwForLine > fftBandwidth)
							break;
						// +0.5 to centerly space bins.
						auto index = cpl::Math::confineTo((std::size_t)(mappedFrequencies[x] * freqToBin + 0.5), 0, numBins - 1);
						auto interpolatedFilter = csf[index];

						wsp[x * 2] = interpolatedFilter.real() * invSize;
						wsp[x * 2 + 1] = interpolatedFilter.imag() * invSize;
					}
					break;
				case BinInterpolation::Linear:
					for (x = 0; x < numPoints - 1; ++x)
					{
						double bwForLine = (mappedFrequencies[x + 1] - mappedFrequencies[x]) / topFrequency;
						if (bwForLine > fftBandwidth)
							break;

						auto interpolatedFilter = cpl::linearFilter<std::complex<ftype>>(csf, N, mappedFrequencies[x] * freqToBin);

						wsp[x * 2] = interpolatedFilter.real() * invSize;
						wsp[x * 2 + 1] = interpolatedFilter.imag() * invSize;
					}
					break;
				case BinInterpolation::Lanczos:
					for (x = 0; x < numPoints - 1; ++x)
					{
						double bwForLine = (mappedFrequencies[x + 1] - mappedFrequencies[x]) / topFrequency;
						if (bwForLine > fftBandwidth)
							break;

						auto interpolatedFilter = cpl::lanczosFilter<std::complex<ftype>, true>(csf, N, mappedFrequencies[x] * freqToBin, 5);
						wsp[x * 2] = interpolatedFilter.real() * invSize;
						wsp[x * 2 + 1] = interpolatedFilter.imag() * invSize;
					}
					break;
				default:
					break;
				}

				oldBin = mappedFrequencies[x] * freqToBin;

				for (; x < numPoints; ++x)
				{
					maxLMag = maxRMag = newLMag = newRMag = 0;

					bin = static_cast<std::size_t>(mappedFrequencies[x] * freqToBin);
#ifdef DEBUG
					if (bin > getNumAudioElements < std::complex < ftype >> ())
						CPL_RUNTIME_EXCEPTION("Corrupt frequency mapping!");
#endif
					maxRBin = maxLBin = bin;

					signed diff = bin - oldBin;
					auto counter = diff ? 1 : 0;
					// here we loop over all the bins that is mapped for a single coordinate
					do
					{
						auto offset = oldBin + counter;
						offset <<= 1;
						auto real = tsf[offset];
						auto imag = tsf[offset + 1];
						newLMag = (real * real + imag * imag);
						// select highest number in this chunk for display. Not exactly correct, though.
						if (newLMag > maxLMag)
						{
							maxLBin = bin + counter;
							maxLMag = newLMag;
							leftMax = std::complex<ftype>(real, imag);
						}
						counter++;
						diff--;
					} while (diff > 0);


					wsp[x * 2] = leftMax.real() * invSize;
					wsp[x * 2 + 1] = leftMax.imag() * invSize;


					oldBin = bin;
				}


				mapAndTransformDFTFilters<ChannelConfiguration::Left, float>(filterStates.data(), wsp, filterResults.data(),
					filterStates.size(), dbRange.low, dbRange.high, peakFilter);
				break;
			}
			// for ChannelConfiguration::Separate and phase cases (where we work with 2 channels)
			// we need a workspace buffer sadly, since the spectrum of the channels
			// are contained quite messily in the output of the fft.
			// since in this setting the filterWorkspace buffer is not used, we will use that.
			case ChannelConfiguration::Phase:
			case ChannelConfiguration::Separate:
			{
				for (int x = 0; x < numPoints; ++x)
				{
					maxLMag = maxRMag = newLMag = newRMag = 0;

					bin = static_cast<std::size_t>(mappedFrequencies[x] * freqToBin);
					maxRBin = maxLBin = bin;

					auto diff = bin - oldBin;
					auto counter = diff ? 1 : 0;

					// here we loop over all the bins that is mapped for a single coordinate
					do
					{
						auto Z = getZFromNFFT<ftype>(tsf, bin + counter, N);
						newLMag = (Z.val[0].real() * Z.val[0].real() + Z.val[0].imag() * Z.val[0].imag());
						// select highest number in this chunk for display. Not exactly correct, though.
						if (newLMag > maxLMag)
						{
							maxLBin = bin + counter;
							newLMag = maxLMag;
							leftMax = Z.val[0];
						}
						newRMag = (Z.val[1].real() * Z.val[1].real() + Z.val[1].imag() * Z.val[1].imag());
						if (newRMag > maxRMag)
						{
							maxRBin = bin + counter;
							newRMag = maxRMag;
							rightMax = Z.val[1];
						}

						counter++;
						diff--;
					} while (diff > 0);

					wsp[x * 2] = leftMax.real() * invSize;
					wsp[x * 2 + 1] = leftMax.imag() * invSize;
					wsp[numFilters * 2 + x * 2] = rightMax.real() * invSize;
					wsp[numFilters * 2 + x * 2 + 1] = rightMax.imag() * invSize;

					oldBin = bin;
				}

				if (state.configuration == ChannelConfiguration::Separate)
					mapAndTransformDFTFilters<ChannelConfiguration::Separate, float>(filterStates.data(), wsp, filterResults.data(),
						filterStates.size(), dbRange.low, dbRange.high, peakFilter);
				else
					mapAndTransformDFTFilters<ChannelConfiguration::Phase, float>(filterStates.data(), wsp, filterResults.data(),
						filterStates.size(), dbRange.low, dbRange.high, peakFilter);

			}
			break;
			}
			break;
		}
		/*case Algorithm::MQDFT:
		{
			auto result = transformer.getTransformResult();
			auto totalData = filterStates.size() * numChannels;

			switch (channelConfiguration)
			{
			case ChannelConfiguration::Left:
			case ChannelConfiguration::Right:
			case ChannelConfiguration::Merge:
				mapAndTransformDFTFilters<ChannelConfiguration::Left, float>(filterStates.data(), result.data(), filterResults.data(),
					filterStates.size(), dbRange.low, dbRange.high, flt);
				break;
			case ChannelConfiguration::Separate:
				mapAndTransformDFTFilters<ChannelConfiguration::Separate, float>(filterStates.data(), result.data(), filterResults.data(),
					filterStates.size(), dbRange.low, dbRange.high, flt);
				break;
			case ChannelConfiguration::Phase:
				mapAndTransformDFTFilters<ChannelConfiguration::Phase, float>(filterStates.data(), result.data(), filterResults.data(),
					filterStates.size(), dbRange.low, dbRange.high, flt);
				break;
			}
			break;
		}*/
		case TransformAlgorithm::RSNT:
		{
			if (state.displayMode != DisplayMode::ColourSpectrum)
			{
				std::complex<float> * wsp = getWorkingMemory<std::complex<float>>();

				switch (state.dspWindow)
				{
				case cpl::dsp::WindowTypes::Hann:
					for (int i = 0; i < numFilters; i++)
					{
						wsp[i] = cresonator.getWindowedResonanceAt<cpl::dsp::WindowTypes::Hann, false>(i, 0);
					}
					break;
				case cpl::dsp::WindowTypes::Hamming:
					for (int i = 0; i < numFilters; i++)
					{
						wsp[i] = cresonator.getWindowedResonanceAt<cpl::dsp::WindowTypes::Hamming, false>(i, 0);
					}
					break;
				case cpl::dsp::WindowTypes::Rectangular:
					for (int i = 0; i < numFilters; i++)
					{
						wsp[i] = cresonator.getResonanceAt(i, 0);
					}
					break;
				}

				mapAndTransformDFTFilters<ChannelConfiguration::Left, float>(filterStates.data(), getWorkingMemory<float>(), filterResults.data(),
					filterStates.size(), dbRange.low, dbRange.high, peakFilter);
			}

		}
		break;
		}


	}


	bool CSpectrum::processNextSpectrumFrame()
	{
		SFrameBuffer::FrameVector * next;
		if (sfbuf.frameQueue.popElement(next))
		{
			SFrameBuffer::FrameVector & curFrame(*next);

			const auto numFilters = getNumFilters();
			if (curFrame.size() == numFilters)
			{
				mapAndTransformDFTFilters<ChannelConfiguration::Left, float>(filterStates.data(), reinterpret_cast<fpoint*>(curFrame.data()), filterResults.data(),
					filterStates.size(), state.dynRange.low, state.dynRange.high, peakFilter); 
			}
			else
			{
				// linearly interpolate bins. if we win the cpu-lottery one day, change this to sinc.
				std::complex<fpoint> * wsp = getWorkingMemory<std::complex<fpoint>>();

				// interpolation factor.
				fpoint wspToNext = (curFrame.size() - 1) / fpoint(std::max(1, numFilters));

				for (std::size_t n = 0; n < numFilters; ++n)
				{
					auto y2 = n * wspToNext;
					auto x = static_cast<std::size_t>(y2);
					auto yFrac = y2 - x;
					wsp[n] = curFrame[x] * (fpoint(1) - yFrac) + curFrame[x + 1] * yFrac;
				}

				mapAndTransformDFTFilters<ChannelConfiguration::Left, float>(filterStates.data(), getWorkingMemory<fpoint>(), filterResults.data(),
					filterStates.size(), state.dynRange.low, state.dynRange.high, peakFilter);
			}
#pragma message cwarn("OPERATOR DELETE OMG!!")
			delete next;
			return true;
		}
		return false;
	}

	void CSpectrum::mapFrequencies()
	{

		throw std::runtime_error("Dont do this");

	}
	
	bool CSpectrum::audioCallback(cpl::CAudioSource & source, float ** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		if (isSuspended || state.isFrozen || state.algo != TransformAlgorithm::RSNT)
			return false;

		// in aux mode, we post the audio data onto the fifo instead.
		if (state.iAuxMode)
		{
			std::size_t bufIndex = 0;

			auto numChannels = 1 + (state.configuration > ChannelConfiguration::Right);

			std::int64_t n = numSamples;

			switch (state.configuration)
			{
			case ChannelConfiguration::Right:
				bufIndex = 1;
			case ChannelConfiguration::Left:
				while (n > 0)
				{

					auto const aSamples = std::min(std::int64_t(AudioFrame::capacity), n - std::max(std::int64_t(0), n - std::int64_t(AudioFrame::capacity)));


					if (aSamples > 0)
					{
						AudioFrame af(aSamples);

						std::memcpy(af.buffer, (buffer[bufIndex] + numSamples - n), aSamples * AudioFrame::element_size);

						if (!audioFifo.pushElement(af))
							droppedAudioFrames++;
					}

					n -= aSamples;
				}
			case ChannelConfiguration::Merge:
			case ChannelConfiguration::Phase:
			case ChannelConfiguration::Separate:
				break;
			}





		}
		else
		{
			auto & cpuinfo = cpl::SysStats::CProcessorInfo::instance();

			if (cpuinfo.test(cpuinfo.AVX))
			{
				// size of magic here to support floats and doubles through changing the typename
				audioProcessing<typename cpl::simd::vector_of<fpoint, 8 / (sizeof(fpoint) / 4)>::type>(buffer, numChannels, numSamples);
			}
			else if (cpuinfo.test(cpuinfo.SSE2))
			{
				audioProcessing<typename cpl::simd::vector_of<fpoint, 4 / (sizeof(fpoint) / 4)>::type>(buffer, numChannels, numSamples);
			}
			else
			{
				audioProcessing<typename cpl::simd::vector_of<fpoint, 1 / (sizeof(fpoint) / 4)>::type>(buffer, numChannels, numSamples);
			}
		}


		return false;
	}

	template<typename V>
	void CSpectrum::addAudioFrame()
	{
		// locking, to ensure the amount of resonators doesn't change inbetween.
		cpl::CFastMutex lock(cresonator);

		const auto numResFilters = cresonator.getNumFilters();

		// yeah this needs to be fixed.
		auto & frame = *(new SFrameBuffer::FrameVector(numResFilters));


		switch (state.dspWindow)
		{
		case cpl::dsp::WindowTypes::Hann:
			for (int i = 0; i < numResFilters; i++)
			{
				frame[i] = cresonator.getWindowedResonanceAt<cpl::dsp::WindowTypes::Hann, false>(i, 0);
			}
			break;
		case cpl::dsp::WindowTypes::Hamming:
			for (int i = 0; i < numResFilters; i++)
			{
				frame[i] = cresonator.getWindowedResonanceAt<cpl::dsp::WindowTypes::Hamming, false>(i, 0);
			}
			break;
		case cpl::dsp::WindowTypes::Rectangular:
			for (int i = 0; i < numResFilters; i++)
			{
				frame[i] = cresonator.getResonanceAt(i, 0);
			}
			break;
		}
		sfbuf.frameQueue.pushElement<true>(&frame);
		sfbuf.currentCounter = 0;
	}

	template<typename V>
		void CSpectrum::audioProcessing(float ** buffer, std::size_t numChannels, std::size_t numSamples)
		{

			if (state.displayMode == DisplayMode::ColourSpectrum)
			{
				cpl::CFastMutex lock(sfbuf);

				std::int64_t n = numSamples;
				std::size_t offset = 0;
				while (n > 0)
				{
					std::int64_t numRemainingSamples = sfbuf.sampleBufferSize - sfbuf.currentCounter;

					const auto availableSamples = numRemainingSamples + std::min(std::int64_t(0), n - numRemainingSamples);

					switch (state.configuration)
					{
					case ChannelConfiguration::Left:
					{
						float * revBuf[2] = { buffer[0] + offset, buffer[1] + offset };
						cresonator.wresonate<V>(revBuf, 1, availableSamples); break;
					}

					case ChannelConfiguration::Right:
					{
						float * revBuf[2] = { buffer[1] + offset, buffer[0] + offset};
						cresonator.wresonate<V>(revBuf, 1, availableSamples); break;
					}
					case ChannelConfiguration::Merge:
					{
						// add relay buffer
						CPL_RUNTIME_EXCEPTION("Channelconfiguration not implemented yet");
					}
					default:
					{
						float * revBuf[2] = { buffer[0] + offset, buffer[1] + offset };
						cresonator.wresonate<V>(revBuf, numChannels, availableSamples);
					}
					}

					sfbuf.currentCounter += availableSamples;
					offset += availableSamples;
					if (sfbuf.currentCounter >= (sfbuf.sampleBufferSize))
					{
						addAudioFrame<V>();
						// change this here. oh really?
						sfbuf.sampleBufferSize = getBlobSamples();
					}

					n -= availableSamples;
				}


				sfbuf.sampleCounter += numSamples;
			}
			else
			{
				switch (state.configuration)
				{
				case ChannelConfiguration::Left:
					cresonator.wresonate<V>(buffer, 1, numSamples); break;
				case ChannelConfiguration::Right:
				{
					float * revBuf[2] = { buffer[1], buffer[0] };
					cresonator.wresonate<V>(revBuf, 1, numSamples); break;
				}
				case ChannelConfiguration::Merge:
				{
					// add relay buffer
					CPL_RUNTIME_EXCEPTION("Channelconfiguration not implemented yet");
				}
				default:
					cresonator.wresonate<V>(buffer, numChannels, numSamples);
				}
			}

			return;
		}

	float CSpectrum::getSampleRate() const noexcept
	{
		return static_cast<float>(audioStream[0].sampleRate);
	}

	int CSpectrum::getAxisPoints() const noexcept
	{
		return state.axisPoints;
	}

	double CSpectrum::getOptimalFramesPerUpdate() const noexcept
	{
#pragma message cwarn("collect this somewhere.")
		const double monitorRefreshRate = 60.0;
		auto res = double(isOpenGL() ? (monitorRefreshRate / getSwapInterval())  : refreshRate) / getBlobSamples();
		assert(std::isnormal(res));
		return res;
	}

	std::size_t CSpectrum::getApproximateStoredFrames() const noexcept
	{
#pragma message cwarn("fix this to include channels, other processing methods.. etc.")
		return sfbuf.frameQueue.enqueuededElements();
	}

	int CSpectrum::getNumFilters() const noexcept
	{
		return getAxisPoints();
	}
};
