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
	std::size_t CSpectrum::getFFTSpace() const noexcept
	{
		auto size = audioMemory.size();
		return size ? (size - 1) / sizeof(T) : 0;
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
	{/*
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
	 #ifdef CPL_LLVM_SUPPORTS_AVX
	 if (cpuinfo.test(cpuinfo.AVX))
	 {
	 // size of magic here to support floats and doubles through changing the typename
	 audioProcessing<typename cpl::simd::vector_of<fpoint, 8 / (sizeof(fpoint) / 4)>::type>(buffers, 2, numSamples);
	 }
	 else
	 #endif
	 if (cpuinfo.test(cpuinfo.SSE2))
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
	 */
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
				windowKernel[i] = 1.0 - std::cos(TAU * i / (sampleSize - 1));
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


		size_t fullSize = getFFTSpace<std::complex<double>>();
		// zero-padding
		for (; i < fullSize; ++i)
		{
			windowKernel[i] = 0.0;
		}
		// uncomment if you want to use window for convolution
		//cpl::dsp::CSignalTransform::sfft(reinterpret_cast<double*>(windowKernel.data()), fullSize);
	}



	void CSpectrum::prepareTransform(const AudioStream::AudioBufferAccess & audio)
	{
		// we can 'freeze' the window pretty effectively, by not allowing new input
		// and still be able to zoom around. 
		if (state.isFrozen)
			return;

		auto size = getWindowSize(); // the size of the transform, containing samples
									 // the quantized (to next power of 2) samples of this transform
									 // that is, the size + additional zero-padding
		auto fullSize = getFFTSpace<std::complex<double>>();

		auto const channelConfiguration = kchannelConfiguration.getZeroBasedSelIndex<ChannelConfiguration>();

		//auto numChannels = 1 + (channelConfiguration > ChannelConfiguration::Merge);

		{
			Stream::AudioBufferView views[2] = { audio.getView(0), audio.getView(1) };


			typedef double fftType;

			// note here, that since our audiobuffer is a circular buffer,
			// index[0] is actually the last input, and index[size - 1] is the first
			// while the fft doesn't care whether the input is reversed or not, 
			// the varying Q-transforms do care, since they vary the amount of samples
			// considered.

			switch (state.algo)
			{
			case TransformAlgorithm::FFT:
			{
				auto buffer = getAudioMemory<std::complex<fftType>>();
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

					for (std::size_t indice = 0; indice < Stream::bufferIndices; ++indice)
					{
						std::size_t range = views[channel].getItRange(indice);
						auto it = views[channel].getItIndex(indice);

						while (range--)
						{
							buffer[i] = *it++ * windowKernel[i];
							i++;
						}
					}
					break;
#endif
				}
				case ChannelConfiguration::Merge:
				{
#ifdef _ZERO_PHASE
					for (; i < halfInputSize; ++i)
					{
						buffer[i + offset] = (audioStream[0].singleCheckAccess(i) + audioStream[1].singleCheckAccess(i)) * windowKernel[i] * 0.5;
					}

					offset = halfInputSize;
					for (; i < size; ++i)
					{
						buffer[i - offset] = (audioStream[0].singleCheckAccess(i) + audioStream[1].singleCheckAccess(i)) * windowKernel[i] * 0.5;;
					}
#else
					for (std::size_t indice = 0; indice < Stream::bufferIndices; ++indice)
					{
						std::size_t range = views[0].getItRange(indice);
						auto left = views[0].getItIndex(indice);
						auto right = views[1].getItIndex(indice);

						while (range--)
						{
							buffer[i] = (*left++ + *right++) * windowKernel[i] * 0.5f;
							i++;
						}
					}
					break;
#endif
				}
				case ChannelConfiguration::Side:
				{

					for (std::size_t indice = 0; indice < Stream::bufferIndices; ++indice)
					{
						std::size_t range = views[0].getItRange(indice);
						auto left = views[0].getItIndex(indice);
						auto right = views[1].getItIndex(indice);

						while (range--)
						{
							buffer[i] = (*left++ - *right++) * windowKernel[i] * 0.5f;
							i++;
						}
					}
					break;
				}
				case ChannelConfiguration::MidSide:
				{
					for (std::size_t indice = 0; indice < Stream::bufferIndices; ++indice)
					{
						std::size_t range = views[0].getItRange(indice);
						auto left = views[0].getItIndex(indice);
						auto right = views[1].getItIndex(indice);

						while (range--)
						{
							buffer[i] = std::complex<fftType>
								(
									(*left - *right) * windowKernel[i] * 0.5f,
									(*left + *right) * windowKernel[i] * 0.5f
									);
							left++;
							right++;
							i++;
						}
					}
					break;
				}
				case ChannelConfiguration::Phase:
				case ChannelConfiguration::Separate:
				{
#ifdef _ZERO_PHASE
					for (; i < halfInputSize; ++i)
					{
						buffer[i + offset] = std::complex<fftType>
							(
								audioStream[0].singleCheckAccess(i)  * windowKernel[i],
								audioStream[1].singleCheckAccess(i)  * windowKernel[i]
								);
					}

					offset = halfInputSize;
					for (; i < size; ++i)
					{
						buffer[i - offset] = std::complex<fftType>
							(
								audioStream[0].singleCheckAccess(i) * windowKernel[i],
								audioStream[1].singleCheckAccess(i) * windowKernel[i]
								);
					}
					break;
#else
					for (std::size_t indice = 0; indice < Stream::bufferIndices; ++indice)
					{
						std::size_t range = views[0].getItRange(indice);
						auto left = views[0].getItIndex(indice);
						auto right = views[1].getItIndex(indice);

						while (range--)
						{
							buffer[i] = std::complex<fftType>
							{
								*left++ * windowKernel[i],
								*right++ * windowKernel[i]
							};
							i++;
						}
					}
					break;
#endif
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
			auto const numSamples = getFFTSpace<std::complex<double>>();
			if (numSamples != 0)
				cpl::signaldust::DustFFT_fwdDa(getAudioMemory<double>(), numSamples);
			/*switch (channelConfiguration)
			{
			case ChannelConfiguration::Left:
			case ChannelConfiguration::Right:
			case ChannelConfiguration::Merge:
			cpl::signaldust::DustFFT_fwdDa(getAudioMemory<double>(), numSamples);
			break;
			case ChannelConfiguration::Phase:
			case ChannelConfiguration::Separate:
			cpl::signaldust::DustFFT_fwdDa(getAudioMemory<double>(), numSamples);

			}*/

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
	template<ChannelConfiguration type, class Scalar, class V2>
	void mapAndTransformDFTFilters(UComplexFilter<Scalar> * CPL_RESTRICT oldVals,
		const V2 & newVals, UComplexFilter<Scalar> * CPL_RESTRICT output, std::size_t size,
		float lowDbs, float highDbs, cpl::CPeakFilter<Scalar> filter)
	{

		double lowerFraction = cpl::Math::dbToFraction<double>(lowDbs);
		double upperFraction = cpl::Math::dbToFraction<double>(highDbs);
		auto deltaYRecip = static_cast<Scalar>(1.0 / log(upperFraction / lowerFraction));
		auto minFracRecip = static_cast<Scalar>(1.0 / lowerFraction);
		auto halfRecip = Scalar(0.5);

		Scalar lowerClip = (Scalar)lowerFraction;

		switch (type)
		{
		case ChannelConfiguration::Left:
		case ChannelConfiguration::Merge:
		case ChannelConfiguration::Right:
		case ChannelConfiguration::Side:
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
				auto result = deltaX > 0 ? log(deltaX) * deltaYRecip : lowerClip;
				output[i].magnitude = result;
				output[i].phase = 0;
			}
			break;
		}
		case ChannelConfiguration::Separate:
		case ChannelConfiguration::MidSide:
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
				auto lResult = deltaLX > 0 ? log(deltaLX) * deltaYRecip : lowerClip;
				auto rResult = deltaRX > 0 ? log(deltaRX) * deltaYRecip : lowerClip;
				output[i].leftMagnitude = lResult;
				output[i].rightMagnitude = rResult;
			}
			break;
		}
		case ChannelConfiguration::Phase:
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
				phase *= mag;

				oldVals[i].phase = phase + phaseFilter * (oldVals[i].phase - phase);


				// log10(y / _min) / log10(_max / _min);
				// deltaX mostly zero here - add simd check
				auto deltaX = oldVals[i].magnitude * minFracRecip;
				auto deltaY = oldVals[i].phase * minFracRecip;
				// deltaX mostly zero here - add simd check
				auto result = deltaX > 0 ? log(deltaX) * deltaYRecip : lowerClip;
				output[i].magnitude = result;
				output[i].phase = (deltaY > 0 ? log(deltaY) * deltaYRecip : lowerClip);
			}
			break;
		}
		};


	}


	void CSpectrum::mapToLinearSpace()
	{
		using namespace cpl;
		auto const numPoints = getAxisPoints();
		auto const & dbRange = getDBs();
		std::size_t numFilters = getNumFilters();

		switch (state.algo)
		{
		case TransformAlgorithm::FFT:
		{
			int bin = 0, oldBin = 0, maxLBin, maxRBin = 0;
			std::size_t N = getFFTSpace<std::complex<double>>();

			// we rely on mapping indexes, so we need N > 2 at least.
			if (N == 0)
				return;

			std::size_t numBins = N >> 1;
			std::size_t halfSize = N >> 1;
			auto const topFrequency = getSampleRate() / 2;
			auto const freqToBin = double(numBins) / topFrequency;
			bool hasZeroPadded = !juce::isPowerOfTwo(getWindowSize());


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

			switch (state.configuration)
			{
			case ChannelConfiguration::Left:
			case ChannelConfiguration::Right:
			case ChannelConfiguration::Merge:
			case ChannelConfiguration::Side:
			{
				oldBin = mappedFrequencies[0] * freqToBin;

				// the DC (0) and nyquist bin are NOT 'halved' due to the symmetric nature of the fft,
				// so halve these:
				csf[0] *= 0.5;
				csf[N >> 1] *= 0.5;


				// following is needed if 'bad' zero-padding is present
				// bad zeropadding (not in the middle) screws with the phase
				// so we replace the complex number with the absolute, so the interpolation
				// works. pretty bad though.

#pragma message cwarn("SIMD")
				//if (hasZeroPadded)
				{

					/*for (std::size_t i = 0; i < numBins; ++i)
					{
					csf[i] = std::abs(csf[i]);
					}*/
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
						auto index = Math::confineTo((std::size_t)(mappedFrequencies[x] * freqToBin + 0.5), 0, numBins - 1);
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

						auto interpolatedFilter = dsp::linearFilter<std::complex<ftype>>(csf, N, mappedFrequencies[x] * freqToBin);

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

						auto interpolatedFilter = dsp::lanczosFilter<std::complex<ftype>, true>(csf, N, mappedFrequencies[x] * freqToBin, 5);
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
					filterStates.size(), dbRange.low, dbRange.high, kMinDBRange, peakFilter);
				break;
			}
			case ChannelConfiguration::Phase:
			case ChannelConfiguration::Separate:
			case ChannelConfiguration::MidSide:
			{
				// two-for-one pass, first channel is 0... N/2 -1, second is N/2 .. N -1
				dsp::separateTransformsIPL(csf, N);

				// fix up DC and nyquist bins (see previous function documentation)
				csf[N] = csf[0].imag() * 0.5;
				csf[0] = csf[0].real() * 0.5;
				csf[N >> 1] *= 0.5;
				csf[(N >> 1) - 1] *= 0.5;

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
						auto index = Math::confineTo((std::size_t)(mappedFrequencies[x] * freqToBin + 0.5), 0, numBins - 1);
						auto interpolatedLeft = csf[index];
						auto interpolatedRight = csf[N - index];
						wsp[x * 2] = interpolatedLeft.real() * invSize;
						wsp[x * 2 + 1] = interpolatedLeft.imag() * invSize;
						wsp[numFilters * 2 + x * 2] = interpolatedRight.real() * invSize;
						wsp[numFilters * 2 + x * 2 + 1] = interpolatedRight.imag() * invSize;
					}
					break;
				case BinInterpolation::Linear:

					// interpolations on phase-mangled buffers are quirky, this just makes everything easier...
					/*for (std::size_t i = 1; i < N; ++i)
					{
					csf[i] = std::abs(csf[i]);
					}*/

					for (x = 0; x < numPoints - 1; ++x)
					{
						double bwForLine = (mappedFrequencies[x + 1] - mappedFrequencies[x]) / topFrequency;
						if (bwForLine > fftBandwidth)
							break;

						auto interpolatedLeft = dsp::linearFilter<std::complex<ftype>>(csf, N + 1, mappedFrequencies[x] * freqToBin);
						auto interpolatedRight = dsp::linearFilter<std::complex<ftype>>(csf, N + 1, N - (mappedFrequencies[x] * freqToBin));

						wsp[x * 2] = interpolatedLeft.real() * invSize;
						wsp[x * 2 + 1] = interpolatedLeft.imag() * invSize;
						wsp[numFilters * 2 + x * 2] = interpolatedRight.real() * invSize;
						wsp[numFilters * 2 + x * 2 + 1] = interpolatedRight.imag() * invSize;
					}
					break;
				case BinInterpolation::Lanczos:

					// interpolations on phase-mangled buffers are quirky, this just makes everything easier...
					/*for (std::size_t i = 1; i < N; ++i)
					{
					csf[i] = std::abs(csf[i]);
					}*/

					for (x = 0; x < numPoints - 1; ++x)
					{

						double bwForLine = (mappedFrequencies[x + 1] - mappedFrequencies[x]) / topFrequency;
						if (bwForLine > fftBandwidth)
							break;

						auto interpolatedLeft = dsp::lanczosFilter<std::complex<ftype>, true>(csf, N + 1, mappedFrequencies[x] * freqToBin, 5);
						auto interpolatedRight = dsp::lanczosFilter<std::complex<ftype>, true>(csf, N + 1, N - (mappedFrequencies[x] * freqToBin), 5);
						wsp[x * 2] = interpolatedLeft.real() * invSize;
						wsp[x * 2 + 1] = interpolatedLeft.imag() * invSize;
						wsp[numFilters * 2 + x * 2] = interpolatedRight.real() * invSize;
						wsp[numFilters * 2 + x * 2 + 1] = interpolatedRight.imag() * invSize;
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
						auto lreal = tsf[offset];
						auto limag = tsf[offset + 1];
						newLMag = (lreal * lreal + limag * limag);

						auto rreal = tsf[N * 2 - offset];
						auto rimag = tsf[N * 2 - offset + 1];
						newRMag = (rreal * rreal + rimag * rimag);
						// select highest number in this chunk for display. Not exactly correct, though.
						if (newLMag > maxLMag)
						{
							maxLBin = bin + counter;
							maxLMag = newLMag;
							leftMax = std::complex<ftype>(lreal, limag);
						}

						if (newRMag > maxRMag)
						{
							maxRBin = bin + counter;
							maxRMag = newRMag;
							rightMax = std::complex<ftype>(rreal, rimag);
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

				// --------------------------------------------------------------------------------------------------------------


				if (state.configuration == ChannelConfiguration::Phase)
					mapAndTransformDFTFilters<ChannelConfiguration::Phase, float>(filterStates.data(), wsp, filterResults.data(),
						filterStates.size(), dbRange.low, dbRange.high, kMinDBRange, peakFilter);
				else
					// handles mid/side also.
					mapAndTransformDFTFilters<ChannelConfiguration::Separate, float>(filterStates.data(), wsp, filterResults.data(),
						filterStates.size(), dbRange.low, dbRange.high, kMinDBRange, peakFilter);

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
					filterStates.size(), dbRange.low, dbRange.high, kMinDBRange, peakFilter);
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

			// the size will be zero for a couple of frames, if there's some messing around with window sizes
			// or we get audio running before anything is actually initiated.
			if (curFrame.size() != 0)
			{
				if (curFrame.size() == numFilters)
				{
					mapAndTransformDFTFilters<ChannelConfiguration::Left, float>(filterStates.data(), reinterpret_cast<fpoint*>(curFrame.data()), filterResults.data(),
						filterStates.size(), state.dynRange.low, state.dynRange.high, kMinDBRange, peakFilter);
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
						filterStates.size(), state.dynRange.low, state.dynRange.high, kMinDBRange, peakFilter);
				}
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

	bool CSpectrum::onAsyncAudio(const AudioStream & source, AudioStream::DataType ** buffer, std::size_t numChannels, std::size_t numSamples)
	{
		if (isSuspended)
			return false;
		auto const factor = std::is_same<AudioStream::DataType, float>::value ? 1 : 2;
		switch (cpl::simd::max_vector_capacity<float>())
		{
		case 8:
			audioProcessing<typename cpl::simd::vector_of<fpoint, 8 * 4 / (sizeof(fpoint))>::type>(buffer, numChannels, numSamples);
			break;
		case 4:
			audioProcessing<typename cpl::simd::vector_of<fpoint, 4 * 4 / (sizeof(fpoint))>::type>(buffer, numChannels, numSamples);
			break;
		case 1:
			audioProcessing<typename cpl::simd::vector_of<fpoint, factor * 4 / (sizeof(fpoint))>::type>(buffer, numChannels, numSamples);
			break;
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
		if (state.algo != TransformAlgorithm::RSNT)
			return;

		// rest only for resonators.
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
					float * revBuf[2] = { buffer[1] + offset, buffer[0] + offset };
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
		return static_cast<float>(audioStream.getInfo().sampleRate);
	}

	int CSpectrum::getAxisPoints() const noexcept
	{
		return state.axisPoints;
	}

	double CSpectrum::getOptimalFramesPerUpdate() const noexcept
	{
#pragma message cwarn("collect this somewhere.")
		const double monitorRefreshRate = 60.0;
		auto res = double(isOpenGL() ? (monitorRefreshRate / getSwapInterval()) : refreshRate) / getBlobSamples();
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
#pragma once
