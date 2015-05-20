#include "CSpectrum.h"

#include <cpl/Mathext.h>

#include <cpl/dsp.h>
//#include <cpl/dsp/CPeakFilter.h>
//#include <cpl/resampling.h>

#include <cpl/simd.h>
#include <cpl/SysStats.h>

namespace Signalizer
{
	void ColorScale(uint8_t * pixel, float intensity)
	{
		intensity *= 3;
		uint8_t red = 0, blue = 0, green = 0;
		// set blue
		if (intensity < 0.333f)
		{
			red = 3 * intensity * 0xFF;
			
		}
		else if (intensity < 0.6666f)
		{
			red = 0xFF;
			green = 3 * (intensity - 0.3333) * 0xFF;
		}
		else if (intensity < 1)
		{
			red = 0xFF;
			green = 0xFF;
			blue = 3 * (intensity - 0.66666) * 0xFF;
		}
		// set green
		pixel[0] = blue;
		pixel[1] = green;
		pixel[2] = red;
		// set red

		// saturate


	}


	template<class ToType, class FromType, class Vector>
	void copyTypeIntoBuffer(Vector & buf, FromType ** vec, std::size_t numChannels, std::size_t numSamples)
	{
		buf.resize(numChannels);
		for (cpl::Types::fint_t c = 0; c < numChannels; ++c)
		{
			buf[c].resize(numSamples);
			for (cpl::Types::fint_t i = 0; i < numSamples; ++i)
			{
				buf[c][i] = static_cast<ToType>(vec[c][i]);
			}
		}
	}


	bool CSpectrum::audioCallback(float ** data, std::size_t numChannels, std::size_t numSamples)
	{
		if (algorithmType == Algorithm::RSNT)
		{
			copyTypeIntoBuffer<double>(relayBuffer, data, numChannels, numSamples);
			auto const & cpuid = cpl::SysStats::CProcessorInfo::instance();
			if (cpuid.test(cpuid.AVX))
				sdft.wresonate<cpl::Types::v4sd>(relayBuffer, numChannels, numSamples);
			else if (cpuid.test(cpuid.SSE2))
				sdft.wresonate<cpl::Types::v2sd>(relayBuffer, numChannels, numSamples);
			else
				sdft.wresonate<double>(relayBuffer, numChannels, numSamples);
		}
		return false;
	}

	float CSpectrum::scale(float in, float min, float max)
	{
		float expScale = cpl::Math::UnityScale::exp<float>(in, 10, audioData[0].sampleRate / 2);
		float unityPoint = (expScale - 10) / (audioData[0].sampleRate / 2 - 10);
		float scaledPoint = min + unityPoint * (max - min);
		return expViewScale * scaledPoint + (1 - expViewScale) * cpl::Math::UnityScale::linear(in, min, max);
	}
	float CSpectrum::invScale(float in, float min, float max)
	{
		return expViewScale * cpl::Math::UnityScale::Inv::exp(in, min, max) + (1 - expViewScale) * cpl::Math::UnityScale::Inv::linear(in, min, max);
	}

	template<typename V>
		void CSpectrum::resonate(float ** data, std::size_t numChannels, std::size_t numSamples)
		{
			using namespace cpl;
			using namespace cpl::simd;

			auto const vfactor = suitable_container<V>::size;
			auto const numFilters = filterResults.size();
			ResonatorSystem<scalar_of<V>::type> & r = resonator;
			if (!r.buffer.size())
				return;
			cpl::CMutex lock(r);

			for (Types::fint_t filter = 0; filter < numFilters; filter += vfactor)
			{
				V pole = load<V>(r.poles + filter);
				V c1 = load<V>(r.c1 + filter);
				V c2 = load<V>(r.c2 + filter);
				V x = load<V>(r.x + filter);
				V y = load<V>(r.y + filter);
				
				V real = load<V>(r.real + filter);
				V imag = load<V>(r.imag + filter);
				V t0;
				V input;

				for (Types::fint_t sample = 0; sample < numSamples; ++sample)
				{
					input = broadcast<V>(data[0] + sample);

					t0 = x * input;
					real = t0 + pole * (real - t0);

					t0 = y * input;
					imag = t0 + pole * (imag - t0);


					t0 = x * c1 - y * c2;
					y = x * c2 + y * c1;
					x = t0;
				}

				store(r.x + filter, x);
				store(r.y + filter, y);
				store(r.real + filter, real);
				store(r.imag + filter, imag);
			}


		}
		template<typename T>
			struct Resonator
			{
				typedef T type;

				std::size_t N;
				T pole[3][2];
				T state[3][2];
			};



		template<typename V>
			void CSpectrum::wresonate(float ** data, std::size_t numChannels, std::size_t numSamples)
			{
				using namespace cpl;
				using namespace cpl::simd;

				auto const vfactor = suitable_container<V>::size;
				auto const numFilters = filterResults.size();

				cpl::CMutex lock(sdft);
				

				auto size = audioData[0].size;
				V t0;

				using namespace cpl;
				using namespace cpl::simd;

				std::ptrdiff_t s2 = -signed(size);
				typename scalar_of<V>::type * ptrs[elements_of<V>::value];

				const auto upperBounds = audioData[0].buffer + size;

				for (Types::fint_t filter = 0; filter < numFilters; filter += vfactor)
				{
					Types::fsint_t ptr = audioData[0].start;

					V
						p_m1_r = load<V>(sdft.realPoles[0] + filter), // pole: e^i*omega-q (real)
						p_m1_i = load<V>(sdft.imagPoles[0] + filter), // pole: e^i*omega-q (imag)
						p_m_r = load<V>(sdft.realPoles[1] + filter), // pole: e^i*omega (real)
						p_m_i = load<V>(sdft.imagPoles[1] + filter), // pole: e^i*omega (imag)
						p_p1_r = load<V>(sdft.realPoles[2] + filter), // pole: e^i*omega+q (real)
						p_p1_i = load<V>(sdft.imagPoles[2] + filter); // pole: e^i*omega+q (imag)

					V
						s_m1_r = load<V>(sdft.realState[0] + filter), // state: e^i*omega-q (real)
						s_m1_i = load<V>(sdft.imagState[0] + filter), // state: e^i*omega-q (imag)
						s_m_r = load<V>(sdft.realState[1] + filter), // state: e^i*omega (real)
						s_m_i = load<V>(sdft.imagState[1] + filter), // state: e^i*omega (imag)
						s_p1_r = load<V>(sdft.realState[2] + filter), // state: e^i*omega+q (real)
						s_p1_i = load<V>(sdft.imagState[2] + filter); // state: e^i*omega+q (imag)

					for (Types::fint_t z = 0; z < vfactor; ++z)
					{
						auto offset = ptr - (signed int)sdft.N[filter + z];
						offset += (offset < 0) * size;
						offset -= (offset >= size) * size;
						ptrs[z] = audioData[0].buffer + offset;

					}



					for (Types::fint_t sample = 0; sample < numSamples; ++sample)
					{

						for (Types::fint_t z = 0; z < vfactor; ++z)
						{

							ptrs[z]++;
							ptrs[z] += (ptrs[z] >= upperBounds) ? s2 : (std::ptrdiff_t)0; // lengths[z] %= size; (wrap around)
						}

						V input = gather<V>(ptrs) + broadcast<V>(data[0] + sample);

						t0 = s_m1_r * p_m1_r - s_m1_i * p_m1_i + input;
						s_m1_i = s_m1_r * p_m1_i + s_m1_i * p_m1_r + input;
						s_m1_r = t0;

						t0 = s_m_r * p_m_r - s_m_i * p_m_i + input;
						s_m_i = s_m_r * p_m_i + s_m_i * p_m_r + input;
						s_m_r = t0;

						t0 = s_p1_r * p_p1_r - s_p1_i * p_p1_i + input;
						s_p1_i = s_p1_r * p_p1_i + s_p1_i * p_p1_r + input;
						s_p1_r = t0;

						//ptr++;
					}


					store(sdft.realPoles[0] + filter, p_m1_r); // pole: e^i*omega-q (real)
					store(sdft.imagPoles[0] + filter, p_m1_i); // pole: e^i*omega-q (imag)
					store(sdft.realPoles[1] + filter, p_m_r); // pole: e^i*omega (real)
					store(sdft.imagPoles[1] + filter, p_m_i); // pole: e^i*omega (imag)
					store(sdft.realPoles[2] + filter, p_p1_r); // pole: e^i*omega+q (real)
					store(sdft.imagPoles[2] + filter, p_p1_i); // pole: e^i*omega+q (imag)

					store(sdft.realState[0] + filter, s_m1_r); // state: e^i*omega-q (real)
					store(sdft.imagState[0] + filter, s_m1_i); // state: e^i*omega-q (imag)
					store(sdft.realState[1] + filter, s_m_r); // state: e^i*omega (real)
					store(sdft.imagState[1] + filter, s_m_i); // state: e^i*omega (imag)
					store(sdft.realState[2] + filter, s_p1_r); // state: e^i*omega+q (real)
					store(sdft.imagState[2] + filter, s_p1_i); // state: e^i*omega+q (imag)

				}


				}
			
	void CSpectrum::mapFrequencies()
	{

		auto const numPoints = numFilters;
		// no filters, dont bother
		if (!numPoints)
			return;
		auto const recipPoints = 1.0 / numPoints;
		auto const low = double(lowestCoord < 0 ? 0 : lowestCoord > numPoints ? numPoints : lowestCoord) / numPoints;
		auto const top = double(highestCoord < 0 ? 0 : highestCoord > numPoints ? numPoints : highestCoord) / numPoints;
		
		for(int i = 0; i < numPoints; ++i)
		{
			// the 'world' coordinate for this current x-coordinate
			auto floatCoord = cpl::Math::UnityScale::linear<double>(i * recipPoints, low, top);
			mappedFrequencies[i] = scale(floatCoord, 10, audioData[0].sampleRate / 2);
			
		}
		transformer.setKernelData(mappedFrequencies, numPoints);
		// calculate resonator coefficients
		std::size_t windowSize = getWindowSize();
		std::size_t minWindow = 8;
		sdft.setWindowSize(minWindow, windowSize);

		sdft.mapSystemHz(mappedFrequencies, mappedFrequencies.size(), audioData[0].sampleRate);
	}

	void CSpectrum::mapResonatingSystem()
	{
		// 7 == state size of all members
		using namespace cpl;
		CMutex lock(sdft);

		/*auto numResonators = numFilters + (8 - numFilters & 0x7); // quantize to next multiple of 8, to ensure vectorization
		auto const dataSize = numResonators * 4 * 3 + numResonators;
		auto const sampleRate = audioData[0].sampleRate;

		sdft.buffer.resize(dataSize);

		for (std::size_t i = 0; i < 3; ++i)
		{
			sdft.realPoles[i] = sdft.buffer.data() + numResonators * i * 4;
			sdft.imagPoles[i] = sdft.buffer.data() + numResonators * i * 4 + numResonators;
			sdft.realState[i] = sdft.buffer.data() + numResonators * i * 4 + numResonators * 2;
			sdft.imagState[i] = sdft.buffer.data() + numResonators * i * 4 + numResonators * 3;
		}



		sdft.N = reinterpret_cast<unsigned*>(sdft.buffer.data() + numResonators * 2 * 4 + numResonators * 4); // ewww

		std::size_t windowSize = getWindowSize();
		std::size_t minWindow = 8;

		for (int i = 0; i < filterResults.size(); ++i)
		{

			auto const omega = 2 * mappedFrequencies[i] * M_PI / sampleRate;


			auto const k = i + 1 >= numFilters ? numFilters - 2 : i;
			auto bandWidth = sampleRate / ((double)mappedFrequencies[k + 1] - mappedFrequencies[k]);
			bandWidth = (qSettings == Quality::Bound) ? cpl::Math::confineTo<double>(bandWidth, minWindow, windowSize) : bandWidth;

			auto const Q = 2 * M_PI / bandWidth;

			sdft.realPoles[0][i] = std::cos(omega - Q);
			sdft.realPoles[1][i] = std::cos(omega);
			sdft.realPoles[2][i] = std::cos(omega + Q);
			sdft.imagPoles[1][i] = std::sin(omega - Q);
			sdft.imagPoles[2][i] = std::sin(omega);
			sdft.imagPoles[3][i] = std::sin(omega + Q);
			sdft.N[i] = cpl::Math::round<unsigned>(bandWidth);
		}
		// all sine phases and accumulators are zero initially
		//std::memset(resonator.y, 0, numFilters * sizeof(float) * 3);
		//setBandwidth();*/
	}


	/*void CSpectrum::mapResonatingSystem()
	{
		// 7 == state size of all members
		using namespace cpl;
		CMutex lock(resonator);

		auto numResonators = numFilters + (8 - numFilters & 0x7); // quantize to next multiple of 8, to ensure vectorization
		auto const dataSize = numResonators * 7;
		auto const sampleRate = audioData[0].sampleRate;

		resonator.buffer.resize(dataSize);
		resonator.poles = resonator.buffer.data();
		resonator.c1 = resonator.poles + numResonators;
		resonator.c2 = resonator.c1 + numResonators;
		resonator.x = resonator.c2 + numResonators;
		resonator.y = resonator.x + numResonators;
		resonator.real = resonator.y + numResonators;
		resonator.imag = resonator.real + numResonators;

		for (int i = 0; i < filterResults.size(); ++i)
		{

			auto const omega = std::tan(mappedFrequencies[i] * M_PI / sampleRate);
			auto const z = 2.0 / (1.0 + omega * omega);

			resonator.c1[i] = static_cast<float>(z - 1.0);
			resonator.c2[i] = static_cast<float>(omega * z);
			resonator.x[i] = 1;
		}
		// all sine phases and accumulators are zero initially
		std::memset(resonator.y, 0, numFilters * sizeof(float)* 3);
		setBandwidth();
	}*/

	void CSpectrum::setWindowSize(std::size_t smps)
	{

		std::size_t n = std::min(audioData[0].maxSize(), smps);
		n -= (n & 0x7); // must be a multiple of 7, due to vectorization
		if (n < 16)
			n = 16;
		auto const bufSize = cpl::Math::nextPow2(n);
		for (auto & buf : audioData)
		{
			buf.setSize(n);
		}
		memory.resize(bufSize * sizeof(std::complex<double>));
		windowKernel.resize(bufSize);
		computeWindowKernel();
		std::memset(memory.data(), 0, memory.size());
		setBandwidth();
	}

	void CSpectrum::setBandwidth()
	{
		auto const sampleRate = audioData[0].sampleRate;
		auto numFilters = filterResults.size();
		std::size_t windowSize = getWindowSize();
		std::size_t minWindow = 8;
		sdft.setWindowSize(minWindow,windowSize);

		for (int i = 0; i < filterResults.size(); ++i)
		{
			auto const k = i + 1 >= numFilters ? numFilters - 2 : i;
			auto bandWidth = sampleRate / (mappedFrequencies[k + 1] - mappedFrequencies[k]);
			bandWidth = (qSettings == Quality::Bound) ? cpl::Math::confineTo<double>(bandWidth, minWindow, windowSize) : bandWidth;
			//resonator.poles[i] = cpl::Math::expDecay(bandWidth);
		}

		sdft.mapSystemHz(mappedFrequencies, mappedFrequencies.size(), audioData[0].sampleRate);
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
	template<CSpectrum::ChannelConfiguration::Enum type, class Scalar, class V2>
		void mapAndTransformDFTFilters(FTFilterResult<Scalar> * oldVals, 
			const V2 & newVals, FTFilterResult<Scalar> * output, std::size_t size,
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
		auto const channelConfiguration = cpl::distribute<ChannelConfiguration>(kchannelConf->bGetValue());
		auto const numPoints = getWidth();
		auto const & dbRange = getDBs();
		std::size_t numFilters = filterStates.size();
		switch(algorithmType)
		{
			case Algorithm::FFT:
			{
				int bin = 0, oldBin = 0, maxLBin, maxRBin = 0;
				std::size_t numBins = getNumAudioElements<std::complex<double>>() / 2;

				std::size_t N = numBins << 1;
				auto const topFrequency = audioData[0].sampleRate / 2;
				auto const freqToBin = double(numBins) / topFrequency;
				auto const logLoft = std::log10(topFrequency);
				auto const invSize = 1.0 / (getWindowSize()); // this will make scaling correct regardless of amount of zero-padding

				typedef double ftype;

				std::complex<ftype> leftMax, rightMax;

				ftype maxLMag, maxRMag, newLMag, newRMag;
				ftype * tsf = getAudioMemory<ftype>();
				std::complex<ftype> * csf = getAudioMemory<std::complex<ftype>>();
				ftype * wsp = getWorkingMemory<ftype>();
				switch (channelConfiguration)
				{
				case ChannelConfiguration::Left:
				case ChannelConfiguration::Right:
				case ChannelConfiguration::Merge:
				{
					oldBin = mappedFrequencies[0] * freqToBin;
					/*
						// following is needed if 'bad' zero-padding is present
						// bad zeropadding (not in the middle) screws with the phase
						// so we replace the complex number with the absolute, so the interpolation
						// works. pretty bad though.

						auto vp = getAudioMemory<std::complex<ftype>>();
						auto vn = getNumAudioElements<std::complex<ftype>>();
						for (int i = 0; i < vn; ++i)
						{
							vp[i] = std::abs(vp[i]);
						}
					*/

					double fftBandwidth = 1.0 / numBins;
					double pxlBandwidth = 1.0 / numPoints;
					cpl::Types::fint_t x;

					for (x = 0; x < numPoints - 1; ++x)
					{
						// bandwidth of the filter for this 'line', or point
						double bwForLine = (mappedFrequencies[x + 1] - mappedFrequencies[x]) / topFrequency;
						// as long as the bandwidth is smaller than our fft resolution, we interpolate the points
						// otherwise, break out and sample the max values of the bins inside the bandwidth
						if (bwForLine > fftBandwidth)
							break;

						auto interpolatedFilter = cpl::lfilter<std::complex<ftype>, true>(csf, N, mappedFrequencies[x] * freqToBin, 5);

						wsp[x * 2] = interpolatedFilter.real() * invSize;
						wsp[x * 2 + 1] = interpolatedFilter.imag() * invSize;
					}
					oldBin = mappedFrequencies[x] * freqToBin;

					for (; x < numPoints; ++x)
					{
						maxLMag = maxRMag = newLMag = newRMag = 0;

						bin = static_cast<std::size_t>(mappedFrequencies[x] * freqToBin);
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
						filterStates.size(), dbRange.low, dbRange.high, flt);
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

						if (channelConfiguration == ChannelConfiguration::Separate)
							mapAndTransformDFTFilters<ChannelConfiguration::Separate, float>(filterStates.data(), wsp, filterResults.data(),
								filterStates.size(), dbRange.low, dbRange.high, flt);
						else
							mapAndTransformDFTFilters<ChannelConfiguration::Phase, float>(filterStates.data(), wsp, filterResults.data(),
								filterStates.size(), dbRange.low, dbRange.high, flt);

					}
					break;
				}
				break;
			}
			case Algorithm::MQDFT:
			{
				auto result = transformer.getTransformResult();
				auto totalData = filterStates.size() * numChannels;
				
				switch(channelConfiguration)
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
			}
			case Algorithm::RSNT:
			{
				/*float * wsp = getAudioMemory<float>();
				for (int i = 1; i < numFilters - 1; i++)
				{
					std::complex<float> zn1(resonator.real[i - 1], resonator.imag[i  - 1 ]);
					std::complex<float> z0(resonator.real[i], resonator.imag[i]);
					std::complex<float> z1(resonator.real[i + 1], resonator.imag[i + 1]);
					z0 = zn1 * -0.25f + z0 * 0.5f + z1 * -0.25f;

					wsp[i * 2] = z0.real();
					wsp[i * 2 + 1] = z0.imag();
				}*/

				float * wsp = getAudioMemory<float>();
				for (int i = 0; i < numFilters; i++)
				{
/*					std::complex<float> n1(sdft.real[0][i], sdft.imag[0][i]);
					std::complex<float> n2(sdft.real[1][i], sdft.imag[1][i]);
					std::complex<float> n3(sdft.real[2][i], sdft.imag[2][i]);


					std::complex<float> m1(-0.25, -0.25);
					std::complex<float> m(0.5, 0.5);
					
					auto fin = -0.25f * n1 + 0.5f * n2 + -0.25f * n3;*/

				//	wsp[i * 2] = (-0.25 * sdft.realState[1][i - 1] + 0.5 * sdft.realState[1][i] + -0.25 * sdft.realState[1][i + 1]) /*/ sdft.getBandwidth(i)*/;
				//	wsp[i * 2 + 1] = (-0.25 * sdft.imagState[1][i - 1] + 0.5 * sdft.imagState[1][i] + -0.25 * sdft.imagState[1][i + 1]) /*/ sdft.getBandwidth(i)*/;
				//	wsp[i * 2] = (-0.25 * sdft.realState[0][i] + 0.5 * sdft.realState[1][i] + -0.25 * sdft.realState[2][i]) / sdft.getBandwidth(i);
				//	wsp[i * 2 + 1] = (-0.25 * sdft.imagState[0][i] + 0.5 * sdft.imagState[1][i] + -0.25 * sdft.imagState[2][i]) / sdft.getBandwidth(i);
					//wsp[i * 2] = (-0.25 * sdft.real[1][i - 1] + 0.5 * sdft.real[1][i] + -0.25 * sdft.real[1][i + 1]);
					//wsp[i * 2 + 1] = (-0.25 * sdft.imag[1][i - 1] + 0.5 * sdft.imag[1][i] + -0.25 * sdft.imag[1][i + 1]);
				//wsp[i * 2] = sdft.realState[1][i] / sdft.getBandwidth(i);
				//wsp[i * 2 + 1] = sdft.imagState[1][i] / sdft.getBandwidth(i);


					auto const & c = sdft.getWindowedResonanceAt<cpl::dsp::WindowTypes::Hann, false>(i,0);
					wsp[i * 2] = c.real();
					wsp[i * 2 + 1] = c.imag();
				}
				mapAndTransformDFTFilters<ChannelConfiguration::Left, float>(filterStates.data(), wsp, filterResults.data(),
																   filterStates.size(), dbRange.low, dbRange.high, flt);
			}
			break;
		}
		
		
	}
	
	void CSpectrum::doTransform()
	{
		auto const channelConfiguration = cpl::Math::round<int>(kchannelConf->bGetValue() * 4);
		switch (algorithmType)
		{
			case Algorithm::FFT:
			{
				auto const numSamples = getNumAudioElements<std::complex<double>>();
				switch (channelConfiguration)
				{
				case ChannelConfiguration::Left:
				case ChannelConfiguration::Right:
				case ChannelConfiguration::Merge:
					transformer.fft<1>(getAudioMemory<double>(), numSamples);
					break;
				case ChannelConfiguration::Phase:
				case ChannelConfiguration::Separate:
					transformer.fft<2>(getAudioMemory<double>(), numSamples);

				}

				break;
			}
			case Algorithm::MQDFT:
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
			/* Resonators dont need preparation as they are continious and real time*/
			case Algorithm::RSNT:
				break;
		}
	}
	
	void CSpectrum::computeWindowKernel()
	{
		size_t sampleSize = getWindowSize();

		int i = 0;
		for (; i < sampleSize; ++i)
		{
			windowKernel[i] = 1.0 - std::cos(TAU * i / (sampleSize - 1));
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
		if (isFrozen)
			return;
		const bool locking = shouldSynchronize();
		auto size = getWindowSize(); // the size of the transform, containing samples
		// the quantized (to next power of 2) samples of this transform
		// that is, the size + additional zero-padding
		auto fullSize = getNumAudioElements<std::complex<double>>();

		auto const channelConfiguration = cpl::Math::round<int>(kchannelConf->bGetValue() * 4);

		numChannels = 1 + (channelConfiguration > ChannelConfiguration::Merge);

		{
			// acquire locks on audio buffegrs
			std::vector<cpl::CMutex> locks;
			if (locking)
			for (auto & buf : audioData)
				locks.emplace_back(buf);

			// note here, that since our audiobuffer is a circular buffer,
			// index[0] is actually the last input, and index[size - 1] is the first
			// while the fft doesn't care whether the input is reversed or not, 
			// the varying Q-transforms do care, since they vary the amount of samples
			// considered.

			// generate an oscillator, which describes the hanning window
			cpl::CFastOscillator<float> winOsc;
			winOsc.reset(size, 1, M_PI / 2); // offset phase by 90 degrees to create a cosine instead.
			switch (algorithmType)
			{
				case Algorithm::FFT:
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
						for (; i < halfInputSize; ++i)
						{
							buffer[i + offset] = audioData[channel].singleCheckAccess(i) * windowKernel[i];
						}

						offset = halfInputSize;
						for (; i < size; ++i)
						{
							buffer[i - offset] = audioData[channel].singleCheckAccess(i) * windowKernel[i];
						}
						break;
					}
					case ChannelConfiguration::Phase:
					case ChannelConfiguration::Separate:
					{
						for (; i < halfInputSize; ++i)
						{
							buffer[i + offset] = std::complex<double>
							(
							audioData[0].singleCheckAccess(i)  * windowKernel[i] ,
							audioData[1].singleCheckAccess(i)  * windowKernel[i]
							);
						}

						offset = halfInputSize;
						for (; i < size; ++i)
						{
							buffer[i - offset] = std::complex<double>
							(
								audioData[0].singleCheckAccess(i) * windowKernel[i],
								audioData[1].singleCheckAccess(i) * windowKernel[i]
							);
						}
						break;
					}
					case ChannelConfiguration::Merge:
						for (; i < halfInputSize; ++i)
						{
							buffer[i + offset] = (audioData[0].singleCheckAccess(i) + audioData[1].singleCheckAccess(i)) 
								* windowKernel[i] * 0.5;
						}

						offset = halfInputSize;
						for (; i < size; ++i)
						{
							buffer[i - offset] = (audioData[0].singleCheckAccess(i) + audioData[1].singleCheckAccess(i))
								* windowKernel[i] * 0.5;;
						}
					}

					/*
						zero-pad until buffer is filled
					*/

					offset = halfSize + halfPadding;
					for (size_t pad = halfSize - halfPadding; pad < offset; ++pad)
					{
						buffer[pad] = 0;
					}


					break;
				}
				// this case is different from CDFT, since the input musnt be windowed
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
				}
			}
		}

	}
};