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

	file:TransformPair.h

		Instance data and processing functions for a pair of signals

*************************************************************************************/

#ifndef SIGNALIZER_TRANSFORM_PAIR_H
#define SIGNALIZER_TRANSFORM_PAIR_H

#include "../Signalizer.h"
#include <cpl/Common.h>
#include <cpl/dsp/CComplexResonator.h>
#include <memory>
#include "SpectrumParameters.h"
#include <cpl/simd.h>
#include <cpl/dsp/CComplexResonator.h>

namespace Signalizer
{
	template<typename T>
	class TransformConstant;

	template<typename T>
	class TransformPair
	{
	public:

		typedef AudioStream::DataType AFloat;
		typedef UComplexFilter<AFloat> UComplex;
		typedef cpl::aligned_vector<UComplex, 32> FrameVector;
		typedef TransformConstant<T> Constant;

		/// <summary>
		/// The complex resonator used for iir spectrums
		/// </summary>
		std::vector<FrameVector> sfbuf;
		/// <summary>
		/// Copies the state from the complex resonator into the output buffer.
		/// The output vector is assumed to accept index assigning of std::complex of fpoints.
		/// It is assumed the output vector can hold at least numChannels (of configuration) times
		/// numFilters.
		/// Output of channels are stored at numFilters offsets.
		///
		/// Returns the total number of complex samples copied into the output
		/// </summary>
		template<typename ISA, class Vector>
		std::size_t copyResonatorStateInto(const Constant& constant, cpl::dsp::WindowTypes windowType, Vector& output, std::size_t outChannels);

		template<typename ISA>
		void resonatingDispatch(const Constant& constant, AFloat** buffer, std::size_t numChannels, std::size_t numSamples);

		template<typename ISA>
		void audioEntryPoint(const Constant& constant, AudioStream::ListenerContext& ctx, AFloat** buffer, std::size_t numChannels, std::size_t numSamples);

		/// <summary>
		/// Maps the current resonating system according to the current model (linear/logarithmic) and the current
		/// subsection of the complete spectrum such that a linear array of output data matches pixels 1:1, as well as
		/// formats the data into the filterResults array according to the channel mode (SpectrumChannels).
		///
		/// Call prepareTransform(), then doTransform(), then mapToLinearSpace().
		/// After the call to mapToLinearSpace, the results are written to getTransformResults().
		/// Returns the complex amount of filters processed.
		/// </summary>
		void mapToLinearSpace(const Constant& constant);

		/// <summary>
		/// For some transform algorithms, it may be a no-op, but for others (like FFTs) that may need zero-padding
		/// or windowing, this is done here.
		///
		/// Call prepareTransform(), then doTransform(), then mapToLinearSpace()
		/// Needs exclusive access to audioResource.
		/// </summary>
		bool prepareTransform(const Constant& constant, const AudioStream::AudioBufferAccess& audio);

		/// <summary>
		/// For some transform algorithms, it may be a no-op, but for others (like FFTs) that may need zero-padding
		/// or windowing, this is done here.
		///
		/// This functions considers the additional arguments as more recent audio than the audio buffers (and as of such, considers numSamples less audio from
		/// the first argument).
		/// Needs exclusive access to audioResource.
		/// </summary>
		bool prepareTransform(const Constant& constant, const AudioStream::AudioBufferAccess& audio, AFloat** preliminaryAudio, std::size_t numChannels, std::size_t numSamples);

		/// <summary>
		/// Again, some algorithms may not need this, but this ensures the transform is done after this call.
		///
		/// Call prepareTransform(), then doTransform(), then mapToLinearSpace()
		/// Needs exclusive access to audioResource.
		/// </summary>
		void doTransform(const Constant& constant);

		/// <summary>
		/// Returns an array of "axis points" size. 
		/// </summary>
		template<typename T>
		const T* getTransformResult() const noexcept
		{
			return getWorkingMemory<T>();
		}

		/// <summary>
		/// </summary>
		const std::complex<T>* getRawFFT() const noexcept
		{
			return getAudioMemory<std::complex<T>>();
		}

		void setStorage(const Constant& constant)
		{
			// some cases it is nice to have an extra entry (see handling of
			// separating real and imaginary transforms)
			// TODO: type these in terms of complex, not god damn chars.
			audioMemory.resize(constant.transformSize + 1);
			workingMemory.resize(constant.axisPoints * 2 * sizeof(std::complex<double>));
		}

		void clearAudioState()
		{
			std::fill(workingMemory.begin(), workingMemory.end(), 0);
			std::fill(audioMemory.begin(), audioMemory.end(), 0);
			cresonator.resetState();
		}

		void remapResonator(Constant& constant)
		{
			cresonator.match(constant.Resonator);
		}

	private:

		std::size_t getStateConfigurationChannels() const noexcept;
		void checkInvariants();

		template<typename T>
		std::size_t getNumWorkingElements() const noexcept;

		template<typename ISA>
		void addAudioFrame(const Constant& constant);


		template<typename T>
		T* getWorkingMemory() noexcept
		{
			return reinterpret_cast<T*>(workingMemory.data());
		}

		template<typename T>
		const T* getWorkingMemory() const noexcept
		{
			return reinterpret_cast<const T*>(workingMemory.data());
		}

		template<typename T>
		T* getAudioMemory() noexcept
		{
			return reinterpret_cast<T*>(audioMemory.data());
		}

		template<typename T>
		const T* getAudioMemory() const noexcept
		{
			return reinterpret_cast<const T*>(audioMemory.data());
		}


		/// <summary>
		/// used for 'real-time' audio processing, when the incoming buffer needs to be rearranged without changing it.
		/// </summary>
		struct RelayBuffer
		{
			cpl::aligned_vector<AFloat, 32> buffer;
			std::size_t samples{}, channels{};
		} relay;

		/// <summary>
		/// Gets the relay buffer for the channel.
		/// Note: you should call ensureRelayBufferSize before.
		/// </summary>
		AFloat* getRelayBufferChannel(std::size_t channel);
		/// <summary>
		/// Ensures the storage for the relay buffer.
		/// Not suited for real-time usage (allocates memory)
		/// </summary>
		void ensureRelayBufferSize(std::size_t channels, std::size_t numSamples);

		/// <summary>
		/// Temporary memory buffer for other applications.
		/// </summary>
		cpl::aligned_vector<char, 32> workingMemory;

		/// <summary>
		/// Temporary memory buffer for audio applications. Resized in setWindowSize (since the size is a function of the window size)
		/// </summary>
		cpl::aligned_vector<std::complex<T>, 32> audioMemory;
		// TODO: Change to T
		cpl::dsp::CComplexResonator<AFloat, 2> cresonator;
		std::size_t currentCounter{};
	};
}


#endif
