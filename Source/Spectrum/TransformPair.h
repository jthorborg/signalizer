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
#include <cpl/ffts.h>
#include <cpl/lib/uarray.h>
#include <array>
#include <optional>

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
		typedef T ProcessingType;
		typedef cpl::simd::consts<T> consts;
		typedef std::array<AudioStream::AudioBufferView, 2> AudioPair;

		struct LineGraphDesc
		{
			friend class TransformPair;

			/// <summary>
			/// The decay/peak-filtered and scaled outputs of the transforms,
			/// with each element corrosponding to a complex output pixel of getAxisPoints() size.
			/// Resized in displayReordered
			/// </summary>
			cpl::uarray<const UComplex> getResults(std::size_t size) const
			{
				results.resize(size);
				return cpl::as_uarray(results);
			}

			void resize(std::size_t n)
			{
				// TODO: This is only called from mapAndTransformDFTFilters. Size can go out of sync with drawing code that uses getResults().
				states.resize(n); results.resize(n);
			}

			void zero() {
				std::memset(states.data(), 0, states.size() * sizeof(UComplex));
				std::memset(results.data(), 0, results.size() * sizeof(UComplex));
			}

		private:
			/// <summary>
			/// The'raw' formatted state output of the mapped transform algorithms.
			/// </summary>
			cpl::aligned_vector<UComplex, 32> states;
			mutable cpl::aligned_vector<UComplex, 32> results;
		};

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
		void resonatingDispatch(const Constant& constant, std::array<AFloat*, 2> buffer, std::size_t numSamples);

		template<typename ISA>
		void audioEntryPoint(const Constant& constant, const std::optional<AudioPair>& pairs, std::array<AFloat*, 2> buffer, std::size_t numSamples);

		/// <summary>
		/// Maps the current resonating system according to the current model (linear/logarithmic) and the current
		/// subsection of the complete spectrum such that a linear array of output data matches pixels 1:1, as well as
		/// formats the data into the filterResults array according to the channel mode (SpectrumChannels).
		///
		/// Call prepareTransform(), then doTransform(), then mapToLinearSpace().
		/// After the call to mapToLinearSpace, the results are written to getTransformResults().
		/// </summary>
		void mapToLinearSpace(const Constant& constant);

		/// <summary>
		/// For some transform algorithms, it may be a no-op, but for others (like FFTs) that may need zero-padding
		/// or windowing, this is done here.
		///
		/// Call prepareTransform(), then doTransform(), then mapToLinearSpace()
		/// </summary>
		bool prepareTransform(const Constant& constant, const AudioPair& audio);

		/// <summary>
		/// For some transform algorithms, it may be a no-op, but for others (like FFTs) that may need zero-padding
		/// or windowing, this is done here.
		///
		/// This functions considers the additional arguments as more recent audio than the audio buffers (and as of such, considers numSamples less audio from
		/// the first argument).
		/// </summary>
		bool prepareTransform(const Constant& constant, const AudioPair& audio, std::array<AFloat*, 2> preliminaryAudio, std::size_t numSamples);

		/// <summary>
		/// Again, some algorithms may not need this, but this ensures the transform is done after this call.
		///
		/// Call prepareTransform(), then doTransform(), then mapToLinearSpace()
		/// Needs exclusive access to audioResource.
		/// </summary>
		void doTransform(const Constant& constant);

		/// <summary>
		/// Runs the transform (of any kind) results through potential post filters and other features, before displaying it.
		/// The transform will be rendered into filterResults after this.
		/// </summary>
		template<class InVector>
		void postProcessTransform(const Constant& constant, const InVector& transform);

		/// <summary>
		/// Post processes the transform that will be interpreted according to what's selected.
		/// </summary>
		void postProcessStdTransform(const Constant& constant);

		/// <summary>
		/// Returns an array of "axis points" size. 
		/// </summary>
		cpl::uarray<const T> getTransformResult(const Constant& constant);

		cpl::uarray<const std::complex<T>> getRawFFT(const Constant& constant)
		{
			return getAudioMemory<std::complex<T>>(constant.transformSize);
		}

		void clearLineGraphStates()
		{
			for (std::size_t i = 0; i < lineGraphs.size(); ++i)
			{
				lineGraphs[i].zero();
			}
		}

		void clearAudioState()
		{
			clearLineGraphStates();

			std::fill(workingMemory.begin(), workingMemory.end(), 0);
			std::fill(audioMemory.begin(), audioMemory.end(), 0);
			cresonator.resetState();
		}

		// dsp objects -- TODO: Make private?
		std::array<LineGraphDesc, SpectrumContent::LineGraphs::LineEnd> lineGraphs;

	private:

		/// <summary>
		/// All inputs must be normalized. Scales the input to the display decibels, and runs it through peak filters.
		/// newVals = current vector of floats / doubles * 2 (complex), output from CSignalTransform::**dft() of size * 2
		/// for mode = left / merge / mid / side / right
		/// 	newVals is a complex vector of floats of size
		/// for mode = separate, mid&side
		/// 	newVals is a complex vector of floats of size * 2
		/// 	newVals[n * 2 + 0] = lreal
		/// 	newVals[n * 2 + 1] = limag
		/// 	newVals[n * 2 + size + 0] = rreal
		/// 	newVals[n * 2 + size + 1] = rimag
		/// for mode = phase
		/// 	newVals[n * 2 + 0] = mag
		/// 	newVals[n * 2 + 1] = phase cancellation(with 1 being totally cancelled)
		/// </summary>
		template<class V2>
		void mapAndTransformDFTFilters(const Constant& constant, const V2& newVals, std::size_t size);

		template<typename ISA>
		void addAudioFrame(const Constant& constant);

		template<typename Y>
		cpl::uarray<Y> getWork(std::size_t size)
		{
			static_assert(sizeof(Y) <= sizeof(std::complex<T>));

			if(workingMemory.size() < size)
				workingMemory.resize(size);

			return cpl::as_uarray(workingMemory).reinterpret<Y>().slice(0, size);
		}

		template<typename Y>
		cpl::uarray<Y> getAudioMemory(std::size_t size) 
		{
			static_assert(sizeof(Y) <= sizeof(std::complex<T>));

			if (audioMemory.size() < size)
				audioMemory.resize(size);

			return cpl::as_uarray(audioMemory).reinterpret<Y>().slice(0, size);
		}

		/// <summary>
		/// Temporary memory buffer for other applications.
		/// </summary>
		cpl::aligned_vector<std::complex<T>, 32> workingMemory;

		/// <summary>
		/// Temporary memory buffer for audio applications. Resized in setWindowSize (since the size is a function of the window size)
		/// </summary>
		cpl::aligned_vector<std::complex<T>, 32> audioMemory;
		// TODO: Change to T
		cpl::dsp::CComplexResonator<T, 2> cresonator;
		std::size_t currentCounter{};

	};
}


#endif
