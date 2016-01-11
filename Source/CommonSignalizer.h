#ifndef COMMON_SIGNALIZER_H
	#define COMMON_SIGNALIZER_H

	#include <cpl/Common.h>
	#include <cpl/CViews.h>
	#include <cpl/CAudioStream.h>
	#include <complex>

	namespace Signalizer
	{
		typedef cpl::CAudioStream<float, 32> AudioStream;

		enum class ChannelConfiguration
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
	};

#endif