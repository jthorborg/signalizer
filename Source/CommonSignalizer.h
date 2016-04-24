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
	#include <cpl/gui/gui.h>
	#include <cpl/CAudioStream.h>
	#include <complex>

	namespace Signalizer
	{
		// TODO: Figure out why sizes around 256 causes buffer overruns
		typedef cpl::CAudioStream<float, 64> AudioStream;
		typedef std::pair<cpl::CBaseControl *, cpl::iCtrlPrec_t> CtrlUpdate;
		
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