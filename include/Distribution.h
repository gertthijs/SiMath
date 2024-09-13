/*
 *  Distribution.h
 *
 *  Created by Gert Thijs on 17/05/06.
 *  Copyright 2006 Silicos NV. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     + Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     + Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY SILICOS NV AND CONTRIBUTORS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL SILICOS NV AND CONTRIBUTORS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
		* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 *
 */


#ifndef __SIMATH_DISTRIBUTION__
#define __SIMATH_DISTRIBUTION__

#include "Vector.h"
#include "Utilities.h"

namespace SiMath
{
	/**
		\class DistributionError
	 \brief  Exception class
	 
	 This exception class is derived from the std::runtime_error class and is used to 
	 indicate that something went with a Distribution object.
	 \ingroup _error_handling_
	 */	
	class SIMATH_EXCEPTION DistributionError : public SiMathError
	{	
		public:
			DistributionError(const std::string & e) : SiMathError("DistributionError : " + e ) {};
	};


	/**
		\class Distribution
		\brief Base class for the representation of a continuous distribution
	 
		A continuous distribution is represented by three parameters:
			- scale
			- offset
			- shape
	 
	 \todo Add sampling method to distributions
	 */
	class Distribution
	{
		protected:
			double _scale;     ///< scale of the distribution
			double _offset;    ///< offset of the distribution
			double _shape;     ///< shape of the distribution
			
		public:
				///\name Structors
				//@{
				/**
					\brief Default constructor
				 */
				Distribution() : _scale(1.0), _offset(0), _shape(1.0) {};
			
				/**
					\brief Constructor
					\param s Scale parameter of the distribution
					\param o Location or offset of the distribution
					\param sh Shape parameter of the distribution
				 */
				Distribution(double s, double o = 0.0, double sh = 1.0) : _scale(s), _offset(o), _shape(sh) {};
				
				/**
					\brief Default destructor
				 */
				virtual ~Distribution() {};
				//@}
			
				///\name Update distribution parameters
				//@{
				/**
					\brief Set the scale parameter of the distribution
					\param v Value to be set
				 */
				virtual void setScale(double v) { _scale = v; };

				/**
					\brief Set the offset parameter of the distribution
					\param v Value to be set
				 */
				virtual void setOffset(double v) { _offset = v;};

				/**
					\brief Set the scale parameter of the distribution
					\param v Value to be set
				 */
				virtual void setShape(double v) { _shape = v;};
				
				/**
					\brief Main method to learn the parameters from data
					\param data   Vector of real valued data.
					\return Status variable to indicate whether ot not the training was succesful.
				 */
				virtual bool train(SiMath::Vector & data) = 0;
				//@}
				
				///\name Distribution statistics
				//@{
				virtual double getScale() const { return _scale;};
				virtual double getOffset() const { return _offset;};
				virtual double getShape() const {return _shape;};
				virtual double mean() const = 0;
				virtual double stDev() const = 0;
				virtual double skewness() const = 0;
				virtual double kurtosis() const = 0;
				//@}
				
				///\name Evaluators
				//@{
				/**
					\brief Evaluate the probability density function at x
					\param x input value
				 */
				virtual double pdf(double x) = 0;

				/**
					\brief Evaluate the cumulative density function at x
					\param x input value
				 */
				virtual double cdf(double x) = 0;

				/**
					\brief Evaluate the survival function at x
					\param x input value
				 */
				virtual double sf(double x) = 0;

				/**
					\brief Draw a sample from the distribution
				 */
				virtual double sample() = 0;

				/**
					\brief Draw n samples from the distribution
				 */
				virtual Vector sample(unsigned int n=1) = 0;
				//@}
	};

	
	
	/**
		\class NormalDist
		\brief Implementation of the normal distribtuion
	 */
	class NormalDist : Distribution 
	{
		public:	
			///\name Structors
			//@{
			/**
				\brief Default constructor sets offset to 0 and scale to 1.0
			 */
			NormalDist() : _c1(0), _c2(0), Distribution() {};
		
			/**
				\brief Constructor
				\param s Scale parameter (standard deviation of the normal)
				\param m Offset parameter (is the mean of the normal)
				\param sh Shape parameter (defaults to 1.0 and is actually neglected)
			 */
			NormalDist(double s, double m, double sh=1.0);
			
			/**
				\brief Destructor
			 */
			~NormalDist() {};
			
			///\name Parameter update
			//@{
			/**
				\brief Function to learn the parameters from data
				\param data Vector of real valued data.
			 */
			virtual bool train(SiMath::Vector & data);

			/**
				\brief Set the scale of the distribution 
				\param v Scale value
				\throw DistributionError if value is smaller than 0.0
			 */
			virtual void setScale(double v);
			//@}
			
			///\name Distribution statistics
			//@{
			virtual double mean() const {return _offset;};     ///< Get the mean of the distribution ( = offset)
			virtual double stDev() const {return _scale;};     ///< Get the standard deviation of the distribution ( = scale)
			virtual double skewness() const {return 0;};       ///< Get the skewness of the distribution ( = 0)
			virtual double kurtosis() const {return 3;};       ///< Get the curtosis of the distribution ( = 3)
			//@}
			
			///\name Evaluators
			//@{
			/**
				\brief Evaluate the probability density function at x
			 \param x input value
			 */
			virtual double pdf(double x);
			
			/**
				\brief Evaluate the cumulative density function at x
			 \param x input value
			 */
			virtual double cdf(double x);
			
			/**
				\brief Evaluate the survival function at x
			 \param x input value
			 */
			virtual double sf(double x);
			
			/**
				\brief Draw a sample from the distribution
			 
				Normal random number generation is based on the Box-Muller transform
				
				\return Vector of sampled values
			 */
			virtual double sample();
			
			/**
				\brief Draw n samples from the normal distribution
			 
				Normal random number generation is based on the Box-Muller transform
			 
				\param n Number of sample to draw
				\return Vector of sampled values
			 */
			virtual Vector sample(unsigned int n);
			//@}
			
		private:
			double _c1;    ///< Local helper variable ( = 1/(sqrt(2PI)* _scale))
			double _c2;    ///< Local helper variable ( = 1/(2*_scale^2))
	};
	

		
	class GammaDist : Distribution 
	{
		public:	
			///\name Structors
			//@{
			/**
				\brief Default constructor
			 */
			GammaDist() : Distribution() {};
		
			/**
				\brief Constructor
				\param s Scale parameter  (default is 1.0)
				\param o Offset parameter (default is 0.0)
				\param sh Shape parameter (default is 1.0)
			 */
			GammaDist(double s=1.0, double o=0.0, double sh=1.0);
			
			/**
				\brief Destructor
			 */
			~GammaDist() {};
			//@}
			
			///\name Update parameters
			//@{
			virtual void setScale(double v);
			virtual void setShape(double v);
			virtual bool train(SiMath::Vector & data);
			//@}
			
			///\name Distribution statistics
			//@{
			virtual double mean() const { return _shape; };                 ///< Get the mean of the distribution ( = shape)
			virtual double stDev() const { return sqrt(_shape); };          ///< Get the standard deviation of the distribution ( = sqrt(shape))
			virtual double skewness() const { return 2.0/sqrt(_shape); };   ///< Get the skewness of the distribution ( = 2/sqrt(shape))
			virtual double kurtosis() const { return 3.0 + 6.0/_shape; };   ///< Get the kurtosis of the distribution ( = 3 + 6/shape))
			//@}
			
			
			///\name Evaluators
			//@{
			/**
				\brief Evaluate the probability density function at x
			 \param x input value
			 */
			virtual double pdf(double x);
			
			/**
				\brief Evaluate the cumulative density function at x
			 \param x input value
			 */
			virtual double cdf(double x);
			
			/**
				\brief Evaluate the survival function at x
			 \param x input value
			 */
			virtual double sf(double x);
			
			/**
				\brief Draw a sample from the distribution
			 */
			virtual double sample() { throw(DistributionError("sampling not yet implemented.")); };
			
			/**
				\brief Draw n samples from the distribution
			 \param n Number of sample to draw
			 \return Vector of sampled values
			 */
			virtual Vector sample(unsigned int n) { throw(DistributionError("sampling not yet implemented.")); };
			//@}
			
	};
	
	
	class MaxExtremeValueDist : Distribution 
	{
		public:	
			///\name Structors
			//@{
			/**
				\brief Default constructor
			 */
			MaxExtremeValueDist() : Distribution() {};
		
			/**
				\brief Constructor
				\param s Scale parameter  (default is 1.0)
				\param o Offset parameter (default is 0.0)
			  \param sh Shape parameter (default is 1.0)
			 */
			MaxExtremeValueDist(double s=1.0, double o=0.0, double sh=1.0) : Distribution(s,o,sh) {};
			
			/**
				\brief Destructor
			 */
			~MaxExtremeValueDist() {};
			//@}
			
			///\name Update parameters
			//@{
			virtual bool train(SiMath::Vector & data);
			//@}
			
			///\name Distribution statistics
			//@{
			virtual double mean() const { return _offset + 0.5772*_scale; };    ///< Get the mean of the distribution
			virtual double stDev() const { return PI * _scale / sqrt(6.0); };   ///< Get the standard deviation of the distribution
			virtual double skewness() const { return 1.13955; };                ///< Get the skewness of the distribution ( = 1.13955)
			virtual double kurtosis() const { return 5.4; };                    ///< Get the kurtosis of the distribution ( = 5.4)
			//@}
			
			///\name Evaluators
			//@{
			/**
				\brief Evaluate the probability density function at x
			 \param x input value
			 */
			virtual double pdf(double x);
			
			/**
				\brief Evaluate the cumulative density function at x
			 \param x input value
			 */
			virtual double cdf(double x);
			
			/**
				\brief Evaluate the survival function at x
			 \param x input value
			 */
			virtual double sf(double x);
			
			/**
				\brief Draw a random sample from the distribution
			 */
			virtual double sample();
			
			/**
				\brief Draw n random samples from the distribution
			 \param n Number of sample to draw
			 \return Vector of sampled values
			 */
			virtual Vector sample(unsigned int n);
			//@}
			
	};
	
	
	
	
	
};

#endif  __SIMATH_DISTRIBUTION__
