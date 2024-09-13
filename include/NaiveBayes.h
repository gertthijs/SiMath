/*
 *  NaiveBayes.h
 *
 *  Created by Gert Thijs on 27/03/06.
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


#ifndef   __SIMATH_NAIVE_BAYES__
#define   __SIMATH_NAIVE_BAYES__

#include "Definitions.h"
#include "Matrix.h"

namespace SiMath
{
	/**
		\class  NaiveBayesError
		\brief  Exception class
	 
		This exception class is derived from the std::runtime_error class and is used to 
		indicate that something went wrong during the naive bayesian classifier algorithm.
	 \ingroup _error_handling_
	 */	
	class SIMATH_EXCEPTION NaiveBayesError : public SiMathError
	{
		public:
			NaiveBayesError(const std::string & e) : SiMathError("NaiveBayesError : " + e) {};
	};
	
	/**
		\struct NaiveBayesParameters
		\brief Basic parameters of a naive bayesian classifier
		\ingroup _model_parameters_
	 
		Usage:
		\code
		SiMath::NaiveBayesParameters p;
		p.dimension = A.nbrColumns();
		p.nbrClasses = 2;
		\endcode
	 */
	struct NaiveBayesParameters {
		unsigned int dimension;    ///< Dimension of the data space
		unsigned int nbrClasses;   ///< Number of classes in the model (should be greater than 2)
		
		/**
			\brief Default parameters:
			- dimension = 0,
			- nbrClasses = 2
		 */
		NaiveBayesParameters() : dimension(0), nbrClasses(2) {};
	};
	
	/**
		\brief Class to represent a naive Bayesian classifier
	 
		A naive Bayesian classifier is based on the assumptions that all variables
	  are independent and the probability of a instance belonging to a class is equal 
		to the product of the individual variables belonging to that class. This is can be computed as
	 
		\f[
			p(C|A) = 1/B * p(C) * \prod_{i=1}^{n} p(C|A[i])
		\f]
		Here C is the class label, A is a data vector of length n and B is a normalisation factor.
		This implementation provides currently a gaussian distribution as model 
		for all variables.
		
		\todo Make it possible to change the type of distribution
	 */
	class NaiveBayes
	{
		public:
			///\name Structors
			//@{
			/** 
				\brief Default constructor from predefined set of parameters
				\param p Set of parameters
			*/
			NaiveBayes(const NaiveBayesParameters & p);
			
			/**
				\brief Destructor
			 */
			~NaiveBayes();
			//@}
			
			///\name Classification
			//@{
			/**
				\brief Function to classify a single data vector
				\param data Reference to single data Vector
				\return Vector of log-probabilities for each class
				\exception NaiveBayesError if the dimension of the data vector do not match the model
			 */
			Vector classify(const Vector & data);

			/**
				\brief Function to classify a single data vector
				\param data Reference to data Matrix
				\return Matrix of log-probabilities for each class
			  \exception NaiveBayesError if the dimension of the data Matrix do not match the model
			 */
			Matrix classify(const Matrix & data);
			//@}
			
			///\name Model definition
			//@{
			/**
				\brief Main method train a specific class in the classifier.
			 
			 \param c Class number
			 \param data reference to Matrix that holds the data set
			 \return nothing
			 */
			void train(unsigned int c, const Matrix & data);
			
			/**
				\brief Define the parameters of a specific class 
				\param c   Class number
				\param m   Mean of the class variables
				\param sd  Std. deviation of the class variables
			 */
			void setModel(unsigned int c, const Vector & m, const Vector & sd);
			
			/**
				\brief Update the class prior
				\param c class number
				\param d class prior
			 */
			void setClassPrior(unsigned int c, double d);
			//@}
			
			/**
				\brief Define the number of classes in the model
				\param n number of classes 
			 */
			void setNbrOfClasses(unsigned int n);

			///\name Model parameter inspector
			//@{
			/**
				\brief Get prior of class c
				\param c class number
			 */
			double getClassPrior(unsigned int c);

			/**
				\brief Get mean of variables of class c
			 \param c class number
			 */
			Vector getClassMean(unsigned int c);
			
			/**
				\brief Get st.dev of variables of class c
				\param c class number
			 */
			Vector getClassStDev(unsigned int c);

			/**
				\brief Check if class c is trained
				\param c class number
			 */
			bool isClassTrained(unsigned int c);
			//@}
			
		private:
			/**
				\struct ClassDefinition
			 */
			struct ClassDefinition {
				Vector   mean;         ///< Vector to hold the mean of each variable
				Vector   stdev;        ///< Vector to hold the standard deviation of each variable
				double   prior;        ///< Class prior (should be between 0 and 1)
				bool     trained;      ///< Status flag
				
				ClassDefinition() : mean(0), stdev(0), prior(1.0), trained(false) {}; 
			};
			
			unsigned int _nbrClasses;  ///< Number of classes in the model
			unsigned int _nbrVar;      ///< Number of variables
			
			std::vector<ClassDefinition> _classes;  ///< Vector of class representations
			
	};
	
};

#endif    __SIMATH_NAIVE_BAYES__
