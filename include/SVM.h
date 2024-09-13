/*
 *  SVM.h
 *
 *  Created by Gert Thijs.
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


#ifndef _SIMATH_SVM__
#define _SIMATH_SVM__

#include <cmath>
#include <string>
#include <vector>
#include <map>
#include "Matrix.h"
#include "Vector.h"


namespace SiMath
{

	/**
		\class SVMError
		\brief  Exception class
	 
		This exception class is derived from the std::runtime_error class and is used to 
		indicate that something went wrong when working with an SVM.
	 \ingroup _error_handling_
	 */	
	class SIMATH_EXCEPTION SVMError : public SiMathError
	{
		public:
		SVMError(const std::string & e) : SiMathError("SVMError : " + e) {};
	};

	/**
		\enum SVMType
		\brief Definition of the type of SVM
	 */
	enum SVMType { 
		C_SVC, 
		NU_SVC, 
		ONE_CLASS, 
		EPSILON_SVR, 
		NU_SVR 
	};

	/**
		\enum KernelType
		\brief Definition of the type of kernel used. Currently, fives types are available.
	 */
	enum KernelType { 
		LINEAR, 
		POLY, 
		RBF, 
		SIGMOID,
		SPLINE
	};

	/** 
		\enum AlphaStatus
		\brief Defintion of different types of alpha 
	 */
	enum AlphaStatus { 
		LOWER_BOUND,
		UPPER_BOUND,
		FREE 
	};

	const std::string SVMTypeTable[] =
	{
		"c_svc",          ///< C-based SVM for classification
		"nu_svc",         ///< nu-based SVM for classification
		"one_class",      ///< One class SVM for novelty detection
		"epsilon_svr",    ///< Epsilon SVM  for regression
		"nu_svr"          ///< nu SVM for regression
	};
	
	const std::string KernelTypeTable[] =
	{
		"linear",
		"polynomial",
		"rbf",
		"sigmoid",
		"spline"
	};
	
	// structure to hold the data set
	struct SVMProblem
	{
		int nbrData;    ///< Number of data points in the problem
		int nbrVar;     ///< Dimension of the data space
		Matrix x;       ///< Input data
		Vector y;       ///< Target data, might be values for regression or class labels
	};
			
	/**
		\struct SVMParameters 
		\brief  Parameters to setup a SVM
		\ingroup _model_parameters_
	 
	 Usage:
	 \code
	 SiMath::SVMParameters p;
	 p.KernelType = SiMath::RBF;
	 p.gamma = 0.3;
	 p.cacheSize = 100; // 100 Mb cache
	 p.C = 5;
	 \endcode
	 
	 */
	struct SVMParameters
	{
		SVMType      svmType;     ///< Type of the SVM (C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR)
		
		///\name Kernel parameters
		//@{
		KernelType   kernelType;  ///< Type of the kernel function to be used
		unsigned int degree;	    ///< Degree of the polynomial kernel
		double gamma;	            ///< Gamma parameter for polynomial, rbf, and sigmoidal kernel
		double coef0;	            ///< Offset for polynomial and sigmoidal kernel
		double cacheSize;         ///< size of the kernel cache in MB
		//@}
		
		// training parameters
		double eps;	            ///< stopping criteria
		double C;               ///< Cost of constraint violation for C_SVC, EPSILON_SVR and NU_SVR
		double nu;	            ///< for NU_SVC, ONE_CLASS, and NU_SVR
		double p;	              ///< for EPSILON_SVR
		int shrinking;	        ///< Define use of the shrinking heuristics
		bool probability;       ///< Define whether or not to use probability estimates
		std::map<int, double> weights;  ///< Storage of class label and corresponding weight for C_SVC
		
		/**
			\brief Default parameters:
		 - svmType = C_SVC,
		 - kernelType = LINEAR,
		 - degree = 0,
		 - gamma = 0.0,
		 - coef0 = 0.0,
		 - cacheSize = 100,
		 - eps = 0.001,
		 - C = 1,
		 - weights (empty),
		 - nu = 0.5,
		 - p = 0.1,
		 - shrinking = 1,
		 - probability = 0
		 */
		SVMParameters() : svmType(C_SVC), kernelType(LINEAR), degree(0), gamma(0.0), coef0(0.0), 
			cacheSize(100), eps(0.001), C(1), weights(), nu(0.5), p(0.1), shrinking(1), probability(0) {}; 
	};
	
	/**
		\struct DecisionFunction
		\brief Struct to hold the decision function of the trained SVM
		*/
	struct DecisionFunction
	{
		std::vector<double> alpha; 	 ///< Holds the alpha values
		double rho;	     ///< ...
	};
	
	/**
		\struct SolutionInfo
		\brief Struct to store information about the solution
	 */
	struct SolutionInfo {
		double obj;          ///< objective 
		double rho;          ///< rho 
		double pUpperBound;  ///< upper bound on the
		double nUpperBound;  ///< upper bound on the 
		double r;	           ///< Value specific for nuSolver
	};
		


	/**
		\class SVM
		\brief Implementation of Support Vector Machines
	 	 
		The SVM class and its helper functions and classes are derived from <b>libsvm v2.81</b>.
	  libsvm was used as a basis to generate a SVM framework that is compatibel with the
	  rest of SiMath. Therefore several of the classes within libsvm have been reorganised
	  and some changes have been made to the API.
	  The main changes and adaptations in this implementation are:
	   - renaming of all classes, methods, functions and variables to be consistent with other naming schemes
	   - the use of the basic Matrix and Vector classes 
	   - the use of std::vector for storage of data instead of arrays and pointers 
	   - using a full Matrix instead of a special sparse %matrix representation
	   - to accomdate the change from arrays of pointers to Matrix for the data, a special indexing scheme had to be established
	   - where possible pointers were substituted by references 
	 
	 More information about libsvm can be found at 
	 http://www.csie.ntu.edu.tw/~cjlin/libsvm/ .
	 
	 The original copyright notice of libsvm reads:
	 \code
	 Copyright (c) 2000-2008 Chih-Chung Chang and Chih-Jen Lin
	 All rights reserved.
	 
	 Redistribution and use in source and binary forms, with or without
	 modification, are permitted provided that the following conditions
	 are met:
	 
	 1. Redistributions of source code must retain the above copyright
	 notice, this list of conditions and the following disclaimer.
	 
	 2. Redistributions in binary form must reproduce the above copyright
	 notice, this list of conditions and the following disclaimer in the
	 documentation and/or other materials provided with the distribution.
	 
	 3. Neither name of copyright holders nor the names of its contributors
	 may be used to endorse or promote products derived from this software
	 without specific prior written permission.
	 
	 
	 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
	 ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
	 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
	 A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR
	 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
	 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
	 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
	 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
	 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
	 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	 \endcode
	 
	 */
	class SVM
	{
		private:
			SVMParameters _param;	         ///< Local copy of the svm parameters

			unsigned int _nbrSV;		       ///< Total number of support vectors (number of rows in _supVecs)
			unsigned int _nbrVar;          ///< Number of variables in the data matrix
			
			Matrix _supVecs;               ///< Matrix to hold the support vectors [_nbrSV x _nbrVar];
			Matrix _svCoef;	               ///< Coefficients for support vectors in decision functions [_nbrClasses-1 x _nbrSV];

			Vector _probA;                 ///< Holds pairwise probability information
			Vector _probB;                 ///< Holds pairwise probability information
		 
			std::vector<DecisionFunction> _decFunctions;   ///< Holds all trained decision functions
			
			// for classification only
			unsigned int _nbrClasses;	                    ///< Number of classes. (is set to 2 in regression and one class svm)
																					
			std::vector<int> _nbrClassSV;                 ///< Vector with the number of support vectors for each class
			std::map<int, unsigned int> _classLabels;     ///< Mapping between original class labels and internal class labels
			std::vector<unsigned int> _dataLabels;        ///< Local storage of data indices ordered by class labels 
			

			/// \name Local helper functions
			//@{
			/**
				\brief Helper function to group classes together on a multi-class problem
			 */
			void _groupClasses(SVMProblem & prob, std::vector<unsigned int> & classStart, std::vector<unsigned int> & classCount);
			
			/**
				\brief Train one specific svm problem
			 */
			void _trainOne(SVMProblem & prob, const SVMParameters & param, double Cp, double Cn, DecisionFunction & f);
			void _binarySVCProbability( SVMProblem & prob, const SVMParameters & param, double Cp, double Cn, double & probA, double & probB);
			Vector _predictValues(const double * xRow);
			double _predictValue(const double * xRow);
			Vector _multiclassProbability(Matrix & r);
			double _svrProbability(SVMProblem & prob, const SVMParameters & param);
			double _predictProbability(const double *xRow);
			double _predict(const double * xRow);

			//@}
			
		public:
			
			///\name Structors
			//@{
			/**
				\brief Default constructor initialises SVM object from predefined parameters
				
				\param param Reference to predefined set of SVM parameters
			 */
			SVM(SVMParameters & param);
			
			/**
				\brief Copy constructor is not yet implemented. Calling this function throws an error.
				\param src Reference to svm object to be copied
				\exception Throws an SVMError
			 */
			SVM(const SVM & src);
			
			/**
				\brief Destructor
			 */
			~SVM();
			//@}

			/// \name Model inspectors
			//@{
			/**
				\brief Get the type of SVM
			 */
			SVMType getSVMType() const { return _param.svmType; };
			
			/**
				\brief Get the number of classes in the problem
			 */
			unsigned int getNbrOfClasses() const { return _nbrClasses; };
			
			/**
				\brief Get the number of support vector of the current problem
			 */
			unsigned int getNbrOfSV() const { return _nbrSV; };

			/**
				\brief Get probability of SV Regression model
	 		 */
			double getSVRProbability();
			
			/**
				\brief Get the class labels
				\return std::vector of integer class labels with as many elements as there are data points in the corresponding SVMProblem
			 */
			std::vector<int> getClassLabels();
			
			/**
				\brief Get the number of support vectors for each class
				*/
			std::vector<int> getClassCount();
			
			/**
				\brief Get the representation of the i-th decision function
				\return A DecisionFunction struct
			 */
			DecisionFunction getDecisionFunction(unsigned int i) const ;
			
			/**
				\brief Get all support vectors
				\return Matrix with coordinates of all support vectors
			 */
			Matrix getSupportVectors() const { return _supVecs; };
			
			/**
				\brief Get the support vector coefficients
			 */
			Matrix getSupportVectorCoefficients() const { return _svCoef; };
			
			/**
				\brief Check if the probability model is set
			 */
			bool checkProbabilityModel();					
			//@}
			
			///\name Training methods
			//@{
			/** 
				\brief Main function to train the SVM
				
				The type and algorithm used to train the SVM is defined using the 
				SVMParameters.
				
				\param data Struct that contains the training data
				*/
			void train(SVMProblem & data);

			/** 
				\brief Cross-validation of SVM
				
				The n-fold crossvalidation procedure, splits the data set in n equal-sized groups
				and then trains n different SVMs.
				
				\param data  Struct that contains the training data
				\param n     Number of groups to split the data in
			 */
			Vector crossValidation(SVMProblem & data, int n);
			//@}
			
			///\name Model parameter adaptors
			//@{
			/**
				\brief Define the number of classes
			 */
			void setNbrOfClasses(unsigned int n);
			
			/**
				\brief Define the mapping of the class labels
			 */
			void setClassLabels(std::vector<int> & cl);
			
			/**
				\brief Define the class labels and number of support vectors per class in one step
				\param cl std::vector of class labels  
				\param cn std::vector of support vector counts per class
			 */
			void setClassLabels(const std::vector<int> & cl, const std::vector<int> cn);

			/**
				\brief Define the number of support vectors per class
			 */
			void setClassCount(std::vector<int> & cn);
			
			/**
				\brief Load a set of stored support vectors and corresponding coefficients
				\param sv Matrix of support vectors
				\param svc Matrix of support vector coefficients
			 */
			void setSupportVectors(const Matrix & sv, const Matrix & svc);

			/**
				\brief Load the i-th decision function
			 */
		void setDecisionFunction(unsigned int i, double rho, std::vector<double> & alpha);
			
			/**
				\brief Load pairwise probabilities
			 */
			void setPairWiseProbabilities(Vector & pA, Vector & pB);
			//@}
			
			
			///\name Prediction methods
			//@{
			/**
				\brief Compute predicted value of single data vector
			 */
			double predict(const Vector & xRow) { return _predict(xRow.getArrayPointer()); };
			
			/**
				\brief Compute predicted values of complete data matrix
			 */
			Vector predict(const Matrix & data);
			
			/**
				\brief Compute predicted probability of single data vector
			 */
			double predictProbability(const Vector & xRow) { return _predictProbability(xRow.getArrayPointer()); };
	
			/**
				\brief Compute predicted probabilities of complete data matrix
			 */
			Vector predictProbability(const Matrix & data);
			//@}

			///\name Helper function to solve specific formulations
			//@{
			static SolutionInfo solveCSVC(SVMProblem & prob, const SVMParameters & param, std::vector<double> & alpha,  double Cp, double Cn);
			static SolutionInfo solveNuSVC(SVMProblem & prob, const SVMParameters & param, std::vector<double> & alpha);
			static SolutionInfo solveOneClass(SVMProblem & prob, const SVMParameters & param, std::vector<double> & alpha);
			static SolutionInfo solveEpsilonSVR(SVMProblem & prob, const SVMParameters & param, std::vector<double> & alpha);
			static SolutionInfo solveNuSVR(SVMProblem & prob, const SVMParameters & param, std::vector<double> & alpha);
			//@}
	};
		
};

#endif __SIMATH_SVM__

