/*
 *  SiMath.h
 *
 *  Created by Gert Thijs on 28/03/2006.
 *  Copyright 2006 Silicos NV, All rights reserved.
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


#ifndef __SIMATH_HEADER_H__
#define __SIMATH_HEADER_H__

// Version constants
#define SIMATH_VERSION                "1"
#define SIMATH_RELEASE                "6"
#define SIMATH_DATE                   "20090320"


// basic functions
#include   "SiMath/Definitions.h"
#include   "SiMath/Matrix.h"      
#include   "SiMath/MatrixOperations.h"
#include   "SiMath/Utilities.h"
#include   "SiMath/Vector.h"

// helper functions
#include   "SiMath/Cache.h"
#include   "SiMath/DiscreteSampling.h"
#include   "SiMath/Distribution.h"
#include   "SiMath/DiscreteDistribution.h"
#include   "SiMath/Kernel.h"
#include   "SiMath/MatrixBinning.h"
#include   "SiMath/PrimeNumber.h" 
#include   "SiMath/Similarity.h"
#include   "SiMath/SVMSolver.h"
#include   "SiMath/SVD.h"

// modeling functions
#include   "SiMath/DivisiveKMeans.h"
#include   "SiMath/InformationGain.h"
#include   "SiMath/KMeans.h"      
#include   "SiMath/LeaderClustering.h"
#include   "SiMath/NaiveBayes.h"  
#include   "SiMath/SOM.h"  
#include   "SiMath/SVM.h"     

#endif  __SIMATH_HEADER_H__

