/*
 *  Cache.h
 *
 *  Created by Gert Thijs on 04/04/06.
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

#ifndef __SIMATH_CACHE_H__
#define __SIMATH_CACHE_H__

#include "Definitions.h"

namespace SiMath
{
	
	/**
		\class  CacheError
		\brief  Exception class
	 
		This exception class is derived from the std::runtime_error class and is used to 
		indicate that something went wrong during the creation and manipulation of the cache.
	 \ingroup _error_handling_
	 */	
	class SIMATH_EXCEPTION CacheError : public SiMathError
	{
		public:
			CacheError(const std::string & e) : SiMathError("CacheError : " + e ) {};
	};

	/**
		\class Cache
		\brief Implementation of Least Recently Used (LRU) Cache to store data as a pointer to double.
	 
		The LRU cache holds a circular list of data blocks in which a pointer to the data is stored.
		When new data is added to the cache, first it checks if there is room to store the data. If not,
		the least recently used blocks are deleted untill there is enough room to store the new data.
	 */
	class Cache
	{	
		public:
			///\name Structors
			//@{
			/**
				\brief Cache constructor of n blocks stored in maximal given Megabytes.
			 
				\param nBlocks Defines the number of data blocks in the cache
				\param size Defines the size of the cache in Mb.
			 */
			Cache(unsigned int nBlocks, unsigned int size);

			/**
				\brief Copy constructor is not yet implemented and throws an error
				\exception throws a CacheError
			 */
			Cache(const Cache & src) {throw(CacheError("Cache copy constructor not implemented."));};
			
			/**
				\brief Destructor
			 */
			~Cache();
			//@}
			
			
			/**
				\brief Assignment operator is not yet implemented and throws an error
			  \exception throws a CacheError
			 */
			Cache & operator= (const Cache & src) {throw(CacheError("Cache copy constructor not yet implemented."));};

			///\name Cache Access
			//@{
			/** 
				\brief Method to access the data in a data block
				
				This method will try to fetch the data from in the range [0,len) from data block index
				On completion it returns a value p. If p < len, then the positions 
				[p,len) need to be filled in the cache. If p >= len, nothing has to be filled 
				and data can be accessed through the pointer (**data) to the internal array.

				\code
				double * d = NULL;
				int l = 100;
				int p = cache.getData(i, &d, l );
				if ( p < l )
				{
					for ( int i=p; i<l; ++i )
						d[i] = ... // compute data to be stored
				}
				\endcode
				\param index   Index of the data block in the cache
				\param data    Pointer to data array
				\param len     Number of data to extract from the cache

				\return int    Position up to where cache is filled. 
				*/
			int getData(const unsigned int index, double **data, unsigned int len);

			/**
				\brief Function to swap the data blocks i and j in the cache
			 */
			void swapIndex(unsigned int i, unsigned int j);	
			//@}
			
		private:
			unsigned int _nbrData;          ///< total number of possible data blocks in this cache 
			unsigned int _size;             ///< total size of the cache
			
			/**
				\brief A datablock is a circular list which holds a pointer to the data
			 */
			struct _dataBlock
			{
				_dataBlock() : prev(NULL), next(NULL), data(NULL), length(0) {};
				
				_dataBlock * prev;      ///< Pointer to the previous element in the list
				_dataBlock * next;	    ///< Pointer to the next element
				double * data;          ///< Pointer to the data itself
				unsigned int length;	  ///< number of data points cached in this entry
			};
			
			_dataBlock * _head;      ///< Pointer to head of the cache
			_dataBlock _lruBlock;    ///< fixed block with no data, which holds pointer to the least recently used block (prev) and pointer to last used (next) 
			
			
			/**
				\brief Helper function to delete a data block
			 
				\param h Pointer to data block to be deleted
			 */
			void _deleteLRU(_dataBlock *h);
			
			/**
				\brief Helper function to insert a new data block
				
				\param h Pointer to a data block
			 */
			void _insertNewBlock(_dataBlock *h);
	};
	
	
	
};

#endif __SIMATH_CACHE_H__

