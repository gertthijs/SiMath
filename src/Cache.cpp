/*
 *  Cache.cpp
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

#include "Cache.h"
#include <iostream>

// constructor
SiMath::Cache::Cache(const unsigned int nBlocks, const unsigned int size)
	: _nbrData(nBlocks), _size(size), _head(NULL), _lruBlock()
{
	_head = new _dataBlock[_nbrData];
	std::cerr << "Cache size: " << _nbrData << "," << _size << std::endl;
	_size /= sizeof(double);
	_size -= _nbrData * sizeof(_dataBlock) / sizeof(double);
	_size = max(_size, 2*_nbrData);	// cache must be large enough for at least two data blocks
	_lruBlock.next = _lruBlock.prev = &_lruBlock;
}

// destructor
SiMath::Cache::~Cache()
{
	// clear each individual data block 
	for( _dataBlock * h = _lruBlock.next; h != &_lruBlock; h=h->next)
	{	
		if ( h->data != NULL )
			delete[] h->data;
	}
	delete[] _head;
}

// delete least-recently used block
void 
SiMath::Cache::_deleteLRU(_dataBlock *h)
{
	// delete from current location
	if ( h->prev == NULL || h->next == NULL )
		throw(CacheError("Unable to delete least recently used data block."));
	
	h->prev->next = h->next;
	h->next->prev = h->prev;
}

// insert new data block
void 
SiMath::Cache::_insertNewBlock(_dataBlock *h)
{
	// insert to last position
	h->next = &_lruBlock;
	h->prev = _lruBlock.prev;
	h->prev->next = h;
	h->next->prev = h;
}


// get the data itself
int 
SiMath::Cache::getData(const unsigned int index, double **data, unsigned int len)
{
	if ( index < 0 || index >= _nbrData )
		throw(SiMath::CacheError("Index outside range of cache"));
	
	_dataBlock *h = &_head[index];
	if( h->length ) _deleteLRU(h);
	
	
	int more = len - h->length;
	if(more > 0)
	{
		// free old space
		while(_size < more)
		{
			_dataBlock *old = _lruBlock.next;
			_deleteLRU(old);
			// free(old->data);
			delete[] old->data;
			_size += old->length;
			old->data = NULL;
			old->length = 0;
		}
		
		// allocate new space
		// h->data = (double *)realloc(h->data,sizeof(double)*len);
		if ( h->data != NULL )
			delete[] h->data;
		h->data = new double[len];
		
		if ( h->data == NULL )
			throw(CacheError("Unable to allocate memory in cache"));
		
		_size -= more;
		std::swap(h->length,len);  
	}
	
	_insertNewBlock(h);
	// store the pointer to the local data
	*data = h->data;
	
	return len;
}



// swap to data blocks in the cache
void 
SiMath::Cache::swapIndex(unsigned int i, unsigned int j)
{
	if( i==j ) return;
	
	if(_head[i].length) _deleteLRU(&_head[i]);
	if(_head[j].length) _deleteLRU(&_head[j]);
	std::swap(_head[i].data,_head[j].data);
	std::swap(_head[i].length,_head[j].length);
	if(_head[i].length) _insertNewBlock(&_head[i]);
	if(_head[j].length) _insertNewBlock(&_head[j]);
	
	if(i>j) std::swap(i,j);
	
	for(_dataBlock *h = _lruBlock.next; h!=&_lruBlock; h=h->next)
	{
		if(h->length > i)
		{
			if(h->length > j)
			{	
				std::swap(h->data[i],h->data[j]); 
			}
			else
			{
				// give up
				_deleteLRU(h);
				if ( h->data != NULL )
					delete[] h->data;
				_size += h->length;
				h->data = NULL;
				h->length = 0;
			}
		}
	}
}

