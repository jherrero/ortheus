/*
 * This is adapted from GNU code as stated below
 */
//
// nazghul - an old-school RPG engine
// Copyright (C) 2002, 2003 Gordon McNutt
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, write to the Free Foundation, Inc., 59 Temple Place,
// Suite 330, Boston, MA 02111-1307 USA
//
// Gordon McNutt  
// gmcnutt@users.sourceforge.net
//    
#include "heapC.h"
#include "fastCMaths.h"  
#include "commonC.h"
  
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h> 
#include <limits.h>

struct heap *heap_create(UNSIGNED_INT_32 max_entries)
{
    struct heap *heap;
	UNSIGNED_INT_32 i;
    
    heap = (struct heap *) mallocLocal(sizeof(struct heap));
	if (!heap)
		return 0;
	memset(heap, 0, sizeof(struct heap));
	heap->max_entries = max_entries;
	heap->num_entries = 0;
	heap->entries = (LONG_64 *) mallocLocal(sizeof(LONG_64) * max_entries);
	if (!heap->entries) {
		free(heap);
		assert(0);
		return 0;
	}
	//memset(heap->entries, LONG_64_MIN, sizeof(LONG_64) * max_entries);
    for(i=0; i<max_entries; i++) {
        heap->entries[i] = LONG_64_MIN;
    }
	return heap;
}

void heap_destroy(struct heap *heap)
{   
    if (heap->entries)
		free(heap->entries);
	free(heap);
}

void heapify(struct heap *heap, LONG_64 i)
{
	UNSIGNED_LONG_64 left;
	UNSIGNED_LONG_64 right;
	LONG_64 largest;

	left = 2 * i;
	right = left + 1;

	if (left < heap->num_entries &&
	    heap->entries[left] > heap->entries[i])
		largest = left;
	else
		largest = i;

	if (right < heap->num_entries &&
	    heap->entries[right] > heap->entries[largest])
		largest = right;

	if (largest != i) {
		LONG_64 tmp = heap->entries[i];
		heap->entries[i] = heap->entries[largest];
		heap->entries[largest] = tmp;
		heapify(heap, largest);
	}
}

INT_32 heap_expand(struct heap *heap)
{
	LONG_64 *tmp;
	tmp =
	    (LONG_64 *) realloc(heap->entries,
			     heap->max_entries * 2 * sizeof(LONG_64));
	if (!tmp) {
		assert(0);
		return -1;
	}
	heap->entries = tmp;
	heap->max_entries *= 2;
	return 0;
}

INT_32 heap_insert(struct heap *heap, LONG_64 entry)
{
	INT_32 i;

	/* Expand the heap if necessary */
	if (heap->num_entries == heap->max_entries) {
		if (heap_expand(heap) < 0) {
			assert(0);
			return -1;
		}
	}

	/* Put the new entry at the bottom of the heap */
	i = heap->num_entries++;

	/* Percolate the new entry up to where it belongs */
	while (i > 0 && heap->entries[i / 2] < entry) {
		heap->entries[i] = heap->entries[i / 2];
		i /= 2;
	}

	heap->entries[i] = entry; 

	return 0;
}

LONG_64 heap_extract(struct heap *heap)
{
	LONG_64 max;

	if (!heap->num_entries)
		return LONG_64_MIN;

	max = heap->entries[0];
	heap->entries[0] = heap->entries[heap->num_entries - 1];
	heap->num_entries--;
	heapify(heap, 0);

	return max;
}

LONG_64 heap_peek(struct heap *heap) {
	if (!heap->num_entries)
		return LONG_64_MIN;
		
	return heap->entries[0];
}

INT_32 heap_empty(struct heap *heap) { 
	return !(heap->num_entries);
}


void heap_clean(struct heap *heap) 
{
	UNSIGNED_INT_32 i;
    
    //memset(heap->entries, LONG_64_MIN, heap->num_entries * sizeof(LONG_64));
	for(i=0; i<heap->num_entries; i++) {
        heap->entries[i] = LONG_64_MIN;
    }
    heap->num_entries = 0;
}

