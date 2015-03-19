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
  
#ifndef HEAPC_H_
#define HEAPC_H_
  
#include "fastCMaths.h"

//#include "macros.h"
  
//BEGIN_DECL

//#define heap_empty(h) (!(h)->num_entries)

struct heap {
        UNSIGNED_INT_32 max_entries;
        UNSIGNED_INT_32 num_entries;
        LONG_64 *entries; 
};

extern INT_32 heap_empty(struct heap *heap);
extern struct heap *heap_create(UNSIGNED_INT_32 max_entries);
extern void heap_destroy(struct heap *heap);
extern void heapify(struct heap *heap, LONG_64 i);
extern INT_32 heap_expand(struct heap *heap);
extern INT_32 heap_insert(struct heap *heap, LONG_64 entry);
extern LONG_64 heap_extract(struct heap *heap);
extern LONG_64 heap_peek(struct heap *heap); 
extern void heap_clean(struct heap *heap);

#define HEADROOM 100

//END_DECL
 
#endif
