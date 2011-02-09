/**
 *
 *
 *                              PIntron
 *
 * A novel pipeline for computational gene-structure prediction based on
 * spliced alignment of expressed sequences (ESTs and mRNAs).
 *
 * Copyright (C) 2010  Yuri Pirola
 *
 * Distributed under the terms of the GNU Affero General Public License (AGPL)
 *
 *
 * This file is part of PIntron.
 *
 * PIntron is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PIntron is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with PIntron.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
#include "ext_array.h"

#include "util.h"
#include <stdlib.h>
#include <string.h>

#define INITIAL_SIZE 10

struct _ext_array {
  item* data;

  int nel;
  int capacity;
};

struct _ext_arrayit {
  item* cur_data;
  int el_restanti;
  int el_visti;
};



pext_array
EA_create(void) {
  pext_array ris= PALLOC(struct _ext_array);
  ris->data= NPALLOC(item, INITIAL_SIZE);
  ris->capacity= INITIAL_SIZE;
  ris->nel= 0;
  return ris;
}


void
EA_destroy(pext_array arr, delete_function del) {
  my_assert(arr!= NULL);
  my_assert(del!=NULL);
  EA_clear(arr, del);
  pfree(arr->data);
  pfree(arr);
}


static void
move(item* newdata, pext_array arr) {
  memcpy(newdata, arr->data, arr->nel*sizeof(item));
  pfree(arr->data);
  arr->data= newdata;
}

/**
 * Inserisce un elemento espandendo, eventualmente, l'array.
 **/
void
EA_insert(pext_array arr, item it) {
  my_assert(arr!=NULL);
  my_assert(arr->nel < arr->capacity);
  arr->data[arr->nel]= it;
  ++arr->nel;
// Controllo il fattore di carico e raddoppio la dimensione se
// supero la soglia dello 0.5
  if ((arr->capacity >> 1)<arr->nel) {
// Espansione!
	 item* newdata= NPALLOC(item, (arr->capacity << 1));
	 move(newdata, arr);
	 arr->capacity= arr->capacity << 1;
  }
  my_assert(arr->nel < arr->capacity);
}

void
EA_clear(pext_array arr, delete_function del) {
  my_assert(arr!= NULL);
  my_assert(del!= NULL);
  for (int i= 0; i<arr->nel; ++i) {
	 del(arr->data[i]);
  }
  arr->nel= 0;
}

bool
EA_is_empty(pext_array arr) {
  my_assert(arr!= NULL);
  return arr->nel==0;
}

item
EA_get(pext_array arr, int pos) {
  my_assert(arr!=NULL);
  my_assert(0<=pos && pos<arr->nel);
  return (arr->data[pos]);
}

void
EA_set(pext_array arr, int pos, item el) {
  my_assert(arr!=NULL);
  my_assert(0<=pos && pos<arr->nel);
  my_assert(el != NULL);
  arr->data[pos]= el;
}

size_t
EA_size(pext_array arr) {
  my_assert(arr!=NULL);
  return (arr->nel);
}


void
EA_sort(pext_array arr, comparator cmp) {
  my_assert(arr != NULL);
  my_assert(cmp != NULL);
  if (arr->nel > 0) {
	 qsort(arr->data, arr->nel, sizeof(item), cmp);
  }
}

int
EA_binary_search(pext_array arr, item val, comparator cmp) {
  my_assert(arr != NULL);
  my_assert(cmp != NULL);
  int l= 0;
  int r= arr->nel;
  int last_ris= -1;
  int ris= -1;
  while (last_ris != 0 && l < r) {
	 int mid= (l+r)>>1;
	 last_ris= cmp( arr->data[mid], val );
	 if (last_ris < 0) {
		r= mid;
	 } else if (last_ris > 0) {
		l= mid+1;
	 } else {
		ris= mid;
	 }
  }
  my_assert(ris>=-1 || ris <arr->nel);
  return ris;
}




pext_arrayit
EA_begin(pext_array arr) {
  my_assert(arr!= NULL);
  pext_arrayit it= PALLOC(struct _ext_arrayit);
  it->cur_data= arr->data;
  it->el_restanti= arr->nel;
  it->el_visti= 0;
  return it;
}

pext_arrayit
ext_arrayit_clone(pext_arrayit eait) {
  my_assert(eait!= NULL);
  pext_arrayit it= PALLOC(struct _ext_arrayit);
  it->cur_data= eait->cur_data;
  it->el_restanti= eait->el_restanti;
  it->el_visti= eait->el_visti;
  return it;
}

void
ext_arrayit_copy(pext_arrayit target, pext_arrayit source) {
  my_assert(target!= NULL);
  my_assert(source!= NULL);
  target->cur_data= source->cur_data;
  target->el_restanti= source->el_restanti;
  target->el_visti= source->el_visti;
}

void
ext_arrayit_destroy(pext_arrayit it) {
  my_assert(it!= NULL);
  pfree(it);
}

bool
ext_arrayit_has_next(pext_arrayit it) {
  my_assert(it!= NULL);
  return it->el_restanti>0;
}

item
ext_arrayit_next(pext_arrayit it) {
  my_assert(it!= NULL);
  if (it->el_restanti<=0) {
	 return NULL;
  }
  item ris= *(it->cur_data);
  --(it->el_restanti);
  ++(it->el_visti);
  if (it->el_restanti>0) {
	 ++(it->cur_data);
  }
  return ris;
}

item
ext_arrayit_prev(pext_arrayit it) {
  my_assert(it!= NULL);
  if (it->el_visti==0) {
	 return NULL;
  }
  ++(it->el_restanti);
  if (it->el_restanti!=1)
	 --(it->cur_data);
  --(it->el_visti);
  item ris= *(it->cur_data);
  return ris;
}

#undef INITIAL_SIZE
