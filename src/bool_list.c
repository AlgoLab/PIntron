/**
 *
 *
 *                              PIntron
 *
 * A novel pipeline for computational gene-structure prediction based on
 * spliced alignment of expressed sequences (ESTs and mRNAs).
 *
 * Copyright (C) 2010  Raffaella Rizzi
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
#include "bool_list.h"
#include "util.h"
#include <stdlib.h>
#include "log.h"

typedef struct _boolnode* _pboolnode;

struct _boolnode {
	_pboolnode next;
	_pboolnode prev;
	BTYPE element;
};

struct _boollist {
  _pboolnode sentinel;
  size_t size;
};


struct _boollistit {
  _pboolnode next;
  _pboolnode prev;
  _pboolnode sentinel;
};


pboollist boollist_create(void) {
  TRACE("Creazione lista");
  pboollist ris= PALLOC(struct _boollist);
  ris->size= 0;
  ris->sentinel= PALLOC(struct _boolnode);
  ris->sentinel->prev= ris->sentinel->next= ris->sentinel;
  ris->sentinel->element= 0;
  return ris;
}


void boollist_clear(pboollist l) {
  my_assert(l!=NULL);
  TRACE("Distruzione di una lista con %zd elementi", boollist_size(l));
  _pboolnode n= l->sentinel;
  while (n->next != l->sentinel) {
	 _pboolnode tmp= n->next->next;
	 pfree(n->next);
	 n->next= tmp;
  }
  l->sentinel->next= l->sentinel->prev= l->sentinel;
  l->size= 0;
}

void boollist_destroy(pboollist l) {
  my_assert(l!=NULL);
  TRACE("Distruzione di una lista con %zd elementi", boollist_size(l));
  _pboolnode n= l->sentinel;
  while (n->next != l->sentinel) {
	 _pboolnode tmp= n->next->next;
	 pfree(n->next);
	 n->next= tmp;
  }
  pfree(l->sentinel);
  pfree(l);
}

void boollist_add_to_head(pboollist l, BTYPE p) {
  my_assert(l!=NULL);
  _pboolnode n= PALLOC(struct _boolnode);
  n->element= p;
  n->next= l->sentinel->next;
  n->prev= l->sentinel;
  n->next->prev= n;
  l->sentinel->next= n;
  l->size= l->size+1;
}

void boollist_add_to_tail(pboollist l, BTYPE p) {
  my_assert(l!=NULL);
  _pboolnode n= PALLOC(struct _boolnode);
  n->element= p;
  n->prev= l->sentinel->prev;
  n->next= l->sentinel;
  n->prev->next= n;
  l->sentinel->prev= n;
  l->size= l->size+1;
}

BTYPE boollist_remove_from_head(pboollist l) {
  my_assert(l!=NULL);
  // Non deve essere vuota
  my_assert(l->sentinel->next != l->sentinel);
  _pboolnode n= l->sentinel->next;
  BTYPE p= n->element;
  l->sentinel->next= n->next;
  n->next->prev= l->sentinel;
  pfree(n);
  l->size= l->size-1;
  return p;
}

BTYPE boollist_remove_from_tail(pboollist l) {
  my_assert(l!=NULL);
  // Non deve essere vuota
  my_assert(l->sentinel->prev != l->sentinel);
  _pboolnode n= l->sentinel->prev;
  BTYPE p= n->element;
  l->sentinel->prev= n->prev;
  n->prev->next= l->sentinel;
  pfree(n);
  l->size= l->size-1;
  return p;
}

BTYPE boollist_head(pboollist l) {
  my_assert(l!=NULL);
  return l->sentinel->next->element;
}

BTYPE boollist_tail(pboollist l) {
  my_assert(l!=NULL);
  return l->sentinel->prev->element;
}

bool boollist_is_empty(pboollist l) {
  my_assert(l!=NULL);
  my_assert( (l->size!=0) || (l->sentinel->next == l->sentinel));
  return (l->sentinel->next == l->sentinel);
}

size_t boollist_size(pboollist l) {
  my_assert(l!=NULL);
  return l->size;
}

void boollist_sort(pboollist l, comparator cmp){
  my_assert(l != NULL);
  my_assert(cmp != NULL);
  bool size = boollist_size(l);
  BTYPE* base= NPALLOC(BTYPE, size);
  bool i = 0;
  pboollistit it = boollist_first(l);
  while(boollistit_has_next(it)){
	 base[i] = boollistit_next(it);
	 ++i;
  }
  boollistit_destroy(it);

  qsort(base, size, sizeof(BTYPE), cmp);

  _pboolnode tmp = l->sentinel->next;
  i= 0;
  while(tmp != l->sentinel){
	 tmp->element = base[i];
	 ++i;
	 tmp = tmp->next;
  }
  pfree(base);
}

pboollist boollist_merge_new(pboollist l1, pboollist l2){
  my_assert(l1!=NULL);
  my_assert(l2!=NULL);
  pboollistit it;
  pboollist mergedList = boollist_create();
  it = boollist_first(l1);
  while(boollistit_has_next(it)) {
	 boollist_add_to_tail(mergedList, boollistit_next(it));
  }
  boollistit_destroy(it);
  it = boollist_first(l2);
  while(boollistit_has_next(it)) {
	 boollist_add_to_tail(mergedList, boollistit_next(it));
  }
  boollistit_destroy(it);
  return mergedList;
}

void boollist_merge(pboollist l1, pboollist l2){
     my_assert(l1!=NULL);
     my_assert(l2!=NULL);
     l2->sentinel->next->prev = l1->sentinel->prev;
     l1->sentinel->prev->next = l2->sentinel->next;
     l1->sentinel->prev = l2->sentinel->prev;
     l2->sentinel->prev->next = l1->sentinel;

     pfree(l2->sentinel);
     pfree(l2);
     }

pboollist boollist_copy(pboollist l){
  my_assert(l != NULL);
  pboollist res = boollist_create();
  pboollistit it = boollist_first(l);

  while (boollistit_has_next(it)) {
	 boollist_add_to_tail(res, boollistit_next(it));
  }
  boollistit_destroy(it);
  return res;

}

void boollist_remove_at_iterator(pboollistit it){
  my_assert(it != NULL);

  it->next->prev = it->prev->prev;
  it->prev->prev->next = it->next;

  pfree(it->prev);
}

void boollist_difference(pboollist l1, pboollist l2){
  pboollistit it1 = boollist_first(l1);
  pboollistit it2 = boollist_first(l2);

  while( boollistit_has_next(it1) &&
			boollistit_has_next(it2) ){
	 bool ris = it1->next->element - it2->next->element;
	 if( ris == 0 ){
		boollistit_next(it1);
		boollist_remove_at_iterator(it1);
	 } else if (ris > 0) {
		boollistit_next(it2);
	 } else {
		boollistit_next(it1);
	 } //endelse
  }//endwhile
  boollistit_destroy(it1);
  boollistit_destroy(it2);
}//endlistdiff


/*
 * boollist iterator definitions
 */

pboollistit boollist_first(pboollist l) {
  my_assert(l!=NULL);
  pboollistit li= PALLOC(struct _boollistit);
  li->next= l->sentinel->next;
  li->prev= l->sentinel;
  li->sentinel= l->sentinel;
  return li;
}

pboollistit boollist_last(pboollist l) {
  my_assert(l!=NULL);
  pboollistit li= PALLOC(struct _boollistit);
  li->prev= l->sentinel->prev;
  li->next= l->sentinel;
  li->sentinel= l->sentinel;
  return li;
}

void boollistit_destroy(pboollistit li) {
  my_assert(li!=NULL);
  pfree(li);
}

bool boollistit_has_next(pboollistit li) {
  my_assert(li!=NULL);
  return li->next != li->sentinel;
}

BTYPE boollistit_next(pboollistit li) {
  my_assert(li!=NULL);
  my_assert(boollistit_has_next(li));
  BTYPE p= li->next->element;
  li->prev= li->next;
  li->next= li->next->next;
  return p;
}

bool boollistit_has_prev(pboollistit li) {
  my_assert(li!=NULL);
  return li->prev != li->sentinel;
}

BTYPE boollistit_prev(pboollistit li) {
  my_assert(li!=NULL);
  my_assert(boollistit_has_prev(li));
  BTYPE p= li->prev->element;
  li->next= li->prev;
  li->prev= li->prev->prev;
  return p;
}

