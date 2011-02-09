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
#include "int_list.h"
#include "util.h"
#include <stdlib.h>
#include "log.h"

typedef struct _intnode* _pintnode;

struct _intnode {
	_pintnode next;
	_pintnode prev;
	ITYPE element;
};

struct _intlist {
  _pintnode sentinel;
  size_t size;
};


struct _intlistit {
  _pintnode next;
  _pintnode prev;
  _pintnode sentinel;
};


pintlist intlist_create(void) {
  TRACE("Creazione lista");
  pintlist ris= PALLOC(struct _intlist);
  ris->size= 0;
  ris->sentinel= PALLOC(struct _intnode);
  ris->sentinel->prev= ris->sentinel->next= ris->sentinel;
  ris->sentinel->element= 0;
  return ris;
}


void intlist_clear(pintlist l) {
  my_assert(l!=NULL);
  TRACE("Distruzione di una lista con %zd elementi", intlist_size(l));
  _pintnode n= l->sentinel;
  while (n->next != l->sentinel) {
	 _pintnode tmp= n->next->next;
	 pfree(n->next);
	 n->next= tmp;
  }
  l->sentinel->next= l->sentinel->prev= l->sentinel;
  l->size= 0;
}

void intlist_destroy(pintlist l) {
  my_assert(l!=NULL);
  TRACE("Distruzione di una lista con %zd elementi", intlist_size(l));
  _pintnode n= l->sentinel;
  while (n->next != l->sentinel) {
	 _pintnode tmp= n->next->next;
	 pfree(n->next);
	 n->next= tmp;
  }
  pfree(l->sentinel);
  pfree(l);
}

void intlist_add_to_head(pintlist l, ITYPE p) {
  my_assert(l!=NULL);
  _pintnode n= PALLOC(struct _intnode);
  n->element= p;
  n->next= l->sentinel->next;
  n->prev= l->sentinel;
  n->next->prev= n;
  l->sentinel->next= n;
  l->size= l->size+1;
}

void intlist_add_to_tail(pintlist l, ITYPE p) {
  my_assert(l!=NULL);
  _pintnode n= PALLOC(struct _intnode);
  n->element= p;
  n->prev= l->sentinel->prev;
  n->next= l->sentinel;
  n->prev->next= n;
  l->sentinel->prev= n;
  l->size= l->size+1;
}

ITYPE intlist_remove_from_head(pintlist l) {
  my_assert(l!=NULL);
  // Non deve essere vuota
  my_assert(l->sentinel->next != l->sentinel);
  _pintnode n= l->sentinel->next;
  ITYPE p= n->element;
  l->sentinel->next= n->next;
  n->next->prev= l->sentinel;
  pfree(n);
  l->size= l->size-1;
  return p;
}

ITYPE intlist_remove_from_tail(pintlist l) {
  my_assert(l!=NULL);
  // Non deve essere vuota
  my_assert(l->sentinel->prev != l->sentinel);
  _pintnode n= l->sentinel->prev;
  ITYPE p= n->element;
  l->sentinel->prev= n->prev;
  n->prev->next= l->sentinel;
  pfree(n);
  l->size= l->size-1;
  return p;
}

ITYPE intlist_head(pintlist l) {
  my_assert(l!=NULL);
  return l->sentinel->next->element;
}

ITYPE intlist_tail(pintlist l) {
  my_assert(l!=NULL);
  return l->sentinel->prev->element;
}

bool intlist_is_empty(pintlist l) {
  my_assert(l!=NULL);
  my_assert( (l->size!=0) || (l->sentinel->next == l->sentinel));
  return (l->sentinel->next == l->sentinel);
}

size_t intlist_size(pintlist l) {
  my_assert(l!=NULL);
  return l->size;
}

void intlist_sort(pintlist l, comparator cmp){
  my_assert(l != NULL);
  my_assert(cmp != NULL);
  int size = intlist_size(l);
  ITYPE* base= NPALLOC(ITYPE, size);
  int i = 0;
  pintlistit it = intlist_first(l);
  while(intlistit_has_next(it)){
	 base[i] = intlistit_next(it);
	 ++i;
  }
  intlistit_destroy(it);

  qsort(base, size, sizeof(ITYPE), cmp);

  _pintnode tmp = l->sentinel->next;
  i= 0;
  while(tmp != l->sentinel){
	 tmp->element = base[i];
	 ++i;
	 tmp = tmp->next;
  }
  pfree(base);
}

pintlist intlist_merge_new(pintlist l1, pintlist l2){
  my_assert(l1!=NULL);
  my_assert(l2!=NULL);
  pintlistit it;
  pintlist mergedList = intlist_create();
  it = intlist_first(l1);
  while(intlistit_has_next(it)) {
	 intlist_add_to_tail(mergedList, intlistit_next(it));
  }
  intlistit_destroy(it);
  it = intlist_first(l2);
  while(intlistit_has_next(it)) {
	 intlist_add_to_tail(mergedList, intlistit_next(it));
  }
  intlistit_destroy(it);
  return mergedList;
}

void intlist_merge(pintlist l1, pintlist l2){
     my_assert(l1!=NULL);
     my_assert(l2!=NULL);
     l2->sentinel->next->prev = l1->sentinel->prev;
     l1->sentinel->prev->next = l2->sentinel->next;
     l1->sentinel->prev = l2->sentinel->prev;
     l2->sentinel->prev->next = l1->sentinel;

     pfree(l2->sentinel);
     pfree(l2);
     }

pintlist intlist_copy(pintlist l){
  my_assert(l != NULL);
  pintlist res = intlist_create();
  pintlistit it = intlist_first(l);

  while (intlistit_has_next(it)) {
	 intlist_add_to_tail(res, intlistit_next(it));
  }
  intlistit_destroy(it);
  return res;

}

void intlist_remove_at_iterator(pintlistit it){
  my_assert(it != NULL);

  it->next->prev = it->prev->prev;
  it->prev->prev->next = it->next;

  pfree(it->prev);
}

void intlist_difference(pintlist l1, pintlist l2){
  pintlistit it1 = intlist_first(l1);
  pintlistit it2 = intlist_first(l2);

  while( intlistit_has_next(it1) &&
			intlistit_has_next(it2) ){
	 int ris = it1->next->element - it2->next->element;
	 if( ris == 0 ){
		intlistit_next(it1);
		intlist_remove_at_iterator(it1);
	 } else if (ris > 0) {
		intlistit_next(it2);
	 } else {
		intlistit_next(it1);
	 } //endelse
  }//endwhile
  intlistit_destroy(it1);
  intlistit_destroy(it2);
}//endlistdiff


/*
 * intlist iterator definitions
 */

pintlistit intlist_first(pintlist l) {
  my_assert(l!=NULL);
  pintlistit li= PALLOC(struct _intlistit);
  li->next= l->sentinel->next;
  li->prev= l->sentinel;
  li->sentinel= l->sentinel;
  return li;
}

pintlistit intlist_last(pintlist l) {
  my_assert(l!=NULL);
  pintlistit li= PALLOC(struct _intlistit);
  li->prev= l->sentinel->prev;
  li->next= l->sentinel;
  li->sentinel= l->sentinel;
  return li;
}

void intlistit_destroy(pintlistit li) {
  my_assert(li!=NULL);
  pfree(li);
}

bool intlistit_has_next(pintlistit li) {
  my_assert(li!=NULL);
  return li->next != li->sentinel;
}

ITYPE intlistit_next(pintlistit li) {
  my_assert(li!=NULL);
  my_assert(intlistit_has_next(li));
  ITYPE p= li->next->element;
  li->prev= li->next;
  li->next= li->next->next;
  return p;
}

bool intlistit_has_prev(pintlistit li) {
  my_assert(li!=NULL);
  return li->prev != li->sentinel;
}

ITYPE intlistit_prev(pintlistit li) {
  my_assert(li!=NULL);
  my_assert(intlistit_has_prev(li));
  ITYPE p= li->prev->element;
  li->next= li->prev;
  li->prev= li->prev->prev;
  return p;
}

