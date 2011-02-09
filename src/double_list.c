/**
 *
 *
 *                              PIntron
 *
 * A novel pipeline for computational gene-structure prediction based on
 * spliced alignment of expressed sequences (ESTs and mRNAs).
 *
 * Copyright (C) 2010  Yuri Pirola, Raffaella Rizzi
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
#include "double_list.h"
#include "util.h"
#include <stdlib.h>
#include "log.h"

typedef struct _doublenode* _pdoublenode;

struct _doublenode {
	_pdoublenode next;
	_pdoublenode prev;
	DTYPE element;
};

struct _doublelist {
  _pdoublenode sentinel;
  size_t size;
};


struct _doublelistit {
  _pdoublenode next;
  _pdoublenode prev;
  _pdoublenode sentinel;
};


pdoublelist doublelist_create(void) {
  TRACE("Creazione lista");
  pdoublelist ris= PALLOC(struct _doublelist);
  ris->size= 0;
  ris->sentinel= PALLOC(struct _doublenode);
  ris->sentinel->prev= ris->sentinel->next= ris->sentinel;
  ris->sentinel->element= 0;
  return ris;
}


void doublelist_clear(pdoublelist l) {
  my_assert(l!=NULL);
  TRACE("Distruzione di una lista con %zd elementi", doublelist_size(l));
  _pdoublenode n= l->sentinel;
  while (n->next != l->sentinel) {
	 _pdoublenode tmp= n->next->next;
	 pfree(n->next);
	 n->next= tmp;
  }
  l->sentinel->next= l->sentinel->prev= l->sentinel;
  l->size= 0;
}

void doublelist_destroy(pdoublelist l) {
  my_assert(l!=NULL);
  TRACE("Distruzione di una lista con %zd elementi", doublelist_size(l));
  _pdoublenode n= l->sentinel;
  while (n->next != l->sentinel) {
	 _pdoublenode tmp= n->next->next;
	 pfree(n->next);
	 n->next= tmp;
  }
  pfree(l->sentinel);
  pfree(l);
}

void doublelist_add_to_head(pdoublelist l, DTYPE p) {
  my_assert(l!=NULL);
  _pdoublenode n= PALLOC(struct _doublenode);
  n->element= p;
  n->next= l->sentinel->next;
  n->prev= l->sentinel;
  n->next->prev= n;
  l->sentinel->next= n;
  l->size= l->size+1;
}

void doublelist_add_to_tail(pdoublelist l, DTYPE p) {
  my_assert(l!=NULL);
  _pdoublenode n= PALLOC(struct _doublenode);
  n->element= p;
  n->prev= l->sentinel->prev;
  n->next= l->sentinel;
  n->prev->next= n;
  l->sentinel->prev= n;
  l->size= l->size+1;
}

DTYPE doublelist_remove_from_head(pdoublelist l) {
  my_assert(l!=NULL);
  // Non deve essere vuota
  my_assert(l->sentinel->next != l->sentinel);
  _pdoublenode n= l->sentinel->next;
  DTYPE p= n->element;
  l->sentinel->next= n->next;
  n->next->prev= l->sentinel;
  pfree(n);
  l->size= l->size-1;
  return p;
}

DTYPE doublelist_remove_from_tail(pdoublelist l) {
  my_assert(l!=NULL);
  // Non deve essere vuota
  my_assert(l->sentinel->prev != l->sentinel);
  _pdoublenode n= l->sentinel->prev;
  DTYPE p= n->element;
  l->sentinel->prev= n->prev;
  n->prev->next= l->sentinel;
  pfree(n);
  l->size= l->size-1;
  return p;
}

DTYPE doublelist_head(pdoublelist l) {
  my_assert(l!=NULL);
  return l->sentinel->next->element;
}

DTYPE doublelist_tail(pdoublelist l) {
  my_assert(l!=NULL);
  return l->sentinel->prev->element;
}

bool doublelist_is_empty(pdoublelist l) {
  my_assert(l!=NULL);
  my_assert( (l->size!=0) || (l->sentinel->next == l->sentinel));
  return (l->sentinel->next == l->sentinel);
}

size_t doublelist_size(pdoublelist l) {
  my_assert(l!=NULL);
  return l->size;
}

void doublelist_sort(pdoublelist l, comparator cmp){
  my_assert(l != NULL);
  my_assert(cmp != NULL);
  double size = doublelist_size(l);
  DTYPE* base= NPALLOC(DTYPE, size);
  int i = 0;
  pdoublelistit it = doublelist_first(l);
  while(doublelistit_has_next(it)){
	 base[i] = doublelistit_next(it);
	 ++i;
  }
  doublelistit_destroy(it);

  qsort(base, size, sizeof(DTYPE), cmp);

  _pdoublenode tmp = l->sentinel->next;
  i= 0;
  while(tmp != l->sentinel){
	 tmp->element = base[i];
	 ++i;
	 tmp = tmp->next;
  }
  pfree(base);
}

pdoublelist doublelist_merge_new(pdoublelist l1, pdoublelist l2){
  my_assert(l1!=NULL);
  my_assert(l2!=NULL);
  pdoublelistit it;
  pdoublelist mergedList = doublelist_create();
  it = doublelist_first(l1);
  while(doublelistit_has_next(it)) {
	 doublelist_add_to_tail(mergedList, doublelistit_next(it));
  }
  doublelistit_destroy(it);
  it = doublelist_first(l2);
  while(doublelistit_has_next(it)) {
	 doublelist_add_to_tail(mergedList, doublelistit_next(it));
  }
  doublelistit_destroy(it);
  return mergedList;
}

void doublelist_merge(pdoublelist l1, pdoublelist l2){
     my_assert(l1!=NULL);
     my_assert(l2!=NULL);
     l2->sentinel->next->prev = l1->sentinel->prev;
     l1->sentinel->prev->next = l2->sentinel->next;
     l1->sentinel->prev = l2->sentinel->prev;
     l2->sentinel->prev->next = l1->sentinel;

     pfree(l2->sentinel);
     pfree(l2);
     }

pdoublelist doublelist_copy(pdoublelist l){
  my_assert(l != NULL);
  pdoublelist res = doublelist_create();
  pdoublelistit it = doublelist_first(l);

  while (doublelistit_has_next(it)) {
	 doublelist_add_to_tail(res, doublelistit_next(it));
  }
  doublelistit_destroy(it);
  return res;

}

void doublelist_remove_at_iterator(pdoublelistit it){
  my_assert(it != NULL);

  it->next->prev = it->prev->prev;
  it->prev->prev->next = it->next;

  pfree(it->prev);
}

void doublelist_difference(pdoublelist l1, pdoublelist l2){
  pdoublelistit it1 = doublelist_first(l1);
  pdoublelistit it2 = doublelist_first(l2);

  while( doublelistit_has_next(it1) &&
			doublelistit_has_next(it2) ){
	 double ris = it1->next->element - it2->next->element;
	 if( ris == 0 ){
		doublelistit_next(it1);
		doublelist_remove_at_iterator(it1);
	 } else if (ris > 0) {
		doublelistit_next(it2);
	 } else {
		doublelistit_next(it1);
	 } //endelse
  }//endwhile
  doublelistit_destroy(it1);
  doublelistit_destroy(it2);
}//endlistdiff


/*
 * doublelist iterator definitions
 */

pdoublelistit doublelist_first(pdoublelist l) {
  my_assert(l!=NULL);
  pdoublelistit li= PALLOC(struct _doublelistit);
  li->next= l->sentinel->next;
  li->prev= l->sentinel;
  li->sentinel= l->sentinel;
  return li;
}

pdoublelistit doublelist_last(pdoublelist l) {
  my_assert(l!=NULL);
  pdoublelistit li= PALLOC(struct _doublelistit);
  li->prev= l->sentinel->prev;
  li->next= l->sentinel;
  li->sentinel= l->sentinel;
  return li;
}

void doublelistit_destroy(pdoublelistit li) {
  my_assert(li!=NULL);
  pfree(li);
}

bool doublelistit_has_next(pdoublelistit li) {
  my_assert(li!=NULL);
  return li->next != li->sentinel;
}

DTYPE doublelistit_next(pdoublelistit li) {
  my_assert(li!=NULL);
  my_assert(doublelistit_has_next(li));
  DTYPE p= li->next->element;
  li->prev= li->next;
  li->next= li->next->next;
  return p;
}

bool doublelistit_has_prev(pdoublelistit li) {
  my_assert(li!=NULL);
  return li->prev != li->sentinel;
}

DTYPE doublelistit_prev(pdoublelistit li) {
  my_assert(li!=NULL);
  my_assert(doublelistit_has_prev(li));
  DTYPE p= li->prev->element;
  li->next= li->prev;
  li->prev= li->prev->prev;
  return p;
}
