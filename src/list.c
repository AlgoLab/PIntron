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
#include "list.h"
#include "log.h"
#include "util.h"
#include <stdlib.h>

plist list_create(void) {
  FINETRACE("New list creation");
  plist ris= PALLOC(struct _list);
  ris->size= 0;
  ris->sentinel= PALLOC(struct _node);
  ris->sentinel->prev= ris->sentinel->next= ris->sentinel;
  ris->sentinel->element= NULL;
  return ris;
}


void list_destroy(plist l, delete_function myfree) {
  NOT_NULL(l);
  FINETRACE("Destroying an %zu element list.", list_size(l));
  _pnode n= l->sentinel;
  while (n->next != l->sentinel) {
	 _pnode tmp= n->next->next;
	 if (myfree!=NULL && n->next->element!=NULL)
		myfree(n->next->element);
	 pfree(n->next);
	 n->next= tmp;
  }
  pfree(l->sentinel);
  pfree(l);
}

void list_add_to_head(plist l, item p) {
  NOT_NULL(l);
  NOT_NULL(p);
  _pnode n= PALLOC(struct _node);
  n->element= p;
  n->next= l->sentinel->next;
  n->prev= l->sentinel;
  n->next->prev= n;
  l->sentinel->next= n;
  l->size= l->size+1;
}

void list_add_to_tail(plist l, item p) {
  NOT_NULL(l);
  NOT_NULL(p);
  _pnode n= PALLOC(struct _node);
  n->element= p;
  n->prev= l->sentinel->prev;
  n->next= l->sentinel;
  n->prev->next= n;
  l->sentinel->prev= n;
  l->size= l->size+1;
}

/*
 * Procedura per aggiungere un item appena prima un plistit
 */
void list_add_before_iterator(plistit it, plist l, item p) {
  NOT_NULL(l);
  NOT_NULL(p);
  my_assert(it != NULL && it->prev != it->l->sentinel);

//Se e' il primo elemento della lista
  if(it->prev->prev == it->l->sentinel) {
	 list_add_to_head(l, p);
  } else {
	 _pnode n= PALLOC(struct _node);
	 n->element= p;

	 n->prev=it->prev->prev;
	 n->next=it->prev;

	 it->prev->prev->next=n;
	 it->prev->prev=n;

	 l->size= l->size+1;
  }
}

item list_remove_from_head(plist l) {
  NOT_NULL(l);
  // Non deve essere vuota
  my_assert(l->sentinel->next != l->sentinel);
  _pnode n= l->sentinel->next;
  item p= n->element;
  l->sentinel->next= n->next;
  n->next->prev= l->sentinel;
  pfree(n);
  l->size= l->size-1;
  return p;
}

item list_remove_from_tail(plist l) {
  NOT_NULL(l);
  // Non deve essere vuota
  my_assert(l->sentinel->prev != l->sentinel);
  _pnode n= l->sentinel->prev;
  item p= n->element;
  l->sentinel->prev= n->prev;
  n->prev->next= l->sentinel;
  pfree(n);
  l->size= l->size-1;
  return p;
}

item list_head(plist l) {
  NOT_NULL(l);
  return l->sentinel->next->element;
}

item list_tail(plist l) {
  NOT_NULL(l);
  return l->sentinel->prev->element;
}

bool list_is_empty(plist l) {
  NOT_NULL(l);
  return (l->sentinel->next == l->sentinel);
}

bool
list_remove_element(plist l, item p, delete_function del) {
  NOT_NULL(l);
  NOT_NULL(p);
  bool removed= false;
  plistit lit= list_first(l);
  while (!removed && listit_has_next(lit)) {
	 if (listit_next(lit)==p) {
		list_remove_at_iterator(lit, del);
		removed= true;
	 }
  }
  listit_destroy(lit);
  return removed;
}

void list_sort(plist l, comparator cmp){
  NOT_NULL(l);
  NOT_NULL(cmp);
  int size = list_size(l);
  item* base= NPALLOC(item, size);
  int i = 0;
  plistit it = list_first(l);
  while(listit_has_next(it)){
	 base[i] = listit_next(it);
	 ++i;
  }
  listit_destroy(it);

  qsort(base, size, sizeof(item), cmp);

  _pnode tmp = l->sentinel->next;
  i= 0;
  while(tmp != l->sentinel){
	 tmp->element = base[i];
	 ++i;
	 tmp = tmp->next;
  }
  pfree(base);
}

plist list_merge_new(plist l1, plist l2){
  NOT_NULL(l1);
  NOT_NULL(l2);
  plistit it;
  plist mergedList = list_create();
  it = list_first(l1);
  while(listit_has_next(it)) {
	 list_add_to_tail(mergedList,listit_next(it));
  }
  listit_destroy(it);
  it = list_first(l2);
  while(listit_has_next(it)) {
	 list_add_to_tail(mergedList,listit_next(it));
  }
  listit_destroy(it);
  return mergedList;
}

void list_merge(plist l1, plist l2){
  NOT_NULL(l1);
  NOT_NULL(l2);
  l2->sentinel->next->prev = l1->sentinel->prev;
  l1->sentinel->prev->next = l2->sentinel->next;
  l1->sentinel->prev = l2->sentinel->prev;
  l2->sentinel->prev->next = l1->sentinel;

  pfree(l2->sentinel);
  pfree(l2);
}

plist list_copy(plist l, copy_item copy){
  NOT_NULL(l);
  NOT_NULL(copy);
  plist res = list_create();
  plistit it = list_first(l);

  while(listit_has_next(it)) {
	 list_add_to_tail(res, copy(listit_next(it)));
  }
  listit_destroy(it);
  return res;
}

void list_remove_at_iterator(plistit it,
									  delete_function delete_item){
  NOT_NULL(it);
  NOT_NULL(delete_item);

  it->next->prev = it->prev->prev;
  it->prev->prev->next = it->next;

  delete_item(it->prev->element);

  --it->l->size;

  _pnode new_prev= it->prev->prev;
  pfree(it->prev);
  it->prev= new_prev;
}

void list_difference(plist l1, plist l2, comparator cmp, delete_function delete_node){

  NOT_NULL(l1);
  NOT_NULL(l2);
  NOT_NULL(cmp);
  NOT_NULL(delete_node);
  plistit it1 = list_first(l1);
  plistit it2 = list_first(l2);

  while( (listit_has_next(it1)) && (listit_has_next(it2)) ){
	 int ris = cmp(&it1->next->element,&it2->next->element);
	 if( ris == 0 ){
		listit_next(it1);
		list_remove_at_iterator(it1, delete_node);
	 } else{
		if( ris > 0 )
		  listit_next(it2);
		else
		  listit_next(it1);
	 }//endelse
  }//endwhile
  listit_destroy(it1);
  listit_destroy(it2);
}//endlistdiff

/****RAFFA****/
bool list_compare(plist l1, plist l2, comparator cmp){
	int ris;
	bool stop;

  if(list_size(l1) != list_size(l2))
	  return false;

  plistit it1 = list_first(l1);
  plistit it2 = list_first(l2);

  stop=false;

  while( (listit_has_next(it1)) && (listit_has_next(it2)) && !stop ){
	 ris = cmp(&it1->next->element,&it2->next->element);
	 if( ris == 0 ){
		listit_next(it1);
		listit_next(it2);
	 } else{
		stop=true;
	 }//endelse
  }//endwhile
  listit_destroy(it1);
  listit_destroy(it2);
  if(stop)
	  return false;
  else
	  return true;
}//endlistcompare
/****RAFFA****/

/****RAFFA****/
/*
 * If allowed_diff is -1, the containment must be perfect (one factorization must be a perfect exon
 * subsequence of the other one)
 * If allowed_diff is 0, the containment must be perfect only on confirmed ss
 * If allowed_diff is > 0 the containment is relaxed (limited by allowed_diff) also on confirmed ss
 * The difference on unconfirmed ss must be less than a threshold
 * If return is 0, then l1 and l2 are different
 * If return is -1, then l1 is contained in l2
 * If return is 1, then l2 is contained in l1
 * If return is -2, then l1 is equal to l2
 */
int relaxed_list_contained(plist l1, plist l2, relaxed_comparator cmp, int allowed_diff){
	int ris;
	bool found, stop;
	int cfr_type;
	int actual_allowed_diff;

	my_assert(allowed_diff >= -1);

	if(list_size(l1) == list_size(l2))
		return relaxed_list_compare(l1, l2, (relaxed_comparator)cmp, allowed_diff);

	if(list_size(l1) == 1 || list_size(l2) == 1)
		return 0;

	actual_allowed_diff=(allowed_diff == -1)?(0):(allowed_diff);

	plistit it1=list_first(l1);
	plistit it2=list_first(l2);
	found=false;
	cfr_type=(allowed_diff == -1)?(0):(-2);

	unsigned int count_long=1;

	if(list_size(l1) > list_size(l2)){
		//Cerco l'elemento di l1 che matcha con il primo di l2
		while(listit_has_next(it1) && !found){
			ris = cmp(&it1->next->element,&it2->next->element, cfr_type, actual_allowed_diff, l1);
			if(ris == 0){
				found=true;
				listit_next(it2);
			}
			else
				count_long++;

			listit_next(it1);

			if(cfr_type == -2)
				cfr_type=2;
		}
	}
	else{
		//Cerco l'elemento di l2 che matcha con il primo di l1
		while(listit_has_next(it2) && !found){
			ris = cmp(&it2->next->element,&it1->next->element, cfr_type, actual_allowed_diff, l2);
			if(ris == 0){
				found=true;
				listit_next(it1);
			}
			else
				count_long++;

			listit_next(it2);

			if(cfr_type == -2)
				cfr_type=2;
		}
	}

	if(!found){
		listit_destroy(it1);
		listit_destroy(it2);
		return 0;
	}

	unsigned int count_factors=1;
	stop=false;

	if(list_size(l1) > list_size(l2)){
		//Controllo se l2 e' contenuta in l1
		while( (listit_has_next(it1)) && (listit_has_next(it2)) && !stop ){

			cfr_type=(allowed_diff == -1)?(0):((count_factors+1 == list_size(l2))?((count_long+1 == list_size(l1))?(-1):(1)):(0));

			ris = cmp(&it1->next->element,&it2->next->element, cfr_type, actual_allowed_diff, l1);

			if( ris == 0 ){
				listit_next(it1);
				listit_next(it2);
			}else{
				stop=true;
			}//endelse
			count_factors++;
			count_long++;
		}//endwhile
	}
	else{
		//Controllo se l1 e' contenuta in l2
		while( (listit_has_next(it1)) && (listit_has_next(it2)) && !stop ){

			cfr_type=(allowed_diff == -1)?(0):((count_factors+1 == list_size(l1))?((count_long+1 == list_size(l2))?(-1):(1)):(0));

			ris = cmp(&it2->next->element,&it1->next->element, cfr_type, actual_allowed_diff, l2);

			if( ris == 0 ){
				listit_next(it1);
				listit_next(it2);
			}else{
				stop=true;
			}//endelse
			count_factors++;
			count_long++;
		}//endwhile
	}

	listit_destroy(it1);
	listit_destroy(it2);

	if(stop)
		return 0;
	else{
		if(list_size(l1) >= list_size(l2)){
			if(count_factors == list_size(l2))
				return 1;
			else
				return 0;
		}
		else{
			if(count_factors == list_size(l1))
				return -1;
			else
				return 0;
		}
	}
}//endlistcompare
/****RAFFA****/

/****RAFFA****/
/*
 * If allowed_diff is -1, the comparison is perfect
 * If allowed_diff is 0, the comparison is perfect only on confirmed ss
 * If allowed_diff is > 0 the comparison is relaxed (limited by allowed_diff)on confirmed ss
 * If return -2, l1 is equal to l2, otherwise the return value is 0
 */
int relaxed_list_compare(plist l1, plist l2, relaxed_comparator cmp, int allowed_diff){
	int ris;
	bool stop;
	int cfr_type;
	int count_factors;
	int actual_allowed_diff;

	my_assert(allowed_diff >= -1);

	if(list_size(l1) != list_size(l2) || list_size(l1) == 1)
		return 0;

	  plistit it1 = list_first(l1);
	  plistit it2 = list_first(l2);

  stop=false;
  count_factors=1;

  actual_allowed_diff=(allowed_diff == -1)?(0):(allowed_diff);

  while( (listit_has_next(it1)) && (listit_has_next(it2)) && !stop ){

	 cfr_type=(allowed_diff == -1)?(0):((count_factors == 1)?(-2):((count_factors == (int)list_size(l1))?(-1):(0)));

	 ris = cmp(&it1->next->element,&it2->next->element, cfr_type, actual_allowed_diff, l1);

	 if( ris == 0 ){
		listit_next(it1);
		listit_next(it2);
	 } else{
		stop=true;
	 }//endelse
	 count_factors++;
  }//endwhile
  listit_destroy(it1);
  listit_destroy(it2);
  if(stop)
	  return 0;
  else
	  return -2;
}//endlistcompare
/****RAFFA****/

void list_complete_difference(plist l1, plist l2, comparator cmp, delete_function delete_node){

  plistit it1 = list_first(l1);
  while (listit_has_next(it1)) {
	 item i1= listit_next(it1);
	 plistit it2 = list_first(l2);
	 int ris= 1;
	 while (ris!=0 && listit_has_next(it2)) {
		item i2= listit_next(it2);
		ris = cmp(&i1, &i2);
	 }
	 if (ris==0) {
		list_remove_at_iterator(it1, delete_node);
	 }
	 listit_destroy(it2);
  }//endwhile
  listit_destroy(it1);
}//endlistdiff


/*
 * List iterator definitions
 */

plistit list_last(plist const l) {
  my_assert(l!=NULL);
  plistit li= PALLOC(struct _listit);
  li->prev= l->sentinel->prev;
  li->next= l->sentinel;
  li->l= l;
  li->sentinel= l->sentinel;
  return li;
}

void list_last_reuse(plist const l, plistit* pli) {
  my_assert(l!=NULL);
  plistit li= NULL;
  if (*pli == NULL) {
	 li= PALLOC(struct _listit);
	 *pli= li;
  } else {
	 li= *pli;
  }
  li->prev= l->sentinel->prev;
  li->next= l->sentinel;
  li->l= l;
  li->sentinel= l->sentinel;
}

plistit listit_copy(plistit li) {
  my_assert(li!=NULL);
  plistit rli= PALLOC(struct _listit);
  rli->prev= li->prev;
  rli->next= li->next;
  rli->l= li->l;
  rli->sentinel= li->sentinel;
  return rli;
}

void
listit_copy_reuse(plistit const li, plistit* prli) {
  my_assert(li!=NULL);
  plistit rli= NULL;
  if (*prli == NULL) {
	 rli= PALLOC(struct _listit);
	 *prli= rli;
  } else {
	 rli= *prli;
  }
  rli->prev= li->prev;
  rli->next= li->next;
  rli->l= li->l;
  rli->sentinel= li->sentinel;
}

item listit_get(plistit li) {
  my_assert(li!=NULL);
  item p= li->prev->element;
  return p;
}

bool listit_has_prev(plistit li) {
  my_assert(li!=NULL);
  return li->prev != li->sentinel;
}

item listit_prev(plistit li) {
  my_assert(li!=NULL);
  my_assert(listit_has_prev(li));
  item p= li->prev->element;
  li->next= li->prev;
  li->prev= li->prev->prev;
  return p;
}

