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
/**
 *
 * @file list.h
 *
 * A data structure that represents a linked list of generic items.
 *
 **/

#ifndef _LIST_H_
#define _LIST_H_

#include "generic.h"
#include <stdbool.h>
#include <stddef.h>

typedef struct _list* plist;
typedef struct _listit* plistit;

plist list_create(void);

void list_destroy(plist, delete_function);

void list_add_to_head(plist, item);

void list_add_to_tail(plist, item);

void list_add_before_iterator(plistit, plist, item);

item list_remove_from_head(plist);

item list_remove_from_tail(plist);

item list_head(plist);

item list_tail(plist);

bool list_is_empty(plist);

static inline
size_t list_size(const plist const);


/*
*Ordina la lista in base alla funzione di ordinamento definita sugli elementi
*che contiene passata come secondo argomento.
*/
void list_sort(plist, comparator);

/*
*Crea una nuova lista che e' la concatenazione della prima con la seconda.
*I dati non vengono copiati, la nuova lista punta alle locazioni di memoria
*contenenti i dati delle due liste passate come parametri in ingresso.
*/

plist list_merge_new(plist, plist);

/*
*Concatena la prima lista passata come argomento con la seconda.
*La locazione di memoria contenente il puntatore alla seconda lista viene deallocata.
*/

void list_merge(plist, plist);

/*
*Crea una nuova lista che e' la copia della lista passata come primo argomento;
*nel creare la nuova lista, i dati che contiene quella passata in ingresso vengono
*copiati e memorizzati in nuove locazioni di memoria mediante la funzione passata
*come secondo argomento: essa infatti si occupa di allocare la memoria per il dato
*da copiare, quindi copia il dato nella nuova locazione di memoria. Ritorna
*il puntatore alla nuova locazione di memoria contenente la copia del dato originale.
*/

plist list_copy(plist, copy_item);

/*
*Rimuove l'elemento prev puntato dall'iteratore. Il secondo argomento
*della funzione consiste in una funzione
*il cui compito e' quello di liberare la memoria dagli item contenuti nella lista
*a cui l'iteratore si riferisce.
*/

void list_remove_at_iterator(plistit, delete_function);

/**
 * Rimuove il nodo della lista @p l che punta all'elemento @p p impiegando la funzione
 * @p del.
 * Nota: in caso piu' nodi puntano a @p p, viene rimosso solo il primo nodo.
 * Ritorna true se la rimozione e' stata effettuata.
 * Complessita': O(list_size)
 **/
bool list_remove_element(plist l, item p, delete_function del);

/*
 * Removes the elements of the first list that appear into the second one.
 * It uses the comparator function to compare the data and the delete_function
 * to deallocate the duplicated data.
 * WARNING!! It assumes that both the lists are sorted!!
 * */
void list_difference(plist, plist, comparator, delete_function);


void list_complete_difference(plist, plist, comparator, delete_function);


/****RAFFA****/
/*
 * Compares two lists and return true if they contain the same elements in the same order.
 */
bool list_compare(plist, plist, comparator);
/****RAFFA***/

/****RAFFA****/
/*
 * Compares two lists and return 0 if they are disjoint, and -1|+1 if one is contained in the other, -2
 * if the first is equal to the second.
 */
int relaxed_list_contained(plist, plist, relaxed_comparator, int);

/****RAFFA****/
/*
 * Compares two lists and return -2 if they contain similar elements in the same order, otherwise 0.
 */
int relaxed_list_compare(plist, plist, relaxed_comparator, int);

/**
 * List iterator definitions
 **/

static inline
plistit list_first(const plist const l);

void list_first_reuse(plist const l, plistit* pli);

plistit list_last(plist l);

void list_last_reuse(plist const l, plistit* pli);

plistit listit_copy(plistit li);

void listit_copy_reuse(plistit const li, plistit* prli);

static inline
void listit_destroy(plistit li);

static inline
bool listit_has_next(const plistit const li);

static inline
item listit_next(plistit li);

item listit_get(plistit li);

bool listit_has_prev(plistit li);

item listit_prev(plistit li);




#include "util.h"

typedef struct _node* _pnode;

struct _node {
	_pnode next;
	_pnode prev;
	item element;
};

struct _list {
  _pnode sentinel;
  size_t size;
};

struct _listit {
  _pnode next;
  _pnode prev;
  _pnode sentinel;
  struct _list* l;
};

static inline
size_t list_size(plist l) {
  NOT_NULL(l);
  return l->size;
}

static inline
plistit list_first(const plist const l) {
  my_assert(l!=NULL);
  plistit li= PALLOC(struct _listit);
  li->l= l;
  const _pnode const sentinel= l->sentinel;
  li->next= sentinel->next;
  li->prev= sentinel;
  li->sentinel= sentinel;
  return li;
}

static inline
void listit_destroy(plistit li) {
  if (li!=NULL)
	 pfree(li);
}

static inline
bool listit_has_next(const plistit const li) {
  my_assert(li!=NULL);
  return li->next != li->sentinel;
}

static inline
item listit_next(plistit li) {
  my_assert(li!=NULL);
  my_assert(listit_has_next(li));
  item p= li->next->element;
  li->prev= li->next;
  li->next= li->next->next;
  return p;
}

#endif
