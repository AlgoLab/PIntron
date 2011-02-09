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
/**
 *
 * @file int_list.h
 *
 * A data structure that represents a linked list of integers.
 *
 **/

#ifndef _INT_LIST_H_
#define _INT_LIST_H_

#include "generic.h"
#include <stdbool.h>
#include <stddef.h>

typedef struct _intlist* pintlist;
typedef struct _intlistit* pintlistit;

#define ITYPE int

pintlist intlist_create(void);

void intlist_destroy(pintlist);

void intlist_clear(pintlist);

void intlist_add_to_head(pintlist, ITYPE);

void intlist_add_to_tail(pintlist, ITYPE);

ITYPE intlist_remove_from_head(pintlist);

ITYPE intlist_remove_from_tail(pintlist);

ITYPE intlist_head(pintlist);

ITYPE intlist_tail(pintlist);

bool intlist_is_empty(pintlist);

size_t intlist_size(pintlist);


/*
*Ordina la lista in base alla funzione di ordinamento definita sugli elementi
*che contiene passata come secondo argomento.
*/
void intlist_sort(pintlist, comparator);

/*
*Crea una nuova lista che e' la concatenazione della prima con la seconda.
*I dati non vengono copiati, la nuova lista punta alle locazioni di memoria
*contenenti i dati delle due liste passate come parametri in ingresso.
*/

pintlist intlist_merge_new(pintlist, pintlist);

/*
*Concatena la prima lista passata come argomento con la seconda.
*La locazione di memoria contenente il puntatore alla seconda lista viene deallocata.
*/

void intlist_merge(pintlist, pintlist);

/*
*Crea una nuova lista che e' la copia della lista passata come primo argomento;
*nel creare la nuova lista, i dati che contiene quella passata in ingresso vengono
*copiati e memorizzati in nuove locazioni di memoria mediante la funzione passata
*come secondo argomento: essa infatti si occupa di allocare la memoria per il dato
*da copiare, quindi copia il dato nella nuova locazione di memoria. Ritorna
*il puntatore alla nuova locazione di memoria contenente la copia del dato originale.
*/

pintlist intlist_copy(pintlist);

/*
*Rimuove l'elemento prev puntato dall'iteratore inizializzato sulla lista passata
*come primo argomento. Il terzo argomento della funzione consiste in una funzione
*il cui compito e' quello di liberare la memoria dagli item contenuti nella lista
*passata come primo argomento.
*/

void intlist_remove_at_iterator(pintlistit);


/*
 * Removes the elements of the first list that appear into the second one.
 * It uses the comparator function to compare the data and the delete_function
 * to deallocate the duplicated data.
 * WARNING!! It assumes that both the lists are sorted!!
 * */
void intlist_difference(pintlist, pintlist);


/**
 * List iterator definitions
 **/

pintlistit intlist_first(pintlist);

pintlistit intlist_last(pintlist);

void intlistit_destroy(pintlistit);

bool intlistit_has_next(pintlistit);

ITYPE intlistit_next(pintlistit);

bool intlistit_has_prev(pintlistit);

ITYPE intlistit_prev(pintlistit);


#endif
