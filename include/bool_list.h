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
/**
 *
 * @file bool_list.h
 *
 * A data structure that represents a linked list of booleans.
 *
 **/

#ifndef _BOOL_LIST_H_
#define _BOOL_LIST_H_

#include "generic.h"
#include <stdbool.h>
#include <stddef.h>

typedef struct _boollist* pboollist;
typedef struct _boollistit* pboollistit;

#define BTYPE bool

pboollist boollist_create(void);

void boollist_destroy(pboollist);

void boollist_clear(pboollist);

void boollist_add_to_head(pboollist, BTYPE);

void boollist_add_to_tail(pboollist, BTYPE);

BTYPE boollist_remove_from_head(pboollist);

BTYPE boollist_remove_from_tail(pboollist);

BTYPE boollist_head(pboollist);

BTYPE boollist_tail(pboollist);

bool boollist_is_empty(pboollist);

size_t boollist_size(pboollist);


/*
*Ordina la lista in base alla funzione di ordinamento definita sugli elementi
*che contiene passata come secondo argomento.
*/
void boollist_sort(pboollist, comparator);

/*
*Crea una nuova lista che e' la concatenazione della prima con la seconda.
*I dati non vengono copiati, la nuova lista punta alle locazioni di memoria
*contenenti i dati delle due liste passate come parametri in ingresso.
*/

pboollist boollist_merge_new(pboollist, pboollist);

/*
*Concatena la prima lista passata come argomento con la seconda.
*La locazione di memoria contenente il puntatore alla seconda lista viene deallocata.
*/

void boollist_merge(pboollist, pboollist);

/*
*Crea una nuova lista che e' la copia della lista passata come primo argomento;
*nel creare la nuova lista, i dati che contiene quella passata in ingresso vengono
*copiati e memorizzati in nuove locazioni di memoria mediante la funzione passata
*come secondo argomento: essa infatti si occupa di allocare la memoria per il dato
*da copiare, quindi copia il dato nella nuova locazione di memoria. Ritorna
*il puntatore alla nuova locazione di memoria contenente la copia del dato originale.
*/

pboollist boollist_copy(pboollist);

/*
*Rimuove l'elemento prev puntato dall'iteratore inizializzato sulla lista passata
*come primo argomento. Il terzo argomento della funzione consiste in una funzione
*il cui compito e' quello di liberare la memoria dagli item contenuti nella lista
*passata come primo argomento.
*/

void boollist_remove_at_iterator(pboollistit);


/*
 * Removes the elements of the first list that appear into the second one.
 * It uses the comparator function to compare the data and the delete_function
 * to deallocate the duplicated data.
 * WARNING!! It assumes that both the lists are sorted!!
 * */
void iboolist_difference(pboollist, pboollist);


/**
 * List iterator definitions
 **/

pboollistit boollist_first(pboollist);

pboollistit boollist_last(pboollist);

void boollistit_destroy(pboollistit);

bool boollistit_has_next(pboollistit);

BTYPE boollistit_next(pboollistit);

bool boollistit_has_prev(pboollistit);

BTYPE boollistit_prev(pboollistit);


#endif
