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
 * @file double_list.h
 *
 * A data structure that represents a linked list of doubles.
 *
 **/

#ifndef _DOUBLE_LIST_H_
#define _DOUBLE_LIST_H_

#include "generic.h"
#include <stdbool.h>
#include <stddef.h>

typedef struct _doublelist* pdoublelist;
typedef struct _doublelistit* pdoublelistit;

#define DTYPE double

pdoublelist doublelist_create(void);

void doublelist_destroy(pdoublelist);

void doublelist_clear(pdoublelist);

void doublelist_add_to_head(pdoublelist, DTYPE);

void doublelist_add_to_tail(pdoublelist, DTYPE);

DTYPE doublelist_remove_from_head(pdoublelist);

DTYPE doublelist_remove_from_tail(pdoublelist);

DTYPE doublelist_head(pdoublelist);

DTYPE doublelist_tail(pdoublelist);

bool doublelist_is_empty(pdoublelist);

size_t doublelist_size(pdoublelist);


/*
*Ordina la lista in base alla funzione di ordinamento definita sugli elementi
*che contiene passata come secondo argomento.
*/
//void doublelist_sort(pdoublelist, comparator);

/*
*Crea una nuova lista che e' la concatenazione della prima con la seconda.
*I dati non vengono copiati, la nuova lista punta alle locazioni di memoria
*contenenti i dati delle due liste passate come parametri in ingresso.
*/

pdoublelist doublelist_merge_new(pdoublelist, pdoublelist);

/*
*Concatena la prima lista passata come argomento con la seconda.
*La locazione di memoria contenente il puntatore alla seconda lista viene deallocata.
*/

void doublelist_merge(pdoublelist, pdoublelist);

/*
*Crea una nuova lista che e' la copia della lista passata come primo argomento;
*nel creare la nuova lista, i dati che contiene quella passata in ingresso vengono
*copiati e memorizzati in nuove locazioni di memoria mediante la funzione passata
*come secondo argomento: essa infatti si occupa di allocare la memoria per il dato
*da copiare, quindi copia il dato nella nuova locazione di memoria. Ritorna
*il puntatore alla nuova locazione di memoria contenente la copia del dato originale.
*/

pdoublelist doublelist_copy(pdoublelist);

/*
*Rimuove l'elemento prev puntato dall'iteratore inizializzato sulla lista passata
*come primo argomento. Il terzo argomento della funzione consiste in una funzione
*il cui compito e' quello di liberare la memoria dagli item contenuti nella lista
*passata come primo argomento.
*/

void doublelist_remove_at_iterator(pdoublelistit);


/*
 * Removes the elements of the first list that appear into the second one.
 * It uses the comparator function to compare the data and the delete_function
 * to deallocate the duplicated data.
 * WARNING!! It assumes that both the lists are sorted!!
 * */
void doublelist_difference(pdoublelist, pdoublelist);


/**
 * List iterator definitions
 **/

pdoublelistit doublelist_first(pdoublelist);

pdoublelistit doublelist_last(pdoublelist);

void doublelistit_destroy(pdoublelistit);

bool doublelistit_has_next(pdoublelistit);

DTYPE doublelistit_next(pdoublelistit);

bool doublelistit_has_prev(pdoublelistit);

DTYPE doublelistit_prev(pdoublelistit);


#endif
