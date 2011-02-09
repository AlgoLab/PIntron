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
** ext_array.h
*/


#ifndef _EXT_ARRAY_H_
#define _EXT_ARRAY_H_

#include <stdbool.h>
#include <stdlib.h>

#include "generic.h"

/**
 * Puntatore a un array estendibile
 **/
typedef struct _ext_array* pext_array;

/**
 * Puntatore a un iteratore di un array estendibile
 **/
typedef struct _ext_arrayit* pext_arrayit;

/**
 * Creazione di un array estendibile.
 * La dimensione iniziale e' INITIAL_SIZE.
 * Tempo costante.
 *
 * @return Un array estendibile.
 **/
pext_array
EA_create(void);

/**
 * Svuotamento di un array estendibile.
 * I singoli elementi vengono liberati con una chiamata a @p del_fn(...)
 *
 * @warning Non libera le singole chiavi
 * perche' sono generalmente un campo dell'item.
 *
 * Tempo lineare nella dimensione dell'array.
 *
 * @param arr Array da svuotare.
 * @param del_fn Funzione per la deallocazione dell'elemento.
 *
 **/
void
EA_clear(pext_array arr, delete_function del_fn);


/**
 * Distruzione di un array estendibile.
 * Libera i singoli elementi con la
 * funzione @p del_fn specificata.
 *
 * @param arr Array da distruggere.
 * @param del_fn Funzione per la deallocazione dell'elemento.
 *
 * Tempo lineare nella dimensione dell'array.
 **/
void
EA_destroy(pext_array arr, delete_function del_fn);

/**
 * Inserimento di un nuovo elemento in coda quindi non
 * preserva l'ordinamento.
 * Se il fattore di carico supera 0.5, la capacita'
 * viene raddoppiata.
 * Tempo costante (ammortizzato).
 *
 * @param arr Array in cui inserire.
 * @param data Dati del nuovo elemento.
 **/
void
EA_insert(pext_array arr, item data);

/**
 * Ritorna true se l'array e' vuoto,
 * false altrimenti.
 * Tempo costante.
 *
 * @param arr Array da interrogare.
 * @return true se @p arr e' vuoto, false altrimenti.
 *
 **/
bool
EA_is_empty(pext_array arr);


/**
 * Ritorna l'elemento in posizione @p pos.
 * La posizione e' 0-based e deve essere valida.
 * Ovvero 0<= @p pos < EA_size(@p arr).
 * Tempo costante.
 *
 * @param arr Array da interrogare.
 * @param pos Posizione dell'elemento che si vuole ottenere.
 * @return item in posizione @p pos in @p arr.
 *
 **/
item
EA_get(pext_array arr, int pos);

/**
 * Assegna l'elemento in posizione @p pos.
 * La posizione e' 0-based e deve essere valida.
 * Ovvero 0<= @p pos < EA_size(@p arr).
 * Tempo costante.
 * L'elemento che viene sostituito viene 'perso'.
 *
 * @param arr Array da interrogare.
 * @param pos Posizione dell'elemento che si vuole ottenere.
 * @param el elemento da assegnare in posizione @p pos in @p arr.
 *
 **/
void
EA_set(pext_array arr, int pos, item el);


/**
 * Ritorna il numero di elementi contenuti nell'array @p arr.
 * Tempo costante.
 *
 * @param arr Array da interrogare.
 * @return numero di elementi in @p arr.
 *
 **/
size_t
EA_size(pext_array arr);


/**
 * Ordina l'array secondo il comparatore @p cmp specificato.
 * Impiega la qsort fornita da C.
 * @param arr Array da ordinare.
 * @param cmp Puntatore a funzione che specifica l'ordinamento.
 **/
void
EA_sort(pext_array arr, comparator cmp);

/**
 * Ricerca e restituisce la posizione della chiave @p val nell'array estendibile
 * ORDINATO @p arr secondo il comparatore specificato @p cmp.
 * Tempo logaritmico nella dimensione (se confronto in tempo costante).
 * @param arr Array estendibile ORDINATO.
 * @param val Elemento da ricercare.
 * @param cmp Puntatore a funzione che specifica l'ordinamento.
 * @return posizione dell'elemento @p val in @p arr, o -1 se non trovato.
 **/
int
EA_binary_search(pext_array arr, item val, comparator cmp);

/**
 * Ritorna un iteratore per l'array estendibile specificato.
 * @remark L'iteratore e' da distruggere
 *
 * Tempo costante.
 *
 * @param arr Array di cui si vuole l'iteratore.
 *
 * @return Un iteratore sull'array.
 **/
pext_arrayit
EA_begin(pext_array arr);

/**
 * Distrugge l'iteratore specificato.
 * Tempo costante.
 *
 * @param ea_it Iteratore da distruggere.
 *
 **/
void
ext_arrayit_destroy(pext_arrayit ea_it);

/**
 * Ritorna true se l'iteratore @p ea_it ha un elemento successivo
 * Tempo costante.
 *
 * @param ea_it Iteratore di un array estendibile.
 *
 * @return true se e' possibile spostare in avanti l'iteratore.
 **/
bool
ext_arrayit_has_next(pext_arrayit ea_it);

/**
 * Ritorna l'item corrente e passa al successivo.
 * Se sono alla fine dell'array restituisco NULL.
 * @remark L'item non e' da distruggere perche' e' ancora di
 * proprieta' dell'array.
 *
 * Tempo costante.
 *
 * @param ea_it Iteratore usato per la navigazione dell'array.
 *
 * @return Il successivo elemento della lista.
 *
 **/
item
ext_arrayit_next(pext_arrayit ea_it);

/**
 * Torna all'item precedente e lo ritorna.
 * Se sono all'inizio dell'array restituisco NULL.
 * @remark L'item non e' da distruggere perche' e' ancora di
 * proprieta' dell'array.
 * @warning Chiamare next e poi prev restituisce lo stesso elemento.
 *
 * Tempo costante.
 *
 * @param ea_it Iteratore usato per la navigazione dell'array.
 *
 * @return L'elemento precedente della lista.
 *
 **/
item
ext_arrayit_prev(pext_arrayit ea_it);

/**
 * Crea un nuovo iteratore e lo inizializza alla posizione di @p ea_it.
 *
 * @param ea_it Iteratore da clonare.
 * @return Un nuovo iteratore allo stesso stato di @p ea_it.
 **/
pext_arrayit
ext_arrayit_clone(pext_arrayit ea_it);

/**
 * Copia lo stato dell'iteratore @p source in @p target.
 *
 * @param source Iteratore sorgente.
 * @param target Iteratore destinazione.
 **/
void
ext_arrayit_copy(pext_arrayit target, pext_arrayit source);


#endif
