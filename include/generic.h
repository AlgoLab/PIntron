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
/*
** generic.h
*/
#ifndef _GENERIC_H_
# define _GENERIC_H_

#include <stdio.h>

/// Tipo di un generico item.
typedef void* item;

/// Tipo di una funzione per la deallocazione di un item.
typedef void (*delete_function)(item);

/// Tipo di una funzione per la stampa di un item.
typedef void (*print_function)(FILE*, item);

/// Tipo di una funzione per l'ordinamento di un item
typedef int (*comparator)(const void*,const void*);

/// RAFFA Tipo di una funzione per l'ordinamento relaxed di un item
typedef int (*relaxed_comparator)(const void*,const void*,int,int, item);

/// Tipo di una funzione per la copia di un item
typedef item* (*copy_item) (item);



#endif /* _GENERIC_H_ */
