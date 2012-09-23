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
#ifndef __BIT_VECTOR_H__
#define __BIT_VECTOR_H__

#include <stdbool.h>
#include <stddef.h>

#define _BTYPE unsigned int
#define _LBTYPE (sizeof(_BTYPE)<<3)

struct _bit_vect
{
  _BTYPE* arr;
  unsigned int n;
  size_t ncells;
};

typedef struct _bit_vect * pbit_vect;


pbit_vect BV_create(const unsigned int n);

void BV_clear(pbit_vect bv);

void BV_or(pbit_vect bv1, pbit_vect bv2);

void BV_destroy(pbit_vect bv);

void BV_set(pbit_vect bv, const unsigned int i, bool value);

void BV_set_block(pbit_vect bv, const unsigned int i, _BTYPE block);

bool BV_get(pbit_vect bv, const unsigned int i);

_BTYPE BV_get_block(pbit_vect bv, const unsigned int i);

_BTYPE BV_get_unaligned_block(pbit_vect v, const unsigned int i, const unsigned int l);

void BV_print(pbit_vect bv);

pbit_vect BV_clone(pbit_vect bv);

void BV_copy(pbit_vect ris, pbit_vect bv);

bool BV_contained(pbit_vect bv1, pbit_vect bv2);

bool BV_all_true(pbit_vect bv);

#ifdef LOG_MSG
char* BV_to_string(pbit_vect);
#else
// Generate a sintax error if the macro is expanded
#define BV_to_string( v )	___________sintax error
#endif

#endif
