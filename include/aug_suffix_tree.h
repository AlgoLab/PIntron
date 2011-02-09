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
** aug_suffix_tree.h
*/

#ifndef _AUG_SUFFIX_TREE_H_
#define _AUG_SUFFIX_TREE_H_

#include "lst_structs.h"
#include "lst_stree.h"
#include "lst_string.h"
#include "configuration.h"
#include "types.h"

#include <stdio.h>

struct _preproc_gen {
  pEST_info gen;
  char* alph;
  size_t* alph_occ;
  size_t alph_size;
  size_t gen_len;
  size_t* keys;
};

typedef struct _preproc_gen* ppreproc_gen;

ppreproc_gen
PGen_create(void);

void
PGen_destroy(ppreproc_gen);

#define get_key( pg, symbol ) ((pg)->keys[(unsigned int)(symbol)])

void
preprocess_text(pEST_info gen, ppreproc_gen const pg);

void stree_preprocess(LST_STree* stree, const ppreproc_gen const pg, const pconfiguration const conf);

void stree_info_destroy(LST_STree* stree);

void
stree_occurrences(LST_STree* stree, pconfiguration config, const int id);

void
stree_occurrences_destroy(LST_STree* stree);

void stree_print(LST_STree *stree, FILE* f);

char* string_print_func_with_dollars(LST_StringIndex *index);

#endif /* !_AUG_SUFFIX_TREE_H_ */
