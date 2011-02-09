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
#ifndef _MIN_FACTORIZATION_H_
#define _MIN_FACTORIZATION_H_

#include "list.h"
#include "types.h"
#include <stdio.h>
#include "bit_vector.h"
#include "semplify_matrix.h"
#include "sempl_info.h"
#include "log.h"

bool create_combinations(int,int,pbit_vect,plist);

bool all_false(pbit_vect);

bool all_true(pbit_vect);

void print_factorizations_result(pbit_vect,plist,plist,psempl);

void color_matrix_semplified_destroy(plist);


pbit_vect min_fact(plist);

int count_true(pbit_vect);

int min_number_of_factor(plist);

int max_of_min(plist);

void addFactoriz(pbit_vect,plist,psempl);

void addEST(pEST,plist,psempl);

plist color_matrix_semplified_create(plist,psempl);

bool valuate_list(pbit_vect , plist);

bool valuate_combination(pbit_vect ,plist);

#if defined (LOG_MSG) && (LOG_LEVEL_DEBUG <= LOG_THRESHOLD)

void print_min_fact(pbit_vect, plist);

#else

#define print_min_fact(v, l)							\
  do {														\
  } while (0)

#endif


#endif

