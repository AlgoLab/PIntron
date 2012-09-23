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
#include "simplify_matrix.h"
#include "simpl_info.h"
#include "log.h"

void print_factorizations_result(pbit_vect,plist,plist,psimpl);

pbit_vect min_fact(plist);

plist color_matrix_simplified_create(plist, psimpl);

void color_matrix_simplified_destroy(plist);


#if defined (LOG_MSG) && (LOG_LEVEL_DEBUG <= LOG_THRESHOLD)

void print_min_fact(pbit_vect, plist);

#else

#define print_min_fact(v, l)							\
  do {														\
  } while (0)

#endif


#endif

