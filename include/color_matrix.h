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
#ifndef _COLOR_MATRIX_H_
#define _COLOR_MATRIX_H_

#include "list.h"
#include "types.h"
#include <stdio.h>
#include "log.h"

plist factors_list_create(plist);
plist windows_list_create(plist);
bool equal_factor(pfactor,pfactor,bool);
bool verify_factor(plist,pfactor);
plist update_windows(plist, pfactor);
int get_factor_position(pfactor,plist,bool);
void add_factoriz(pEST ,plist,plist,bool);
plist color_matrix_create(plist, bool);

#if defined (LOG_MSG) && (LOG_LEVEL_DEBUG <= LOG_THRESHOLD)

void color_matrix_print(plist);
void print_factors_list(plist, bool);

#else

#define color_matrix_print( l ) \
  do {								  \
  } while(0)

#define print_factors_list( l ) \
  do {								  \
  } while(0)

#endif

#endif
