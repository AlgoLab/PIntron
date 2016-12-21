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
#ifndef _SEMPL_INFO_H_
#define _SEMPL_INFO_H_

#include "list.h"
#include "types.h"
#include "bit_vector.h"
#include <stdio.h>
#include "log.h"

typedef struct _simpl_info* psimpl;

struct _simpl_info
{
  pbit_vect factors_used;
  pbit_vect ests_ok;
  pbit_vect factors_not_used;
};

psimpl psimpl_create(void);

void psimpl_destroy(psimpl);

int countTrue(pbit_vect);


#if defined (LOG_MSG) && (LOG_LEVEL_INFO <= LOG_THRESHOLD)
void psimpl_print(psimpl);
#else
#define psimpl_print( ps ) \
  do {							\
  } while(0)

#endif

#endif
