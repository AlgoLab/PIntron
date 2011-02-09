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
#ifndef _MYTIME_H_
#define _MYTIME_H_

#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>

#include "log.h"

#define MYTIME_LOG( log_level, timer ) \
  do {														\
	 log_level("@Timer %-22s. "							\
				  "Time elapsed: %15llu microsec",	\
				  MYTIME_getname(timer),				\
				  MYTIME_getinterval(timer));			\
  } while (0)


#define DTYPE unsigned long long

typedef struct _mytime* pmytime;

void
MYTIME_print_interval(FILE*, pmytime);

void
MYTIME_print_current(FILE*, pmytime);

pmytime
MYTIME_create(void);

pmytime
MYTIME_create_with_name(const char*);

void
MYTIME_destroy(pmytime pt);

void
MYTIME_start(pmytime pt);

void
MYTIME_reset(pmytime pt);

void
MYTIME_stop(pmytime pt);

const char*
MYTIME_getname(pmytime pt);

unsigned long long
MYTIME_getinterval(pmytime pt);


#endif
