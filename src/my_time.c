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
#include "my_time.h"
#include "util.h"
#include <stdio.h>


struct _mytime {
  const char* timer_name;
  MYTIME_DTYPE interval;
  struct timeval start;
  struct timeval stop;
  bool active;
};

struct _mytime_parallel {
  pmytime parent;
  struct timeval start;
  struct timeval stop;
};


static const char* default_timer_name= "generic timer";


static MYTIME_DTYPE
diff_usec(struct timeval start, struct timeval stop) {
  MYTIME_DTYPE start_usec= ((MYTIME_DTYPE)start.tv_sec*(MYTIME_DTYPE)1000000)+(MYTIME_DTYPE)start.tv_usec;
  MYTIME_DTYPE stop_usec= ((MYTIME_DTYPE)stop.tv_sec*(MYTIME_DTYPE)1000000)+(MYTIME_DTYPE)stop.tv_usec;
  return stop_usec-start_usec;
}

static void
compute_difference(pmytime pt) {
  pt->interval+= diff_usec(pt->start, pt->stop);
}

void
MYTIME_print_interval(FILE* file, pmytime pt) {
  my_assert(pt!=NULL);
  my_assert(file!=NULL);
  MYTIME_DTYPE diff= pt->interval;
  if (diff>=(MYTIME_DTYPE)1000) {
// in secs
	 fprintf(file, "@Timer %s. Time elapsed: ", pt->timer_name);
	 MYTIME_DTYPE min= diff/60000000;
	 diff= diff%60000000;
	 if (min>0)
		fprintf(file, "%llum ", min);
	 fprintf(file, "%.3lfs\n", ((double)diff/1000000.0));
  } else {
	 fprintf(file, "@Timer %s. Time elapsed: %llumicrosec", pt->timer_name, diff);
  }
}

void
MYTIME_print_current(FILE* file, pmytime pt) {
  MYTIME_stop(pt);
  MYTIME_print_interval(file, pt);
  MYTIME_start(pt);
}



pmytime
MYTIME_create(void) {
  return MYTIME_create_with_name(default_timer_name);
}

pmytime
MYTIME_create_with_name(const char* timer_name) {
  my_assert(timer_name!=NULL);
  pmytime pt= PALLOC(struct _mytime);
  pt->active= false;
  pt->interval= 0;
  pt->timer_name= timer_name;
  return pt;
}

void
MYTIME_destroy(pmytime pt)
{
  my_assert(pt!=NULL);
  pfree(pt);
}

void
MYTIME_start(pmytime pt) {
  my_assert(pt!=NULL);
  pt->active= true;
  gettimeofday(&pt->start, NULL);
}

pmytime_parallel
MYTIME_start_parallel(pmytime pt) {
  my_assert(pt!=NULL);
  pmytime_parallel ppt= PALLOC(struct _mytime_parallel);
  ppt->parent= pt;
  gettimeofday(&ppt->start, NULL);
  return ppt;
}

void
MYTIME_stop_parallel(pmytime_parallel ppt) {
  my_assert(ppt!=NULL);
  gettimeofday(&(ppt->stop), NULL);
  ppt->parent->interval+= diff_usec(ppt->start, ppt->stop);
  pfree(ppt);
}

void
MYTIME_reset(pmytime pt)
{
  my_assert(pt!=NULL);
  pt->active= false;
  pt->interval= 0;
}

void
MYTIME_stop(pmytime pt)
{
  my_assert(pt!=NULL);
  gettimeofday(&(pt->stop), NULL);
  pt->active= false;
  compute_difference(pt);
}

const char*
MYTIME_getname(pmytime pt)
{
  my_assert(pt!=NULL);
  return pt->timer_name;
}

MYTIME_DTYPE
MYTIME_getinterval(pmytime pt)
{
  my_assert(pt!=NULL);
  return pt->interval;
}


struct _mytime_timeout {
  struct timeval start;
  struct timeval current;
  MYTIME_DTYPE time_limit;
  bool expired;
};

pmytime_timeout
MYTIME_timeout_create(MYTIME_DTYPE time_limit)
{
  pmytime_timeout ptt= PALLOC(struct _mytime_timeout);
  ptt->expired= false;
  ptt->time_limit= time_limit;
  gettimeofday(&ptt->start, NULL);
  return ptt;
}


bool
MYTIME_timeout_expired(pmytime_timeout ptt)
{
  my_assert(ptt!=NULL);
  if (ptt->expired)
	 return true;
  gettimeofday(&ptt->current, NULL);
  ptt->expired= diff_usec(ptt->start, ptt->current) > ptt->time_limit;
  return ptt->expired;
}

void
MYTIME_timeout_destroy(pmytime_timeout ptt)
{
  my_assert(ptt!=NULL);
  pfree(ptt);
}

