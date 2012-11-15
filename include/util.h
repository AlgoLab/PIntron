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
/**
 *
 * @file util.h
 *
 * Funzioni di utilita' generale.
 *
 **/

#ifndef _UTIL_H_
#define _UTIL_H_

#include <unistd.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include "log.h"

#define NPALLOC( type, dim ) (type*)palloc((dim)*sizeof(type))
#define PALLOC( type ) (type*)palloc(sizeof(type))

#define MY_SWAP( type, el1, el2 ) \
  {										 \
	 type tmp= (el1);					 \
	 (el1)= el2;						 \
	 (el2)= tmp;						 \
  }

#ifndef MIN
#define MIN( x, y ) (((x)<=(y))?(x):(y))
#endif
#ifndef MAX
#define MAX( x, y ) (((x)>=(y))?(x):(y))
#endif

#ifndef NDEBUG

#define fail()																				\
  do {																						\
	 fprintf(stderr, "Failing at %s:%d.\n", __FILE__, __LINE__);			\
	 char* pc= NULL;																		\
	 *pc= '\65';																			\
	 exit(1);																				\
  } while (0)

#define my_assert( cond )																\
  do {																						\
	 if (!(cond)) {																		\
		fprintf(stderr, "Assertion " #cond " failed at %s:%d.\n", __FILE__, __LINE__); \
		fail();																				\
	 }																							\
  } while (0)

#else

#define fail()																			\
  do {																					\
	 fprintf(stderr, "Failing at %s:%d.\n", __FILE__, __LINE__);		\
	 char* pc= NULL;																	\
	 *pc= '\65';																		\
	 exit(1);																			\
  } while (0)

#define my_assert( cond ) do { } while (0)

#endif

static inline
void* palloc(const size_t size) {

  void* p= malloc(size);
  if (p==NULL) {
	 FATAL("Allocation memory error. Trying to allocate %zu bytes.", size);
	 fail();
  }
  FINETRACE("Pointer %p allocated", p);
  return p;
}

static inline
void pfree(const void* const p) {
  if (p==NULL) {
	 FATAL("Cannot free a NULL pointer.");
	 fail();
  }
  free((void *)p);
}

void pfree_function(const void* const p);

char* c_palloc(size_t dim) __attribute__ ((malloc));

char* alloc_and_copy(const char* const source)  __attribute__ ((malloc));

void noop_free(void*) __attribute__ ((const));

char* substring(const int, const char* ) __attribute__ ((pure));

char* real_substring(int index, int length, const char* const string) __attribute__ ((pure));

char* reverse(const char* const s, const size_t len);

#ifdef __APPLE__

// Rewrite of GNU getline for MacOS compatibility
ssize_t custom_getline(char **lineptr, size_t *n, FILE *stream);

#else // on UNIX/Linux targets

#define custom_getline getline

#endif

// Chiama la getline e rimuove i caratteri non stampabili
// finali (ad es. \n)
ssize_t my_getline(char **lineptr, size_t *n, FILE *stream)
#ifndef __ICC
  __attribute__ ((nonnull))
#endif
  ;

void print_repetitions(FILE* f, const char c, int rep);

void resource_usage_log(void);

FILE* open_statm_file(void);

void log_info(FILE* const logfile, char* description);

void log_info_extended(FILE* const logfile, char* description, void* additional_info);

#define NOT_NULL( v ) my_assert( v != NULL )

#define fail_if( cond ) \
  do {						\
	 if (cond) {			\
		fail();				\
	 }							\
  } while (0)

#endif

#define DOT_100_STRING "...................................................................................................."

#define DOT_1K_STRING DOT_100_STRING DOT_100_STRING DOT_100_STRING DOT_100_STRING DOT_100_STRING DOT_100_STRING DOT_100_STRING DOT_100_STRING DOT_100_STRING DOT_100_STRING

#define DOT_10K_STRING DOT_1K_STRING DOT_1K_STRING DOT_1K_STRING DOT_1K_STRING DOT_1K_STRING DOT_1K_STRING DOT_1K_STRING DOT_1K_STRING DOT_1K_STRING DOT_1K_STRING

#define DOT_100K_STRING DOT_10K_STRING DOT_10K_STRING DOT_10K_STRING DOT_10K_STRING DOT_10K_STRING DOT_10K_STRING DOT_10K_STRING DOT_10K_STRING DOT_10K_STRING DOT_10K_STRING

#define DOT_1M_STRING DOT_100K_STRING DOT_100K_STRING DOT_100K_STRING DOT_100K_STRING DOT_100K_STRING DOT_100K_STRING DOT_100K_STRING DOT_100K_STRING DOT_100K_STRING DOT_100K_STRING

#define DOT_STRING DOT_1M_STRING

#define SPACE_32_STRING "                                "

#define SPACE_64_STRING SPACE_32_STRING SPACE_32_STRING

#define SPACE_256_STRING SPACE_64_STRING SPACE_64_STRING SPACE_64_STRING SPACE_64_STRING

#define SPACE_1024_STRING SPACE_256_STRING SPACE_256_STRING SPACE_256_STRING SPACE_256_STRING

#define SPACE_4095_STRING											\
  SPACE_1024_STRING SPACE_1024_STRING SPACE_1024_STRING	\
  SPACE_256_STRING SPACE_256_STRING SPACE_256_STRING		\
  SPACE_64_STRING SPACE_64_STRING SPACE_64_STRING			\
  SPACE_32_STRING "                               "

#define SPACE_STRING SPACE_4095_STRING
