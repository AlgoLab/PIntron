/**
 *
 *
 *                              PIntron
 *
 * A novel pipeline for computational gene-structure prediction based on
 * spliced alignment of expressed sequences (ESTs and mRNAs).
 *
 * Copyright (C) 2012  Yuri Pirola
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
 * @file factorization-util.h
 *
 * Utility procedures for factorizations
 *
 **/

#ifndef _FACTORIZATION_UTIL_H_
#define _FACTORIZATION_UTIL_H_

#include "list.h"

#ifndef LOG_THRESHOLD
#define LOG_THRESHOLD LOG_LEVEL_INFO
#endif


void _impl_print_factorization_on_log_full(const int log_level,
												 plist factorization, const char* const gen_seq);

void _impl_print_factorizations_on_log_full(const int log_level,
												  plist factorization_list, const char* const gen_seq);

void _impl_print_factorization_on_log(const int log_level, plist factorization);

void _impl_print_factorizations_on_log(const int log_level, plist factorization_list);


#define print_factorization_on_log_full(log_level, factorization, gen_seq) \
  do {																						\
	 if ((log_level)<=LOG_THRESHOLD) {												\
		_impl_print_factorization_on_log_full(log_level,						\
														  factorization,					\
														  gen_seq);							\
	 }																							\
  } while (0)

#define print_factorizations_on_log_full(log_level, factorization_list, gen_seq) \
  do {																						\
	 if ((log_level)<=LOG_THRESHOLD) {												\
		_impl_print_factorizations_on_log_full(log_level,						\
															factorization_list,			\
															gen_seq);						\
	 }																							\
  } while (0)


#define print_factorization_on_log(log_level, factorization)				\
  do {																						\
	 if ((log_level)<=LOG_THRESHOLD) {												\
		_impl_print_factorization_on_log_full(log_level,						\
														  factorization,					\
														  NULL);								\
	 }																							\
  } while (0)


#define print_factorizations_on_log(log_level, factorization_list)		\
  do {																						\
	 if ((log_level)<=LOG_THRESHOLD) {												\
		_impl_print_factorizations_on_log_full(log_level,						\
															factorization_list,			\
															NULL);							\
	 }																							\
  } while (0)


#endif //_FACTORIZATION_UTIL_H_
