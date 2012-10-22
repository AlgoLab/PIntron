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
/**
 *
 * @file log.h
 *
 * Funzioni per la stampa di messaggi di log.
 *
 **/

#ifndef _LOG_H_
#define _LOG_H_

#define LOG_LEVEL_FATAL (0)
#define LOG_LEVEL_ERROR (1)
#define LOG_LEVEL_WARN (2)
#define LOG_LEVEL_INFO (3)
#define LOG_LEVEL_DEBUG (4)
#define LOG_LEVEL_TRACE (5)
#define LOG_LEVEL_FINETRACE (6)

#ifndef LOG_THRESHOLD
#define LOG_THRESHOLD LOG_LEVEL_INFO
#endif

#ifndef LOG_PREFIX
#define LOG_PREFIX "* "
#endif


#endif // _LOG_H_


#ifdef LOG
#undef __INTERNAL_LOG
#undef LOG
#undef FATAL
#undef ERROR
#undef WARN
#undef INFO
#undef DEBUG
#undef TRACE
#undef FINETRACE
#endif

#ifdef LOG_MSG

#include <stdio.h>
#include <string.h>

#define MAX_LEN_FUNC_NAME 16
#define MAX_LEN_FILE_NAME 20

extern const char* const __LOG_PREFIXES__[];

#define LOG(level, ...) __INTERNAL_LOG(level, __LOG_PREFIXES__[level], __VA_ARGS__, "")

#define __INTERNAL_LOG(level, prefix, format, ...) do {						\
	 if (level<=LOG_THRESHOLD) {														\
		char __my_internal_funz__[MAX_LEN_FUNC_NAME+1];							\
		const int __my_internal_log_len__= strlen(__func__);					\
		strncpy(__my_internal_funz__, __func__, MAX_LEN_FUNC_NAME);			\
		for (int i= __my_internal_log_len__; i<MAX_LEN_FUNC_NAME; ++i) {	\
		  __my_internal_funz__[i]=' ';												\
		}																						\
		if (__my_internal_log_len__+2>=MAX_LEN_FUNC_NAME)						\
		  __my_internal_funz__[MAX_LEN_FUNC_NAME-1]=								\
			 __my_internal_funz__[MAX_LEN_FUNC_NAME-2]= '.';					\
		__my_internal_funz__[MAX_LEN_FUNC_NAME]= '\0';							\
		char __my_internal_file__[MAX_LEN_FILE_NAME+1];							\
		const int __my_internal_file_len__= strlen(__FILE__);					\
		strncpy(__my_internal_file__, __FILE__+									\
				  (__my_internal_file_len__>MAX_LEN_FILE_NAME ?					\
					__my_internal_file_len__-MAX_LEN_FILE_NAME : 0)				\
				  , MAX_LEN_FILE_NAME);													\
		fprintf (stderr, LOG_PREFIX "%s(%s@%*.*s:%-4d) " format "  %s\n",	\
					prefix, __my_internal_funz__,										\
					MAX_LEN_FILE_NAME, MAX_LEN_FILE_NAME,							\
					__my_internal_file__, __LINE__,									\
					__VA_ARGS__);														\
	 }																							\
  } while (0)

#else

#define LOG(level, prefix, ...) do { } while (0)

#endif


#if (LOG_LEVEL_FATAL <= LOG_THRESHOLD) && defined LOG_MSG
#define LOG_FATAL_ENABLED
#define FATAL(...) LOG(LOG_LEVEL_FATAL, __VA_ARGS__)
#else
#undef LOG_FATAL_ENABLED
#define FATAL(...) do { } while (0)
#endif

#if (LOG_LEVEL_ERROR <= LOG_THRESHOLD) && defined LOG_MSG
#define LOG_ERROR_ENABLED
#define ERROR(...) LOG(LOG_LEVEL_ERROR, __VA_ARGS__)
#else
#undef LOG_ERROR_ENABLED
#define ERROR(...) do { } while (0)
#endif

#if (LOG_LEVEL_WARN <= LOG_THRESHOLD) && defined LOG_MSG
#define LOG_WARN_ENABLED
#define WARN(...) LOG(LOG_LEVEL_WARN, __VA_ARGS__)
#else
#undef LOG_WARN_ENABLED
#define WARN(...) do { } while (0)
#endif

#if (LOG_LEVEL_INFO <= LOG_THRESHOLD) && defined LOG_MSG
#define LOG_INFO_ENABLED
#define INFO(...) LOG(LOG_LEVEL_INFO, __VA_ARGS__)
#else
#undef LOG_INFO_ENABLED
#define INFO(...) do { } while (0)
#endif

#if (LOG_LEVEL_DEBUG <= LOG_THRESHOLD) && defined LOG_MSG
#define LOG_DEBUG_ENABLED
#define DEBUG(...) LOG(LOG_LEVEL_DEBUG, __VA_ARGS__)
#else
#undef LOG_DEBUG_ENABLED
#define DEBUG(...) do { } while (0)
#endif

#if (LOG_LEVEL_TRACE <= LOG_THRESHOLD) && defined LOG_MSG
#define LOG_TRACE_ENABLED
#define TRACE(...) LOG(LOG_LEVEL_TRACE, __VA_ARGS__)
#else
#undef LOG_TRACE_ENABLED
#define TRACE(...) do { } while (0)
#endif

#if (LOG_LEVEL_FINETRACE <= LOG_THRESHOLD) && defined LOG_MSG
#define LOG_FINETRACE_ENABLED
#define FINETRACE(...) LOG(LOG_LEVEL_FINETRACE, __VA_ARGS__)
#else
#undef LOG_FINETRACE_ENABLED
#define FINETRACE(...) do { } while (0)
#endif

