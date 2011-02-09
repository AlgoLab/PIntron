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
** log-build-info.h
**
** Made by Yuri Pirola
** Login   <yuri@yuri>
**
** Started on  Mon Aug 24 00:38:28 2009 Yuri Pirola
** Last update Mon Aug 24 00:38:28 2009 Yuri Pirola
*/

#ifndef __LOG_BUILD_INFO_H__
#define __LOG_BUILD_INFO_H__

#include "log.h"

#ifndef __SRC_DESC
#define __SRC_DESC "not available"
#endif

#ifndef __BUILD_DESC
#define __BUILD_DESC "not available"
#endif

#ifndef __BUILD_DATETIME
#define __BUILD_DATETIME "unknown"
#endif

#ifndef __BUILD_HOST
#define __BUILD_HOST "unknown"
#endif

#ifndef __COMPILER_VER
#ifdef __VERSION__
#define QUOTEME( x ) #x
#define __COMPILER_VER QUOTEME(__VERSION__)
#else
#define __COMPILER_VER "unknown"
#endif
#endif

#define PRINT_SYSTEM_INFORMATION													\
  do {																						\
	 INFO("Program compiled %s @ %s.", __BUILD_DATETIME,						\
			__BUILD_HOST);																	\
	 INFO("Source version >%s<.", __SRC_DESC);									\
	 INFO("Configuration  %s.", __BUILD_DESC);									\
	 INFO("Compiler %s", __COMPILER_VER);											\
  } while (0)

#define PRINT_LICENSE_INFORMATION \
  do {									 \
  INFO("This program is part of PIntron package. Copyright (C) 2010  Paola Bonizzoni, Gianluca Della Vedova, Yuri Pirola, Raffaella Rizzi."); \
  INFO("This program is distributed under the terms of the GNU Affero General Public License (AGPL), either version 3 of the License, or (at your option) any later version."); \
  INFO("This program comes with ABSOLUTELY NO WARRANTY. See the GNU Affero General Public License for more details."); \
  INFO("This is free software, and you are welcome to redistribute it under the conditions specified by the license.");	\
  } while (0)

#endif /* !LOG-BUILD-INFO_H_ */
