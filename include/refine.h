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
/*
** refine.h
**
** Made by Yuri
** Login   <yuri@dottxxii-12.dottorato.disco.unimib.it>
**
** Started on  Wed Jul 22 14:05:39 2009 Yuri
** Last update Wed Jul 22 14:05:39 2009 Yuri
*/

#ifndef _REFINE_H_
#define _REFINE_H_

#include <stdbool.h>
#include <string.h>

/**
 * Compute the minimum cost alignment between p and t
 * considering a gap in the middle of t as 0 cost gap.
 * The maximum number of errors is max_errs.
 *
 * Return true if such alignment exists, and, if so,
 * out_offset_p indicate the position on p corresponding
 * to the gap on t, and out_offset_t1/2 the starting position
 * and the ending position of the gap on t
 * */
bool
refine_borders(const char* const p,
					const size_t len_p,
					const char* const t,
					const size_t len_t,
					const unsigned int max_errs,
					size_t* out_offset_p,
					size_t* out_offset_t1,
					size_t* out_offset_t2,
					unsigned int* out_edit_distance);

unsigned int*
edit_distance(const char* const s1, const size_t ls1,
				  const char* const s2, const size_t ls2);



#endif /* !_REFINE_H_ */
