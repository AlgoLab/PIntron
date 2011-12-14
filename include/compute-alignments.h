/**
 *
 *
 *                              PIntron
 *
 * A novel pipeline for computational gene-structure prediction based on
 * spliced alignment of expressed sequences (ESTs and mRNAs).
 *
 * Copyright (C) 2010  Raffaella Rizzi
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
 * @file compute-alignment.h
 *
 * Procedures for computing alignments.
 *
 **/

#ifndef _COMPUTE_ALIGNMENT_H_
#define _COMPUTE_ALIGNMENT_H_

#include "types.h"

/*
 * Returns a list of alignments
 */
plist compute_alignment(char *, char *, bool);

/*
 * Returns the direction matrix and the score in the score parameter
 */
unsigned int
ComputeAlignMatrix(const char * const EST_seq,
						 const size_t EST_len,
						 const char * const genomic_seq,
						 const size_t genomic_len,
						 char * const Mdir);

void TracebackAlignment(size_t, palignment, char *, char *, char *, int, int);

size_t*
edit_distance_matrix(const char* const s1, const size_t l1,
							const char* const s2, const size_t l2);

size_t
compute_edit_distance(const char* const s1, const size_t l1,
							 const char* const s2, const size_t l2);

bool K_band_edit_distance(char *, char *, unsigned int, unsigned int *);

#endif

