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
 * @file refine-intron.h
 *
 * Procedures for refining the intron splice sites.
 *
 **/

#ifndef _REFINE_INTRON_H_
#define _REFINE_INTRON_H_

#include "types.h"
#include "configuration.h"

//Include

//Refining of 5' and 3' intron splice sites
bool refine_intron(pconfiguration, pEST_info, pEST_info, pfactor, pfactor, bool);

int Try_Burset_after_match(char *, char *, int *, int *, int *, int, int);

int Check_Burset_patterns(char *, int, int);

int getBursetFrequency_adaptor(const char* const t,
										 const size_t cut1, const size_t cut2);

int getBursetFrequency(char *, char *);

//Returns a list of gap alignments
plist compute_gap_alignment(char *, char *, bool, int, int, int);

void ComputeGapAlignMatrix(char *, char *, char **, char **, char **, char *, bool);

void TracebackGapAlignment(size_t, pgap_alignment, char *, char *, char **, char **, char **, int, int, char);

void Find_AG_after_on_the_right(pgap_alignment, int, int *, int *, int *);

void Find_ACCEPTOR_before_on_the_left(pgap_alignment, int, int *, int *, int *, char *);

void Find_ACCEPTOR_after_on_the_left(pgap_alignment, int, int *, char *);

bool Shift_right_to_left_1(char *, char *, int, pgap_alignment, int *, int *, int *, char *);

bool Shift_right_to_left_2(char *, char *, int, pgap_alignment, int *, int *, int *, char *);

bool Shift_left_to_right_1(char *, char *, int, pgap_alignment, int *, int *, int *, char *);

bool Shift_left_to_right_2(char *, char *, int, pgap_alignment, int *, int *, int *, char *);

void Find_AG_before_on_the_right(pgap_alignment, int, int *);

char *To_lower(char *);

char *To_upper(char *);

#endif
