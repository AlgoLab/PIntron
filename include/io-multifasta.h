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
//Header file for IO-MULTIFASTA.c

#include <stdio.h>
#include "list.h"
#include "types.h"

/*
  Reads a file that contains data in MULTIFASTA format
  and returns a plist that
  contains data provided by the input file
*/
plist read_multifasta(FILE* );

/*
  Receives a plist and a file opened in write mode and prints
  data contained in the plist
  whithin the file in the MULTIFASTA format
*/
void write_multifasta(plist , FILE* );

/*
  Receives a pEST_info (genomic), a plist (of pEST) and a file opened in write mode and prints
  data contained in the plist
  whithin the file in the MULTIFASTA format
  >Header
  est_left est_right gen_left gen_right
  if the boolean arg is set to true the esternal factors are deleted
*/
void write_multifasta_output(pEST_info, pEST , FILE* , char);

pEST_info read_single_EST_info(FILE*);

void write_single_EST_info(FILE*, pEST_info);

void set_EST_GB_identification(pEST_info);

void set_Chromosome(pEST_info);

void set_Chromosome_coordinates(pEST_info);

void set_Genomic_Strand(pEST_info);

void set_EST_Strand_and_RC(pEST_info, pEST_info);

void reverse_and_complement(pEST_info);

/*
  Look for a polyA or a polyT at the ends of the sequence and
  substitute them with a string of _POLYA_CHR or _POLYT_CHR
*/

#define _POLYA_CHR ('*')
#define _POLYT_CHR ('#')
#define _POLYA_MIN_LEN (14)
#define _POLYA_MIN_FRACTION (0.72)

void polyAT_substitution(pEST_info);

void Ntails_removal(pEST_info);

