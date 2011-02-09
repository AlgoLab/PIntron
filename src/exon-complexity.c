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
#include <ctype.h>
#include <math.h>

#include "exon-complexity.h"
#include "est-factorizations.h"

#include "log.h"

double dustScoreByLeftAndRight(char *genomic_sequence, int start, int end){
	my_assert(genomic_sequence != NULL);
	size_t gen_length=strlen(genomic_sequence);
	my_assert(start >= 0 && end < (int)gen_length);

	char *sequence=real_substring(start, end-start+1, genomic_sequence);
	double dust=dustScore(sequence);

	pfree(sequence);
	return dust;
}

double dustScore(char *sequence){
	my_assert(sequence != NULL);

	size_t length=strlen(sequence);

	if((int)length  <= 2)
		return 0.0;

	int *dinucleotide_freq=NPALLOC(int, 17);
	my_assert(dinucleotide_freq != NULL);

	int i;
	for(i=0; i<17; i++)
		dinucleotide_freq[i]=0;

	int running_count=0;
	for(i=0; i < (int)length-1; i++){
		int index=getDinucleotideIndex(sequence[i], sequence[i+1]);
		running_count+=dinucleotide_freq[index];
		dinucleotide_freq[index]++;
	}

	double dust=(10.0 * (double)running_count)/((double)(length-2));

	pfree(dinucleotide_freq);

	//Dust score wrt to the sequence length
	return dust/length;
}

int getDinucleotideIndex(char firstChar, char secondChar){

	if((firstChar == 'a' || firstChar == 'A') && (secondChar == 'a' || secondChar == 'A'))
		return 0;

	if((firstChar == 'a' || firstChar == 'A') && (secondChar == 'c' || secondChar == 'C'))
		return 1;

	if((firstChar == 'a' || firstChar == 'A') && (secondChar == 'g' || secondChar == 'G'))
		return 2;

	if((firstChar == 'a' || firstChar == 'A') && (secondChar == 't' || secondChar == 'T'))
		return 3;

	if((firstChar == 'c' || firstChar == 'C') && (secondChar == 'a' || secondChar == 'A'))
		return 4;

	if((firstChar == 'c' || firstChar == 'C') && (secondChar == 'c' || secondChar == 'C'))
		return 5;

	if((firstChar == 'c' || firstChar == 'C') && (secondChar == 'g' || secondChar == 'G'))
		return 6;

	if((firstChar == 'c' || firstChar == 'C') && (secondChar == 't' || secondChar == 'T'))
		return 7;

	if((firstChar == 'g' || firstChar == 'G') && (secondChar == 'a' || secondChar == 'A'))
		return 8;

	if((firstChar == 'g' || firstChar == 'G') && (secondChar == 'c' || secondChar == 'C'))
		return 9;

	if((firstChar == 'g' || firstChar == 'G') && (secondChar == 'g' || secondChar == 'G'))
		return 10;

	if((firstChar == 'g' || firstChar == 'G') && (secondChar == 't' || secondChar == 'T'))
		return 11;

	if((firstChar == 't' || firstChar == 'T') && (secondChar == 'a' || secondChar == 'A'))
		return 12;

	if((firstChar == 't' || firstChar == 'T') && (secondChar == 'c' || secondChar == 'C'))
		return 13;

	if((firstChar == 't' || firstChar == 'T') && (secondChar == 'g' || secondChar == 'G'))
		return 14;

	if((firstChar == 't' || firstChar == 'T') && (secondChar == 't' || secondChar == 'T'))
		return 15;

	return 16;
}
