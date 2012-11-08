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

#include "est-factorizations.h"
#include "list.h"
#include "types.h"
#include "detect-polya.h"

#include "log.h"

//#define LOG_THRESHOLD LOG_LEVEL_TRACE

//Detect the polyA signal at the end of a composition

plist correct_composition_tail(plist factorization, char *genomic_sequence, char *est_sequence){
	my_assert(genomic_sequence != NULL);
	my_assert(est_sequence != NULL);
	my_assert(factorization != NULL);
	my_assert(!list_is_empty(factorization));

	pfactor tail=list_tail(factorization);

	size_t i=tail->EST_end+1;
	size_t j=tail->GEN_end+1;

	const size_t est_length=strlen(est_sequence);
	const size_t gen_length=strlen(genomic_sequence);

	while((i < est_length) &&
			(j < gen_length) &&
			(genomic_sequence[j] == est_sequence[i])) {
		DEBUG("Tail correction: added one base");
		i++;
		j++;
	}

	tail->EST_end=i-1;
	tail->GEN_end=j-1;

	return factorization;
}

//Per il momento il segnale di poliadenilazione non viene restituito (rientrava solo
//nel GTF).
//Viene restituito solo la polyA che per l'assemblaggio dei full-length.
bool detect_polyA_signal(plist factorization, char *genomic_sequence, char *est_sequence, bool *polyadenil){
	my_assert(genomic_sequence != NULL);
	my_assert(est_sequence != NULL);
	my_assert(factorization != NULL);
	my_assert(!list_is_empty(factorization));

	pfactor tail=list_tail(factorization);

	DEBUG("Correct factor %d-%d (%d-%d on EST)", tail->GEN_start, tail->GEN_end, tail->EST_start, tail->EST_end);

	size_t est_length=strlen(est_sequence);
	DEBUG("EST length %zu", est_length);
	char *est_cleavage_seq=real_substring(tail->EST_end+1, est_length-tail->EST_end-1, est_sequence);


	int i=0, matches=0;
	bool stop=false;
	size_t est_cleavage_seq_length=strlen(est_cleavage_seq);
	while(i < (int)est_cleavage_seq_length && stop == false){
		if(est_cleavage_seq[i] == 'a' || est_cleavage_seq[i] == 'A'){
			//XXX
			if(matches >= 8)
				stop=true;
			else{
				matches++;
				i++;
			}
		}
		else{
			//XXX
			if(matches >= 8)
				stop=true;
			else
				i=(int)est_cleavage_seq_length;
		}
	}

	pfree(est_cleavage_seq);

	*polyadenil=false;

	if(stop == true){
		//PAS su genomica
		//XXX
	  i=MAX(0, tail->GEN_end-39);
		while(i <= tail->GEN_end && *polyadenil == false){
			if(genomic_sequence[i] == 'a' || genomic_sequence[i] == 'A'){
				//XXX
				char *pas=real_substring(i, 6, genomic_sequence);
				if((!strcmp(pas, "aataaa") || !strcmp(pas, "AATAAA")) || (!strcmp(pas, "attaaa") || !strcmp(pas, "ATTAAA")))
					*polyadenil=true;
				pfree(pas);
			}
			i++;
		}
	}

	if(stop == true){
		//XXX
	  i=MAX(0, tail->GEN_end-9);
		matches=0;
		while(i <= tail->GEN_end+10 && stop == true && genomic_sequence[i]!='\0'){
			//XXX
			if(matches >= 6)
				stop=false;
			else{
				if(genomic_sequence[i] == 'a' || genomic_sequence[i] == 'A')
					matches++;
				else
					matches=0;

				i++;
			}
		}

		if(stop == true){
			i=tail->GEN_end+1;
			int count=0;
			//XXX
			while(i <= tail->GEN_end+10 && stop == true && genomic_sequence[i]!='\0'){
				//XXX
				if(count >= 7)
					stop=false;
				else{
					if(genomic_sequence[i] == 'a' || genomic_sequence[i] == 'A')
						count++;
					i++;
				}
			}
		}
	}

	return stop;
}
