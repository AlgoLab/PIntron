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

#include "agree-introns.h"
#include "compute-alignments.h"
#include "types.h"
#include "list.h"
#include "est-factorizations.h"
#include "refine.h"
#include "refine-intron.h"
#include "classify-intron.h"

#include "log.h"

//#define LOG_THRESHOLD LOG_LEVEL_TRACE

bool try_agreement_to_intron_list(char *genomic_sequence, pintron intron_from, plist genomic_intron_list, int allowed_error){
	my_assert(genomic_intron_list != NULL);
	my_assert(genomic_sequence != NULL);
	my_assert(intron_from->isReal == true);
	my_assert(intron_from->agree_type != 0);
	my_assert(intron_from->try_agree == true);
	my_assert(intron_from->agreed == false);
	my_assert(intron_from->donor != NULL && intron_from->acceptor != NULL);
	my_assert(intron_from->gen_intron != NULL);

	bool agree_ok=false;
	plistit agree_it=list_first(genomic_intron_list);
	while(agree_ok == false && listit_has_next(agree_it)){
 		ppointer pp=(ppointer)listit_next(agree_it);
		pgenomic_intron gen_intron_to=(pgenomic_intron)pp->pointer;
		if(gen_intron_to->supportingESTs > 0)
			agree_ok=try_agreement(genomic_sequence, intron_from, gen_intron_to, allowed_error);
	}
	listit_destroy(agree_it);

	return agree_ok;
}

bool try_agreement_to_intron_list_on_single_site(char *genomic_sequence, pintron intron_from, plist gen_intron_list_to, plist gen_intron_list){
	my_assert(gen_intron_list_to != NULL);
	my_assert(genomic_sequence != NULL);
	my_assert(intron_from->isReal == true);
	my_assert(intron_from->agree_type != 0);
	my_assert(intron_from->try_agree == true);
	my_assert(intron_from->agreed == false);
	my_assert(intron_from->donor != NULL && intron_from->acceptor != NULL);
	my_assert(intron_from->gen_intron != NULL);

	bool agree_ok=false;
	plistit agree_it=list_first(gen_intron_list_to);
	while(agree_ok == false && listit_has_next(agree_it)){
 		ppointer pp=(ppointer)listit_next(agree_it);
		pgenomic_intron gen_intron_to=(pgenomic_intron)pp->pointer;
		if(gen_intron_to->supportingESTs > 0)
			agree_ok=try_agreement_on_single_site(genomic_sequence, intron_from, gen_intron_to, gen_intron_list);
	}
	listit_destroy(agree_it);

	return agree_ok;
}

bool try_agreement(char *genomic_sequence, pintron intron_from, pgenomic_intron gen_intron_to, int allowed_error){
	my_assert(intron_from->isReal == true);
	my_assert(intron_from->agree_type != 0);
	my_assert(intron_from->try_agree == true);
	my_assert(intron_from->agreed == false);
	my_assert(intron_from->donor != NULL && intron_from->acceptor != NULL);
	my_assert(intron_from->gen_intron != NULL);
	my_assert(gen_intron_to != NULL);

	int reducing_range=12; //From ASPIC parameters

	DEBUG("...to %d-%d (agree type %d)", gen_intron_to->start, gen_intron_to->end, gen_intron_to->agree_type);
	int start_diff=intron_from->gen_intron->start-gen_intron_to->start;
	int end_diff=intron_from->gen_intron->end-gen_intron_to->end;
	start_diff=(start_diff >= 0)?(start_diff):(-start_diff);
	end_diff=(end_diff >= 0)?(end_diff):(-end_diff);
	bool agree_ok=false;
	if(start_diff < reducing_range && end_diff < reducing_range){
		if(intron_from->donor->GEN_start < gen_intron_to->start && intron_from->acceptor->GEN_end > gen_intron_to->end){
			unsigned int error=get_agreement_error(genomic_sequence, intron_from, gen_intron_to);
			TRACE("...with agreeing error ==> %d", error);
			//XXX
			if((int)error <= allowed_error){ //ASPIC parameters ==> 2
				agree_ok=true;
				intron_from->agreed=true;
				intron_from->gen_intron->supportingESTs-=1;
				intron_from->gen_intron=gen_intron_to;
				intron_from->gen_intron->supportingESTs+=1;
				intron_from->donor->GEN_end=intron_from->gen_intron->start-1;
				intron_from->acceptor->GEN_start=intron_from->gen_intron->end+1;
				Correct_est_alignment(genomic_sequence, intron_from);
			}
		}
	}
	else{
		TRACE("...too much distant!");
	}

	 return agree_ok;
}

bool try_agreement_on_single_site(char *genomic_sequence, pintron intron_from, pgenomic_intron gen_intron_to, plist gen_intron_list){
	my_assert(intron_from->isReal == true);
	my_assert(intron_from->agree_type == 2);
	my_assert(intron_from->try_agree == true);
	my_assert(intron_from->agreed == false);
	my_assert(intron_from->donor != NULL && intron_from->acceptor != NULL);
	my_assert(intron_from->gen_intron != NULL);
	my_assert(gen_intron_to != NULL);

	DEBUG("...to %d-%d (agree type %d)", gen_intron_to->start, gen_intron_to->end, gen_intron_to->agree_type);

	int start_diff=intron_from->gen_intron->start-gen_intron_to->start;
	int end_diff=intron_from->gen_intron->end-gen_intron_to->end;
	start_diff=(start_diff >= 0)?(start_diff):(-start_diff);
	end_diff=(end_diff >= 0)?(end_diff):(-end_diff);

	int reducing_range=16; //ASPIC parameter: 12

	bool agree_ok=false;

	if(start_diff < reducing_range){
		DEBUG("...try agreement on donor site");
		agree_ok=try_agreement_on_donor_site(genomic_sequence, intron_from, gen_intron_to, gen_intron_list);
	}

	if(agree_ok == false && end_diff < reducing_range){
		DEBUG("...try agreement on acceptor site");
		agree_ok=try_agreement_on_acceptor_site(genomic_sequence, intron_from, gen_intron_to, gen_intron_list);
	}

	return agree_ok;
}

bool try_agreement_on_donor_site(char *genomic_sequence, pintron intron_from, pgenomic_intron gen_intron_to, plist gen_intron_list){
	my_assert(genomic_sequence != NULL);
	my_assert(intron_from != NULL);
	my_assert(gen_intron_list != NULL);
	my_assert(gen_intron_to != NULL);

	plist candidate_burset_list=list_create();

	int candidate_intron_start=gen_intron_to->start;
	bool eq_intron_start=(candidate_intron_start == (intron_from->gen_intron->start))?(true):(false);

	int reducing_range=16; //From ASPIC parameters

	int candidate_intron_end=intron_from->gen_intron->end-reducing_range;
	my_assert(candidate_intron_end > intron_from->gen_intron->start);
	int k=intron_from->gen_intron->end+reducing_range;
	if(k > intron_from->acceptor->GEN_end){
		k=intron_from->gen_intron->end+(intron_from->acceptor->GEN_end-intron_from->acceptor->GEN_start+1)/2;
		TRACE("...the acceptor exon may be too small!");
	}

	int current_freq=-1;
	if(eq_intron_start){
		my_assert(intron_from->gen_intron->burset_frequency != -1);
		current_freq=intron_from->gen_intron->burset_frequency;
		//current_freq=get_intron_Burset_frequency_start_end(genomic_sequence, candidate_intron_start, intron_from->gen_intron->end);
	}

	while(candidate_intron_end <= k){
		int frequency=get_intron_Burset_frequency_start_end(genomic_sequence, candidate_intron_start, candidate_intron_end);
		if(frequency > current_freq){
			pburset_frequency pbf=burset_frequency_create(candidate_intron_start, candidate_intron_end);
			pbf->frequency=frequency;
			list_add_to_tail(candidate_burset_list, pbf);
		}
		candidate_intron_end++;
	}

	list_sort(candidate_burset_list, (comparator)burset_frequency_compare);

	bool agree_ok=try_agreement_to_a_burset_frequency_list(genomic_sequence, intron_from, candidate_burset_list, gen_intron_list, 2);

	list_destroy(candidate_burset_list, (delete_function)burset_frequency_destroy);

	return agree_ok;
}

bool try_agreement_on_acceptor_site(char *genomic_sequence, pintron intron_from, pgenomic_intron gen_intron_to, plist gen_intron_list){
	my_assert(genomic_sequence != NULL);
	my_assert(intron_from != NULL);
	my_assert(gen_intron_list != NULL);
	my_assert(gen_intron_to != NULL);

	plist candidate_burset_list=list_create();

	int candidate_intron_end=gen_intron_to->end;
	bool eq_intron_end=(candidate_intron_end == (intron_from->gen_intron->end))?(true):(false);

	int reducing_range=16; //From ASPIC parameters

	int candidate_intron_start=intron_from->gen_intron->start-reducing_range;
	if(candidate_intron_start < intron_from->donor->GEN_start){
		candidate_intron_start=intron_from->gen_intron->start-(intron_from->donor->GEN_end-intron_from->donor->GEN_start+1)/2;
		TRACE("...the donor exon may be too small!");
	}
	int k=intron_from->gen_intron->start+reducing_range;
	my_assert(k < intron_from->gen_intron->end);

	int current_freq=-1;
	if(eq_intron_end){
		my_assert(intron_from->gen_intron->burset_frequency != -1);
		current_freq=intron_from->gen_intron->burset_frequency;
		//current_freq=get_intron_Burset_frequency_start_end(genomic_sequence, intron_from->gen_intron->start, candidate_intron_end);
	}

	while(candidate_intron_start <= k){
		int frequency=get_intron_Burset_frequency_start_end(genomic_sequence, candidate_intron_start, candidate_intron_end);
		if(frequency > current_freq){
			pburset_frequency pbf=burset_frequency_create(candidate_intron_start, candidate_intron_end);
			pbf->frequency=frequency;
			list_add_to_tail(candidate_burset_list, pbf);
		}
		candidate_intron_start++;
	}

	list_sort(candidate_burset_list, (comparator)burset_frequency_compare);

	bool agree_ok=try_agreement_to_a_burset_frequency_list(genomic_sequence, intron_from, candidate_burset_list, gen_intron_list, 2);

	list_destroy(candidate_burset_list, (delete_function)burset_frequency_destroy);

	return agree_ok;
}

bool find_better_intron(char *genomic_sequence, pintron intron_from, plist gen_intron_list){
	my_assert(intron_from != NULL);
	my_assert(gen_intron_list != NULL);
	my_assert(genomic_sequence != NULL);

	plist candidate_burset_list=list_create();

	int reducing_range=3;

 	int candidate_intron_start=intron_from->gen_intron->start-reducing_range;
	if(candidate_intron_start < intron_from->donor->GEN_start){
		candidate_intron_start=intron_from->gen_intron->start-(intron_from->donor->GEN_end-intron_from->donor->GEN_start+1)/2;
		TRACE("...the donor exon may be too small!");
	}

	int init_candidate_intron_end=intron_from->gen_intron->end-reducing_range;
	my_assert(init_candidate_intron_end > intron_from->gen_intron->start);

	int k_start=intron_from->gen_intron->start+reducing_range;
	my_assert(k_start < intron_from->gen_intron->end);

	int k_end=intron_from->gen_intron->end+reducing_range;
	if(k_end > intron_from->acceptor->GEN_end){
		k_end=intron_from->gen_intron->end+(intron_from->acceptor->GEN_end-intron_from->acceptor->GEN_start+1)/2;
		TRACE("...the acceptor exon may be too small!");
	}

	my_assert(intron_from->gen_intron->burset_frequency != -1);
	int current_freq=intron_from->gen_intron->burset_frequency;
	//int current_freq=get_intron_Burset_frequency_start_end(genomic_sequence, intron_from->gen_intron->start, intron_from->gen_intron->end);

 	while(candidate_intron_start <= k_start){
 		int candidate_intron_end=init_candidate_intron_end;
		while(candidate_intron_end <= k_end){
			int frequency=get_intron_Burset_frequency_start_end(genomic_sequence, candidate_intron_start, candidate_intron_end);
			if(frequency > current_freq){
				pburset_frequency pbf=burset_frequency_create(candidate_intron_start, candidate_intron_end);
				pbf->frequency=frequency;
				list_add_to_tail(candidate_burset_list, pbf);
			}
			candidate_intron_end++;
		}
		candidate_intron_start++;
	}

	list_sort(candidate_burset_list, (comparator)burset_frequency_compare);

	bool agree_ok=try_agreement_to_a_burset_frequency_list(genomic_sequence, intron_from, candidate_burset_list, gen_intron_list, 0);

	list_destroy(candidate_burset_list, (delete_function)burset_frequency_destroy);

	return agree_ok;
}

/*
 * allowed_error is the max error for canonical pt (for other pt must be 0)
 */
bool try_agreement_to_a_burset_frequency_list(char *genomic_sequence, pintron intron_from, plist burset_frequency_list, plist gen_intron_list, unsigned int allowed_error){

	bool agree_ok=false;
	plistit burset_it=list_first(burset_frequency_list);
	while(agree_ok == false && listit_has_next(burset_it)){
		pburset_frequency pbf=(pburset_frequency)listit_next(burset_it);
		unsigned int error=get_agreement_error_start_end(genomic_sequence, intron_from, pbf->start, pbf->end);
		DEBUG("...try to %d-%d (freq=%d, error=%d)", pbf->start, pbf->end, pbf->frequency, error);
		char *donor_burset_pt=real_substring(pbf->start, 2, genomic_sequence);
		char *acceptor_burset_pt=real_substring(pbf->end-1, 2, genomic_sequence);
		unsigned int max_error=allowed_error;

		if((strcmp(donor_burset_pt, "GT") && strcmp(donor_burset_pt, "gt")) && (strcmp(donor_burset_pt, "GC") && strcmp(donor_burset_pt, "gc"))){
			if(strcmp(donor_burset_pt, "AT") && strcmp(donor_burset_pt, "at")){
				max_error=0;
			}
			else{
				if(strcmp(acceptor_burset_pt, "AC") && strcmp(acceptor_burset_pt, "ac")){
					max_error=0;
				}
			}
		}
		else{
			if(strcmp(acceptor_burset_pt, "AG") && strcmp(acceptor_burset_pt, "ag")){
				max_error=0;
			}
		}
		pfree(donor_burset_pt);
		pfree(acceptor_burset_pt);

		if(intron_from->donor->GEN_start < pbf->start && intron_from->acceptor->GEN_end > pbf->end){
			if(error <= max_error){ //ASPIC parameters ==> 2
				agree_ok=true;
				intron_from->agreed=true;
				pgenomic_intron new_gen_intron;
				gen_intron_list=add_genomic_intron(genomic_sequence, gen_intron_list, pbf->start, pbf->end, &new_gen_intron);
				if(new_gen_intron->classified == false)
					new_gen_intron=classify_genomic_intron(genomic_sequence, new_gen_intron);
				intron_from->gen_intron->supportingESTs-=1;
				intron_from->gen_intron=new_gen_intron;
				intron_from->donor->GEN_end=intron_from->gen_intron->start-1;
				intron_from->acceptor->GEN_start=intron_from->gen_intron->end+1;
				Correct_est_alignment(genomic_sequence, intron_from);
			}
		}
	}
	listit_destroy(burset_it);

	return agree_ok;
}

pintron set_agree_flags(pintron intron){

	  my_assert(intron != NULL);
	  my_assert(intron->gen_intron != NULL);

	  intron->try_agree=true;
	  intron->agreed=false;
	  intron->agree_type=2;

	  if(intron->isReal==true){
		  my_assert(intron->gen_intron != NULL);
		  my_assert(intron->est_info != NULL);
		  my_assert(intron->gen_intron->donor_pt != NULL);
		  my_assert(intron->gen_intron->acceptor_pt != NULL);

		  if(intron->est_info->EST_gb[0] != 'N' || intron->est_info->EST_gb[1] != 'M'){
			  if((strcmp(intron->gen_intron->donor_pt, "gt") && strcmp(intron->gen_intron->donor_pt, "GT")) && (strcmp(intron->gen_intron->donor_pt, "gc") && strcmp(intron->gen_intron->donor_pt, "GC"))){
				  if(!strcmp(intron->gen_intron->donor_pt, "at") || !strcmp(intron->gen_intron->donor_pt, "AT")){
					  if(!strcmp(intron->gen_intron->acceptor_pt, "ac") || !strcmp(intron->gen_intron->acceptor_pt, "AC")){
						  if(intron->gen_intron->type != 2)
							  intron->agree_type=1;
					  }
				  }
			  }
			  else{
				  if(!strcmp(intron->gen_intron->acceptor_pt, "ag") || !strcmp(intron->gen_intron->acceptor_pt, "AG"))
					  intron->agree_type=1;
			  }
		  }
		  else{
			  intron->try_agree=false;
			  intron->agree_type=0;
		  }
	  }

	  return intron;
}

plist get_exon_composition_from_an_intron_composition(plist intron_composition){
	my_assert(intron_composition != NULL);
	my_assert(!list_is_empty(intron_composition));

	plist exon_composition=list_create();

	pintron head=list_remove_from_head(intron_composition);
	my_assert(head->donor == NULL);

	plistit intron_it=list_first(intron_composition);
	while(listit_has_next(intron_it)){
		pintron intron=(pintron)listit_next(intron_it);
		my_assert(intron->donor != NULL);
		list_add_to_tail(exon_composition, intron->donor);
	}
	listit_destroy(intron_it);

	return exon_composition;
}

plist get_intron_composition_from_an_exon_composition(pEST_info info, int gen_length, char *genomic_sequence, plist exon_composition, plist gen_intron_list, bool exon_ends_from_one){

	my_assert(genomic_sequence != NULL);
	my_assert(info != NULL);
	my_assert(info->EST_seq != NULL);
	my_assert(exon_composition != NULL);
	my_assert(list_size(exon_composition) >= 1);
	my_assert(gen_intron_list != NULL);

	plist intron_composition=list_create();

	size_t est_length=strlen(info->EST_seq);

	plistit exon_it=list_first(exon_composition);
	pfactor donor=NULL, acceptor=NULL;
	int start=-1, end;

	while(listit_has_next(exon_it)){
		acceptor=(pfactor)listit_next(exon_it);
		if(exon_ends_from_one){
			acceptor->EST_start=acceptor->EST_start-1;
			acceptor->EST_end=acceptor->EST_end-1;
			acceptor->GEN_start=acceptor->GEN_start-1;
			acceptor->GEN_end=acceptor->GEN_end-1;
		}

		end=acceptor->GEN_start-1;

		TRACE("Adding the genomic intron from %d to %d", start, end);

		pintron intron=intron_create();

		intron->donor=donor;
		intron->acceptor=acceptor;

		pgenomic_intron p=NULL;
		if(start >= 0 && end < (int)gen_length){
			TRACE("\tTRUE!");
			gen_intron_list=add_genomic_intron(genomic_sequence, gen_intron_list, start, end, &p);
			intron->isReal=true;
			my_assert(p != NULL);
		  }
		else{
			TRACE("\tFALSE!");
			p=genomic_intron_create(start, end);
			p->type=2;
			intron->isReal=false;
		}

		intron->gen_intron=p;
		intron->est_info=info;
		intron->try_agree=false;
		intron->agreed=false;

		TRACE("%s intron created!", (intron->isReal)?("TRUE"):("FALSE"));
		TRACE("\tfrom %d to %d on the genomic sequence", intron->gen_intron->start, intron->gen_intron->end);
		TRACE("\twith EST cut in %d", intron->acceptor->EST_start-1);
		TRACE("\twith exon start in donor in %d", (intron->donor==NULL)?(-1):(intron->donor->GEN_start));
		TRACE("\twith exon end in acceptor in %d", intron->acceptor->GEN_end);
		TRACE("\twith EST start in donor in %d", (intron->donor==NULL)?(-1):(intron->donor->EST_start));
		TRACE("\twith EST end in acceptor in %d", intron->acceptor->EST_end);

		TRACE("\twith donor pattern %s", intron->gen_intron->donor_pt);
		TRACE("\twith acceptor pattern %s", intron->gen_intron->acceptor_pt);

		list_add_to_tail(intron_composition, intron);

		start=acceptor->GEN_end+1;

		donor=acceptor;
 	  }

	  end=gen_length;

	  TRACE("Adding the genomic intron from %d to %d", start, end);
	  TRACE("\tFALSE!");

	  pintron last_intron=intron_create();
	  pgenomic_intron lastp=genomic_intron_create(start, end);
	  lastp->type=2;

	  last_intron->isReal=false;

	  last_intron->gen_intron=lastp;
	  last_intron->est_info=info;
	  last_intron->try_agree=false;
	  last_intron->agreed=false;

	  last_intron->donor=acceptor;
	  last_intron->acceptor=NULL;

	  TRACE("%s intron created!", (last_intron->isReal)?("TRUE"):("FALSE"));
	  TRACE("\tfrom %d to %d on the genomic sequence", last_intron->gen_intron->start, last_intron->gen_intron->end);
	  TRACE("\twith EST cut in %d", last_intron->donor->EST_end);
	  TRACE("\twith exon start in donor in %d", last_intron->donor->GEN_start);
	  TRACE("\twith exon end in acceptor in %d", gen_length);
	  TRACE("\twith EST start in donor in %d", last_intron->donor->EST_start);
	  TRACE("\twith EST end in acceptor in %zu", est_length);

	  TRACE("\twith donor pattern %s", last_intron->gen_intron->donor_pt);
	  TRACE("\twith acceptor pattern %s", last_intron->gen_intron->acceptor_pt);

	  list_add_to_tail(intron_composition, last_intron);

	  listit_destroy(exon_it);

	  return intron_composition;
}

plist add_genomic_intron(char *genomic_sequence, plist gen_intron_list, int start, int end, pgenomic_intron *p){

	my_assert(gen_intron_list != NULL);

	my_assert(genomic_sequence != NULL);
	size_t gen_length=strlen(genomic_sequence);

	my_assert(start >= 0 && start < (int)gen_length);
	my_assert(end >= 0 && end < (int)gen_length);

	plistit list_it_for_gen_intron=list_first(gen_intron_list);
	bool found=false;
	while(found == false && listit_has_next(list_it_for_gen_intron)){
		 pgenomic_intron current_gen_intron=(pgenomic_intron)listit_next(list_it_for_gen_intron);
		 if(current_gen_intron->start == start && current_gen_intron->end == end){
			 found=true;
			 current_gen_intron->supportingESTs=current_gen_intron->supportingESTs+1;
			 *p=current_gen_intron;
		 }
	}
	listit_destroy(list_it_for_gen_intron);

	if(found == false){
		pgenomic_intron gen_intron=genomic_intron_create(start, end);

		TRACE("\tNEW!");
		TRACE("\tfrom %d to %d", gen_intron->start, gen_intron->end);

		gen_intron=set_pattern(genomic_sequence, gen_intron);
		gen_intron=set_intron_Burset_frequency(gen_intron);
		my_assert(gen_intron->burset_frequency != -1);

		TRACE("\tdonor pattern ==> %s", gen_intron->donor_pt);
		TRACE("\tacceptor pattern ==> %s", gen_intron->acceptor_pt);

		gen_intron->supportingESTs=1;

		list_add_to_head(gen_intron_list, gen_intron);
		*p=gen_intron;
	}

	return gen_intron_list;
}

unsigned int get_agreement_error(char *genomic_sequence, pintron intron_from, pgenomic_intron gen_intron_to){

	my_assert(genomic_sequence != NULL);
	my_assert(intron_from != NULL);
	my_assert(gen_intron_to != NULL);
	my_assert(intron_from->isReal == true);
	my_assert(intron_from->donor != NULL && intron_from->acceptor != NULL);

	return get_agreement_error_start_end(genomic_sequence, intron_from, gen_intron_to->start, gen_intron_to->end);
}

unsigned int get_agreement_error_start_end(char *genomic_sequence, pintron intron_from, int gen_start, int gen_end){

	my_assert(genomic_sequence != NULL);
	my_assert(intron_from != NULL);
	my_assert(intron_from->isReal == true);
	my_assert(intron_from->donor != NULL && intron_from->acceptor != NULL);

	char *donor_seq_reduced;
	char *acceptor_seq_reduced;

	//If the start on the genomic sequence of the intron to be reduced is
	//greater than the start on the genomic sequence of the reduction intron
	if(intron_from->gen_intron->start > gen_start){
		int diff=intron_from->gen_intron->start-gen_start;

		int donor_EST_end=intron_from->donor->EST_end;
		int donor_EST_suffix_start=donor_EST_end-3*diff;
		if(donor_EST_suffix_start < intron_from->donor->EST_start)
			donor_EST_suffix_start=intron_from->donor->EST_start;

		char *donor_EST_suffix=real_substring(donor_EST_suffix_start, donor_EST_end-donor_EST_suffix_start+1, intron_from->est_info->EST_seq);

		int donor_GEN_end=intron_from->gen_intron->start-1;
		int donor_GEN_suffix_start=donor_GEN_end-3*diff;
		if(donor_GEN_suffix_start < intron_from->donor->GEN_start)
			donor_GEN_suffix_start=intron_from->donor->GEN_start;

		char *donor_GEN_suffix=real_substring(donor_GEN_suffix_start, donor_GEN_end-donor_GEN_suffix_start+1, genomic_sequence);

		plist alignments=compute_alignment(donor_EST_suffix, donor_GEN_suffix, true);

		palignment alignment=list_head(alignments);

		donor_seq_reduced=NPALLOC(char, alignment->alignment_dim+1);

		int i=0, j=0, k=1;

		my_assert(alignment->alignment_dim<=strlen(alignment->EST_alignment));
		my_assert(alignment->alignment_dim<=strlen(alignment->GEN_alignment));
		while ((i < alignment->alignment_dim) &&
				 (k <= diff)) {
		  if (alignment->EST_alignment[alignment->alignment_dim-i-1] != '-') {
			 donor_seq_reduced[j]=alignment->EST_alignment[alignment->alignment_dim-i-1];
			 j++;
		  }
		  if (alignment->GEN_alignment[alignment->alignment_dim-i-1] != '-') {
			 k++;
		  }
		  i++;
		}

		donor_seq_reduced[j]='\0';

		i=0;
		size_t dsr_length=strlen(donor_seq_reduced);

		while(i < (int)dsr_length/2){
			char help=donor_seq_reduced[dsr_length-i-1];
			donor_seq_reduced[dsr_length-i-1]=donor_seq_reduced[i];
			donor_seq_reduced[i]=help;
			i++;
		}

		alignments_destroy(alignments);

		pfree(donor_EST_suffix);
		pfree(donor_GEN_suffix);
	}
	else{
		donor_seq_reduced=NPALLOC(char, 1);
		my_assert(donor_seq_reduced != NULL);
		donor_seq_reduced[0]='\0';
	}

	char *donor_seq_reducing=real_substring(intron_from->gen_intron->start, (gen_start > intron_from->gen_intron->start)?(gen_start-intron_from->gen_intron->start):(0), genomic_sequence);

	if(intron_from->gen_intron->end < gen_end){
		int diff=gen_end-intron_from->gen_intron->end;

		int acceptor_EST_start=intron_from->acceptor->EST_start;
		int acceptor_EST_prefix_end=acceptor_EST_start+3*diff;
		if(acceptor_EST_prefix_end > intron_from->acceptor->EST_end)
			acceptor_EST_prefix_end=intron_from->acceptor->EST_end;

		char *acceptor_EST_prefix=real_substring(acceptor_EST_start, acceptor_EST_prefix_end-acceptor_EST_start+1, intron_from->est_info->EST_seq);

		int acceptor_GEN_start=intron_from->gen_intron->end+1;
		int acceptor_GEN_prefix_end=acceptor_GEN_start+3*diff;
		if(acceptor_GEN_prefix_end > intron_from->acceptor->GEN_end)
			acceptor_GEN_prefix_end=intron_from->acceptor->GEN_end;

		char *acceptor_GEN_prefix=real_substring(acceptor_GEN_start, acceptor_GEN_prefix_end-acceptor_GEN_start+1, genomic_sequence);

		plist alignments=compute_alignment(acceptor_EST_prefix, acceptor_GEN_prefix, true);

		palignment alignment=list_head(alignments);

		acceptor_seq_reduced=NPALLOC(char, alignment->alignment_dim+1);

		int i=0, j=0, k=1;

		my_assert(alignment->alignment_dim<=strlen(alignment->EST_alignment));
		my_assert(alignment->alignment_dim<=strlen(alignment->GEN_alignment));
		while ((i < alignment->alignment_dim) &&
				 (k <= diff)) {
		  if (alignment->EST_alignment[i] != '-') {
			 acceptor_seq_reduced[j]=alignment->EST_alignment[i];
			 j++;
		  }
		  if (alignment->GEN_alignment[i] != '-') {
			 k++;
		  }
		  i++;
		}
		acceptor_seq_reduced[j]='\0';

		alignments_destroy(alignments);

		pfree(acceptor_EST_prefix);
		pfree(acceptor_GEN_prefix);
	}
	else{
		acceptor_seq_reduced=NPALLOC(char, 1);
		my_assert(acceptor_seq_reduced != NULL);
		acceptor_seq_reduced[0]='\0';
	}

	DEBUG("Acceptor seq red >%s<", acceptor_seq_reduced);

	char *acceptor_seq_reducing=real_substring(gen_end+1, (intron_from->gen_intron->end > gen_end)?(intron_from->gen_intron->end-gen_end):(0), genomic_sequence);

	size_t reduced_donor_length=strlen(donor_seq_reduced);
	size_t reduced_acceptor_length=strlen(acceptor_seq_reduced);

	char *seq_reduced=NPALLOC(char, reduced_donor_length+reduced_acceptor_length+1);
	my_assert(seq_reduced != NULL);

	strcpy(seq_reduced, donor_seq_reduced);
	strcat(seq_reduced, acceptor_seq_reduced);

	pfree(donor_seq_reduced);
	pfree(acceptor_seq_reduced);

	size_t reducing_donor_length=strlen(donor_seq_reducing);
	size_t reducing_acceptor_length=strlen(acceptor_seq_reducing);

	char *seq_reducing=NPALLOC(char, reducing_donor_length+reducing_acceptor_length+1);
	my_assert(seq_reducing != NULL);

	strcpy(seq_reducing, donor_seq_reducing);
	strcat(seq_reducing, acceptor_seq_reducing);

	free(donor_seq_reducing);
	free(acceptor_seq_reducing);

	size_t l1=strlen(seq_reduced);
	size_t l2=strlen(seq_reducing);

	unsigned int* M=edit_distance(seq_reduced, l1, seq_reducing, l2);

	unsigned int error=M[(l1+1)*(l2+1)-1];
	pfree(M);
	return error;
}

void Correct_est_alignment(char *genomic_sequence, pintron intron){

	my_assert(genomic_sequence != NULL);
	my_assert(intron != NULL);
	my_assert(intron->gen_intron != NULL);
	my_assert(intron->est_info != NULL);
	my_assert(intron->donor != NULL && intron->acceptor);

	int est_suffix_dim=15, est_prefix_dim=15; //From ASPIC parameters
	int gen_suffix_dim=20, gen_prefix_dim=20; //From ASPIC parameters

	int donor_EST_start=intron->donor->EST_start;
	int donor_EST_end=intron->donor->EST_end;

	int donor_EST_suffix_start=donor_EST_end-est_suffix_dim;
	if(donor_EST_suffix_start < donor_EST_start)
		donor_EST_suffix_start=donor_EST_start;
	int donor_suffix_dim=donor_EST_end-donor_EST_suffix_start+1;
	char *donor_EST_factor=real_substring(donor_EST_suffix_start, donor_EST_end-donor_EST_suffix_start+1, intron->est_info->EST_seq);

	int acceptor_EST_start=intron->acceptor->EST_start;
	int acceptor_EST_end=intron->acceptor->EST_end;

	int acceptor_EST_prefix_end=acceptor_EST_start+est_prefix_dim;
	if(acceptor_EST_prefix_end > acceptor_EST_end)
		acceptor_EST_prefix_end=acceptor_EST_end;
	char *acceptor_EST_factor=real_substring(acceptor_EST_start, acceptor_EST_prefix_end-acceptor_EST_start+1, intron->est_info->EST_seq);

	int donor_GEN_start=intron->donor->GEN_start;
	int donor_GEN_end=intron->donor->GEN_end;

	int donor_GEN_suffix_start=donor_GEN_end-gen_suffix_dim;
	if(donor_GEN_suffix_start < donor_GEN_start)
		donor_GEN_suffix_start=donor_GEN_start;
	char *donor_GEN_factor=real_substring(donor_GEN_suffix_start, donor_GEN_end-donor_GEN_suffix_start+1, genomic_sequence);

	int acceptor_GEN_start=intron->acceptor->GEN_start;
	int acceptor_GEN_end=intron->acceptor->GEN_end;

	int acceptor_GEN_prefix_end=acceptor_GEN_start+gen_prefix_dim;
	if(acceptor_GEN_prefix_end > acceptor_GEN_end)
		acceptor_GEN_prefix_end=acceptor_GEN_end;
	char *acceptor_GEN_factor=real_substring(acceptor_GEN_start, acceptor_GEN_prefix_end-acceptor_GEN_start+1, genomic_sequence);

	char *intron_seq=NPALLOC(char, 30);
	strcpy(intron_seq, "xxxxxxxxxxxxxxxxxxxx");
	my_assert(intron_seq != NULL);

	size_t donor_GEN_length=strlen(donor_GEN_factor);
	size_t acceptor_GEN_length=strlen(acceptor_GEN_factor);
	size_t intron_length=strlen(intron_seq);

	char *gen_seq=NPALLOC(char, (donor_GEN_length+acceptor_GEN_length+intron_length+10));
	my_assert(gen_seq != NULL);

	strcpy(gen_seq, donor_GEN_factor);
	strcat(gen_seq, intron_seq);
	strcat(gen_seq, acceptor_GEN_factor);

	size_t donor_EST_length=strlen(donor_EST_factor);
	size_t acceptor_EST_length=strlen(acceptor_EST_factor);

	char *est_seq=NPALLOC(char, (donor_EST_length+acceptor_EST_length+10));
	my_assert(est_seq != NULL);

	strcpy(est_seq, donor_EST_factor);
	strcat(est_seq, acceptor_EST_factor);

	plist alignments=compute_gap_alignment(est_seq, gen_seq, 1, 0, 0, 0);
	pgap_alignment alignment=(pgap_alignment)list_head(alignments);

	int new_donor_EST_end=donor_EST_end-donor_suffix_dim+alignment->factor_cut;

	intron->donor->EST_end=new_donor_EST_end;
	intron->acceptor->EST_start=new_donor_EST_end+1;
	my_assert(intron->donor->EST_end > intron->donor->EST_start);
	my_assert(intron->acceptor->EST_end > intron->acceptor->EST_start);

	gap_alignments_destroy(alignments);

	pfree(est_seq);
	pfree(gen_seq);
	pfree(intron_seq);
	pfree(donor_EST_factor);
	pfree(acceptor_EST_factor);
	pfree(donor_GEN_factor);
	pfree(acceptor_GEN_factor);
}

int get_intron_Burset_frequency(char *genomic_sequence, pgenomic_intron gen_intron){
	my_assert(gen_intron != NULL);
	my_assert(genomic_sequence != NULL);

	return get_intron_Burset_frequency_start_end(genomic_sequence, gen_intron->start, gen_intron->end);
}

pgenomic_intron set_intron_Burset_frequency(pgenomic_intron gen_intron){
	my_assert(gen_intron != NULL);
	my_assert(gen_intron->donor_pt != NULL && gen_intron->acceptor_pt != NULL);

	gen_intron->burset_frequency=getBursetFrequency(gen_intron->donor_pt, gen_intron->acceptor_pt);

	return gen_intron;
}

int get_intron_Burset_frequency_start_end(char *genomic_sequence, int start, int end){
	my_assert(genomic_sequence != NULL);
	size_t gen_length=strlen(genomic_sequence);
	my_assert(start >= 0 && end < (int)gen_length);

	char *donor_pt=real_substring(start, 2, genomic_sequence);
	char *acceptor_pt=real_substring(end-1, 2, genomic_sequence);
	my_assert(donor_pt != NULL && acceptor_pt != NULL);

	int frequency=getBursetFrequency(donor_pt, acceptor_pt);

	pfree(donor_pt);
	pfree(acceptor_pt);

	return frequency;
}
