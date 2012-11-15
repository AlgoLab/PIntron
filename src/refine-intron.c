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

#include "refine-intron.h"
#include "refine.h"
#include "est-factorizations.h"
#include "list.h"
#include "types.h"

#include "log.h"

//#define LOG_THRESHOLD LOG_LEVEL_TRACE

//Refining of 5' and 3' intron splice sites (every alignment is considered)
//bool multiple_refine_intron(pconfiguration config, pEST_info gen_info, pEST_info est_info, pfactor donor, pfactor acceptor, bool first_intron){
//}

//Refining of 5' and 3' intron splice sites (just one alignment is considered)
bool refine_intron(pconfiguration config, pEST_info gen_info, pEST_info est_info, pfactor donor, pfactor acceptor, bool first_intron){
	int suffpref_length_on_est=config->suffpref_length_on_est;
	int suffpref_length_for_intron=config->suffpref_length_for_intron;
	int suffpref_length_on_gen=config->suffpref_length_on_gen;

	my_assert(donor->EST_end < acceptor->EST_start);
	my_assert(donor->GEN_end < acceptor->GEN_start);

	int donor_suffix_left_on_gen=donor->GEN_start;
	if(donor->GEN_end-suffpref_length_on_gen+1 >= donor_suffix_left_on_gen)
		donor_suffix_left_on_gen=donor->GEN_end-suffpref_length_on_gen+1;
	char *donor_suffix_on_gen=real_substring(donor_suffix_left_on_gen, donor->GEN_end-donor_suffix_left_on_gen+1, gen_info->EST_seq);
	TRACE("\t...donor exon suffix: %s (%zd)", donor_suffix_on_gen, strlen(donor_suffix_on_gen));

	int donor_suffix_left_on_est=donor->EST_start;
	if(donor->EST_end-suffpref_length_on_est+1 >= donor_suffix_left_on_est)
		donor_suffix_left_on_est=donor->EST_end-suffpref_length_on_est+1;
	char *donor_suffix_on_est=real_substring(donor_suffix_left_on_est, donor->EST_end-donor_suffix_left_on_est+1, est_info->EST_seq);
	TRACE("\t...donor est factor suffix: %s (%zd)", donor_suffix_on_est, strlen(donor_suffix_on_est));

	int acceptor_prefix_right_on_gen=acceptor->GEN_end;
	if(acceptor->GEN_start+suffpref_length_on_gen-1 <= acceptor_prefix_right_on_gen)
		acceptor_prefix_right_on_gen=acceptor->GEN_start+suffpref_length_on_gen-1;
	char *acceptor_prefix_on_gen=real_substring(acceptor->GEN_start, acceptor_prefix_right_on_gen-acceptor->GEN_start+1, gen_info->EST_seq);
	TRACE("\t...acceptor exon prefix: %s (%zd)", acceptor_prefix_on_gen, strlen(acceptor_prefix_on_gen));

	int acceptor_prefix_right_on_est=acceptor->EST_end;
	if(acceptor->EST_start+suffpref_length_on_est-1 <= acceptor_prefix_right_on_est)
		acceptor_prefix_right_on_est=acceptor->EST_start+suffpref_length_on_est-1;
	char *acceptor_prefix_on_est=real_substring(acceptor->EST_start, acceptor_prefix_right_on_est-acceptor->EST_start+1, est_info->EST_seq);
	TRACE("\t...acceptor est factor prefix: %s (%zd)", acceptor_prefix_on_est, strlen(acceptor_prefix_on_est));

	char *gap_on_est=NULL;
	if(donor->EST_end != acceptor->EST_start-1){
		gap_on_est=real_substring(donor->EST_end+1, acceptor->EST_start-donor->EST_end-1, est_info->EST_seq);
	}

	char *sequence_on_est=NPALLOC(char, strlen(donor_suffix_on_est)+strlen(acceptor_prefix_on_est)+1);
	strcpy(sequence_on_est, donor_suffix_on_est);
	if(gap_on_est != NULL)
		strcat(sequence_on_est, gap_on_est);
	strcat(sequence_on_est, acceptor_prefix_on_est);
	pfree(donor_suffix_on_est);
	pfree(acceptor_prefix_on_est);
	if(gap_on_est != NULL)
		pfree(gap_on_est);

	char *intron_prefix=real_substring(donor->GEN_end+1, suffpref_length_for_intron, gen_info->EST_seq);
	char *intron_suffix=real_substring(acceptor->GEN_start-suffpref_length_for_intron, suffpref_length_for_intron, gen_info->EST_seq);

	TRACE("\t...intron prefix: %s (%zd)", intron_prefix, strlen(intron_prefix));
	TRACE("\t...intron suffix: %s (%zd)", intron_suffix, strlen(intron_suffix));

	char *sequence_on_gen=NPALLOC(char, strlen(donor_suffix_on_gen)+strlen(intron_prefix)+strlen(intron_suffix)+strlen(acceptor_prefix_on_gen)+1);
	strcpy(sequence_on_gen, donor_suffix_on_gen);
	strcat(sequence_on_gen, intron_prefix);
	strcat(sequence_on_gen, intron_suffix);
	strcat(sequence_on_gen, acceptor_prefix_on_gen);
	pfree(donor_suffix_on_gen);
	pfree(acceptor_prefix_on_gen);
	pfree(intron_suffix);
	pfree(intron_prefix);

	int deleted_intron_dim=acceptor->GEN_start-donor->GEN_end-1-2*suffpref_length_for_intron;

	//Prima si tenta con un solo allineamento
	int gen_cut=0;
	int gen_cut_left=0;
	int gen_cut_right=0;

	plist alignments=compute_gap_alignment(sequence_on_est, sequence_on_gen, true, gen_cut, gen_cut_left, gen_cut_right);

 	pgap_alignment alignment=(pgap_alignment)list_head(alignments);

  	TRACE("\tAlignment matrix on gen: %s", alignment->GEN_gap_alignment);
 	TRACE("\tAlignment matrix on est: %s", alignment->EST_gap_alignment);

	alignment->new_acceptor_factor_left=donor_suffix_left_on_est+alignment->factor_cut;
	alignment->new_donor_right_on_gen=donor_suffix_left_on_gen+alignment->intron_start-1;
	alignment->new_acceptor_left_on_gen=donor_suffix_left_on_gen+alignment->intron_end+deleted_intron_dim+1;

	if(alignment->new_acceptor_factor_left == donor->EST_start){
		DEBUG("     The donor exon should be attached to the acceptor exon in the gap alignment!");
		if(first_intron){
			DEBUG("     ...and the donor exon has been attached to the acceptor exon in the gap alignment, since it is not the first one!");
			acceptor->EST_start=alignment->new_acceptor_factor_left;
			acceptor->GEN_start=alignment->new_acceptor_left_on_gen;
			gap_alignments_destroy(alignments);
			return true;
		}
		else{
			DEBUG("     ...but the intron is not the first one. The intron will not be refined!");
			gap_alignments_destroy(alignments);
			return false;
		}
	}

	if(alignment->new_acceptor_left_on_gen-alignment->new_donor_right_on_gen < config->min_intron_length){
		DEBUG("     The intron may be too small and the intron will not be refined!");
		gap_alignments_destroy(alignments);
		return false;
	}

	int donor_right_shift=(alignment->new_donor_right_on_gen > donor->GEN_end)?(alignment->new_donor_right_on_gen-donor->GEN_end):(donor->GEN_end-alignment->new_donor_right_on_gen);
	int acceptor_left_shift=(alignment->new_acceptor_left_on_gen > acceptor->GEN_start)?(alignment->new_acceptor_left_on_gen-acceptor->GEN_start):(acceptor->GEN_start-alignment->new_acceptor_left_on_gen);
	if(donor_right_shift > 20 || acceptor_left_shift > 20){
		DEBUG("     The donor/acceptor shift may be too big and the intron will not be refined!");
		gap_alignments_destroy(alignments);
		return false;
	}

	TRACE("\tNew acceptor est factor left: %d", alignment->new_acceptor_factor_left);
	TRACE("\tNew donor exon right: %d", alignment->new_donor_right_on_gen);
	TRACE("\tNew acceptor exon left: %d", alignment->new_acceptor_left_on_gen);

	//ATTENZIONE!!! Mettere un controllo sui valori trovati!!!

	int left_cut_on_align=0, left_genomic_cut_dim=0, left_EST_cut_dim=0;
	int right_cut_on_align=0, right_genomic_cut_dim=0, right_EST_cut_dim=0;

	Find_ACCEPTOR_before_on_the_left(alignment, alignment->intron_start_on_align-1, &left_cut_on_align, &left_genomic_cut_dim, &left_EST_cut_dim, "GT");
	Find_AG_after_on_the_right(alignment, alignment->intron_end_on_align+1, &right_cut_on_align, &right_genomic_cut_dim, &right_EST_cut_dim);

	//Type intron flag:
	// 0: gt-ag con al piu' 2 errori
	// 1: gt-ag con al piu' 4 errori
	// 2: Burset a 0 errori
	// 3: non classificato
	// 4: gc-ag con al piu' 2 errori
	// 5: gc-ag con al piu' 4 errori

	//ATTENZIONE!!!
	//Sostituire l_exon_right con donor_right_on_gen
	//Sostituire r_exon_left con acceptor_left_on_gen
	//Sostituire factor_left con acceptor_factor_left

	int final_new_donor_right_on_gen=0, final_new_acceptor_left_on_gen=0, final_new_acceptor_factor_left=0;
	int shifted_donor_right_on_gen=0, shifted_acceptor_left_on_gen=0, shifted_acceptor_factor_left=0;
	int shifted_donor_right_on_gen_gt=0, shifted_acceptor_left_on_gen_gt=0, shifted_acceptor_factor_left_gt=0;
	int shiftedlr_donor_right_on_gen_gt=0, shiftedlr_acceptor_left_on_gen_gt=0, shiftedlr_acceptor_factor_left_gt=0;
	int shifted_donor_right_on_gen_gc=0, shifted_acceptor_left_on_gen_gc=0, shifted_acceptor_factor_left_gc=0;
	int shiftedlr_donor_right_on_gen_gc=0, shiftedlr_acceptor_left_on_gen_gc=0, shiftedlr_acceptor_factor_left_gc=0;

	if(left_genomic_cut_dim == 0 && right_genomic_cut_dim == 0){
		DEBUG("     The intron is already canonical");
		//Il taglio e' rimasto invariato perche' e' gia' gt-ag
		final_new_donor_right_on_gen=alignment->new_donor_right_on_gen;
		final_new_acceptor_left_on_gen=alignment->new_acceptor_left_on_gen;
		final_new_acceptor_factor_left=alignment->new_acceptor_factor_left;
	}
	else{
		bool error_shift_gt=Shift_right_to_left_1(est_info->EST_seq, gen_info->EST_seq, 2, alignment, &shifted_donor_right_on_gen_gt, &shifted_acceptor_left_on_gen_gt, &shifted_acceptor_factor_left_gt, "GT");
		if(error_shift_gt){
			DEBUG("     ...successful 3' to 5' shifting for GT-AG!");
			shifted_donor_right_on_gen=shifted_donor_right_on_gen_gt;
			shifted_acceptor_left_on_gen=shifted_acceptor_left_on_gen_gt;
			shifted_acceptor_factor_left=shifted_acceptor_factor_left_gt;
		}
		else{
			error_shift_gt=Shift_left_to_right_1(est_info->EST_seq, gen_info->EST_seq, 2, alignment, &shiftedlr_donor_right_on_gen_gt, &shiftedlr_acceptor_left_on_gen_gt, &shiftedlr_acceptor_factor_left_gt, "GT");
			if(error_shift_gt){
				DEBUG("     ...successful 5' to 3' shifting for GT-AG!");
				shifted_donor_right_on_gen=shiftedlr_donor_right_on_gen_gt;
				shifted_acceptor_left_on_gen=shiftedlr_acceptor_left_on_gen_gt;
				shifted_acceptor_factor_left=shiftedlr_acceptor_factor_left_gt;
			}
			else{
				bool error_shift_gc=Shift_right_to_left_2(est_info->EST_seq, gen_info->EST_seq, 2, alignment, &shifted_donor_right_on_gen_gc, &shifted_acceptor_left_on_gen_gc, &shifted_acceptor_factor_left_gc, "GC");
				if(error_shift_gc){
					DEBUG("     ...successfull 3' to 5' shifting for GC-AG!");
					shifted_donor_right_on_gen=shifted_donor_right_on_gen_gc;
					shifted_acceptor_left_on_gen=shifted_acceptor_left_on_gen_gc;
					shifted_acceptor_factor_left=shifted_acceptor_factor_left_gc;
				}
				else{
					error_shift_gc=Shift_left_to_right_2(est_info->EST_seq, gen_info->EST_seq, 2, alignment, &shiftedlr_donor_right_on_gen_gc, &shiftedlr_acceptor_left_on_gen_gc, &shiftedlr_acceptor_factor_left_gc, "GC");
					if(error_shift_gc){
						DEBUG("     ...successfull 5' to 3' shifting for GC-AG!");
						shifted_donor_right_on_gen=shiftedlr_donor_right_on_gen_gc;
						shifted_acceptor_left_on_gen=shiftedlr_acceptor_left_on_gen_gc;
						shifted_acceptor_factor_left=shiftedlr_acceptor_factor_left_gc;
					}
					else{
						int burset_acceptor_factor_left=alignment->new_acceptor_factor_left;
						int burset_donor_right_on_gen=alignment->new_donor_right_on_gen;
						int burset_acceptor_left_on_gen=alignment->new_acceptor_left_on_gen;
						Try_Burset_after_match(est_info->EST_seq, gen_info->EST_seq, &burset_acceptor_factor_left, &burset_donor_right_on_gen, &burset_acceptor_left_on_gen, donor->EST_start, acceptor->EST_end);
						shifted_donor_right_on_gen=burset_donor_right_on_gen;
						shifted_acceptor_left_on_gen=burset_acceptor_left_on_gen;
						shifted_acceptor_factor_left=burset_acceptor_factor_left;
					}
				}
			}
		}

		final_new_donor_right_on_gen=shifted_donor_right_on_gen;
		final_new_acceptor_left_on_gen=shifted_acceptor_left_on_gen;
		final_new_acceptor_factor_left=shifted_acceptor_factor_left;

		if(final_new_acceptor_left_on_gen > acceptor->GEN_end || final_new_donor_right_on_gen < donor->GEN_start){
			DEBUG("     Cannot refine intron!");
			pfree(sequence_on_est);
			pfree(sequence_on_gen);
			gap_alignments_destroy(alignments);
			return false;
		}
	}

	donor->GEN_end=final_new_donor_right_on_gen;
	acceptor->GEN_start=final_new_acceptor_left_on_gen;
	acceptor->EST_start=final_new_acceptor_factor_left;
	donor->EST_end=acceptor->EST_start-1;

	pfree(sequence_on_est);
	pfree(sequence_on_gen);

	gap_alignments_destroy(alignments);

	return true;
}

int Try_Burset_after_match(char *est_sequence, char *genomic_sequence, int *acceptor_factor_left, int *donor_right_on_gen, int *acceptor_left_on_gen, int shifting_donor_factor_left, int shifting_acceptor_factor_right){
	int shifting_acceptor_factor_left=*acceptor_factor_left;
	int shifting_acceptor_left_on_gen=*acceptor_left_on_gen;
	int shifting_donor_right_on_gen=*donor_right_on_gen;
	int upd_acceptor_factor_left=0;
	int upd_acceptor_left_on_gen=0;
	int upd_donor_right_on_gen=0;
	int frequency=0, tmp_freq=0;
	bool right_to_left=false;

	DEBUG("\tTry Burset pattern");

	upd_acceptor_factor_left=shifting_acceptor_factor_left;
	upd_acceptor_left_on_gen=shifting_acceptor_left_on_gen;
	upd_donor_right_on_gen=shifting_donor_right_on_gen;

	bool stop=false;
	//Provo da sinistra a destra
	while((stop == false && est_sequence[shifting_acceptor_factor_left] == genomic_sequence[shifting_acceptor_left_on_gen]) && shifting_acceptor_factor_left > shifting_donor_factor_left+1){
		TRACE("\t...character %c (%d) of est in 5' matches to character %c (%d) of genomic in 3'", est_sequence[shifting_acceptor_factor_left], shifting_acceptor_factor_left, genomic_sequence[shifting_acceptor_left_on_gen], shifting_acceptor_left_on_gen);
		if(shifting_acceptor_factor_left == 0 || shifting_donor_right_on_gen == -1){
			TRACE("\t...the first est/genomic character is reached!");
			stop=true;
		}
		else{
			tmp_freq=Check_Burset_patterns(genomic_sequence, shifting_donor_right_on_gen, shifting_acceptor_left_on_gen);
			if(tmp_freq > frequency){
				frequency=tmp_freq;
				upd_acceptor_factor_left=shifting_acceptor_factor_left;
				upd_acceptor_left_on_gen=shifting_acceptor_left_on_gen;
				upd_donor_right_on_gen=shifting_donor_right_on_gen;
			}

			shifting_acceptor_factor_left--;
			shifting_donor_right_on_gen--;
			shifting_acceptor_left_on_gen--;
		}
	}

	shifting_acceptor_factor_left=*acceptor_factor_left;
	shifting_acceptor_left_on_gen=*acceptor_left_on_gen+1;
	shifting_donor_right_on_gen=*donor_right_on_gen+1;

	stop=false;
	//Provo da destra a sinistra
	while((stop == false && est_sequence[shifting_acceptor_factor_left] == genomic_sequence[shifting_donor_right_on_gen]) && shifting_acceptor_factor_left < shifting_acceptor_factor_right){
		TRACE("\t...character %c of est in 3' (%d) matches to character %c (%d) of genomic in 5'", est_sequence[shifting_acceptor_factor_left], shifting_acceptor_factor_left, genomic_sequence[shifting_donor_right_on_gen], shifting_donor_right_on_gen);
		if((unsigned int)shifting_acceptor_factor_left == strlen(est_sequence) || (unsigned int)shifting_acceptor_left_on_gen == strlen(genomic_sequence)){
			TRACE("\t...the last est/genomic character is reached!");
			stop=true;
		}
		else{
			tmp_freq=Check_Burset_patterns(genomic_sequence, shifting_donor_right_on_gen, shifting_acceptor_left_on_gen);
			if(tmp_freq > frequency){
				frequency=tmp_freq;
				upd_acceptor_factor_left=shifting_acceptor_factor_left;
				upd_acceptor_left_on_gen=shifting_acceptor_left_on_gen;
				upd_donor_right_on_gen=shifting_donor_right_on_gen;
				right_to_left=true;
			}

			shifting_acceptor_factor_left++;
			shifting_donor_right_on_gen++;
			shifting_acceptor_left_on_gen++;
		}
	}

	if(right_to_left == true)
		upd_acceptor_factor_left+=1;

	*acceptor_factor_left=upd_acceptor_factor_left;
	*donor_right_on_gen=upd_donor_right_on_gen;
	*acceptor_left_on_gen=upd_acceptor_left_on_gen;

	//Non si e' trovato nessun pattern di Burset e quindi si lascia il match di partenza per frequency=0
	return frequency;
}

//Returns the frequency of the pairs, or 0 it the pairs does not exist
int Check_Burset_patterns(char *genomic_sequence, int donor_left_on_gen, int acceptor_right_on_gen){
	int frequency=0;

	char *donor_pt=real_substring(donor_left_on_gen+1, 2, genomic_sequence);
	char *acceptor_pt=real_substring(acceptor_right_on_gen-2, 2, genomic_sequence);

	TRACE("\t...pattern (%d-%d) %s-%s", donor_left_on_gen+1, acceptor_right_on_gen-2, donor_pt, acceptor_pt);

	frequency=getBursetFrequency(donor_pt, acceptor_pt);

	pfree(donor_pt);
	pfree(acceptor_pt);

	return frequency;
}

int getBursetFrequency_adaptor(const char* const t,
										 const size_t cut1, const size_t cut2) {
  if (cut2<2)
	 return 0;
  char donor[3], acceptor[3];
  donor[2]= acceptor[2]= '\0';
  donor[0]= t[cut1];
  donor[1]= t[cut1+1];
  acceptor[0]= t[cut2-2];
  acceptor[1]= t[cut2-1];
  TRACE("Get Burset frequency of intron %s-%s.", donor, acceptor);
  return getBursetFrequency(donor, acceptor);
}

int getBursetFrequency(char *donor_pt, char *acceptor_pt){

	char *up_donor_pt=To_upper(donor_pt);
	char *up_acceptor_pt=To_upper(acceptor_pt);

	if(!strcmp(up_donor_pt, "AA") && !strcmp(up_acceptor_pt, "AG"))
	 return 1;

	if(!strcmp(up_donor_pt, "AA") && !strcmp(up_acceptor_pt, "AT"))
	 return 1;

	if(!strcmp(up_donor_pt, "AA") && !strcmp(up_acceptor_pt, "GT"))
	 return 1;

	if(!strcmp(up_donor_pt, "AC") && !strcmp(up_acceptor_pt, "CC"))
	 return 1;

	if(!strcmp(up_donor_pt, "AG") && !strcmp(up_acceptor_pt, "AC"))
	 return 1;

	if(!strcmp(up_donor_pt, "AG") && !strcmp(up_acceptor_pt, "AG"))
	 return 5;

	if(!strcmp(up_donor_pt, "AG") && !strcmp(up_acceptor_pt, "CT"))
	 return 2;

	if(!strcmp(up_donor_pt, "AG") && !strcmp(up_acceptor_pt, "GC"))
	 return 1;

	if(!strcmp(up_donor_pt, "AG") && !strcmp(up_acceptor_pt, "TG"))
	 return 2;

	if(!strcmp(up_donor_pt, "AT") && !strcmp(up_acceptor_pt, "AA"))
	 return 1;

	if(!strcmp(up_donor_pt, "AT") && !strcmp(up_acceptor_pt, "AC"))
	 return 8;

	if(!strcmp(up_donor_pt, "AT") && !strcmp(up_acceptor_pt, "AG"))
	 return 7;

	if(!strcmp(up_donor_pt, "AT") && !strcmp(up_acceptor_pt, "AT"))
		return 2;

	if(!strcmp(up_donor_pt, "AT") && !strcmp(up_acceptor_pt, "GC"))
		return 1;

	if(!strcmp(up_donor_pt, "AT") && !strcmp(up_acceptor_pt, "GT"))
		return 1;

	if(!strcmp(up_donor_pt, "CA") && !strcmp(up_acceptor_pt, "AG"))
		return 1;

	if(!strcmp(up_donor_pt, "CA") && !strcmp(up_acceptor_pt, "TT"))
		return 1;

	if(!strcmp(up_donor_pt, "CC") && !strcmp(up_acceptor_pt, "AG"))
		return 2;

	if(!strcmp(up_donor_pt, "CG") && !strcmp(up_acceptor_pt, "AG"))
		return 1;

	if(!strcmp(up_donor_pt, "CG") && !strcmp(up_acceptor_pt, "CA"))
		return 1;

	if(!strcmp(up_donor_pt, "CT") && !strcmp(up_acceptor_pt, "AC"))
		return 2;

	if(!strcmp(up_donor_pt, "CT") && !strcmp(up_acceptor_pt, "CA"))
		return 1;

	if(!strcmp(up_donor_pt, "GA") && !strcmp(up_acceptor_pt, "AG"))
		return 8;

	if(!strcmp(up_donor_pt, "GA") && !strcmp(up_acceptor_pt, "GT"))
		return 1;

	if(!strcmp(up_donor_pt, "GA") && !strcmp(up_acceptor_pt, "TC"))
		return 1;

	if(!strcmp(up_donor_pt, "GA") && !strcmp(up_acceptor_pt, "TG"))
		return 1;

	if(!strcmp(up_donor_pt, "GC") && !strcmp(up_acceptor_pt, "AG"))
		return 126;

	if(!strcmp(up_donor_pt, "GC") && !strcmp(up_acceptor_pt, "GG"))
		return 1;

	if(!strcmp(up_donor_pt, "GC") && !strcmp(up_acceptor_pt, "TA"))
		return 1;

	if(!strcmp(up_donor_pt, "GG") && !strcmp(up_acceptor_pt, "AC"))
		return 1;

	if(!strcmp(up_donor_pt, "GG") && !strcmp(up_acceptor_pt, "AG"))
		return 11;

	if(!strcmp(up_donor_pt, "GG") && !strcmp(up_acceptor_pt, "CA"))
		return 1;

	if(!strcmp(up_donor_pt, "GG") && !strcmp(up_acceptor_pt, "GA"))
		return 2;

	if(!strcmp(up_donor_pt, "GG") && !strcmp(up_acceptor_pt, "TC"))
		return 2;

	if(!strcmp(up_donor_pt, "GT") && !strcmp(up_acceptor_pt, "AG"))
		return 200;

	if(!strcmp(up_donor_pt, "GT") && !strcmp(up_acceptor_pt, "AC"))
		return 4;

	if(!strcmp(up_donor_pt, "GT") && !strcmp(up_acceptor_pt, "AT"))
		return 2;

	if(!strcmp(up_donor_pt, "GT") && !strcmp(up_acceptor_pt, "CA"))
		return 9;

	if(!strcmp(up_donor_pt, "GT") && !strcmp(up_acceptor_pt, "CG"))
		return 4;

	if(!strcmp(up_donor_pt, "GT") && !strcmp(up_acceptor_pt, "CT"))
		return 3;

	if(!strcmp(up_donor_pt, "GT") && !strcmp(up_acceptor_pt, "GC"))
		return 1;

	if(!strcmp(up_donor_pt, "GT") && !strcmp(up_acceptor_pt, "GG"))
		return 10;

	if(!strcmp(up_donor_pt, "GT") && !strcmp(up_acceptor_pt, "GT"))
		return 1;

	if(!strcmp(up_donor_pt, "GT") && !strcmp(up_acceptor_pt, "TA"))
		return 7;

	if(!strcmp(up_donor_pt, "GT") && !strcmp(up_acceptor_pt, "TC"))
		return 2;

	if(!strcmp(up_donor_pt, "GT") && !strcmp(up_acceptor_pt, "TG"))
		return 8;

	if(!strcmp(up_donor_pt, "GT") && !strcmp(up_acceptor_pt, "TT"))
		return 2;

	if(!strcmp(up_donor_pt, "TA") && !strcmp(up_acceptor_pt, "AG"))
		return 6;

	if(!strcmp(up_donor_pt, "TA") && !strcmp(up_acceptor_pt, "CG"))
		return 1;

	if(!strcmp(up_donor_pt, "TA") && !strcmp(up_acceptor_pt, "TC"))
		return 1;

	if(!strcmp(up_donor_pt, "TC") && !strcmp(up_acceptor_pt, "AG"))
		return 1;

	if(!strcmp(up_donor_pt, "TC") && !strcmp(up_acceptor_pt, "GG"))
		return 1;

	if(!strcmp(up_donor_pt, "TG") && !strcmp(up_acceptor_pt, "AC"))
		return 1;

	if(!strcmp(up_donor_pt, "TG") && !strcmp(up_acceptor_pt, "AG"))
		return 7;

	if(!strcmp(up_donor_pt, "TG") && !strcmp(up_acceptor_pt, "GG"))
		return 2;

	if(!strcmp(up_donor_pt, "TT") && !strcmp(up_acceptor_pt, "AG"))
		return 5;

	if(!strcmp(up_donor_pt, "TT") && !strcmp(up_acceptor_pt, "AT"))
		return 1;

	if(!strcmp(up_donor_pt, "TT") && !strcmp(up_acceptor_pt, "GG"))
		return 1;

	return 0;
}

//Returns a list of gap alignments
//Arguments EST_seq and genomic_seq are the sequence regions to be aligned
plist compute_gap_alignment(char *EST_seq, char *genomic_seq, bool only_one_align, int gen_cut, int gen_cut_left, int gen_cut_right){
	size_t n=strlen(EST_seq);
	size_t m=strlen(genomic_seq);

	//Per ora computa sempre un solo allineamento
	only_one_align=true;
	gen_cut=0;
	gen_cut_left=0;
	gen_cut_right=0;
	//Per ora computa sempre un solo allineamento (DOPO TOGLIERE QUESTA PARTE!)

	plist alignments=list_create();

	char **Ldir=NPALLOC(char*, 4);
	char **Gdir=NPALLOC(char*, 4);
	char **Rdir=NPALLOC(char*, 4);

	my_assert(Ldir != NULL && Gdir != NULL && Rdir != NULL);

	unsigned int i;

	for(i=0; i<4; i++){
		Ldir[i]=NPALLOC(char, (n+1)*(m+1));
		Gdir[i]=NPALLOC(char, (n+1)*(m+1));
		Rdir[i]=NPALLOC(char, (n+1)*(m+1));
		memset(Ldir[i], 0, (n+1)*(m+1)*sizeof(char));
		memset(Gdir[i], 0, (n+1)*(m+1)*sizeof(char));
		memset(Rdir[i], 0, (n+1)*(m+1)*sizeof(char));
	}


	char start_matrix;

	//Per ora fare calcolare un solo allineamento
	if(only_one_align){
		ComputeGapAlignMatrix(EST_seq, genomic_seq, Ldir, Gdir, Rdir, &start_matrix, only_one_align);
		pgap_alignment alignment=gap_alignment_create(n+m+10);
		TracebackGapAlignment(m, alignment, EST_seq, genomic_seq, Ldir, Gdir, Rdir, (int)n, (int)m, start_matrix);

		alignment->EST_gap_alignment[alignment->gap_alignment_dim]='\0';
		alignment->GEN_gap_alignment[alignment->gap_alignment_dim]='\0';

		list_add_to_tail(alignments, alignment);
	}
	else{
		//Computazione di piu' allineamenti ==> DA FARE
	}

	for(i=0; i<4; i++){
		pfree(Ldir[i]);
		pfree(Gdir[i]);
		pfree(Rdir[i]);
	}

	pfree(Ldir);
	pfree(Gdir);
	pfree(Rdir);

	return alignments;
}

//n: EST_seq length
//m: genomic_seq length
void ComputeGapAlignMatrix(char *EST_seq, char *genomic_seq, char **Ldir, char **Gdir, char **Rdir, char *start_matrix, bool only_one_align){
	int *L=NULL;
	int *G=NULL;
	int *R=NULL;

	unsigned int i, j;

	int i_del;

	my_assert(EST_seq != NULL);
	my_assert(genomic_seq != NULL);
	my_assert(Ldir != NULL);
	my_assert(Gdir != NULL);
	my_assert(Rdir != NULL);

	size_t n=strlen(EST_seq);
	size_t m=strlen(genomic_seq);

	L=NPALLOC(int, (n+1)*(m+1));
	G=NPALLOC(int, (n+1)*(m+1));
	R=NPALLOC(int, (n+1)*(m+1));
	my_assert(L != NULL && G != NULL && R != NULL);

	for(i=0; i<n+1; i++){
		for(j=0; j<m+1; j++){
			L[i*(m+1)+j]=0;
			G[i*(m+1)+j]=0;
			R[i*(m+1)+j]=0;
		}
	}

	//Casi base
	for(i=0; i<n+1; i++){
		Rdir[2][i*(m+1)]=1;
		Gdir[2][i*(m+1)]=1;
		Ldir[2][i*(m+1)]=1;
	}
	for(i=1; i<m+1; i++){
		Rdir[3][i]=1;
		Gdir[3][i]=1;
		Ldir[3][i]=1;
	}

	//Costruzione della matrice L
	for(i=1; i<n+1; i++){
		for(j=1; j<m+1; j++){

			L[i*(m+1)+j]=L[(i-1)*(m+1)+j-1];

			if(EST_seq[i-1] == genomic_seq[j-1] || EST_seq[i-1] == 'n' || genomic_seq[j-1] == 'n' || EST_seq[i-1] == 'N' || genomic_seq[j-1] == 'N')
				L[i*(m+1)+j]=L[i*(m+1)+j]+1;	//Costo match a +1
			else
				L[i*(m+1)+j]=L[i*(m+1)+j]-1;	//Costo mismatch a -1

			if(only_one_align)
				Ldir[0][i*(m+1)+j]=0;		//0 per allineamento caratteri in i-1 e j-1
			else
				Ldir[1][i*(m+1)+j]=1;		//0 per allineamento caratteri in i-1 e j-1


			if(L[i*(m+1)+j] < L[(i-1)*(m+1)+j]-1){
				L[i*(m+1)+j]= L[(i-1)*(m+1)+j]-1;	//Costo spazio a -1
				if(only_one_align)
					Ldir[0][i*(m+1)+j]=1;	//1 per cancellazione in genomic_seq; il carattere in i-1 matcha con -
				else{
					Ldir[2][i*(m+1)+j]=1;	//1 per cancellazione in genomic_seq; il carattere in i-1 matcha con -
					Ldir[1][i*(m+1)+j]=0;
				}
			}
			else{
				if(L[i*(m+1)+j] == L[(i-1)*(m+1)+j]-1 && !only_one_align)
					Ldir[2][i*(m+1)+j]=1;
			}

			if(L[i*(m+1)+j] < L[i*(m+1)+j-1]-1){
				L[i*(m+1)+j]= L[i*(m+1)+j-1]-1;	//Costo spazio a -1
				if(only_one_align)
					Ldir[0][i*(m+1)+j]=2;	//2 per cancellazione in EST_seq; il carattere in j-1 matcha con -
				else{
					Ldir[3][i*(m+1)+j]=1;	//2 per cancellazione in EST_seq; il carattere in j-1 matcha con -
					Ldir[1][i*(m+1)+j]=0;
					Ldir[2][i*(m+1)+j]=0;
				}
			}
			else{
				if(L[i*(m+1)+j] == L[i*(m+1)+j-1]-1 && !only_one_align)
					Ldir[3][i*(m+1)+j]=1;
			}
		}
	}

	//Costruzione della matrice G
	for(i=1; i<n+1; i++){
		for(j=1; j<m+1; j++){
			if(only_one_align)
				Gdir[0][i*(m+1)+j]=2;		//2 per cancellazione in EST_seq; il carattere in j-1 matcha con -
			else
				Gdir[3][i*(m+1)+j]=1;		//2 per cancellazione in EST_seq; il carattere in j-1 matcha con -

			G[i*(m+1)+j]=G[i*(m+1)+j-1];	//Costo gap a 0

			if(G[i*(m+1)+j] < L[i*(m+1)+j-1]){
				G[i*(m+1)+j]=L[i*(m+1)+j-1];	//Costo gap a 0
				if(only_one_align)
					Gdir[0][i*(m+1)+j]=-2;		//Salto alla L
				else{
					Gdir[0][i*(m+1)+j]=1;		//Salto alla L
					Gdir[3][i*(m+1)+j]=0;
				}
			}
			else{
				if(G[i*(m+1)+j] == L[i*(m+1)+j-1] && !only_one_align)
					Gdir[0][i*(m+1)+j]=1;		//Salto alla L
			}
		}
	}

	//Costruzione della matrice R
	for(i=1; i<n+1; i++){
		for(j=1; j<m+1; j++){

			R[i*(m+1)+j]=R[(i-1)*(m+1)+j-1];

			if(EST_seq[i-1] == genomic_seq[j-1] || EST_seq[i-1] == 'n' || genomic_seq[j-1] == 'n' || EST_seq[i-1] == 'N' || genomic_seq[j-1] == 'N')
				R[i*(m+1)+j]=R[i*(m+1)+j]+1;	//Costo match a +1
			else
				R[i*(m+1)+j]=R[i*(m+1)+j]-1;	//Costo mismatch a -1

			if(only_one_align)
				Rdir[0][i*(m+1)+j]=0;		//0 per allineamento caratteri in i-1 e j-1
			else
				Rdir[1][i*(m+1)+j]=1;		//0 per allineamento caratteri in i-1 e j-1

			if(i != n)
				i_del=R[i*(m+1)+j-1]-1;	//Costo spazio a -1
			else
				i_del=R[i*(m+1)+j-1];		//Costo gap a 0

			if(R[i*(m+1)+j] < i_del){
				R[i*(m+1)+j]=i_del;
				if(only_one_align)
					Rdir[0][i*(m+1)+j]=2;	//2 per cancellazione in EST_seq; il carattere in j-1 matcha con -
				else{
					Rdir[3][i*(m+1)+j]=1;	//2 per cancellazione in EST_seq; il carattere in j-1 matcha con -
					Rdir[1][i*(m+1)+j]=0;
				}
			}
			else{
				if(R[i*(m+1)+j] == i_del && !only_one_align)
					Rdir[3][i*(m+1)+j]=1;
			}

			if(R[i*(m+1)+j] < G[i*(m+1)+j-1]){
				R[i*(m+1)+j]=G[i*(m+1)+j-1];
				if(only_one_align)
					Rdir[0][i*(m+1)+j]=-2;	//2 per cancellazione in EST_seq; il carattere in j-1 matcha con -
				else{
					Rdir[0][i*(m+1)+j]=1;	//2 per cancellazione in EST_seq; il carattere in j-1 matcha con -
					Rdir[1][i*(m+1)+j]=0;
					Rdir[3][i*(m+1)+j]=0;
				}
			}
			else{
				if(R[i*(m+1)+j] == G[i*(m+1)+j-1] && !only_one_align)
					Rdir[0][i*(m+1)+j]=1;
			}


			if(R[i*(m+1)+j] < R[(i-1)*(m+1)+j]-1){
				R[i*(m+1)+j]=R[(i-1)*(m+1)+j]-1;	//Costo spazio a -1
				if(only_one_align)
					Rdir[0][i*(m+1)+j]=1;	//1 per cancellazione in genomic_seq; il carattere in i-1 matcha con -
				else{
					Rdir[2][i*(m+1)+j]=1;	//1 per cancellazione in genomic_seq; il carattere in i-1 matcha con -
					Rdir[1][i*(m+1)+j]=0;
					Rdir[3][i*(m+1)+j]=0;
				}
			}
			else{
				if(R[i*(m+1)+j] == R[(i-1)*(m+1)+j]-1 && !only_one_align)
					Rdir[2][i*(m+1)+j]=1;
			}
		}
	}

	if(R[n*(m+1)+m] >= G[n*(m+1)+m]){
		if(R[n*(m+1)+m] >= L[n*(m+1)+m])
			*start_matrix=2;
		else
			*start_matrix=0;
	}
	else{
		if(G[n*(m+1)+m] >= L[n*(m+1)+m])
			*start_matrix=1;
		else
			*start_matrix=0;
	}

	pfree(L);
	pfree(G);
	pfree(R);
}

//Argument alignment must be created before calling this procedure
//Arguments EST_seq and genomic_seq are the sequence regions to be aligned
void TracebackGapAlignment(size_t m, pgap_alignment alignment, char *EST_seq, char *genomic_seq, char **Ldir, char **Gdir, char **Rdir, int i, int j, char start_matrix){

	char direction=0;

	if(i > 0 && j > 0){
		direction=(start_matrix == 2)?(Rdir[0][i*(m+1)+j]):((start_matrix == 1)?(Gdir[0][i*(m+1)+j]):(Ldir[0][i*(m+1)+j]));

		if (direction == 0){
			TracebackGapAlignment(m, alignment, EST_seq, genomic_seq, Ldir, Gdir, Rdir, i-1, j-1, start_matrix);
			alignment->EST_gap_alignment[alignment->gap_alignment_dim]=EST_seq[i-1];
			alignment->GEN_gap_alignment[alignment->gap_alignment_dim]=genomic_seq[j-1];
			alignment->gap_alignment_dim=alignment->gap_alignment_dim+1;
		}
		else{
			if (direction == 1){
				TracebackGapAlignment(m, alignment, EST_seq, genomic_seq, Ldir, Gdir, Rdir, i-1, j, start_matrix);
				alignment->EST_gap_alignment[alignment->gap_alignment_dim]=EST_seq[i-1];
				alignment->GEN_gap_alignment[alignment->gap_alignment_dim]='-';
				alignment->gap_alignment_dim=alignment->gap_alignment_dim+1;
			}
			else{
				if(direction == -2){
					if(start_matrix == 2){
						alignment->intron_end=j-1;
						alignment->factor_cut=i;
					}
					else{
						alignment->intron_start=j-1;
					}
					start_matrix=start_matrix-1;
				}

				TracebackGapAlignment(m, alignment, EST_seq, genomic_seq, Ldir, Gdir, Rdir, i, j-1, start_matrix);
				alignment->EST_gap_alignment[alignment->gap_alignment_dim]='-';

				if(direction == -2){
					if(start_matrix == 1)
						alignment->intron_end_on_align=alignment->gap_alignment_dim;
					else
						alignment->intron_start_on_align=alignment->gap_alignment_dim;
				}
				alignment->GEN_gap_alignment[alignment->gap_alignment_dim]=genomic_seq[j-1];
				alignment->gap_alignment_dim=alignment->gap_alignment_dim+1;
			}
		}
	}
	else{
		if(i > 0){
			TracebackGapAlignment(m, alignment, EST_seq, genomic_seq, Ldir, Gdir, Rdir, i-1, j, start_matrix);
			alignment->EST_gap_alignment[alignment->gap_alignment_dim]=EST_seq[i-1];
			alignment->GEN_gap_alignment[alignment->gap_alignment_dim]='-';
			alignment->gap_alignment_dim=alignment->gap_alignment_dim+1;
		}
		else{
			if(j > 0){
				TracebackGapAlignment(m, alignment, EST_seq, genomic_seq, Ldir, Gdir, Rdir, i, j-1, start_matrix);
				alignment->EST_gap_alignment[alignment->gap_alignment_dim]='-';
				alignment->GEN_gap_alignment[alignment->gap_alignment_dim]=genomic_seq[j-1];
				alignment->gap_alignment_dim=alignment->gap_alignment_dim+1;
			}
		}
	}
}

void Find_AG_after_on_the_right(pgap_alignment alignment, int init, int *cut_on_align, int *genomic_cut_dim, int *EST_cut_dim){
	*cut_on_align=-1;
	*genomic_cut_dim=-1;
	*EST_cut_dim=-1;

	my_assert(init >= 2);
	my_assert(init >= alignment->intron_end_on_align+1);

	size_t index= init-2;

	bool stop=false;
	while ((!stop) &&
			 (index < strlen(alignment->GEN_gap_alignment)-1)) {
		while(alignment->GEN_gap_alignment[index] == '-')
			index++;
		char acceptor_pt[3];
		acceptor_pt[0]=alignment->GEN_gap_alignment[index];
		index++;
		while(alignment->GEN_gap_alignment[index] == '-')
			index++;
		acceptor_pt[1]=alignment->GEN_gap_alignment[index];
		acceptor_pt[2]='\0';
		const int acceptor_cmp=strcmp(acceptor_pt, "AG");
		stop= (acceptor_cmp == 0);
	}

	if(!stop)
		return;

	int cut_dim_on_gen=0;
	int cut_dim_on_EST=0;

	//Il taglio e' da init a index (compreso) sulla matrice di allineamento
	*cut_on_align= index+1;
	my_assert(alignment->intron_end_on_align >= -1);
	size_t i= alignment->intron_end_on_align + 1;
	while(i <= index){
		if(alignment->GEN_gap_alignment[i] != '-'){
			cut_dim_on_gen++;
		}
		if(alignment->EST_gap_alignment[i] != '-'){
			cut_dim_on_EST++;
		}
		i++;
	}

	*genomic_cut_dim=cut_dim_on_gen;
	*EST_cut_dim=cut_dim_on_EST;
}

void Find_ACCEPTOR_before_on_the_left(pgap_alignment alignment, int init, int *cut_on_align, int *genomic_cut_dim, int *EST_cut_dim, char *acceptor_str){
	*cut_on_align=-1;
	*genomic_cut_dim=-1;
	*EST_cut_dim=-1;

	int index=init+2;

	my_assert(init <= alignment->intron_start_on_align-1);

	bool stop=false;
	while(!stop && index > 0){
		while(alignment->GEN_gap_alignment[index] == '-')
			index--;
		char acceptor_pt[3];
		acceptor_pt[1]=alignment->GEN_gap_alignment[index];
		index--;
		while(index>=0 && alignment->GEN_gap_alignment[index] == '-')
			index--;
		if (index < 0) {
		  acceptor_pt[0]= '\0';
		} else {
		  acceptor_pt[0]=alignment->GEN_gap_alignment[index];
		}
		acceptor_pt[2]='\0';
		char acceptor_cmp=strcmp(acceptor_pt, acceptor_str);
		if(acceptor_cmp == 0)
			stop=true;
	}

	if(!stop)
		return;

	int cut_dim_on_gen=0;
	int cut_dim_on_EST=0;

	//Il taglio e' da init a index (compreso) sulla matrice di allineamento
	*cut_on_align=index-1;
	int i=alignment->intron_start_on_align-1;
	while(i >= index){
		if(alignment->GEN_gap_alignment[i] != '-')
			cut_dim_on_gen++;
		if(alignment->EST_gap_alignment[i] != '-')
			cut_dim_on_EST++;
		i--;
	}

	*genomic_cut_dim=cut_dim_on_gen;
	*EST_cut_dim=cut_dim_on_EST;
}

bool Shift_right_to_left_1(char *estseq, char *genseq, int cycle, pgap_alignment alignment, int *shifted_donor_right_on_gen, int *shifted_acceptor_left_on_gen, int *shifted_acceptor_factor_left, char *acceptor_str){
	int init_right=alignment->intron_end_on_align+1;
	int init_left=alignment->intron_start_on_align;

	int cut_on_align=0;
	int *genomic_cut_dim=NULL, *EST_cut_dim=NULL, *genomic_substr_dim=NULL;
	char **cut_factor=NULL, **match_str=NULL, **prev_match_str=NULL;

	genomic_cut_dim=NPALLOC(int, cycle);
	my_assert(genomic_cut_dim != NULL);

	EST_cut_dim=NPALLOC(int, cycle);
	my_assert(EST_cut_dim != NULL);

	genomic_substr_dim=NPALLOC(int, cycle);
	my_assert(genomic_substr_dim != NULL);

	cut_factor=NPALLOC(char*, cycle);
	my_assert(cut_factor != NULL);

	match_str=NPALLOC(char*, cycle);
	my_assert(match_str != NULL);

	prev_match_str=NPALLOC(char*, cycle);
	my_assert(prev_match_str != NULL);

	DEBUG("     Try shifting from 3' to 5' for searching %s-AG pattern (%d cycle(s))", acceptor_str, cycle);

	int i;
	for(i=0; i < cycle; i++){
		Find_AG_after_on_the_right(alignment, init_right, &cut_on_align, genomic_cut_dim+i, EST_cut_dim+i);

		if(EST_cut_dim[i] > -1){
			prev_match_str[i]=real_substring(alignment->new_acceptor_left_on_gen, genomic_cut_dim[i], genseq);
			cut_factor[i]=real_substring(alignment->new_acceptor_factor_left, EST_cut_dim[i], estseq);
			init_right=cut_on_align+1;
		}
		else{
			prev_match_str[i]=NULL;
			cut_factor[i]=NULL;
		}

		Find_ACCEPTOR_after_on_the_left(alignment, init_left, genomic_substr_dim+i, acceptor_str);

		if(genomic_substr_dim[i] > -1){
			match_str[i]=real_substring(alignment->new_donor_right_on_gen+1, genomic_substr_dim[i], genseq);
			init_left=alignment->intron_start_on_align+genomic_substr_dim[i]+1;
		}
		else
			match_str[i]=NULL;

		TRACE("\tcycle %d:", i+1);
		TRACE("\t\t3' genomic cut: %s", prev_match_str[i]);
		TRACE("\t\t3' est cut: %s", cut_factor[i]);
		TRACE("\t\t5' genomic matching: %s", match_str[i]);
	}

	unsigned int error=1000;
	unsigned int edit_prev=1000;
	int j;
	bool stop=false;
	i=0;
	while(i < cycle && !stop){
		j=0;
		while(j < cycle && !stop){
			if(cut_factor[i] != NULL && match_str[j] != NULL){
				size_t l1=strlen(cut_factor[i]);
				size_t l2=strlen(prev_match_str[i]);
				unsigned int* M=edit_distance(cut_factor[i], l1, prev_match_str[i], l2);
				edit_prev=M[(l1+1)*(l2+1)-1];
				pfree(M);
				if(edit_prev <= 5){
					l1=strlen(cut_factor[i]);
					l2=strlen(match_str[j]);
					M=edit_distance(cut_factor[i], l1, match_str[j], l2);
					error=M[(l1+1)*(l2+1)-1]-edit_prev;
					pfree(M);
				}
			}
			//XXX - ridiscutere questo parametro
			//if(error == 0){
			if(error <= 1){
				TRACE("\t3' good est cut: %d", EST_cut_dim[i]);
				TRACE("\t3' good genomic cut: %d", genomic_cut_dim[i]);
				TRACE("\t5' good genomic matching: %d", genomic_substr_dim[i]);

				*shifted_acceptor_factor_left=alignment->new_acceptor_factor_left+EST_cut_dim[i];
				*shifted_donor_right_on_gen=alignment->new_donor_right_on_gen+genomic_substr_dim[j];
				*shifted_acceptor_left_on_gen=alignment->new_acceptor_left_on_gen+genomic_cut_dim[i];
				stop=true;
			}
			j++;
		}
		i++;
	}

	for(i=0; i < cycle; i++){
		if(cut_factor[i] != NULL)
			pfree(cut_factor[i]);
		if(match_str[i] != NULL)
			pfree(match_str[i]);
		if(prev_match_str[i] != NULL)
			pfree(prev_match_str[i]);
	}

	pfree(cut_factor);
	pfree(match_str);
	pfree(prev_match_str);
	pfree(genomic_cut_dim);
	pfree(EST_cut_dim);
	pfree(genomic_substr_dim);

	//Solo se error e' = 0 si mantiene inalterato l'errore di allineamento
	return stop;
}

//Usata in ASPIC per GC-AG
bool Shift_right_to_left_2(char *estseq, char *genseq, int cycle, pgap_alignment alignment, int *shifted_donor_right_on_gen, int *shifted_acceptor_left_on_gen, int *shifted_acceptor_factor_left, char *acceptor_str){
	int init_right=alignment->intron_end_on_align+1;
	int init_left=alignment->intron_start_on_align;

	int cut_on_align=0;
	int *genomic_cut_dim=NULL, *EST_cut_dim=NULL, *genomic_substr_dim=NULL;
	char **cut_factor=NULL, **match_str=NULL;

	genomic_cut_dim=NPALLOC(int, cycle);
	EST_cut_dim=NPALLOC(int, cycle);
	genomic_substr_dim=NPALLOC(int, cycle);
	my_assert(genomic_cut_dim != NULL && EST_cut_dim != NULL && genomic_substr_dim != NULL);

	cut_factor=NPALLOC(char*, cycle);
	match_str=NPALLOC(char*, cycle);
	my_assert(cut_factor != NULL && match_str != NULL);

	DEBUG("     Try shifting from 3' to 5' for searching %s-AG pattern (%d cycle(s))", acceptor_str, cycle);

	int i;
	for(i=0; i < cycle; i++){
		Find_AG_after_on_the_right(alignment, init_right, &cut_on_align, genomic_cut_dim+i, EST_cut_dim+i);

		if(EST_cut_dim[i] > -1){
			cut_factor[i]=real_substring(alignment->new_acceptor_factor_left, EST_cut_dim[i], estseq);
			init_right=cut_on_align+1;
		}
		else
			cut_factor[i]=NULL;

		Find_ACCEPTOR_after_on_the_left(alignment, init_left, genomic_substr_dim+i, acceptor_str);

		if(genomic_substr_dim[i] > -1){
			match_str[i]=real_substring(alignment->new_donor_right_on_gen+1, genomic_substr_dim[i], genseq);
			init_left=alignment->intron_start_on_align+genomic_substr_dim[i]+1;
		}
		else
			match_str[i]=NULL;

		TRACE("\tcycle %d:", i+1);
		TRACE("\t\t3' est cut: %s", cut_factor[i]);
		TRACE("\t\t5' genomic matching: %s", match_str[i]);
	}

	i=0;
	int j;
	int error=1000, edit=1000;
	bool stop=false;
	while(i < cycle && !stop){
		j=0;
		while(j < cycle && !stop){
			if(cut_factor[i] != NULL && match_str[j] != NULL){
				size_t l1=strlen(cut_factor[i]);
				size_t l2=strlen(match_str[j]);
				unsigned int* M=edit_distance(cut_factor[i], l1, match_str[j], l2);
				edit=M[(l1+1)*(l2+1)-1];
				pfree(M);
			}
			else
				edit=1000;

			if(edit < error){
				TRACE("\t3' good est cut: %d", EST_cut_dim[i]);
				TRACE("\t3' good genomic cut: %d", genomic_cut_dim[i]);
				TRACE("\t5' good genomic matching: %d", genomic_substr_dim[i]);

				error=edit;
				*shifted_acceptor_factor_left=alignment->new_acceptor_factor_left+EST_cut_dim[i];
				*shifted_donor_right_on_gen=alignment->new_donor_right_on_gen+genomic_substr_dim[j];
				*shifted_acceptor_left_on_gen=alignment->new_acceptor_left_on_gen+genomic_cut_dim[i];
			}

			//XXX - ridiscutere questo parametro
			if(error == 0)
				stop=true;

			j++;
		}
		i++;
	}

	for(i=0; i < cycle; i++){
		if(cut_factor[i] != NULL)
			pfree(cut_factor[i]);
		if(match_str[i] != NULL)
			pfree(match_str[i]);
	}

	pfree(cut_factor);
	pfree(match_str);
	pfree(genomic_cut_dim);
	pfree(EST_cut_dim);
	pfree(genomic_substr_dim);

	return stop;
}

bool Shift_left_to_right_1(char *estseq, char *genseq, int cycle, pgap_alignment alignment, int *shifted_donor_right_on_gen, int *shifted_acceptor_left_on_gen, int *shifted_acceptor_factor_left, char *acceptor_str){
	int init_right=alignment->intron_end_on_align;
	int init_left=alignment->intron_start_on_align-1;

	int cut_on_align=0;
	int *genomic_cut_dim=NULL, *EST_cut_dim=NULL, *genomic_substr_dim=NULL;
	char **cut_factor=NULL, **match_str=NULL, **prev_match_str=NULL;

	genomic_cut_dim=NPALLOC(int, cycle);
	my_assert(genomic_cut_dim != NULL);

	EST_cut_dim=NPALLOC(int, cycle);
	my_assert(EST_cut_dim != NULL);

	genomic_substr_dim=NPALLOC(int, cycle);
	my_assert(genomic_substr_dim != NULL);

	cut_factor=NPALLOC(char*, cycle);
	my_assert(cut_factor != NULL);

	match_str=NPALLOC(char*, cycle);
	my_assert(match_str != NULL);

	prev_match_str=NPALLOC(char*, cycle);
	my_assert(prev_match_str != NULL);

	DEBUG("     Try shifting from 5' to 3' for searching %s-AG pattern (%d cycle(s))", acceptor_str, cycle);

	int i;
	for(i=0; i < cycle; i++){
		Find_ACCEPTOR_before_on_the_left(alignment, init_left, &cut_on_align, genomic_cut_dim+i, EST_cut_dim+i, acceptor_str);

		if(EST_cut_dim[i] > -1){
			prev_match_str[i]=real_substring(alignment->new_donor_right_on_gen-genomic_cut_dim[i]+1, genomic_cut_dim[i], genseq);
			cut_factor[i]=real_substring(alignment->new_acceptor_factor_left-EST_cut_dim[i], EST_cut_dim[i], estseq);
			init_left=cut_on_align-1;
		}
		else{
			prev_match_str[i]=NULL;
			cut_factor[i]=NULL;
		}

		Find_AG_before_on_the_right(alignment, init_right, genomic_substr_dim+i);

		if(genomic_substr_dim[i] > -1){
			match_str[i]=real_substring(alignment->new_acceptor_left_on_gen-genomic_substr_dim[i], genomic_substr_dim[i], genseq);
			init_right=alignment->intron_end_on_align-genomic_substr_dim[i]-1;
		}
		else
			match_str[i]=NULL;

		TRACE("\tcycle %d:", i+1);
		TRACE("\t\t5' genomic cut: %s", prev_match_str[i]);
		TRACE("\t\t5' est cut: %s", cut_factor[i]);
		TRACE("\t\t3' genomic matching: %s", match_str[i]);
	}

	int j;
	i=0;
	unsigned int error=1000;
	unsigned int edit_prev=1000;
	bool stop=false;
	while(i < cycle && !stop){
		j=0;
		while(j < cycle && !stop){
			if(cut_factor[i] != NULL && match_str[j] != NULL){
				size_t l1=strlen(cut_factor[i]);
				size_t l2=strlen(prev_match_str[i]);
				unsigned int *M=edit_distance(cut_factor[i], l1, prev_match_str[i], l2);
				edit_prev=M[(l1+1)*(l2+1)-1];
				pfree(M);
				if(edit_prev <= 5){
					l1=strlen(cut_factor[i]);
					l2=strlen(match_str[j]);
					M=edit_distance(cut_factor[i], l1, match_str[j], l2);
					error=M[(l1+1)*(l2+1)-1]-edit_prev;
					pfree(M);
				}
			}
			//XXX - ridiscutere questo parametro
			//if(error == 0){
			if(error <= 1){
				TRACE("\t5' good est cut: %d", EST_cut_dim[i]);
				TRACE("\t5' good genomic cut: %d", genomic_cut_dim[i]);
				TRACE("\t3' good genomic matching: %d", genomic_substr_dim[i]);

				*shifted_acceptor_factor_left=alignment->new_acceptor_factor_left-EST_cut_dim[i];
				*shifted_donor_right_on_gen=alignment->new_donor_right_on_gen-genomic_cut_dim[i];
				*shifted_acceptor_left_on_gen=alignment->new_acceptor_left_on_gen-genomic_substr_dim[j];
				stop=true;
			}
			j++;
		}
		i++;
	}

	for(i=0; i < cycle; i++){
		if(cut_factor[i] != NULL)
			pfree(cut_factor[i]);
		if(match_str[i] != NULL)
			pfree(match_str[i]);
		if(prev_match_str[i] != NULL)
			pfree(prev_match_str[i]);
	}
	pfree(cut_factor);
	pfree(match_str);
	pfree(prev_match_str);
	pfree(genomic_cut_dim);
	pfree(EST_cut_dim);
	pfree(genomic_substr_dim);

	//Solo se error e' = 0 si mantiene inalterato l'errore di allineamento
	return stop;
}

//Usata in ASPIC per GC-AG
bool Shift_left_to_right_2(char *estseq, char *genseq, int cycle, pgap_alignment alignment, int *shifted_donor_right_on_gen, int *shifted_acceptor_left_on_gen, int *shifted_acceptor_factor_left, char *acceptor_str){
	int init_right=alignment->intron_end_on_align;
	int init_left=alignment->intron_start_on_align-1;

	int cut_on_align=0;
	int *genomic_cut_dim=NULL, *EST_cut_dim=NULL, *genomic_substr_dim=NULL;
	char **cut_factor=NULL, **match_str=NULL;

	genomic_cut_dim=NPALLOC(int, cycle);
	EST_cut_dim=NPALLOC(int, cycle);
	genomic_substr_dim=NPALLOC(int, cycle);
	my_assert(genomic_cut_dim != NULL && EST_cut_dim != NULL && genomic_substr_dim != NULL);

	cut_factor=NPALLOC(char*, cycle);
	match_str=NPALLOC(char*, cycle);
	my_assert(cut_factor != NULL && match_str != NULL);

	DEBUG("     Try shifting from 5' to 3' for searching %s-AG pattern (%d cycle(s))", acceptor_str, cycle);

	int i;
	for(i=0; i < cycle; i++){
		Find_ACCEPTOR_before_on_the_left(alignment, init_left, &cut_on_align, genomic_cut_dim+i, EST_cut_dim+i, acceptor_str);

		if(EST_cut_dim[i] > -1){
			cut_factor[i]=real_substring(alignment->new_acceptor_factor_left-EST_cut_dim[i], EST_cut_dim[i], estseq);
			init_left=cut_on_align-1;
		}
		else
			cut_factor[i]=NULL;

		Find_AG_before_on_the_right(alignment, init_right, genomic_substr_dim+i);

		if(genomic_substr_dim[i] > -1){
			match_str[i]=real_substring(alignment->new_acceptor_left_on_gen-genomic_substr_dim[i], genomic_substr_dim[i], genseq);
			init_right=alignment->intron_end_on_align-genomic_substr_dim[i]-1;
		}
		else
			match_str[i]=NULL;

		TRACE("\tcycle %d:", i+1);
		TRACE("\t\t5' est cut: %s", cut_factor[i]);
		TRACE("\t\t3' genomic matching: %s", match_str[i]);
	}

	i=0;
	int j;
	bool stop=false;
	int error=1000, edit=1000;
	while(i < cycle && !stop){
		j=0;
		while(j < cycle && !stop){
			if(cut_factor[i] != NULL && match_str[j] != NULL){
				size_t l1=strlen(cut_factor[i]);
				size_t l2=strlen(match_str[j]);
				unsigned int *M=edit_distance(cut_factor[i], l1, match_str[j], l2);
				edit=M[(l1+1)*(l2+1)-1];
				pfree(M);
			}
			else
				edit=1000;

			if(edit < error){
				TRACE("\t5' good est cut: %d", EST_cut_dim[i]);
				TRACE("\t5' good genomic cut: %d", genomic_cut_dim[i]);
				TRACE("\t3' good genomic matching: %d", genomic_substr_dim[i]);

				error=edit;
				*shifted_acceptor_factor_left=alignment->new_acceptor_factor_left-EST_cut_dim[i];
				*shifted_donor_right_on_gen=alignment->new_donor_right_on_gen-genomic_cut_dim[i];
				*shifted_acceptor_left_on_gen=alignment->new_acceptor_left_on_gen-genomic_substr_dim[j];
			}

			//XXX - ridiscutere questo parametro
			if(error == 0)
				stop=true;

			j++;
		}
		i++;
	}

	for(i=0; i < cycle; i++){
		if(cut_factor[i] != NULL)
			pfree(cut_factor[i]);
		if(match_str[i] != NULL)
			pfree(match_str[i]);
	}

	pfree(cut_factor);
	pfree(match_str);
	pfree(genomic_cut_dim);
	pfree(EST_cut_dim);
	pfree(genomic_substr_dim);

	return stop;
}

void Find_ACCEPTOR_after_on_the_left(pgap_alignment alignment, int init, int *genomic_substr_dim, char *acceptor_str){
	*genomic_substr_dim=-1;

	int index=init;
	my_assert(init >= alignment->intron_start_on_align);

	bool stop=false;
	while(!stop && index < alignment->intron_end_on_align){
		char acceptor_pt[3];
		acceptor_pt[0]=alignment->GEN_gap_alignment[index];
		index++;
		acceptor_pt[1]=alignment->GEN_gap_alignment[index];
		acceptor_pt[2]='\0';
		char acceptor_cmp=strcmp(acceptor_pt, acceptor_str);
		if(acceptor_cmp == 0)
			stop=true;
	}

	if(!stop)
		return;

	*genomic_substr_dim=index-alignment->intron_start_on_align-1;
}

void Find_AG_before_on_the_right(pgap_alignment alignment, int init, int *genomic_substr_dim){

	*genomic_substr_dim=-1;

	int index=init;
	my_assert(init <= alignment->intron_end_on_align);

	bool stop=false;
	while(!stop && index > alignment->intron_start_on_align){
		char donor_pt[3];
		donor_pt[1]=alignment->GEN_gap_alignment[index];
		index--;
		donor_pt[0]=alignment->GEN_gap_alignment[index];
		donor_pt[2]='\0';
		char donor_cmp=strcmp(donor_pt, "AG");
		if(donor_cmp == 0)
			stop=true;
	}

	if(!stop)
		return;

	*genomic_substr_dim=alignment->intron_end_on_align-index-1;
}

char *To_lower(char *str){
	size_t length=strlen(str);
	unsigned int i=0;

	for(i=0; i<length; i++){
			str[i]=tolower(str[i]);
	}

	return str;
}

char *To_upper(char *str){
	size_t length=strlen(str);
	unsigned int i=0;

	for(i=0; i<length; i++){
			str[i]=toupper(str[i]);
	}

	return str;
}
