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
 * @file main.c
 *
 * Corpo principale del programma.
 *
 **/

#include <stdio.h>

#include "types.h"
#include "util.h"
#include "list.h"
#include "configuration.h"

#include "agree-introns.h"
#include "classify-intron.h"
#include "refine.h"

#include "io-multifasta.h"
#include "io-factorizations.h"
#include "est-factorizations.h"
#include "conversions.h"

#include "my_time.h"
#include "log.h"
#include "log-build-info.h"

int main(void) {
  INFO("INTRON-AGREEMENT v1");
  PRINT_LICENSE_INFORMATION;
  PRINT_SYSTEM_INFORMATION;
  DEBUG("Initialization");
  pmytime pt_tot= MYTIME_create_with_name("Total");
  pmytime pt_pre= MYTIME_create_with_name("Preprocessing");	//Per preparazione
  pmytime pt_alg= MYTIME_create_with_name("Algorithm");		//Per agreement
  pmytime pt_io= MYTIME_create_with_name("IO");

  MYTIME_start(pt_tot);
//Per ora i parametri sono tutti interni
//pconfiguration config= config_create(argc, argv);
  MYTIME_start(pt_io);

  char buf[1000];
  snprintf(buf, 1000, "info-pid-%u.log", (unsigned)getpid());

  FILE* floginfo= fopen(buf, "w");
  if (!floginfo) {
	 FATAL("Cannot create file info.log! Terminating");
	 fail();
  }

// Log resource utilization
  log_info(floginfo, "start");

  DEBUG("Reading genomic sequence");
  FILE* fgen= fopen("genomic.txt", "r");
  if (!fgen) {
	 FATAL("File genomic.txt not found! Terminating");
	 fail();
  }
  plist gen_list= read_multifasta(fgen);
  fclose(fgen);
  my_assert(list_size(gen_list)==1);
  pEST_info gen= (pEST_info)list_head(gen_list);
  parse_genomic_header(gen);

  list_destroy(gen_list, noop_free);

//In this file ESTs must be compatible with the genomic strand
  DEBUG("Reading EST sequences");
  FILE* fests= fopen("processed-ests.txt", "r");
  if (!fests) {
	 FATAL("File processed-ests.txt not found! Terminating");
	 fail();
  }
  plist estinfo_list= read_multifasta(fests);
  fclose(fests);
  DEBUG("EST sequences read");

/*plistit debug_it=list_first(estinfo_list);
  while(listit_has_next(debug_it)){
  pEST_info est_prova=(pEST_info)listit_next(debug_it);
  DEBUG(">%s", est_prova->EST_id);
  DEBUG("%s", est_prova->EST_seq);
  }
  listit_destroy(debug_it);*/

  DEBUG("Reading EST exon compositions");
  FILE* fcomp= fopen("out-agree.txt", "r");
  if (!fcomp) {
	 FATAL("File out-agree.txt not found! Terminating");
	 fail();
  }
//pEST list
  plist est_with_intron_list=read_factorizations(fcomp);

  fclose(fcomp);
  DEBUG("EST exon compositions read");

/*plistit debug_it0=list_first(est_with_intron_list);
  while(listit_has_next(debug_it0)){
  pEST est_prova=(pEST)listit_next(debug_it0);
  DEBUG(">%s", est_prova->info->EST_id);

  my_assert(est_prova->factorizations != NULL);

  if(list_size(est_prova->factorizations) != 1){
  DEBUG("...there are more than one composition and only the first one will be considered!");
  }
  plist comp_prova=list_head(est_prova->factorizations);
  plistit debug_it2=list_first(comp_prova);
  while(listit_has_next(debug_it2)){
  pfactor e=(pfactor)listit_next(debug_it2);
  DEBUG("Exon %d-%d (%d-%d)", e->GEN_start, e->GEN_end, e->EST_start, e->EST_end);
  }
  listit_destroy(debug_it2);
  }
  listit_destroy(debug_it0);*/

  FILE* f_multif_out= fopen("out-after-intron-agree.txt", "w");
  if (!f_multif_out) {
	 FATAL("Cannot create file out-after-intron-agree.txt! Terminating");
	 fail();
  }

  FILE* gtf_out= fopen("predicted-introns.txt", "w");
   if (!gtf_out) {
 	 FATAL("Cannot create file predicted-introns.gtf! Terminating");
 	 fail();
   }

// Log resource utilization
  log_info(floginfo, "data-io-end");

  MYTIME_stop(pt_io);

  MYTIME_start(pt_pre);

  size_t gen_length=strlen(gen->EST_seq);
  plist gen_intron_list=list_create();

  DEBUG("Building the intron compositions for each EST and the list of the genomic introns");

//Create a list of pEST objects for storing intron compositions
  plistit est_list_it=list_first(est_with_intron_list);
  while(listit_has_next(est_list_it)){
	 pEST est=(pEST)listit_next(est_list_it);

	 DEBUG("EST ==> %s", est->info->EST_id);

	 char *ID=est->info->EST_id;

	 bool found=false;

	 plistit estinfo_list_it=list_first(estinfo_list);
	 while(found == false && listit_has_next(estinfo_list_it)){
		pEST_info estinfo=(pEST_info)listit_next(estinfo_list_it);
		if(!strcmp(estinfo->EST_id, ID)){
		  EST_info_destroy(est->info);
		  set_EST_GB_identification(estinfo);
		  est->info=estinfo;
		  found=true;
		}
	 }

	 listit_destroy(estinfo_list_it);
	 my_assert(found == true);

	 my_assert(est->factorizations != NULL);
	 if(list_size(est->factorizations) != 1){
		DEBUG("...there are more than one composition and only the first one will be considered!");
	 }

//Transform the exon compositions (list of pfactor) into an intron composition
//(list of pintron) and build the list of genomic introns
	 plist exon_composition=(plist)list_head(est->factorizations);

	 my_assert(!list_is_empty(exon_composition));

	 plist intron_composition_list=list_create();

	 plist intron_composition=get_intron_composition_from_an_exon_composition(est->info, gen_length, gen->EST_seq, exon_composition, gen_intron_list, true);

	 DEBUG("Intron composition retrieved!");

/*plistit debug_it=list_first(intron_composition);
  while(listit_has_next(debug_it)){
  pintron i=(pintron)listit_next(debug_it);
  DEBUG("%s intron created!", (i->isReal)?("TRUE"):("FALSE"));
  DEBUG("\tfrom %d to %d on the genomic sequence", i->gen_intron->start, i->gen_intron->end);
  DEBUG("\twith EST cut in %d", i->EST_cut);
  DEBUG("\twith exon start in donor in %d", i->donor_GEN_start);
  DEBUG("\twith exon end in acceptor in %d", i->acceptor_GEN_end);
  DEBUG("\twith EST start in donor in %d", i->donor_EST_start);
  DEBUG("\twith EST end in acceptor in %d", i->acceptor_EST_end);
  DEBUG("\twith donor pattern %s", i->gen_intron->donor_pt);
  DEBUG("\twith acceptor pattern %s", i->gen_intron->acceptor_pt);
  }
  listit_destroy(debug_it);*/

	 DEBUG("Genomic intron list updated!");

/*plistit debug_gen_intron_it=list_first(gen_intron_list);
  while(listit_has_next(debug_gen_intron_it)){
  pgenomic_intron gi=(pgenomic_intron)listit_next(debug_gen_intron_it);
  DEBUG("\tfrom %d to %d on the genomic sequence", gi->start, gi->end);
  DEBUG("\twith donor pattern %s", gi->donor_pt);
  DEBUG("\twith acceptor pattern %s", gi->acceptor_pt);

  }
  listit_destroy(debug_gen_intron_it);*/

	 list_add_to_head(intron_composition_list, intron_composition);

	 plistit destroy_it=list_first(est->factorizations);
	 bool first_comp=true;
	 while(listit_has_next(destroy_it)){
		 plist exon_comp_to_be_destroyed=(plist)listit_next(destroy_it);
		 if(first_comp){
			 list_destroy(exon_comp_to_be_destroyed,(delete_function)noop_free);
			 first_comp=false;
		 }
		 else{
			 list_destroy(exon_comp_to_be_destroyed,(delete_function)factor_destroy);
		}
	 }
	 listit_destroy(destroy_it);

	 est->factorizations=intron_composition_list;

	 /*size_t est_length=strlen(est->info->EST_seq);
	 plist first_list=list_head(est->factorizations);
	 plistit debug_it=list_first(first_list);
	 while(listit_has_next(debug_it)){
		  pintron i=(pintron)listit_next(debug_it);
		  my_assert(i->donor != NULL || i->acceptor != NULL);

		  DEBUG("%s intron created!", (i->isReal)?("TRUE"):("FALSE"));
		  DEBUG("\tfrom %d to %d on the genomic sequence", i->gen_intron->start, i->gen_intron->end);

		  DEBUG("\twith EST cut in %d", (i->donor != NULL)?(i->donor->EST_end):(i->acceptor->EST_start-1));
		  DEBUG("\twith exon start in donor in %d", (i->donor == NULL)?(-1):(i->donor->GEN_start));
		  DEBUG("\twith exon end in acceptor in %d", (i->acceptor == NULL)?(gen_length):(i->acceptor->GEN_end));
		  DEBUG("\twith EST start in donor in %d", (i->donor == NULL)?(-1):(i->donor->EST_start));
		  DEBUG("\twith EST end in acceptor in %d", (i->acceptor == NULL)?(est_length):(i->acceptor->EST_end));

		  DEBUG("\twith donor pattern %s", i->gen_intron->donor_pt);
		  DEBUG("\twith acceptor pattern %s", i->gen_intron->acceptor_pt);
	 }
	 listit_destroy(debug_it);*/
  }
  listit_destroy(est_list_it);

  list_destroy(estinfo_list, (delete_function)noop_free);

  //Create a list of pEST objects for storing intron compositions
  gen_intron_list=classify_genomic_intron_list(gen->EST_seq, gen_intron_list);

  DEBUG("At total of %zu genomic introns was classified!", list_size(gen_intron_list));

  /*plistit debug_gen_intron_it=list_first(gen_intron_list);
  while(listit_has_next(debug_gen_intron_it)){
	  pgenomic_intron gi=(pgenomic_intron)listit_next(debug_gen_intron_it);
	  DEBUG("\tIntron %d-%d", gi->start+1, gi->end+1);
	  DEBUG("\t\twith pattern %s-%s", gi->donor_pt, gi->acceptor_pt);
	  DEBUG("\t\twith type %s", (gi->type == 0)?("U12"):((gi->type == 1)?("U2"):("UNCLASSIFIED")));
	  DEBUG("\t\twith donor-acceptor scores %f-%f", gi->score5, gi->score3);
	  DEBUG("\t\twith BPS %f", gi->BPS_score);
	  DEBUG("\t\twith BPS position %d", gi->BPS_position);
  }
  listit_destroy(debug_gen_intron_it);*/

//Set the try_agree flag for the introns and construction of the agreement list
  plist refseq_list=list_create();
  plist canonical_list=list_create();
  plist agreement_list=list_create();

  est_list_it=list_first(est_with_intron_list);
  while(listit_has_next(est_list_it)){
	 pEST est=(pEST)listit_next(est_list_it);

	 my_assert(est->factorizations != NULL);
	 my_assert(list_size(est->factorizations) == 1);

	 plist intron_composition=(plist)list_head(est->factorizations);
	 plistit intron_it=list_first(intron_composition);
	 while(listit_has_next(intron_it)){
		pintron intron=(pintron)listit_next(intron_it);
		set_agree_flags(intron);

		if(intron->agree_type <= intron->gen_intron->agree_type)
			intron->gen_intron->agree_type=intron->agree_type;

		DEBUG("The intron %d-%d on the genomic sequence linked to the EST %s", intron->gen_intron->start, intron->gen_intron->end, intron->est_info->EST_gb);
		if(intron->isReal== true){
			 ppointer pp=pointer_create();
			 pp->pointer=(pintron)intron;
			 if(intron->agree_type == 0){
				 my_assert(intron->try_agree == false);
				 //Add to the refseq list or to the canonical list
				 if(intron->agree_type == 0){
					 DEBUG("\twill be added to the REFSEQ list!");
					 list_add_to_tail(refseq_list, pp);
				 }
			 }
			 else{
				 my_assert(intron->try_agree == true);
				 //Add to the agree list
				 if(intron->agree_type == 1){
					 DEBUG("\twill be added to the canonical list!");
					 list_add_to_tail(canonical_list, pp);
				 }
				 else{
					 DEBUG("\twill be added to the agreement list!");
					 list_add_to_tail(agreement_list, pp);
				 }
			 }
		}
		else{
			DEBUG("\tis FALSE and will not be taken into account!");
		}
	 }
	 listit_destroy(intron_it);
  }
  listit_destroy(est_list_it);

  plist genomic_refseq_list=list_create();
  plist genomic_canonical_list=list_create();
  plist genomic_agreement_list=list_create();
  plistit build_gen_list_it=list_first(gen_intron_list);

  while(listit_has_next(build_gen_list_it)){
	  pgenomic_intron gi=(pgenomic_intron)listit_next(build_gen_list_it);
	  ppointer pp=pointer_create();
	  pp->pointer=(pgenomic_intron)gi;
	  if(gi->agree_type == 0)
		 list_add_to_tail(genomic_refseq_list, pp);
	  else{
		  if(gi->agree_type == 1)
			 list_add_to_tail(genomic_canonical_list, pp);
		  else
			 list_add_to_tail(genomic_agreement_list, pp);
	  }
  }
  listit_destroy(build_gen_list_it);

  /*plistit debug_gen_intron_it=list_first(gen_intron_list);
  while(listit_has_next(debug_gen_intron_it)){
	  pgenomic_intron gi=(pgenomic_intron)listit_next(debug_gen_intron_it);
	  DEBUG("\tIntron %d-%d", gi->start+1, gi->end+1);
	  DEBUG("\t\twith pattern %s-%s", gi->donor_pt, gi->acceptor_pt);
	  DEBUG("\t\twith type %s", (gi->type == 0)?("U12"):((gi->type == 1)?("U2"):("UNCLASSIFIED")));
	  DEBUG("\t\twith donor-acceptor scores %f-%f", gi->score5, gi->score3);
	  DEBUG("\t\twith BPS %f", gi->BPS_score);
	  DEBUG("\t\twith BPS position %d", gi->BPS_position);
	  DEBUG("\t\tagree type %d", gi->agree_type);
	  my_assert(gi->supportingESTs > 0);
	  DEBUG("\t\tsupporting ESTs %d", gi->supportingESTs);
  }
  listit_destroy(debug_gen_intron_it);*/

  DEBUG("All the agree flags set!");

   /*plistit debug_it1=list_first(est_with_intron_list);
  while(listit_has_next(debug_it1)){
		pEST est_prova=(pEST)listit_next(debug_it1);
		DEBUG(">%s", est_prova->info->EST_id);

		my_assert(est_prova->factorizations != NULL);
		my_assert(list_size(est_prova->factorizations) == 1);

		plist comp_prova=list_head(est_prova->factorizations);
		plistit debug_it2=list_first(comp_prova);
		while(listit_has_next(debug_it2)){
			pintron i=(pintron)listit_next(debug_it2);
			my_assert(i->donor != NULL || i->acceptor != NULL);

			DEBUG("Intron %d-%d (EST cut %d)", i->gen_intron->start, i->gen_intron->end, (i->donor == NULL)?(i->acceptor->EST_start-1):(i->donor->EST_end));
			DEBUG("\t%s!", (i->isReal)?("TRUE"):("FALSE"));
			DEBUG("\tpattern %s-%s", i->gen_intron->donor_pt, i->gen_intron->acceptor_pt);
			DEBUG("\t\t%s", (i->gen_intron->type == 0)?("U12"):((i->gen_intron->type == 1)?("U2"):("UNCLASSIFIED")));
			DEBUG("\ttry agreement: %s!", (i->try_agree)?("TRUE"):("FALSE"));
			DEBUG("\tagree type: %d!", i->agree_type);
			DEBUG("Link to EST %s", i->est_info->EST_id);
			DEBUG("Link to EST gb %s", i->est_info->EST_gb);
			DEBUG("Link to seq %s", i->est_info->EST_seq);
		}
		listit_destroy(debug_it2);
  }
  listit_destroy(debug_it1);*/

  /*plistit debug_it1=list_first(refseq_list);
  DEBUG("RefSeq introns:");
  while(listit_has_next(debug_it1)){
		ppointer pp=(ppointer)listit_next(debug_it1);
		pintron i=(pintron)pp->pointer;
		my_assert(i->isReal == true);
		my_assert(i->agree_type == 0);
		my_assert(i->try_agree == false);
		my_assert(i->donor != NULL && i->acceptor != NULL);

		DEBUG("Intron %d-%d (EST cut %d)", i->gen_intron->start, i->gen_intron->end, i->donor->EST_end);
		DEBUG("\tpattern %s-%s", i->gen_intron->donor_pt, i->gen_intron->acceptor_pt);
		DEBUG("\t\t%s", (i->gen_intron->type == 0)?("U12"):((i->gen_intron->type == 1)?("U2"):("UNCLASSIFIED")));
		DEBUG("Link to EST %s", i->est_info->EST_id);
		DEBUG("Link to EST gb %s", i->est_info->EST_gb);
		DEBUG("Link to seq %s", i->est_info->EST_seq);

  }
  listit_destroy(debug_it1);*/

  /*plistit debug_it1=list_first(canonical_list);
  DEBUG("Canonical introns (no RefSeq):");
  while(listit_has_next(debug_it1)){
		ppointer pp=(ppointer)listit_next(debug_it1);
		pintron i=(pintron)pp->pointer;
		my_assert(i->isReal == true);
		my_assert(i->agree_type == 1);
		my_assert(i->try_agree == false);
		my_assert(i->donor != NULL && i->acceptor != NULL);

		DEBUG("Intron %d-%d (EST cut %d)", i->gen_intron->start, i->gen_intron->end, i->donor->EST_end);
		DEBUG("\tpattern %s-%s", i->gen_intron->donor_pt, i->gen_intron->acceptor_pt);
		DEBUG("\t\t%s", (i->gen_intron->type == 0)?("U12"):((i->gen_intron->type == 1)?("U2"):("UNCLASSIFIED")));
		DEBUG("Link to EST %s", i->est_info->EST_id);
		DEBUG("Link to EST gb %s", i->est_info->EST_gb);
		DEBUG("Link to seq %s", i->est_info->EST_seq);

  }
  listit_destroy(debug_it1);*/

  /*plistit debug_it1=list_first(agreement_list);
  DEBUG("Introns to be agreed:");
  while(listit_has_next(debug_it1)){
		ppointer pp=(ppointer)listit_next(debug_it1);
		pintron i=(pintron)pp->pointer;
		my_assert(i->isReal == true);
		my_assert(i->agree_type == 2);
		my_assert(i->try_agree == true);
		my_assert(i->donor != NULL && i->acceptor != NULL);

		DEBUG("Intron %d-%d (EST cut %d)", i->gen_intron->start, i->gen_intron->end, i->donor->EST_end);
		DEBUG("\tpattern %s-%s", i->gen_intron->donor_pt, i->gen_intron->acceptor_pt);
		DEBUG("\t\t%s", (i->gen_intron->type == 0)?("U12"):((i->gen_intron->type == 1)?("U2"):("UNCLASSIFIED")));
		DEBUG("Link to EST %s", i->est_info->EST_id);
		DEBUG("Link to EST gb %s", i->est_info->EST_gb);
		DEBUG("Link to seq %s", i->est_info->EST_seq);

  }
  listit_destroy(debug_it1);*/

  MYTIME_stop(pt_pre);

  log_info(floginfo, "preprocessing-end");

  MYTIME_start(pt_alg);

//  unsigned int get_agreement_error(char *genomic_sequence, pintron intron_from, pintron intron_to){

  DEBUG("Try agreement canonical (%zu) -> refseq (%zu):", list_size(canonical_list), list_size(genomic_refseq_list));
  //Try to agree canonical introns to some refseq intron
  plistit can_to_ref_agree_it=list_first(canonical_list);
  while(listit_has_next(can_to_ref_agree_it)){
		ppointer pp=(ppointer)listit_next(can_to_ref_agree_it);
		pintron intron_from=(pintron)pp->pointer;
 		DEBUG("Try agree canonical intron %d-%d (EST %s)", intron_from->gen_intron->start, intron_from->gen_intron->end, intron_from->est_info->EST_gb);

 		DEBUG("...to a RefSeq intron:");
 		bool agree_ok=try_agreement_to_intron_list(gen->EST_seq, intron_from, genomic_refseq_list, 0);
 		if(agree_ok == true){
	 		DEBUG("Agree to a RefSeq intron!");
 		}
   }
  listit_destroy(can_to_ref_agree_it);

   //printf("Try agreement canonical (%zu) -> refseq (%zu)\n", list_size(canonical_list), list_size(genomic_refseq_list));

  //printf("Try agreement canonical (%zu) -> canonical (%zu)\n", list_size(canonical_list), list_size(genomic_canonical_list));

  DEBUG("Try agreement canonical -> canonical:");
   //Try to agree canonical introns to some other canonical intron
  can_to_ref_agree_it=list_first(canonical_list);
  while(listit_has_next(can_to_ref_agree_it)){
		ppointer pp=(ppointer)listit_next(can_to_ref_agree_it);
		pintron intron_from=(pintron)pp->pointer;
		if(intron_from->agreed == false){
			my_assert(intron_from->gen_intron->burset_frequency != -1);
			int freq_from=intron_from->gen_intron->burset_frequency;
			//int freq_from=get_intron_Burset_frequency(gen->EST_seq, intron_from->gen_intron);
			DEBUG("Try agree canonical intron %d-%d (EST %s, Burset frequency=%d)", intron_from->gen_intron->start, intron_from->gen_intron->end, intron_from->est_info->EST_gb, freq_from);
			DEBUG("...to a better (Burset) intron:");
		    plistit to_agree_it=list_first(genomic_canonical_list);
		    bool agree_ok=false;
		    while(agree_ok == false && listit_has_next(to_agree_it)){
				ppointer pp_to=(ppointer)listit_next(to_agree_it);
		 		pgenomic_intron gen_intron_to=(pgenomic_intron)pp_to->pointer;
		 		if(gen_intron_to->start != intron_from->gen_intron->start || gen_intron_to->end != intron_from->gen_intron->end){
					my_assert(gen_intron_to->burset_frequency != -1);
					int freq_to=gen_intron_to->burset_frequency;
		 			//int freq_to=get_intron_Burset_frequency(gen->EST_seq, gen_intron_to);
		 			if(freq_to > freq_from){
		 				DEBUG("...%d-%d (Burset frequency=%d)", gen_intron_to->start, gen_intron_to->end, freq_to);
		 				agree_ok=try_agreement(gen->EST_seq, intron_from, gen_intron_to, 0);
		 			}
		 		}
		    }
		    listit_destroy(to_agree_it);
			if(agree_ok == true){
			 	DEBUG("Agree to a better Burset intron!");
		 	}
			else{
		 		DEBUG("NO AGREE!");
			}
		}
  }
  listit_destroy(can_to_ref_agree_it);

  //Agree of others introns
  plist agreed_list=list_create();
  plist not_agreed_list=list_create();
  plistit agree_it=list_first(agreement_list);

  DEBUG("Try agreement others -> RefSeq+canonical:");
  while(listit_has_next(agree_it)){
 		ppointer pp=(ppointer)listit_next(agree_it);
 		pintron intron_from=(pintron)pp->pointer;

 		//First to a RefSeq intron
  		DEBUG("Try agree intron %d-%d (EST %s)", intron_from->gen_intron->start, intron_from->gen_intron->end, intron_from->est_info->EST_gb);

 		DEBUG("...to a RefSeq intron:");
 		bool agree_ok=try_agreement_to_intron_list(gen->EST_seq, intron_from, genomic_refseq_list, 4);
 		if(agree_ok == false){
 	 		DEBUG("...to a canonical intron:");
 	 		agree_ok=try_agreement_to_intron_list(gen->EST_seq, intron_from, genomic_canonical_list, 4);
 			if(agree_ok == true){
 		 		DEBUG("Agree to a canonical intron!");
				ppointer pp=pointer_create();
 				pp->pointer=(pintron)intron_from;
				list_add_to_tail(agreed_list, pp);
	 		}
 			else{
 				DEBUG("...to a RefSeq intron on a single site:");
		 		bool agree_ok=try_agreement_to_intron_list_on_single_site(gen->EST_seq, intron_from, genomic_refseq_list, gen_intron_list);
  				if(agree_ok == false){
 			 		DEBUG("...to a canonical intron on a single site:");
 		 	 		agree_ok=try_agreement_to_intron_list_on_single_site(gen->EST_seq, intron_from, genomic_canonical_list, gen_intron_list);
 		 			if(agree_ok == true){
 		 		 		DEBUG("Agree to a canonical intron on a single site!");
 						ppointer pp=pointer_create();
 		 				pp->pointer=(pintron)intron_from;
 						list_add_to_tail(agreed_list, pp);
 		 	 		}
 		 			else{
 				 		DEBUG("NO AGREE!");
 						ppointer pp=pointer_create();
 		 				pp->pointer=(pintron)intron_from;
 						list_add_to_tail(not_agreed_list, pp);
 		 			}
 				}
 				else{
 		 	 		DEBUG("Agree to a RefSeq intron on a single site!");
 					ppointer pp=pointer_create();
 					pp->pointer=(pintron)intron_from;
 					list_add_to_tail(agreed_list, pp);
 				}
 			}
 		}
 		else{
 	 		DEBUG("Agree to a RefSeq intron!");
			ppointer pp=pointer_create();
			pp->pointer=(pintron)intron_from;
			list_add_to_tail(agreed_list, pp);
		}
   }
   listit_destroy(agree_it);

   list_destroy(agreement_list, (delete_function)pointer_destroy);

   /*plistit debug_it1=list_first(agreed_list);
   DEBUG("Introns agreed before Burset:");
   while(listit_has_next(debug_it1)){
 		ppointer pp=(ppointer)listit_next(debug_it1);
 		pintron i=(pintron)pp->pointer;
 		my_assert(i->isReal == true);
 		my_assert(i->agreed == true);
 		my_assert(i->donor != NULL && i->acceptor != NULL);

 		DEBUG("EST %s ==> intron %d-%d (EST cut %d)", i->est_info->EST_gb, i->gen_intron->start, i->gen_intron->end, i->donor->EST_end);
 		DEBUG("\tpattern %s-%s", i->gen_intron->donor_pt, i->gen_intron->acceptor_pt);
    }
   listit_destroy(debug_it1);*/

   /*plistit debug_it1=list_first(not_agreed_list);
    DEBUG("Introns not agreed before Burset:");
    while(listit_has_next(debug_it1)){
  		ppointer pp=(ppointer)listit_next(debug_it1);
  		pintron i=(pintron)pp->pointer;
  		my_assert(i->isReal == true);
  		my_assert(i->agreed == false);
  		my_assert(i->donor != NULL && i->acceptor != NULL);

  		DEBUG("EST %s ==> intron %d-%d (EST cut %d)", i->est_info->EST_gb, i->gen_intron->start, i->gen_intron->end, i->donor->EST_end);
  		DEBUG("\tpattern %s-%s", i->gen_intron->donor_pt, i->gen_intron->acceptor_pt);
     }
    listit_destroy(debug_it1);*/

   DEBUG("Try agreement others -> others:");
  //Try agree to an intron with a better Burset pattern
   plist final_not_agreed_list=list_create();
   agree_it=list_first(not_agreed_list);
   while(listit_has_next(agree_it)){
 		ppointer pp=(ppointer)listit_next(agree_it);
 		pintron intron_from=(pintron)pp->pointer;
 		my_assert(intron_from->gen_intron->burset_frequency != -1);
 		int freq_from=intron_from->gen_intron->burset_frequency;
 		//int freq_from=get_intron_Burset_frequency(gen->EST_seq, intron_from->gen_intron);
		DEBUG("Try agree intron %d-%d (EST %s, Burset frequency=%d)", intron_from->gen_intron->start, intron_from->gen_intron->end, intron_from->est_info->EST_gb, freq_from);
 		DEBUG("...to a better (Burset) intron:");
	    plistit to_agree_it=list_first(genomic_agreement_list);
	    bool agree_ok=false;
	    while(agree_ok == false && listit_has_next(to_agree_it)){
			ppointer pp_to=(ppointer)listit_next(to_agree_it);
	 		pgenomic_intron gen_intron_to=(pgenomic_intron)pp_to->pointer;
	 		if(gen_intron_to->start != intron_from->gen_intron->start || gen_intron_to->end != intron_from->gen_intron->end){
	 			my_assert(gen_intron_to->burset_frequency != -1);
	 			int freq_to=gen_intron_to->burset_frequency;
	 			//int freq_to=get_intron_Burset_frequency(gen->EST_seq, gen_intron_to);
	 			if(freq_to > freq_from){
	 				DEBUG("...%d-%d (Burset frequency=%d)", gen_intron_to->start, gen_intron_to->end, freq_to);
	 				if(gen_intron_to->supportingESTs > 0)
	 					agree_ok=try_agreement(gen->EST_seq, intron_from, gen_intron_to, 4);
	 			}
	 		}
	    }
	    listit_destroy(to_agree_it);
		if(agree_ok == true){
		 	DEBUG("Agree to a better Burset intron!");
			ppointer pp=pointer_create();
			pp->pointer=(pintron)intron_from;
			list_add_to_tail(agreed_list, pp);
	 	}
		else{
	 		DEBUG("NO AGREE!");
			ppointer pp=pointer_create();
			pp->pointer=(pintron)intron_from;
			list_add_to_tail(final_not_agreed_list, pp);
		}
   }
   listit_destroy(agree_it);

   list_destroy(not_agreed_list, (delete_function)pointer_destroy);

   plistit debug_it1=list_first(agreed_list);
   DEBUG("Introns agreed:");
   while(listit_has_next(debug_it1)){
 		ppointer pp=(ppointer)listit_next(debug_it1);
 		pintron i=(pintron)pp->pointer;
 		my_assert(i->isReal == true);
 		my_assert(i->agreed == true);
 		my_assert(i->donor != NULL && i->acceptor != NULL);

 		DEBUG("EST %s ==> intron %d-%d (EST cut %d)", i->est_info->EST_gb, i->gen_intron->start, i->gen_intron->end, i->donor->EST_end);
 		DEBUG("\tpattern %s-%s", i->gen_intron->donor_pt, i->gen_intron->acceptor_pt);
    }
   listit_destroy(debug_it1);

   /*plistit debug_it1=list_first(final_not_agreed_list);
     DEBUG("Introns not agreed:");
     while(listit_has_next(debug_it1)){
   		ppointer pp=(ppointer)listit_next(debug_it1);
   		pintron i=(pintron)pp->pointer;
   		my_assert(i->isReal == true);
   		my_assert(i->agreed == false);
   		my_assert(i->donor != NULL && i->acceptor != NULL);

   		DEBUG("EST %s ==> intron %d-%d (EST cut %d)", i->est_info->EST_gb, i->gen_intron->start, i->gen_intron->end, i->donor->EST_end);
   		DEBUG("\tpattern %s-%s", i->gen_intron->donor_pt, i->gen_intron->acceptor_pt);
      }
     listit_destroy(debug_it1);*/

   //ATTENZIONE: in questa fase la lista agreed_list non viene aggiornata aggiungendo gli introni modificati
   //e la lista final_not_agreed_list non viene aggiornata togliendo gli introni modificati
	DEBUG("Search a better intron:");
    plistit not_agree_it=list_first(final_not_agreed_list);
    while(listit_has_next(not_agree_it)){
  		ppointer pp=(ppointer)listit_next(not_agree_it);
  		pintron intron_from=(pintron)pp->pointer;
		DEBUG("Try agree intron %d-%d (EST %s)", intron_from->gen_intron->start, intron_from->gen_intron->end, intron_from->est_info->EST_gb);
		bool agree_ok=find_better_intron(gen->EST_seq, intron_from, gen_intron_list);
		if(agree_ok == true){
 		 	DEBUG("A better intron was found!");
 		}
 		else{
 	 		DEBUG("A better intron was not found!");
 		}
    }
    listit_destroy(not_agree_it);

    MYTIME_stop(pt_alg);

  log_info(floginfo, "intron-agreement-end");

  MYTIME_start(pt_io);

  est_list_it=list_first(est_with_intron_list);
  while(listit_has_next(est_list_it)){
 	 pEST est=(pEST)listit_next(est_list_it);
 	 DEBUG("EST ==> %s", est->info->EST_id);

 	 my_assert(est->factorizations != NULL);
 	 my_assert(list_size(est->factorizations) == 1);

	 plist intron_composition=(plist)list_head(est->factorizations);

	 my_assert(!list_is_empty(intron_composition));

 	 plist exon_composition_list=list_create();
 	 plist exon_composition=list_create();

 	 pintron head=list_remove_from_head(intron_composition);
 	 my_assert(head->donor == NULL);

 	 plistit intron_it=list_first(intron_composition);
 	 while(listit_has_next(intron_it)){
 		pintron intron=(pintron)listit_next(intron_it);
 		my_assert(intron->donor != NULL);
 		list_add_to_tail(exon_composition, intron->donor);
 		if(intron->isReal == true){
 			my_assert(intron->gen_intron != NULL);
 			if(intron->gen_intron->info == NULL)
 				intron->gen_intron->info=list_create();
 			pgenomic_intron_info gii=genomic_intron_info_create(est->info, intron->donor->EST_end);
 			list_add_to_tail(intron->gen_intron->info, gii);
 		}
 	 }
 	 listit_destroy(intron_it);

	 list_destroy(intron_composition,(delete_function)intron_destroy_free_if_false);

	 DEBUG("Exon composition retrieved!");

	 /* plistit debug_it3=list_first(exon_composition);
	   DEBUG("Exon composition:");
	   while(listit_has_next(debug_it3)){
	 		pfactor factor=(pfactor)listit_next(debug_it3);
	 		DEBUG("%d-%d (%d-%d)", factor->GEN_start, factor->GEN_end, factor->EST_start, factor->EST_end);
	    }
	   listit_destroy(debug_it3);*/

	 list_add_to_head(exon_composition_list, exon_composition);

	 est->factorizations=exon_composition_list;

	 //Write also the external exons
	 write_multifasta_output(gen, est, f_multif_out, 1);
 }
  listit_destroy(est_list_it);

  /*plistit debug_gen_intron_it=list_first(gen_intron_list);
  while(listit_has_next(debug_gen_intron_it)){
	  pgenomic_intron gi=(pgenomic_intron)listit_next(debug_gen_intron_it);
	  DEBUG("\tIntron %d-%d", gi->start+1, gi->end+1);
	  DEBUG("\t\twith pattern %s-%s", gi->donor_pt, gi->acceptor_pt);
	  DEBUG("\t\twith type %s", (gi->type == 0)?("U12"):((gi->type == 1)?("U2"):("UNCLASSIFIED")));
	  DEBUG("\t\twith donor-acceptor scores %f-%f", gi->score5, gi->score3);
	  DEBUG("\t\twith BPS %f", gi->BPS_score);
	  DEBUG("\t\twith BPS position %d", gi->BPS_position);
	  my_assert((int)list_size(gi->info) == gi->supportingESTs);
	  DEBUG("\t\tESTs ==> %zu", list_size(gi->info));
  }
  listit_destroy(debug_gen_intron_it);*/

  int strand=atoi(gen->EST_strand_as_read);
  list_sort(gen_intron_list, (comparator)genomic_intron_compare);
  plistit out_gen_intron_it=list_first(gen_intron_list);
  int gen_intron_index=1;
  bool first_time=true;
   while(listit_has_next(out_gen_intron_it)){
	  pgenomic_intron gi=(pgenomic_intron)listit_next(out_gen_intron_it);
	  if(!list_is_empty(gi->info)){
		  my_assert(gi->classified == true);

		  if(first_time == false){
			  fprintf(gtf_out, "\n");
		  }
		  else{
			  first_time=false;
		  }

	 	  //fprintf(gtf_out, "intron#%d\t", gen_intron_index);
	 	  fprintf(gtf_out, "%d\t%d\t", gi->start+1, gi->end+1);

	 	  int abs_start, abs_end;
	 	  if(gen->abs_start < gen->abs_end)
	 		  get_ABS_region_start_end(gen->abs_start, gen->abs_end, strand, gi->start+1, gi->end+1, &abs_start, &abs_end);
	 	  else
	 		  get_ABS_region_start_end(gen->abs_end, gen->abs_start, strand, gi->start+1, gi->end+1, &abs_start, &abs_end);

	 	  fprintf(gtf_out, "%d\t%d\t", abs_start, abs_end);	//Absolute coordinates
	 	  fprintf(gtf_out, "%d\t", gi->end-gi->start+1);	//Intron length
		  fprintf(gtf_out, "%zu\t", list_size(gi->info));

		  char *repeat=GetRepeatSequence(gen->EST_seq, gi->start, gi->end);

		  char *donor_suffix=get_donor_suffix(gen->EST_seq, gi, 15);
		  char *acceptor_prefix=get_acceptor_prefix(gen->EST_seq, gi, 15);
		  char *intron_prefix=get_intron_prefix(gen->EST_seq, gi, 20);
		  char *intron_suffix=get_intron_suffix(gen->EST_seq, gi, 20);

		  unsigned int tot_donor_edit=0, tot_acceptor_edit=0;

		  plistit info_it=list_first(gi->info);
		   while(listit_has_next(info_it)){
			   pgenomic_intron_info i_info=listit_next(info_it);
			   fprintf(gtf_out, "%s,", i_info->info->EST_gb);

			   char *donor_EST_suffix=get_donor_EST_suffix(i_info->info->EST_seq, i_info->EST_cut+1, 15);
			   char *acceptor_EST_prefix=get_acceptor_EST_prefix(i_info->info->EST_seq, i_info->EST_cut+1, 15);

			   unsigned int error;
			   size_t l1, l2;
			   unsigned int* M;

			   l1=strlen(donor_EST_suffix);
			   l2=strlen(donor_suffix);
			   M=edit_distance(donor_EST_suffix, l1, donor_suffix, l2);
			   error=M[(l1+1)*(l2+1)-1];
			   pfree(M);
			   tot_donor_edit+=error;

			   l1=strlen(acceptor_EST_prefix);
			   l2=strlen(acceptor_prefix);
			   M=edit_distance(acceptor_EST_prefix, l1, acceptor_prefix, l2);
			   error=M[(l1+1)*(l2+1)-1];
			   pfree(M);
			   tot_acceptor_edit+=error;

			   pfree(donor_EST_suffix);
			   pfree(acceptor_EST_prefix);
		   }

		   double mean_donor_edit=(double)tot_donor_edit/(double)list_size(gi->info);
		   double mean_acceptor_edit=(double)tot_acceptor_edit/(double)list_size(gi->info);

		   fprintf(gtf_out, "\t%f\t%f\t", mean_donor_edit, mean_acceptor_edit);	//Mean edit in donor/acceptor
		   fprintf(gtf_out, "%f\t%f\t", gi->score5, gi->score3);
		   fprintf(gtf_out, "%f\t%d\t", gi->BPS_score, gi->BPS_position);
		   fprintf(gtf_out, "%d\t", gi->type);
		   fprintf(gtf_out, "%s%s\t", gi->donor_pt, gi->acceptor_pt);

			if(repeat == NULL)
			  fprintf(gtf_out, ".\t");
			else
			  fprintf(gtf_out, "%s\t", repeat);	//Repeat sequence

		   fprintf(gtf_out, "%s\t", donor_suffix);	//Donor suffix
		   fprintf(gtf_out, "%s\t", intron_prefix);	//Intron prefix
		   fprintf(gtf_out, "%s\t", intron_suffix);	//Intron suffix
		   fprintf(gtf_out, "%s", acceptor_prefix);	//Acceptor prefix

		   pfree(donor_suffix);
		   pfree(intron_prefix);
		   pfree(intron_suffix);
		   pfree(acceptor_prefix);
			if (repeat != NULL)
			  pfree(repeat);

		   listit_destroy(info_it);

		   gen_intron_index++;
	  }
   }
   listit_destroy(out_gen_intron_it);

  MYTIME_stop(pt_io);

  log_info(floginfo, "output-end");

  DEBUG("Finalizing structures");
//config_destroy(config);

  list_destroy(gen_intron_list,(delete_function)genomic_intron_destroy);

  EST_info_destroy(gen);
  est_list_it=list_first(est_with_intron_list);
  while(listit_has_next(est_list_it)){
	 pEST est=(pEST)listit_next(est_list_it);
	 factorization_list_destroy(est->factorizations);
	 est->factorizations=NULL;
  }
  listit_destroy(est_list_it);

  list_destroy(est_with_intron_list, (delete_function)EST_destroy);
  list_destroy(refseq_list, (delete_function)pointer_destroy);
  list_destroy(canonical_list, (delete_function)pointer_destroy);
  list_destroy(agreed_list, (delete_function)pointer_destroy);
  list_destroy(final_not_agreed_list, (delete_function)pointer_destroy);

  list_destroy(genomic_refseq_list, (delete_function)pointer_destroy);
  list_destroy(genomic_canonical_list, (delete_function)pointer_destroy);
  list_destroy(genomic_agreement_list, (delete_function)pointer_destroy);


  fclose(f_multif_out);
  fclose(gtf_out);

  MYTIME_stop(pt_tot);

  MYTIME_LOG(INFO, pt_pre);
  MYTIME_LOG(INFO, pt_alg);
  MYTIME_LOG(INFO, pt_io);
  MYTIME_LOG(INFO, pt_tot);

  MYTIME_destroy(pt_tot);
  MYTIME_destroy(pt_pre);
  MYTIME_destroy(pt_alg);
  MYTIME_destroy(pt_io);

  log_info(floginfo, "end");

  INFO("End");
  resource_usage_log();
  fclose(floginfo);
  return 0;
}
