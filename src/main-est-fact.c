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

#include "max-emb-graph.h"
#include "aug_suffix_tree.h"

#include "io-multifasta.h"
#include "io-meg.h"

#include "meg-simplification.h"
#include "est-factorizations.h"

#include "my_time.h"
#include "log.h"
#include "log-build-info.h"

#ifndef MUMMER_EMULATION
static void
log_meg(pext_array V) {
  DEBUG("Start MEG");
  DEBUG("(   p,    t,    l)");
  for (unsigned int i= 0; i< EA_size(V); ++i) {
	 plist Vi= EA_get(V, i);
	 plistit lit= list_first(Vi);
	 while (listit_has_next(lit)) {
		ppairing p= listit_next(lit);
		DEBUG("(%4d, %4d, %4d)", PAIRING(p));
		plistit ait= list_first(p->adjs);
		while (listit_has_next(ait)) {
		  ppairing a= listit_next(ait);
		  DEBUG("     -> (%4d, %4d, %4d)", PAIRING(a));
		}
		listit_destroy(ait);
	 }
	 listit_destroy(lit);
  }
  DEBUG("End MEG");
}
#endif

#ifdef MUMMER_EMULATION

void
save_vertex_set_to_file(FILE* fvset,
								char* est_id,
								pext_array V,
								bool reversed) {
  if (!reversed)
	 fprintf(fvset, "> %s\n", est_id);
  else
	 fprintf(fvset, "> %s Reverse\n", est_id);
  for (size_t i= 1; i<EA_size(V)-1; ++i) {
	 plist Vi= EA_get(V, i);
	 plistit Vit= list_first(Vi);
	 while (listit_has_next(Vit)) {
		ppairing I= listit_next(Vit);
		fprintf(fvset,"%d %d %d\n", I->t+1, I->p+1, I->l);
	 }
	 listit_destroy(Vit);
  }
  fflush(fvset);
}

#endif

static char*
create_and_copy(const char const* src) {
  if (src==NULL)
	 return NULL;
  char* dest= alloc_and_copy(src);
  return dest;
}

pEST_info copy_and_reverse(pEST_info est) {
  pEST_info rev_est= EST_info_create();
  rev_est->EST_seq= create_and_copy(est->EST_seq);
  rev_est->original_EST_seq= create_and_copy(est->original_EST_seq);
  reverse_and_complement(rev_est);
  rev_est->EST_id= create_and_copy(est->EST_id);
  rev_est->EST_gb= create_and_copy(est->EST_gb);
  rev_est->EST_chr= create_and_copy(est->EST_chr);
  rev_est->EST_strand_as_read=
	 create_and_copy(est->EST_strand_as_read);
  rev_est->EST_strand= - est->EST_strand;
  rev_est->pref_polyA_length= est->suff_polyT_length;
  rev_est->suff_polyA_length= est->pref_polyT_length;
  rev_est->pref_polyT_length= est->suff_polyA_length;
  rev_est->suff_polyT_length= est->pref_polyA_length;
  return rev_est;
}


int main(int argc, char** argv) {
  INFO("EST-FACTORIZATION v2");
  PRINT_LICENSE_INFORMATION;
  PRINT_SYSTEM_INFORMATION;
  DEBUG("Initialization");
  pmytime pt_tot= MYTIME_create_with_name("Total");
  pmytime pt_st= MYTIME_create_with_name("Suffix Tree");
  pmytime pt_alg= MYTIME_create_with_name("Algorithm");
  pmytime pt_comp= MYTIME_create_with_name("Compositions");
  pmytime pt_ccomp= MYTIME_create_with_name("Internal Comp.");
  pmytime pt_io= MYTIME_create_with_name("IO");
  pmytime pt_meg= MYTIME_create_with_name("MEGs");

  MYTIME_start(pt_tot);
  pconfiguration config= config_create(argc, argv);
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
  list_destroy(gen_list, noop_free);

  parse_genomic_header(gen);

  DEBUG("Removing N tails");
  Ntails_removal(gen);

  DEBUG("Reading EST sequences");
  FILE* fests= fopen("ests.txt", "r");
  if (!fests) {
	 FATAL("File ests.txt not found! Terminating");
	 fail();
  }
  plist est_list= read_multifasta(fests);
  fclose(fests);
  DEBUG("EST sequences read");

#ifdef MUMMER_EMULATION

  FILE* fvset= fopen("vertex_set.txt", "w");
  if (!fvset) {
	 FATAL("Cannot create file vertex_set.txt! Terminating");
	 fail();
  }

#else

  FILE* f_multif_out= fopen("raw-multifasta-out.txt", "w");
  if (!f_multif_out) {
	 FATAL("Cannot create file raw-multifasta-out.txt! Terminating");
	 fail();
  }

  FILE* fmeg= fopen("megs.txt", "w");
  if (!fmeg) {
	 FATAL("Cannot create file megs.txt! Terminating");
	 fail();
  }

  FILE* fpmeg= fopen("processed-megs.txt", "w");
  if (!fpmeg) {
	 FATAL("Cannot create file processed-megs.txt! Terminating");
	 fail();
  }

  FILE* ftmeg= fopen("processed-megs-info.txt", "w");
  if (!ftmeg) {
	 FATAL("Cannot create file processed-megs-time.txt! Terminating");
	 fail();
  }

  FILE* est_multif_out= fopen("processed-ests.txt", "w");
  if (!est_multif_out) {
	 FATAL("Cannot create file processed-ests.txt! Terminating");
	 fail();
  }

  FILE* fintronic= fopen("meg-edges.txt", "w");
  if (!fintronic) {
	 FATAL("Cannot create file meg-edges.txt! Terminating");
	 fail();
  }
#endif

// Log resource utilization
  log_info(floginfo, "data-io-end");

  MYTIME_stop(pt_io);

  const size_t n_est= list_size(est_list);

  INFO("Read %zd sequences.", n_est);
  plist new_est_list= list_create();

  plistit estit= list_first(est_list);
  while (listit_has_next(estit)) {
	 pEST_info est= (pEST_info)listit_next(estit);
	 INFO("EST: %s", est->EST_id);

	 DEBUG("Set the GB id from the fasta header");
	 set_EST_GB_identification(est);

	 DEBUG("Set the EST strand and RC");
	 set_EST_Strand_and_RC(est, gen);

	 pEST_info rev_est= copy_and_reverse(est);
	 list_add_to_tail(new_est_list, est);
	 list_add_to_tail(new_est_list, rev_est);

#ifndef MUMMER_EMULATION
	 DEBUG("Replace polyA/T with fake characters");
	 polyAT_substitution(est);
	 polyAT_substitution(rev_est);
#endif

  }
  listit_destroy(estit);

  list_destroy(est_list, (delete_function)noop_free);
  est_list= new_est_list;

  INFO("Creating the suffix tree");

// Log resource utilization
  log_info(floginfo, "gst-construction-begin");

  MYTIME_start(pt_st);
  LST_StringSet *set= lst_stringset_new();
  LST_String * lst= PALLOC(LST_String);
  lst_string_init(lst, gen->EST_seq, sizeof(char),
						strlen(gen->EST_seq));
  lst_stringset_add(set, lst);
  LST_STree* tree = lst_stree_new(set);
  MYTIME_stop(pt_st);

// Log resource utilization
  log_info(floginfo, "gst-preprocessing-begin");

  DEBUG("Preprocessing the GST");
  MYTIME_start(pt_alg);
  ppreproc_gen pg= PGen_create();
  preprocess_text(gen, pg);
  stree_preprocess(tree, pg, config);
  MYTIME_stop(pt_alg);

// Log resource utilization
  log_info(floginfo, "gst-preprocessing-end");

  size_t id_p= 1;
  estit= list_first(est_list);

  //Print work in progress...
  //int counter=0;
  //int input_size=list_size(est_list);

  double processing_percentage;
  while (listit_has_next(estit)) {
	 pEST_info est= (pEST_info)listit_next(estit);
	 INFO("EST: %s", est->EST_id);

	  //Print work in progress...
	 /*counter++;
	 processing_percentage=((double)counter)/((double)input_size);
	 processing_percentage=processing_percentage*100.0;
	 if(1){
		fprintf(stdout, "%f of ESTs processed!\n", processing_percentage);
	 }*/

	 MYTIME_start(pt_alg);
	 DEBUG("Calculating the occurrence set");

	 pext_array V= NULL;
	 bool too_complex= false;
	 size_t inc_pairing_len= 0;
	 log_info(floginfo, "meg-construction-begin");
	 do {
		DEBUG("Building the MEG vertex set");

		config->min_factor_len += inc_pairing_len;
		V= build_vertex_set(est, tree, pg, config);
#ifdef MUMMER_EMULATION
		save_vertex_set_to_file(fvset, est->EST_id, V, id_p % 2 == 0);
#else
		MYTIME_reset(pt_meg);
		MYTIME_start(pt_meg);
		DEBUG("Building the MEG edge set");
		build_edge_set(V, config);
		save_meg_to_filename(V, "meg-1-untouched.dot");
		simplify_meg(V, config);
		save_meg_to_filename(V, "meg-2-after-basic-simplification.dot");
		if (config->trans_red) {
		  pgraph g= meg2graph(V);
		  transitive_reduction(g);
		  graph_destroy(g);
		  save_meg_to_filename(V, "meg-3-after-transitive-reduction.dot");
		}
		too_complex= is_too_complex_for_compaction(V, config);
		if (!too_complex && config->short_edge_comp) {
		  compact_short_edges(V, config);
		  save_meg_to_filename(V, "meg-4-after-short-edge-contraction.dot");
		}
		DEBUG("Analyzing complexity of the MEG vertex set");
		too_complex= too_complex || is_too_complex(V, config);
		config->min_factor_len -= inc_pairing_len;
		if (too_complex) {
		  ++inc_pairing_len;
		  EA_destroy(V, (delete_function)vi_destroy);
		  INFO("MEG too much complex. Re-trying with min-factor-len= %zd.",
				 config->min_factor_len+inc_pairing_len);
		}
		MYTIME_stop(pt_meg);
#endif
	 } while(too_complex);
	 log_info(floginfo, "meg-construction-end");

	 MYTIME_stop(pt_alg);
#ifndef MUMMER_EMULATION
	 MYTIME_start(pt_io);
	 log_meg(V);
#ifndef NDEBUG
	 print_meg(V, stdout);
	 fflush(stdout);
#endif
	 DEBUG("Appending MEG to the file");
	 fprintf(fmeg, "\n\n***********\n\n");
	 write_single_EST_info(fmeg, est);
	 meg_write(fmeg, V);
	 fflush(fmeg);
	 MYTIME_stop(pt_io);

	 //PROVA EST-FACTORIZATIONS

	 log_info(floginfo, "est-factorization-begin");

	 MYTIME_start(pt_comp);
	 MYTIME_reset(pt_ccomp);
	 MYTIME_start(pt_ccomp);
	 pEST factorized_est=get_EST_factorizations(est, V, config, gen);
	 MYTIME_stop(pt_ccomp);
	 MYTIME_stop(pt_comp);

	 if(!list_is_empty(factorized_est->factorizations)) {
		MYTIME_start(pt_io);
		INFO("...EST aligned!");
		write_multifasta_output(gen, factorized_est, f_multif_out, config->retain_externals);
		write_single_EST_info(est_multif_out, factorized_est->info);
		fprintf(fintronic, ">%s\n", est->EST_id);
		add_intronic_edges_to_file(fintronic, V);
		write_single_EST_info(fpmeg, est);
		meg_write(fpmeg, V);
		fprintf(ftmeg, "%llu %llu %u\n",
				  MYTIME_getinterval(pt_meg),
				  MYTIME_getinterval(pt_ccomp),
				  list_size(factorized_est->factorizations));
		MYTIME_stop(pt_io);
		if(id_p % 2 == 1){
		  listit_next(estit);

		  //Print work in progress...
		  /*counter++;
		  processing_percentage=((double)counter)/((double)input_size);
		  processing_percentage=processing_percentage*100.0;
		  if(1){
			 fprintf(stdout, "%f of ESTs processed!\n", processing_percentage);
		  }*/

		  ++id_p;
		}
	 } else {
		 if(id_p % 2 == 1){
			 INFO("...the strand from the input file may be wrong!");
		 } else {
			 INFO("...the EST %s has no alignment!", est->EST_gb);
		 }
	 }

	 EST_destroy_just_factorizations(factorized_est);

#endif

	 ++id_p;

	 log_info(floginfo, "est-factorization-end");

	 DEBUG("Destroying the MEG and the occurrence set");
	 MYTIME_start(pt_alg);
	 EA_destroy(V, (delete_function)vi_destroy);
	 MYTIME_stop(pt_alg);
	 log_info(floginfo, "est-processing-end");
  }

  DEBUG("Destroying the GST additional informations");
  MYTIME_start(pt_alg);
  stree_info_destroy(tree);
  MYTIME_stop(pt_alg);
  DEBUG("Destroying the GST");
  MYTIME_start(pt_st);
  lst_stree_free(tree);
  pfree(set);
  MYTIME_stop(pt_st);

  listit_destroy(estit);

  DEBUG("Finalizing structures");
  config_destroy(config);
  pg->gen= NULL;
  PGen_destroy(pg);
  EST_info_destroy(gen);
  list_destroy(est_list, (delete_function)EST_info_destroy);

#ifdef MUMMER_EMULATION
  fclose(fvset);
#else
  fclose(fmeg);
  fclose(fpmeg);
  fclose(ftmeg);
  fclose(f_multif_out);
  fclose(est_multif_out);
  fclose(fintronic);
#endif

  MYTIME_stop(pt_tot);

  MYTIME_LOG(INFO, pt_st);
  MYTIME_LOG(INFO, pt_alg);
  MYTIME_LOG(INFO, pt_comp);
  MYTIME_LOG(INFO, pt_io);
  MYTIME_LOG(INFO, pt_tot);

  MYTIME_destroy(pt_tot);
  MYTIME_destroy(pt_st);
  MYTIME_destroy(pt_alg);
  MYTIME_destroy(pt_comp);
  MYTIME_destroy(pt_io);

  log_info(floginfo, "end");

  INFO("End");
  resource_usage_log();
  fclose(floginfo);
  return 0;
}
