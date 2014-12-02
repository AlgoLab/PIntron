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

#include "factorization-refinement.h"
#include "compute-est-fact.h"


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
  rev_est->fixed_strand= est->fixed_strand;
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
  pmytime pt_io= MYTIME_create_with_name("IO");

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

	 list_add_to_tail(new_est_list, est);

	 DEBUG("Replace polyA/T with fake characters");
	 polyAT_substitution(est);

        if (!est->fixed_strand) {
          DEBUG("Strand is not fixed. Adding also its reverse and complement.");
          pEST_info rev_est= copy_and_reverse(est);
          list_add_to_tail(new_est_list, rev_est);
          polyAT_substitution(rev_est);
        }
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
  bool reversed= false;

  while (listit_has_next(estit)) {
    pEST_info est= (pEST_info)listit_next(estit);
    pEST factorized_est=
      compute_est_fact(gen, est, tree, pg,
                       floginfo, fmeg, fpmeg, ftmeg,
                       fintronic,
                       pt_alg, pt_comp, pt_io, config);

    if (!list_is_empty(factorized_est->factorizations)) {
      // Found valid factorizations
      INFO("Found valid factorization(s) for EST %s on the %s strand of the genomic.",
           est->EST_gb,
           ((est->EST_strand==1)?"same":"opposite"));
      MYTIME_start(pt_io);
      write_multifasta_output(gen, factorized_est, f_multif_out, config->retain_externals);
      write_single_EST_info(est_multif_out, factorized_est->info);
      MYTIME_stop(pt_io);
      // Skip next sequence if it is the reverse of this one
      if (!est->fixed_strand && !reversed) {
        listit_next(estit);
        ++id_p;
      }
      reversed= false; // Next sequence is "original"
    } else {
      // No valid factorization found
      if (reversed || est->fixed_strand) {
        // It is already the rev&compl sequence or it cannot be rev&complement
        INFO("...the EST %s has no alignment! (Fixed strand? %s)",
             est->EST_gb, (est->fixed_strand ? "true" : "false"));
        reversed= false; // Next sequence is "original"
      } else {
        // !reversed && !fixed_strand
        INFO("...the strand from the input file may be wrong (read: '%s')!", est->EST_strand_as_read);
        reversed= true; // Next sequence is reversed
      }
    }

    EST_destroy_just_factorizations(factorized_est);

    ++id_p;

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

  fclose(fmeg);
  fclose(fpmeg);
  fclose(ftmeg);
  fclose(f_multif_out);
  fclose(est_multif_out);
  fclose(fintronic);

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
