/**
 *
 *
 *                              PIntron
 *
 * A novel pipeline for computational gene-structure prediction based on
 * spliced alignment of expressed sequences (ESTs and mRNAs).
 *
 * Copyright (C) 2010,2012  Yuri Pirola, Raffaella Rizzi
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

#include <stdio.h>

#include "util.h"
#include "list.h"

#include "max-emb-graph.h"

#include "io-multifasta.h"
#include "io-meg.h"

#include "meg-simplification.h"
#include "est-factorizations.h"

#include "log.h"
#include "log-build-info.h"

#include "factorization-refinement.h"

#include "compute-est-fact.h"

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

static void
report_meg(pEST_info est, pmytime pt_io,
			  FILE* fmeg, pext_array V) {
  MYTIME_START_PARALLEL(pt_io);
  log_meg(V);
#ifndef NDEBUG
  print_meg(V, stdout);
  fflush(stdout);
#endif
  DEBUG("Appending MEG to the MEGs file..");
  fprintf(fmeg, "\n\n***********\n\n");
  write_single_EST_info(fmeg, est);
  meg_write(fmeg, V);
  fflush(fmeg);
  MYTIME_STOP_PARALLEL(pt_io);
}

static void
build_meg(pEST_info est,
			 LST_STree* tree,
			 ppreproc_gen pg,
			 FILE* floginfoext,
			 pmytime pt_alg, pmytime pt_meg,
			 pconfiguration shared_config,
			 size_t* pt_inc_pairing_len,
			 pext_array* pV) {

// Create a local copy of configuration parameters
  pconfiguration config= config_clone(shared_config);

  MYTIME_START_PARALLEL(pt_alg);
  DEBUG("Calculating the occurrence set");

  bool too_complex= false;
  log_info_extended(floginfoext, "meg-construction-begin", (void*)est);
  do {
	 DEBUG("Building the MEG vertex set");

	 config->min_factor_len += *pt_inc_pairing_len;
	 *pV= build_vertex_set(est, tree, pg, config);
	 MYTIME_reset(pt_meg);
	 MYTIME_start(pt_meg);
	 DEBUG("Building the MEG edge set");
	 build_edge_set(*pV, config);
	 save_meg_to_filename(*pV, "meg-1-untouched.dot");
	 simplify_meg(*pV, config);
	 save_meg_to_filename(*pV, "meg-2-after-basic-simplification.dot");
	 if (config->trans_red) {
		pgraph g= meg2graph(*pV);
		transitive_reduction(g);
		graph_destroy(g);
		save_meg_to_filename(*pV, "meg-3-after-transitive-reduction.dot");
	 }
	 too_complex= is_too_complex_for_compaction(*pV, config);
	 if (!too_complex && config->short_edge_comp) {
		compact_short_edges(*pV, config);
		save_meg_to_filename(*pV, "meg-4-after-short-edge-contraction.dot");
	 }
	 DEBUG("Analyzing complexity of the MEG vertex set");
	 too_complex= too_complex || is_too_complex(*pV, config);
	 config->min_factor_len -= *pt_inc_pairing_len;
	 if (too_complex) {
		++(*pt_inc_pairing_len);
		EA_destroy(*pV, (delete_function)vi_destroy);
		INFO("MEG too much complex. Re-trying with min-factor-len= %zd.",
			  config->min_factor_len+(*pt_inc_pairing_len));
	 }
	 MYTIME_stop(pt_meg);
  } while(too_complex);
  log_info_extended(floginfoext, "meg-construction-end", (void*)est);

  MYTIME_STOP_PARALLEL(pt_alg);
// Destroy local configuration
  config_destroy(config);
}

static void
internal_get_EST_factorizations(pEST_info gen,
										  pEST_info est,
										  FILE* floginfoext,
										  pmytime pt_comp, pmytime pt_ccomp,
										  pconfiguration shared_config,
										  pext_array V,
										  pEST * pfactorized_est,
										  bool* is_timeout_expired) {

  pmytime_timeout pt_fact_timeout= MYTIME_timeout_create(shared_config->max_single_factorization_time);

  log_info_extended(floginfoext, "est-factorization-begin", (void*)est);
  MYTIME_START_PARALLEL(pt_comp);
  MYTIME_reset(pt_ccomp);
  MYTIME_start(pt_ccomp);
  *pfactorized_est= get_EST_factorizations(est, V, shared_config, gen, pt_fact_timeout);
  *is_timeout_expired= MYTIME_timeout_expired(pt_fact_timeout);
  if (*pfactorized_est != NULL) {
	 DEBUG("Computed %zu factorizations.", list_size((*pfactorized_est)->factorizations));
	 refine_EST_factorizations(gen, *pfactorized_est, shared_config);
	 remove_factorizations_with_very_small_exons((*pfactorized_est)->factorizations);
	 DEBUG("Remained %zu factorizations.", list_size((*pfactorized_est)->factorizations));
	 if(!list_is_empty((*pfactorized_est)->factorizations)) {
		remove_duplicated_factorizations((*pfactorized_est)->factorizations);
		my_assert(!list_is_empty((*pfactorized_est)->factorizations));
	 }
  } else {
	 DEBUG("Timeout expired!");
	 my_assert(*is_timeout_expired);
  }
  MYTIME_stop(pt_ccomp);
  MYTIME_STOP_PARALLEL(pt_comp);

  *is_timeout_expired= (*is_timeout_expired) || MYTIME_timeout_expired(pt_fact_timeout);
  MYTIME_timeout_destroy(pt_fact_timeout);
}

pEST
compute_est_fact(pEST_info gen,
					  pEST_info est,
					  LST_STree* tree,
					  ppreproc_gen pg,
					  FILE* floginfoext,
					  FILE* fmeg, FILE* fpmeg, FILE* ftmeg,
					  FILE* fintronic,
					  pmytime pt_alg, pmytime pt_comp, pmytime pt_io,
					  pconfiguration shared_config) {
  INFO("EST: %s", est->EST_id);

// Create local timers
  pmytime pt_ccomp= MYTIME_create_with_name("Internal Comp.");
  pmytime pt_meg= MYTIME_create_with_name("MEGs");

  size_t inc_pairing_len= 0;
  pext_array V= NULL;

  bool is_timeout_expired= false;

  pEST factorized_est= NULL;
  size_t prev_tot_pairings= 0, prev_tot_edges= 0;
  do {
	 size_t tot_pairings, tot_edges;
	 bool same_MEG_as_before= false;
	 do {
		same_MEG_as_before= false;
		build_meg(est, tree, pg, floginfoext, pt_alg, pt_meg, shared_config,
					 &inc_pairing_len, &V);

		MEG_stats(V, &tot_pairings, &tot_edges);
		same_MEG_as_before= prev_tot_pairings > 2 &&
		  prev_tot_edges > 0 &&
		  (prev_tot_pairings <= tot_pairings || prev_tot_edges <= tot_edges);
		if (same_MEG_as_before) {
		  WARN("The computed MEG is the same computed in the previous round. "
			  "Re-trying with longer with min-factor-len= %zd.",
			  shared_config->min_factor_len+inc_pairing_len+1);
		  ++inc_pairing_len;
		  DEBUG("Destroying the MEG and the occurrence set");
		  MYTIME_START_PARALLEL(pt_alg);
		  EA_destroy(V, (delete_function)vi_destroy);
		  MYTIME_STOP_PARALLEL(pt_alg);
		}
	 } while (same_MEG_as_before);
	 prev_tot_pairings= tot_pairings;
	 prev_tot_edges= tot_edges;

	 DEBUG("A possible MEG has been built. Trying to get the factorizations...");

	 is_timeout_expired= false;
	 internal_get_EST_factorizations(gen, est, floginfoext, pt_comp, pt_ccomp, shared_config,
												V, &factorized_est, &is_timeout_expired);

	 DEBUG("Timeout expired?        %s", (is_timeout_expired)?"YES":"no");
	 DEBUG("Factorization returned? %s", (factorized_est!=NULL &&
													  !list_is_empty(factorized_est->factorizations))?"YES":"no");
	 if (!is_timeout_expired ||
		  (factorized_est!=NULL && !list_is_empty(factorized_est->factorizations))) {
		DEBUG("EST factorization procedure correctly terminated "
				"(possibly without factorizations).");
		report_meg(est, pt_io, fmeg, V);
	 }

	 if (factorized_est!=NULL && !list_is_empty(factorized_est->factorizations)) {
		INFO("...EST aligned!");
		is_timeout_expired= false; // Reset timeout expiration
											// since we computed some factorizations
		fprintf(fintronic, ">%s\n", est->EST_id);
		add_intronic_edges_to_file(fintronic, V);
		write_single_EST_info(fpmeg, est);
		meg_write(fpmeg, V);
		fprintf(ftmeg, "%llu %llu %zu\n",
				  MYTIME_getinterval(pt_meg),
				  MYTIME_getinterval(pt_ccomp),
				  list_size(factorized_est->factorizations));
	 } else if (!is_timeout_expired) {
		INFO("...the EST %s has no alignment!", est->EST_gb);
	 } else {
		my_assert(is_timeout_expired);
		WARN("The timeout has expired while computing the EST factorizations. "
			  "Re-trying with longer with min-factor-len= %zd.",
			  shared_config->min_factor_len+inc_pairing_len+1);
		++inc_pairing_len;
	 }

	 log_info_extended(floginfoext, "est-factorization-end", (void*)est);

	 DEBUG("Destroying the MEG and the occurrence set");
	 MYTIME_START_PARALLEL(pt_alg);
	 EA_destroy(V, (delete_function)vi_destroy);
	 MYTIME_STOP_PARALLEL(pt_alg);

  } while (is_timeout_expired);

// Destroy local timers
  MYTIME_destroy(pt_meg);
  MYTIME_destroy(pt_ccomp);

  return factorized_est;
}
