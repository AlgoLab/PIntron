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

// Create a local copy of configuration parameters
  pconfiguration config= config_clone(shared_config);
// Create local timers
  pmytime pt_ccomp= MYTIME_create_with_name("Internal Comp.");
  pmytime pt_meg= MYTIME_create_with_name("MEGs");

  MYTIME_START_PARALLEL(pt_alg);
  DEBUG("Calculating the occurrence set");

  pext_array V= NULL;
  bool too_complex= false;
  size_t inc_pairing_len= 0;
  log_info_extended(floginfoext, "meg-construction-begin", (void*)est);
  do {
	 DEBUG("Building the MEG vertex set");

	 config->min_factor_len += inc_pairing_len;
	 V= build_vertex_set(est, tree, pg, config);
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
  } while(too_complex);
  log_info_extended(floginfoext, "meg-construction-end", (void*)est);

  MYTIME_STOP_PARALLEL(pt_alg);
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

	 // EST-FACTORIZATIONS

  log_info_extended(floginfoext, "est-factorization-begin", (void*)est);

  MYTIME_START_PARALLEL(pt_comp);
  MYTIME_reset(pt_ccomp);
  MYTIME_start(pt_ccomp);
  pEST factorized_est=get_EST_factorizations(est, V, config, gen);
  DEBUG("Computed %zu factorizations.", list_size(factorized_est->factorizations));
  refine_EST_factorizations(gen, factorized_est, config);
  MYTIME_stop(pt_ccomp);
  MYTIME_STOP_PARALLEL(pt_comp);

  if(!list_is_empty(factorized_est->factorizations)) {
	 INFO("...EST aligned!");
	 remove_duplicated_factorizations(factorized_est->factorizations);
	 my_assert(!list_is_empty(factorized_est->factorizations));
	 fprintf(fintronic, ">%s\n", est->EST_id);
	 add_intronic_edges_to_file(fintronic, V);
	 write_single_EST_info(fpmeg, est);
	 meg_write(fpmeg, V);
	 fprintf(ftmeg, "%llu %llu %zu\n",
				MYTIME_getinterval(pt_meg),
				MYTIME_getinterval(pt_ccomp),
				list_size(factorized_est->factorizations));
  } else {
	 INFO("...the EST %s has no alignment!", est->EST_gb);
  }

  EST_destroy_just_factorizations(factorized_est);

  log_info_extended(floginfoext, "est-factorization-end", (void*)est);

  DEBUG("Destroying the MEG and the occurrence set");
  MYTIME_START_PARALLEL_REUSE(pt_alg);
  EA_destroy(V, (delete_function)vi_destroy);
  MYTIME_STOP_PARALLEL(pt_alg);

// Destroy local configuration
  config_destroy(config);
// Destroy local timers
  MYTIME_destroy(pt_meg);
  MYTIME_destroy(pt_ccomp);

  return factorized_est;
}
