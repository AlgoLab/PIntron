/**
 *
 *
 *                              PIntron
 *
 * A novel pipeline for computational gene-structure prediction based on
 * spliced alignment of expressed sequences (ESTs and mRNAs).
 *
 * Copyright (C) 2011  Yuri Pirola
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
/*
** factorization-refinement.c
*/

#include "factorization-refinement.h"

#include "list.h"
#include "log.h"

// ------------------------------------------------------------------------
//
// Recovering "lost" prefixes and/or suffixes
// Issue: #4
//

#define _MAX_ERROR_RATE_ 0.17   //approx. 1 error every 6 bases

#include "compute-alignments.h"

static
bool
find_longest_affix(char* est, size_t estl,
						 char* genomic, size_t genomicl,
						 size_t* out_est_cut, size_t* out_genomic_cut) {
  size_t* matrix= edit_distance_matrix(est, estl, genomic, genomicl);
// Get the best genomic cut for each EST cut
  bool valid_cut= false;
  size_t best_gcut= 0;
  size_t best_ecut= 0;
  double best_weight= 1.0;
  size_t base= genomicl+1;
  for (size_t ecut= 1; ecut < estl+1; ++ecut) {
	 for (size_t gcut= 1; gcut < genomicl+1; ++gcut) {
		const double cut_weight= 2.0*((double)matrix[base+gcut])/(ecut+gcut);
		if ((est[ecut-1] == genomic[gcut-1]) &&
			 (cut_weight <= _MAX_ERROR_RATE_) &&
			 (cut_weight <= best_weight)) {
		  best_gcut= gcut;
		  best_ecut= ecut;
		  best_weight= cut_weight;
		  valid_cut= true;
		}
	 }
	 base+= (genomicl+1);
  }
  pfree(matrix);
  if (valid_cut) {
	 DEBUG("We recovered the following almost-matching factors: "
			 "EST(%zunt): %.*s   GENOMIC(%zunt): %.*s",
			 best_ecut, (int)best_ecut, est,
			 best_gcut, (int)best_gcut, genomic);
	 *out_est_cut= best_ecut;
	 *out_genomic_cut= best_gcut;
  } else {
	 DEBUG("No previously discarded external factors have been recovered.");
  }
  return valid_cut;
}


static
void
recover_lost_prefixes_and_suffixes(pEST_info genomic,
											  pEST factorized_est,
											  pconfiguration config) {
  DEBUG("Recovering possible lost prefixes and/or suffixes...");
  const size_t totglen= strlen(genomic->EST_seq);
  const size_t totelen= strlen(factorized_est->info->EST_seq);
  plistit pl_f_it= list_first(factorized_est->factorizations);
  while (listit_has_next(pl_f_it)) {
	 DEBUG("Analyzing a factorization...");
	 pfactorization pfact= listit_next(pl_f_it);
	 pfactor pff= list_head(pfact);
	 DEBUG("  ...analyzing the prefix of the first factor %d-%d %d-%d...",
			 pff->EST_start, pff->EST_end,
			 pff->GEN_start, pff->GEN_end);
	 if ((pff->EST_start > 0) && (pff->GEN_start > 0)) {
// Computing the factors
		const size_t flen= MIN(pff->EST_start, pff->GEN_start);
		const size_t elen= MIN(pff->EST_start, (int)((1.0+_MAX_ERROR_RATE_)*flen));
		const size_t glen= MIN(pff->GEN_start, (int)((1.0+_MAX_ERROR_RATE_)*flen));
		char* efact= c_palloc(elen+1);
		strncpy(efact, factorized_est->info->EST_seq+pff->EST_start-elen, elen);
		efact[elen]= '\0';
		for (size_t i= 0; i<elen/2; ++i) MY_SWAP(char, efact[i], efact[elen-i-1]);
		char* gfact= c_palloc(glen+1);
		strncpy(gfact, genomic->EST_seq+pff->GEN_start-glen, glen);
		gfact[glen]= '\0';
		for (size_t i= 0; i<glen/2; ++i) MY_SWAP(char, gfact[i], gfact[glen-i-1]);
// If we find that efact[0] == gfact[0], then the factor was already refined. Skip!
		if (efact[0] != gfact[0]) {
		  size_t ecut, gcut;
		  const bool valid_cut= find_longest_affix(efact, elen, gfact, glen,
																 &ecut, &gcut);
		  if (valid_cut) {
			 INFO("A previously discarded prefix has been recovered:");
			 INFO("OLD: %d-%d %d-%d   NEW: %d-%d %d-%d",
					pff->EST_start, pff->EST_end,
					pff->GEN_start, pff->GEN_end,
					pff->EST_start - (int)ecut, pff->EST_end,
					pff->GEN_start - (int)gcut, pff->GEN_end);
			 pff->EST_start -= ecut;
			 pff->GEN_start -= gcut;
		  }
		}
		pfree(efact);
		pfree(gfact);
	 }
// EST_end and GEN_end are inclusive!
	 pfactor pfl= list_tail(pfact);
	 DEBUG("  ...analyzing the suffix of the last factor %d-%d %d-%d...",
			 pfl->EST_start, pfl->EST_end,
			 pfl->GEN_start, pfl->GEN_end);
	 my_assert(pfl->GEN_end >= 0 && (unsigned int)pfl->GEN_end <= totglen);
	 my_assert(pfl->EST_end >= 0 && (unsigned int)pfl->EST_end <= totelen);
	 if (((totelen - pfl->EST_end) > 1) &&
		  ((totglen - pfl->GEN_end) > 1)) {
// Computing the factors
		const size_t flen= MIN(totelen - pfl->EST_end - 1, totglen - pfl->GEN_end - 1);
		const size_t elen= MIN(totelen - pfl->EST_end - 1, (int)(1.0+_MAX_ERROR_RATE_)*flen);
		const size_t glen= MIN(totglen - pfl->GEN_end - 1, (int)(1.0+_MAX_ERROR_RATE_)*flen);
		my_assert(glen>0);
		my_assert(elen>0);
		char* efact= c_palloc(elen+1);
		strncpy(efact, factorized_est->info->EST_seq+pfl->EST_end, elen);
		efact[elen]= '\0';
		char* gfact= c_palloc(glen+1);
		strncpy(gfact, genomic->EST_seq+pfl->GEN_end, glen);
		gfact[glen]= '\0';
// If we find that efact[0] == gfact[0], then the factor was already refined. Skip!
		if (efact[0] != gfact[0]) {
		  size_t ecut, gcut;
		  const bool valid_cut= find_longest_affix(efact, elen, gfact, glen,
																 &ecut, &gcut);
		  if (valid_cut) {
			 INFO("A previously discarded suffix has been recovered:");
			 INFO("OLD: %d-%d %d-%d   NEW: %d-%d %d-%d",
					pfl->EST_start, pfl->EST_end,
					pfl->GEN_start, pfl->GEN_end,
					pfl->EST_start, pfl->EST_end + (int)ecut,
					pfl->GEN_start, pfl->GEN_end + (int)gcut);
			 pfl->EST_end += ecut;
			 pfl->GEN_end += gcut;
		  }
		}
		pfree(efact);
		pfree(gfact);
	 }
  }
  listit_destroy(pl_f_it);
}



void
refine_EST_factorizations(pEST_info genomic,
								  pEST factorized_est,
								  pconfiguration config) {
  INFO("Further refinement of EST factorizations...");
  recover_lost_prefixes_and_suffixes(genomic, factorized_est, config);
}

