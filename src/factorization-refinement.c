/**
 *
 *
 *                              PIntron
 *
 * A novel pipeline for computational gene-structure prediction based on
 * spliced alignment of expressed sequences (ESTs and mRNAs).
 *
 * Copyright (C) 2011,2012  Yuri Pirola
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

//#define LOG_THRESHOLD LOG_LEVEL_TRACE

#include "list.h"
#include "log.h"

#include "compute-alignments.h"
#include "refine.h"
#include "refine-intron.h"

// ------------------------------------------------------------------------
//
// Removing factorizations with very small exons
// Issue #22: Wrong call of procedure check_small_exons
//
#define _UB_VERY_SMALL_EXON_LENGTH_ 5

void
remove_factorizations_with_very_small_exons(plist factorizations) {
  DEBUG("Removing factorizations with very small exons...");
  plistit pl_f_it= list_first(factorizations);
  size_t fact_id= 0;
  while (listit_has_next(pl_f_it)) {
	 ++fact_id;
	 DEBUG("Analyzing factorization %zu...", fact_id);
	 pfactorization pfact= listit_next(pl_f_it);
	 plistit pl_factor_it= list_first(pfact);
	 bool has_very_small_exons= false;
	 while (!has_very_small_exons && listit_has_next(pl_factor_it)) {
		pfactor pcurr= listit_next(pl_factor_it);
		has_very_small_exons= (pcurr->EST_end + 1 - pcurr->EST_start) < _UB_VERY_SMALL_EXON_LENGTH_;
		if (has_very_small_exons) {
		  DEBUG("Exon (%d-%d) -- (%d-%d) is very small (less than %dbp). Factorization discarded.",
				  pcurr->EST_start, pcurr->EST_end, pcurr->GEN_start, pcurr->GEN_end,
				  _UB_VERY_SMALL_EXON_LENGTH_);
		}
	 }
	 listit_destroy(pl_factor_it);
	 if (has_very_small_exons) {
		INFO("Factorization %zu has very small exons. Deleting...", fact_id);
		list_remove_at_iterator(pl_f_it, (delete_function)factorization_destroy);
	 }
  }
  listit_destroy(pl_f_it);
}



// ------------------------------------------------------------------------
//
// Removing duplicated factorizations
// Issue #14: Duplicated factorizations reported by 'est-fact'
//

void
remove_duplicated_factorizations(plist factorizations) {
  DEBUG("Removing possible duplicated factorizations...");
  bool has_possible_duplicated= false;
  unsigned int members= 0;
  plistit pl_f_it= list_first(factorizations);
  while (!has_possible_duplicated && listit_has_next(pl_f_it)) {
	 DEBUG("Analyzing a new factorization...");
	 pfactorization pfact= listit_next(pl_f_it);
	 plistit pl_factor_it= list_first(pfact);
	 unsigned int hashfact= 1;
	 while (listit_has_next(pl_factor_it)) {
		pfactor pcurr= listit_next(pl_factor_it);
		unsigned int shift= (pcurr->EST_start + pcurr->EST_end +
									pcurr->GEN_start + pcurr->GEN_end) % (sizeof(hashfact)*CHAR_BIT);
		hashfact = (hashfact >> shift) | (hashfact << (sizeof(hashfact)*CHAR_BIT - shift));
	 }
	 has_possible_duplicated = (members & hashfact) != 0u;
	 members |= hashfact;
	 listit_destroy(pl_factor_it);
  }
  listit_destroy(pl_f_it);
  if (has_possible_duplicated) {
	 DEBUG("The computed factorizations could be duplicated. Performing a full-check...");
	 plistit pl_f_it1= list_first(factorizations);
	 size_t fact1_id= 0;
	 while (listit_has_next(pl_f_it1)) {
		++fact1_id;
		DEBUG("Analyzing factorization %zu...", fact1_id);
		pfactorization pfact1= listit_next(pl_f_it1);
		plistit pl_f_it2= list_first(factorizations);
		size_t fact2_id= 0;
		while (listit_has_next(pl_f_it2)) {
		  pfactorization pfact2= listit_next(pl_f_it2);
		  if (pfact1 == pfact2)
			 break;
		  ++fact2_id;
		  DEBUG("   ...comparing with factorization %zu...", fact2_id);
		  if (list_size(pfact1) != list_size(pfact2))
			 continue;
		  plistit pl_factor_it1= list_first(pfact1);
		  plistit pl_factor_it2= list_first(pfact2);
		  bool is_equal= true;
		  while (is_equal && listit_has_next(pl_factor_it1)) {
			 pfactor pcurr1= listit_next(pl_factor_it1);
			 pfactor pcurr2= listit_next(pl_factor_it2);
			 is_equal= (pcurr1->EST_start==pcurr2->EST_start) &&
				(pcurr1->EST_end==pcurr2->EST_end) &&
				(pcurr1->GEN_start==pcurr2->GEN_start) &&
				(pcurr1->GEN_end==pcurr2->GEN_end);
		  }
		  listit_destroy(pl_factor_it1);
		  listit_destroy(pl_factor_it2);
		  if (is_equal) {
			 INFO("Factorizations %zu and %zu are duplicated. Deleting factorization %zu...",
					fact1_id, fact2_id, fact1_id);
			 list_remove_at_iterator(pl_f_it1, (delete_function)factorization_destroy);
			 break;
		  }
		}
		listit_destroy(pl_f_it2);
	 }
	 listit_destroy(pl_f_it1);
  } else {
	 DEBUG("The computed factorizations are not duplicated. Continuing...");
  }
}



// ------------------------------------------------------------------------
//
// Removing false small exons
// Issue #3: Unnecessary exons
//

// Exons longer than this are always considered NOT SMALL
#define _UB_SMALL_EXON_LENGTH_ 100
// Lenght of prefixes/suffixes considered during alignment
#define _AFFIXES_LENGTH_ 5

static
void
analyze_possibly_small_exon(pfactor * ppprev,
									 pfactor * ppcurr,
									 pfactor pnext,
									 pfactorization pfact,
									 pEST_info genomic,
									 pEST factorized_est,
									 pconfiguration config) {
  pfactor pprev= *ppprev;
  pfactor pcurr= *ppcurr;
// Skip external exons
  if (pprev==NULL || pnext==NULL) {
	 DEBUG("  ...skipped external factor %d-%d %d-%d.",
			 pcurr->EST_start, pcurr->EST_end,
			 pcurr->GEN_start, pcurr->GEN_end);
	 return;
  }
  DEBUG("  ...analyzing factor %d-%d %d-%d...",
		  pcurr->EST_start, pcurr->EST_end,
		  pcurr->GEN_start, pcurr->GEN_end);
  my_assert(pcurr->EST_start <= pcurr->EST_end);
  my_assert(pcurr->GEN_start <= pcurr->GEN_end);
  const size_t elen= pcurr->EST_end + 1 - pcurr->EST_start;
  const size_t glen= pcurr->GEN_end + 1 - pcurr->GEN_start;
  if (//(elen <= config->min_factor_len) &&
		(elen <= _UB_SMALL_EXON_LENGTH_)) {
	 DEBUG("     it is a small exon! Trying to remove it...");
	 const char* const efact= factorized_est->info->EST_seq + pcurr->EST_start;
	 const char* const gfact= genomic->EST_seq + pcurr->GEN_start;

	 const size_t orig_ed= compute_edit_distance(efact, elen, gfact, glen);
	 DEBUG("        original edit distance: %zu", orig_ed);

// Take a short suffix of the previous exon and a short prefix of the following exon
	 const size_t estart= MAX(pprev->EST_start, pprev->EST_end + 1 - _AFFIXES_LENGTH_);
	 const size_t eend= MIN(pnext->EST_end + 1, pnext->EST_start + _AFFIXES_LENGTH_);
	 const size_t epreflen= pprev->EST_end + 1 - estart;
	 const size_t esufflen= eend - pnext->EST_start;
	 const size_t allelen= eend - estart;
	 const char* const allefact= factorized_est->info->EST_seq + estart;
	 my_assert(epreflen + esufflen + elen == allelen);
	 const size_t gstart= MAX(pprev->GEN_start, pprev->GEN_end + 1 - _AFFIXES_LENGTH_);
	 const size_t gend= MIN(pnext->GEN_end + 1, pnext->GEN_start + _AFFIXES_LENGTH_);
	 const size_t gpreflen= pprev->GEN_end + 1 - gstart;
	 const size_t gsufflen= gend - pnext->GEN_start;
	 const size_t allglen= gend - gstart;
	 const char* const allgfact= genomic->EST_seq + gstart;
	 my_assert(gpreflen + gsufflen + glen + pcurr->GEN_start - pprev->GEN_end - 1 +
				  pnext->GEN_start - pcurr->GEN_end - 1 == allglen);

	 TRACE("  Considered EST prefix:     %.*s", (int)epreflen, allefact);
	 TRACE("  Considered genomic prefix: %.*s", (int)gpreflen, allgfact);
	 const size_t orig_ed_pref= compute_edit_distance(allefact, epreflen, allgfact, gpreflen);
	 DEBUG("  pref. original edit distance: %zu", orig_ed_pref);

	 TRACE("  Considered EST suffix:     %.*s", (int)esufflen, allefact - esufflen);
	 TRACE("  Considered genomic suffix: %.*s", (int)gsufflen, allgfact - gsufflen);
	 const size_t orig_ed_suff= compute_edit_distance(allefact - esufflen, esufflen,
																				allgfact - gsufflen, gsufflen);
	 DEBUG("  suff. original edit distance: %zu", orig_ed_suff);

	 size_t offset_p, offset_t1, offset_t2;
	 unsigned int new_ed= allglen;
	 TRACE("Factor of the EST sequence [%9zu-%9zu]: %.*s",
			 estart, eend, (int)allelen, allefact);
	 TRACE("Factor of the genomic seq  [%9zu-%9zu]: %.*s",
			 gstart, gend, (int)allglen, allgfact);
	 const bool ref_res= refine_borders(allefact, allelen,
													allgfact, allglen,
													orig_ed + orig_ed_pref + orig_ed_suff,
													&offset_p, &offset_t1, &offset_t2,
													&new_ed);
	 TRACE("Border refinement gave: "
			 "Success? %s  "
			 "offset_p: %zu  offset_t1: %zu  offset_t2: %zu  "
			 "new_edit_distance: %u",
			 ref_res?"OK":"NO", offset_p, offset_t1, offset_t2, new_ed);

	 if (ref_res) {
		INFO("Found a possibly better placement for exon %d-%d %d-%d.",
			  pcurr->EST_start, pcurr->EST_end,
			  pcurr->GEN_start, pcurr->GEN_end);
		my_assert(new_ed <= orig_ed + orig_ed_pref + orig_ed_suff);
		my_assert(offset_p <= allelen);
		my_assert(offset_t1 <= offset_t2);
		my_assert(offset_t2 <= allglen);
		INFO("Checking if the new intron would be not worse than the previous two introns...");
		const double prev_avg_burset_freq=
		  (getBursetFrequency_adaptor(genomic->EST_seq, pprev->GEN_end+1, pcurr->GEN_start) +
			getBursetFrequency_adaptor(genomic->EST_seq, pcurr->GEN_end+1, pnext->GEN_start))/2.0;
		const double new_burset_freq=
		  getBursetFrequency_adaptor(genomic->EST_seq,
											  gstart + offset_t1,
											  gend - allglen + offset_t2);
		DEBUG("Previous average 'Burset score': %.1f", prev_avg_burset_freq);
		DEBUG("New 'Burset score': %.1f", new_burset_freq);
		if (new_burset_freq >= prev_avg_burset_freq) {
		  INFO("The new placement does not decrease intron score, thus the exon will be split "
				 "(old score: %.1f, new score %.1f).", prev_avg_burset_freq, new_burset_freq);
		  INFO("Previous 'local' factorization on EST:     (%8d -%8d) (%8d -%8d) (%8d -%8d).",
				 pprev->EST_start, pprev->EST_end,
				 pcurr->EST_start, pcurr->EST_end,
				 pnext->EST_start, pnext->EST_end);
		  INFO("Previous 'local' factorization on genomic: (%8d -%8d) (%8d -%8d) (%8d -%8d).",
				 pprev->GEN_start, pprev->GEN_end,
				 pcurr->GEN_start, pcurr->GEN_end,
				 pnext->GEN_start, pnext->GEN_end);
		  pprev->EST_end= estart + offset_p - 1;
		  pnext->EST_start= eend + offset_p - allelen;
		  pprev->GEN_end= gstart + offset_t1 - 1;
		  pnext->GEN_start= gend + offset_t2 - allglen;
		  INFO("Modified 'local' factorization on EST:     (%8d -%8d)   ... skipped ...    (%8d -%8d).",
				 pprev->EST_start, pprev->EST_end,
				 pnext->EST_start, pnext->EST_end);
		  INFO("Modified 'local' factorization on genomic: (%8d -%8d)   ... skipped ...    (%8d -%8d).",
				 pprev->GEN_start, pprev->GEN_end,
				 pnext->GEN_start, pnext->GEN_end);
		  list_remove_element(pfact, pcurr, (delete_function)factor_destroy);
		  *ppcurr= pprev;
		  *ppprev= NULL;
		} else {
		  INFO("The new placement decreases intron score, thus the exon will NOT be split.");
		}
	 }
  }
}

static
void
remove_false_small_exons(pEST_info genomic,
								 pEST factorized_est,
								 pconfiguration config) {
  DEBUG("Removing 'false' small exons...");
  plistit pl_f_it= list_first(factorized_est->factorizations);
  while (listit_has_next(pl_f_it)) {
	 DEBUG("Analyzing a factorization...");
	 pfactorization pfact= listit_next(pl_f_it);

	 pfactor pprev= NULL, pcurr= NULL, pnext= NULL;
	 plistit pl_factor_it= list_first(pfact);
	 if (listit_has_next(pl_factor_it))
		pnext= listit_next(pl_factor_it);
	 while (pnext != NULL) {
		pprev= pcurr;
		pcurr= pnext;
		pnext= NULL;
		if (listit_has_next(pl_factor_it))
		  pnext= listit_next(pl_factor_it);
		analyze_possibly_small_exon(&pprev, &pcurr, pnext, pfact,
											 genomic, factorized_est, config);
	 }
	 listit_destroy(pl_factor_it);
  }
  listit_destroy(pl_f_it);
}

// ------------------------------------------------------------------------
//
// Recovering "lost" prefixes and/or suffixes
// Issue: #4
//

#define _MAX_ERROR_RATE_ 0.17   //approx. 1 error every 6 bases

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
  remove_false_small_exons(genomic, factorized_est, config);
}

