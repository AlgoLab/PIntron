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

#ifdef LOG_TRACE_ENABLED

#define DOT_100_STRING "...................................................................................................."

#define DOT_1K_STRING DOT_100_STRING DOT_100_STRING DOT_100_STRING DOT_100_STRING DOT_100_STRING DOT_100_STRING DOT_100_STRING DOT_100_STRING DOT_100_STRING DOT_100_STRING

#define DOT_10K_STRING DOT_1K_STRING DOT_1K_STRING DOT_1K_STRING DOT_1K_STRING DOT_1K_STRING DOT_1K_STRING DOT_1K_STRING DOT_1K_STRING DOT_1K_STRING DOT_1K_STRING

#define DOT_100K_STRING DOT_10K_STRING DOT_10K_STRING DOT_10K_STRING DOT_10K_STRING DOT_10K_STRING DOT_10K_STRING DOT_10K_STRING DOT_10K_STRING DOT_10K_STRING DOT_10K_STRING

#define DOT_1M_STRING DOT_100K_STRING DOT_100K_STRING DOT_100K_STRING DOT_100K_STRING DOT_100K_STRING DOT_100K_STRING DOT_100K_STRING DOT_100K_STRING DOT_100K_STRING DOT_100K_STRING

#define DOT_STRING DOT_1M_STRING

#endif

// ------------------------------------------------------------------------
//
// Common definitions
//
// Exons shorter than that cause the factorization to be discarded
#define _UB_VERY_SMALL_EXON_LENGTH_ 2
// Perfect matches shorter than this are never considered as exons
#define _LB_SMALL_EXON_LENGTH_ 3
// Exons longer than this are always considered NOT SMALL
#define _UB_SMALL_EXON_LENGTH_ 15
#define _UB_MED_EXON_LENGTH_ 100
// Lenght of prefixes/suffixes considered during alignment
#define _AFFIXES_LENGTH_ 5
// Maximum error rate allowed in factors (i.e. "lost" prefixes/suffixes) to be added
#define _MAX_ERROR_RATE_ 0.17   //approx. 1 error every 6 bases



// ------------------------------------------------------------------------
//
// Removing factorizations with very small exons
// Issue #22: Wrong call of procedure check_small_exons
//

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
		has_very_small_exons= (pcurr->EST_end + 1 - pcurr->EST_start) <= _UB_VERY_SMALL_EXON_LENGTH_;
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
// Search for potentially true small exons
// Issue #11: Small exons not found
// Issue #28: Small external exons are not detected
//

#include "util.h"

#include "lst_structs.h"
#include "lst_stree.h"
#include "lst_string.h"

struct struct_lcs_visit_node;

typedef struct struct_lcs_visit_node {
  LST_Node* node;
  SLIST_ENTRY(struct_lcs_visit_node) visit_list;
} lcs_visit_node;

static
void
find_longest_common_factor_dp(const char* const restrict s1, const size_t l1,
										const char* const restrict s2, const size_t l2,
										size_t* const pocc1, size_t* const pocc2,
										size_t* const plen) {
  DEBUG("Searching for a longest common factor using a dynamic programming algorithm...");
  if (l2 > l1) {
	 find_longest_common_factor_dp(s2, l2, s1, l1, pocc2, pocc1, plen);
  }
  size_t* curr= NPALLOC(size_t, l2+1);
  size_t* prev= NPALLOC(size_t, l2+1);
  size_t* swap= NULL;
  *pocc1= 0;
  *pocc2= 0;
  *plen= 0;

  for (size_t i2= 0; i2<=l2; ++i2) {
	 prev[i2]= 0;
  }
  for (size_t i1= 0; i1<l1; ++i1) {
	 curr[0]= 0;
	 for (size_t i2= 0; i2<l2; ++i2) {
		if (s1[i1] != s2[i2]) {
		  curr[i2+1]= 0;
		} else {
		  curr[i2+1]= 1 + prev[i2];
		}
		if (*plen < curr[i2+1]) {
		  *plen= curr[i2+1];
		  *pocc1= i1+1-*plen;
		  *pocc2= i2+1-*plen;
		}
	 }
	 swap= curr;
	 curr= prev;
	 prev= swap;
  }
  DEBUG("Found a %zu-long common substring.", *plen);
  DEBUG("S1 (%8zu): %.2s>%.*s<%.2s", *pocc1,
		  (*pocc1>=2)? s1+*pocc1-2 : "  ",
		  (int)*plen, s1+*pocc1,
		  s1+*pocc1+*plen);
  DEBUG("S2 (%8zu): %.2s>%.*s<%.2s", *pocc2,
		  (*pocc2>=2)? s2+*pocc2-2 : "  ",
		  (int)*plen, s2+*pocc2,
		  s2+*pocc2+*plen);
  pfree(curr);
  pfree(prev);
}


static
void
find_longest_common_factor(const char* const s1, const size_t l1,
									const char* const s2, const size_t l2,
									const size_t cut_position2,
									size_t* const pocc1, size_t* const pocc2,
									size_t* const plen) {
  DEBUG("Searching for the longest common factor...");
  my_assert((cut_position2 > 0) && (cut_position2 < l2));
  *pocc1= *pocc2= *plen= 0;
  LST_StringSet *set= lst_stringset_new();
  LST_String * lst1= PALLOC(LST_String);
  lst_string_init(lst1, s1, sizeof(char), l1);
  lst_stringset_add(set, lst1);
  LST_String * lst2= PALLOC(LST_String);
  lst_string_init(lst2, s2, sizeof(char), l2);
  lst_stringset_add(set, lst2);
  LST_STree* tree = lst_stree_new(set);
// Visit
  lcs_visit_node* vnode;
  LST_Edge *edge;
  SLIST_HEAD(qhead, struct_lcs_visit_node) stack;
  SLIST_INIT(&stack);
  vnode= PALLOC(lcs_visit_node);
  vnode->node= tree->root_node;
  vnode->node->single_char= (char)0;
  vnode->node->string_depth= 0;
  vnode->node->slices= NULL;
  SLIST_INSERT_HEAD(&stack, vnode, visit_list);

  size_t max_comm_substr_len= 0;
  LST_Node* max_comm_substr= NULL;
  while (!SLIST_EMPTY(&stack)) {
	 vnode= SLIST_FIRST(&stack);
	 if ((vnode->node->kids.lh_first==NULL)|| // Is a leaf or it is the second time that I visit the (internal) node
		  (vnode->node->single_char != (char)0)) {
		my_assert(((char)0 <= vnode->node->single_char) &&
					 (vnode->node->single_char < (char)8));
		SLIST_REMOVE_HEAD(&stack, visit_list);
		if (vnode->node->up_edge != NULL) {
		  my_assert(vnode->node->up_edge->range.string == lst1 ||
						vnode->node->up_edge->range.string == lst2);
		  if (vnode->node->kids.lh_first==NULL) {
			 my_assert(vnode->node->single_char == (char)0);
			 if (vnode->node->up_edge->range.string==lst1) {
				vnode->node->single_char= (char)1;
			 } else {
// Check that the common substring crosses the cut on s2
				const size_t this_pos= *(vnode->node->up_edge->range.end_index)+1-vnode->node->string_depth;
				if (this_pos <= cut_position2) {
				  LST_Node* up_node= vnode->node;
				  while ((up_node != NULL) &&
							(this_pos + up_node->string_depth >= cut_position2) &&
							!(up_node->single_char & (char)2)) {
					 up_node->single_char= (char)2;
					 up_node= (up_node->up_edge!=NULL)?(up_node->up_edge->src_node):NULL;
				  }
				}
			 }
		  }
		  const char this_node= vnode->node->single_char;
// Propagate only s1 (since s2 has been already propagated)
		  vnode->node->up_edge->src_node->single_char |= ((char)1 & this_node);
		  if (this_node == (char)3) {
// Common node
			 const size_t new_comm_substr_len= vnode->node->string_depth;
			 if (new_comm_substr_len>max_comm_substr_len) {
				TRACE("Found a longer common substring... (length: %zu)", new_comm_substr_len);
				max_comm_substr_len= new_comm_substr_len;
				max_comm_substr= vnode->node;
			 }
		  }
		}
		pfree(vnode);
	 } else {
		vnode->node->single_char |= (char)4;
		const size_t string_depth= vnode->node->string_depth;
		for (edge = vnode->node->kids.lh_first;
			  edge != NULL;
			  edge = edge->siblings.le_next) {
		  vnode= PALLOC(lcs_visit_node);
		  vnode->node= edge->dst_node;
		  vnode->node->single_char= (char)0;
		  vnode->node->string_depth= string_depth+lst_edge_get_length(edge);
		  vnode->node->slices= NULL;
		  SLIST_INSERT_HEAD(&stack, vnode, visit_list);
		}
	 }
  }
  DEBUG("The longest common substring has length: %zu", max_comm_substr_len);
#ifdef LOG_DEBUG_ENABLED
  {
	 char* common= c_palloc(max_comm_substr_len+1);
	 common[max_comm_substr_len]= '\0';
	 LST_Node* nt= max_comm_substr;
	 size_t cursor= max_comm_substr_len;
	 while (nt != NULL && nt->up_edge != NULL) {
		cursor -= lst_edge_get_length(nt->up_edge);
		strncpy(common + cursor,
				  (char*)(nt->up_edge->range.string->data) + nt->up_edge->range.start_index,
				  lst_edge_get_length(nt->up_edge));
		nt= nt->up_edge->src_node;
	 }
	 DEBUG("The longest common substring is: %s", common);
	 pfree(common);
  }
#endif

  if (max_comm_substr_len==0) {
	 DEBUG("The longest common substring is empty. Returning...");
	 return;
  }

// Searching the first occurrence on the two strings
  size_t pos[2]= { 0, 0};
  SLIST_INIT(&stack);
  vnode= PALLOC(lcs_visit_node);
  vnode->node= max_comm_substr;
  SLIST_INSERT_HEAD(&stack, vnode, visit_list);

  while (!SLIST_EMPTY(&stack)) {
	 vnode= SLIST_FIRST(&stack);
	 SLIST_REMOVE_HEAD(&stack, visit_list);
	 if (vnode->node->kids.lh_first==NULL) {
		NOT_NULL(vnode->node->up_edge);
		NOT_NULL(max_comm_substr->up_edge);
		my_assert(vnode->node->up_edge->range.string == lst1 ||
					 vnode->node->up_edge->range.string == lst2);
		pos[(vnode->node->up_edge->range.string==lst1?0:1)]= *(vnode->node->up_edge->range.end_index)-vnode->node->string_depth+1;
	 } else {
		for (edge = vnode->node->kids.lh_first;
			  edge != NULL;
			  edge = edge->siblings.le_next) {
		  lcs_visit_node* ivnode= PALLOC(lcs_visit_node);
		  ivnode->node= edge->dst_node;
		  SLIST_INSERT_HEAD(&stack, ivnode, visit_list);
		}
	 }
	 pfree(vnode);
  }
  DEBUG("New maximal pairing: ( %zu, %zu, %zu).",
		  pos[1], pos[0], max_comm_substr_len);
  *pocc1= pos[0];
  *pocc2= pos[1];
  *plen= max_comm_substr_len;
  lst_stree_free(tree);
  pfree(set);
}

static
bool
is_canonical_intron(const char* const GEN_seq,
						  int intron_start, int intron_end) {
  return ( (GEN_seq[intron_start  ]=='G' &&
				GEN_seq[intron_start+1]=='T' &&
				GEN_seq[intron_end-1]=='A' &&
				GEN_seq[intron_end  ]=='G') ||
			  (GEN_seq[intron_start  ]=='g' &&
				GEN_seq[intron_start+1]=='t' &&
				GEN_seq[intron_end-1]=='a' &&
				GEN_seq[intron_end  ]=='g') );
}

static
void
search_small_exon_at_prefix(pfactor * pp1,
									 plistit pfactit,
									 pfactorization pfact,
									 pEST_info genomic,
									 pEST factorized_est,
									 pconfiguration config) {
  NOT_NULL(pp1);
  NOT_NULL(*pp1);
  pfactor p1= *pp1;
  DEBUG("  ...analyzing the prefix before [%d-%d, %d-%d]...",
		  p1->EST_start, p1->EST_end,
		  p1->GEN_start, p1->GEN_end);
  my_assert(p1->EST_start <= p1->EST_end);
  my_assert(p1->GEN_start <= p1->GEN_end);
  const size_t e1len= p1->EST_end + 1 - p1->EST_start;
  const size_t g1len= p1->GEN_end + 1 - p1->GEN_start;
// Check if it is possible to decompose the prefix and the first exon into one small exon and one normal exon
  if ((e1len + p1->EST_start) >= (_LB_SMALL_EXON_LENGTH_ + _UB_SMALL_EXON_LENGTH_)) {
	 const size_t eplen= MIN(MIN(p1->EST_start, p1->GEN_start), 2*_UB_SMALL_EXON_LENGTH_);
	 const char* const epfact= factorized_est->info->EST_seq + p1->EST_start - eplen;

	 const char* const EST_seq= factorized_est->info->EST_seq;
	 const char* const GEN_seq= genomic->EST_seq;

	 const size_t e1plen= MIN(MIN(e1len, g1len), _UB_SMALL_EXON_LENGTH_);
	 const char* const e1pfact= EST_seq + p1->EST_start;
	 const size_t g1plen= MIN(MIN(e1len, g1len), _UB_SMALL_EXON_LENGTH_);
	 const char* const g1pfact= GEN_seq + p1->GEN_start;
	 TRACE("EST discarded prefix: %.*s", (int)eplen, epfact);

// Search the longest common factor of epfact and g1sfact
	 size_t pg, pe, cflen;
	 find_longest_common_factor_dp(GEN_seq, p1->GEN_start,
											 epfact, eplen,
											 &pg, &pe, &cflen);
	 if (cflen >= _LB_SMALL_EXON_LENGTH_) {
		DEBUG("Found a possible small exon of length (at most) %zu.", cflen);
		TRACE("First border refinement.");
		TRACE("Factor of the EST sequence     [pos=%7zu]: %.*s %.*s",
				pe, (pg<=3)?(int)pg:3, "   ", (int)cflen, epfact+pe);
		TRACE("Factor of the genomic sequence [pos=%7zu]: %.*s>%.*s<%.3s",
				pg, (pg<=3)?(int)pg:3, GEN_seq+pg-((pg<=3)?pg:3),
				(int)cflen, GEN_seq+pg,
				GEN_seq+pg+cflen);
		TRACE("EST 1st exon prefix:  %.*s", (int)e1plen, e1pfact);
		TRACE("GEN 1st exon prefix:  %.*s", (int)g1plen, g1pfact);
		const unsigned int edp= compute_edit_distance(e1pfact, e1plen,
																	 g1pfact, g1plen);
		TRACE("Edit distance= %u", edp);
		const size_t allelen= MIN(p1->EST_end+1,
										  p1->EST_start+_UB_SMALL_EXON_LENGTH_)-pe;
		const size_t allglen= MIN(p1->GEN_end+1,
										  p1->GEN_start+_UB_SMALL_EXON_LENGTH_)-pg;
		size_t offset_p, offset_t1, offset_t2;
		unsigned int new_ed;
		const bool ref_res=
		  general_refine_borders(EST_seq+pe, allelen,
										 _LB_SMALL_EXON_LENGTH_, allelen-_LB_SMALL_EXON_LENGTH_,
										 GEN_seq+pg, allglen,
										 edp,
										 &offset_p, &offset_t1, &offset_t2,
										 &new_ed);
		TRACE("Border refinement gave: "
				"Success? %s  "
				"offset_p: %zu  offset_t1: %zu  offset_t2: %zu  "
				"new_edit_distance: %u",
				ref_res?"OK":"NO", offset_p, offset_t1, offset_t2, new_ed);

		if (!ref_res) {
		  DEBUG("The border alignments of the new intron is not good enough. Small exon discarded.");
		} else if ((int)offset_t2 - (int)offset_t1 < config->min_intron_length) {
		  DEBUG("The new intron induced by the small exon is too short (%zunt). "
				  "Small exon discarded.", offset_t2 - offset_t1);
		} else if ((int)offset_t2 - (int)offset_t1 < config->min_intron_length) {
		  DEBUG("The new intron induced by the small exon is too short (%zunt). "
				  "Small exon discarded.", offset_t2 - offset_t1);
		} else if (!is_canonical_intron(GEN_seq, pg + offset_t1, pg + offset_t2 - 1)) {
		  DEBUG("The new intron induced by the small exon is not canonical (%.2s..%.2s). "
				  "Small exon discarded.",
				  GEN_seq + pg + offset_t1, GEN_seq + pg + offset_t2 - 2);
		} else if (offset_p - pe < _LB_SMALL_EXON_LENGTH_) {
		  DEBUG("The new small exon is too short (%zunt). "
				  "Small exon discarded.", offset_p);
		} else {
		  INFO("A new external small exon has been found.");
		  DEBUG("Previous 'local' factorization on EST:                           [%8d-%8d]",
				  p1->EST_start, p1->EST_end);
		  DEBUG("Previous 'local' factorization on genomic:                       [%8d-%8d]",
				  p1->GEN_start, p1->GEN_end);
		  pfactor pnew= factor_create();
		  pnew->EST_start= pe;
		  pnew->EST_end  = pe + offset_p - 1;
		  pnew->GEN_start= pg;
		  pnew->GEN_end  = pg + offset_t1 - 1;
		  p1->EST_start= pe + offset_p;
		  p1->GEN_start= pg + offset_t2;
		  INFO("New 'local' factorization on EST:     [%8d-%8d]   --   [%8d-%8d]",
				 pnew->EST_start, pnew->EST_end, p1->EST_start, p1->EST_end);
		  INFO("New 'local' factorization on genomic: [%8d-%8d] %.2s..%.2s [%8d-%8d]",
				 pnew->GEN_start, pnew->GEN_end,
				 GEN_seq + pnew->GEN_end + 1, GEN_seq + p1->GEN_start - 2,
				 p1->GEN_start, p1->GEN_end);
		  list_add_before_iterator(pfactit, pfact, pnew);
		}
	 }
  }
}

static
void
search_small_exon(pfactor * pp1,
						pfactor * pp2,
						plistit pfactit,
						pfactorization pfact,
						pEST_info genomic,
						pEST factorized_est,
						pconfiguration config) {
  NOT_NULL(pp1);
  NOT_NULL(pp2);
  NOT_NULL(*pp1);
  NOT_NULL(*pp2);
  pfactor p1= *pp1;
  pfactor p2= *pp2;
// Skip external exons
  DEBUG("  ...analyzing pair of factors [%d-%d, %d-%d] -- [%d-%d, %d-%d]...",
		  p1->EST_start, p1->EST_end,
		  p1->GEN_start, p1->GEN_end,
		  p2->EST_start, p2->EST_end,
		  p2->GEN_start, p2->GEN_end);
  my_assert(p1->EST_start <= p1->EST_end);
  my_assert(p1->GEN_start <= p1->GEN_end);
  my_assert(p2->EST_start <= p2->EST_end);
  my_assert(p2->GEN_start <= p2->GEN_end);
  const size_t e1len= p1->EST_end + 1 - p1->EST_start;
  const size_t g1len= p1->GEN_end + 1 - p1->GEN_start;
  const size_t e2len= p2->EST_end + 1 - p2->EST_start;
  const size_t g2len= p2->GEN_end + 1 - p2->GEN_start;
// Check if it is possible to decompose the three exons into one small exon and two normal exons
  if ((e1len + e2len) >= (_LB_SMALL_EXON_LENGTH_ + 2*_UB_SMALL_EXON_LENGTH_)) {
	 size_t e1slen= MIN(MIN(e1len, g1len), _UB_SMALL_EXON_LENGTH_);
	 size_t e1sstart= p1->EST_end + 1 - e1slen;
	 char* e1sfact= c_palloc(e1slen+1);
	 char* const orige1sfact= e1sfact;
	 strncpy(e1sfact, factorized_est->info->EST_seq + e1sstart, e1slen);
	 e1sfact[e1slen]= '\0';
	 size_t g1slen= MIN(MIN(e1len, g1len), _UB_SMALL_EXON_LENGTH_);
	 size_t g1sstart= p1->GEN_end + 1 - g1slen;
	 char* g1sfact= c_palloc(g1slen+1);
	 char* const origg1sfact= g1sfact;
	 strncpy(g1sfact, genomic->EST_seq + g1sstart, g1slen);
	 g1sfact[g1slen]= '\0';
	 TRACE("EST suffix: %s", e1sfact);
	 TRACE("GEN suffix: %s", g1sfact);

	 size_t e2plen= MIN(MIN(e2len, g2len), _UB_SMALL_EXON_LENGTH_);
	 char* e2pfact= c_palloc(e2plen+1);
	 strncpy(e2pfact, factorized_est->info->EST_seq + p2->EST_start, e2plen);
	 e2pfact[e2plen]= '\0';
	 size_t g2plen= MIN(MIN(e2len, g2len), _UB_SMALL_EXON_LENGTH_);
	 char* g2pfact= c_palloc(g2plen+1);
	 strncpy(g2pfact, genomic->EST_seq + p2->GEN_start, g2plen);
	 g2pfact[g2plen]= '\0';
	 TRACE("EST prefix: %s", e2pfact);
	 TRACE("GEN prefix: %s", g2pfact);

	 size_t e1scut, g1scut;
	 const size_t sed= compute_best_prefix_cut(e1sfact, e1slen, g1sfact, g1slen,
															 &e1scut, &g1scut);
	 size_t e2pcut, g2pcut;
	 const size_t ped= compute_best_suffix_cut(e2pfact, e2plen, g2pfact, g2plen,
															 &e2pcut, &g2pcut);
	 DEBUG("        original edit distance of the suffix of the first exon:  %zu (cuts= %zu, %zu)",
			 sed, e1scut, g1scut);
	 DEBUG("        original edit distance of the prefix of the second exon: %zu (cuts= %zu, %zu)",
			 ped, e2pcut, g2pcut);

	 const size_t prev_ed= (sed + ped + e1scut + g1scut + e2plen + g2plen - e2pcut - g2pcut);
	 bool perform_search= false;
	 if (prev_ed > (_MAX_ERROR_RATE_*_UB_SMALL_EXON_LENGTH_)) {
		DEBUG("Border alignment is not good. There could be a small exon.");
		perform_search= true;
	 }
	 if ((prev_ed > 0) && !is_canonical_intron(genomic->EST_seq,
															 p1->GEN_end+1, p2->GEN_start-1)) {
		DEBUG("Border alignment is quite good but the intron is not canonical. "
				"There could be a small exon.");
		perform_search= true;
	 }
	 if (perform_search) {
// Make the cut
		e1sfact += e1scut;
		e1sstart += e1scut;
		e1slen = e1slen - e1scut;
		g1sfact += g1scut;
		g1sstart += g1scut;
		g1slen = g1slen - g1scut;
		e1sfact[e1slen]= '\0';
		g1sfact[g1slen]= '\0';
		TRACE("EST new suffix: %s", e1sfact);
		TRACE("GEN new suffix: %s", g1sfact);
		e2plen= e2pcut;
		e2pfact[e2plen]= '\0';
		g2plen= g2pcut;
		g2pfact[g2plen]= '\0';
		TRACE("EST new prefix: %s", e2pfact);
		TRACE("GEN new prefix: %s", g2pfact);
// Search the longest common factor of e1sfact+e2pfact and allgfact
		const size_t elen= e1slen+e2plen;
		char* const efact= c_palloc(elen+1);
		strncpy(efact, e1sfact, e1slen);
		strncpy(efact + e1slen, e2pfact, e2plen);
		efact[elen]= '\0';
		const size_t allglen= p2->GEN_start + g2plen - g1sstart;
		char* allgfact= c_palloc(allglen+1);
		strncpy(allgfact, genomic->EST_seq + g1sstart, allglen);
		allgfact[allglen]= '\0';
		size_t pg, pe, cflen;
		find_longest_common_factor(allgfact+g1slen, allglen-g1slen-g2plen, efact, elen, e1slen,
											&pg, &pe, &cflen);
		if (cflen >= _LB_SMALL_EXON_LENGTH_) {
		  DEBUG("Found a possible small exon of length (at most) %zu.", cflen);
		  TRACE("First border refinement.");
		  TRACE("Factor of the EST sequence: %.*s", (int)(pe+cflen), efact);
		  TRACE("Factor of the genomic seq:  %.*s", (int)(pg+cflen+g1slen), allgfact);
		  size_t offset_p_1, offset_t1_1, offset_t2_1;
		  unsigned int new_ed_1= allglen;
		  const bool ref_res_1= refine_borders(efact, pe+cflen,
															allgfact, pg+cflen+g1slen,
															sed+ped,
															&offset_p_1, &offset_t1_1, &offset_t2_1,
															&new_ed_1);
		  TRACE("1st border refinement gave: "
				  "Success? %s  "
				  "offset_p: %zu  offset_t1: %zu  offset_t2: %zu  "
				  "new_edit_distance: %u",
				  ref_res_1?"OK":"NO", offset_p_1, offset_t1_1, offset_t2_1, new_ed_1);
		  TRACE("Resulting alignment:");
		  TRACE("Factor of the EST sequence: %.*s%.*s%s", (int)(offset_p_1), efact, (int)(offset_t2_1-offset_t1_1), DOT_STRING, efact+offset_p_1);
		  TRACE("Factor of the genomic seq:  %.*s", (int)(pg+cflen+g1slen), allgfact);
		  if (!ref_res_1) {
			 DEBUG("The border alignments of the first intron is not good enough. Small exon discarded.");
		  } else if ((int)offset_t2_1 - (int)offset_t1_1 < config->min_intron_length) {
			 DEBUG("The first intron induced by the small exon is too short (%zunt). "
					 "Small exon discarded.", offset_t2_1 - offset_t1_1);
		  } else {
			 TRACE("Second border refinement.");
			 TRACE("Factor of the EST sequence: %.*s", (int)(elen-offset_p_1), efact+offset_p_1);
			 TRACE("Factor of the genomic seq:  %.*s", (int)(allglen-offset_t2_1), allgfact+offset_t2_1);
			 size_t offset_p_2, offset_t1_2, offset_t2_2;
			 unsigned int new_ed_2= allglen;
			 const bool ref_res_2= refine_borders(efact+offset_p_1, elen-offset_p_1,
															  allgfact+offset_t2_1, allglen-offset_t2_1,
															  sed+ped-new_ed_1,
															  &offset_p_2, &offset_t1_2, &offset_t2_2,
															  &new_ed_2);
			 TRACE("2nd border refinement gave: "
					 "Success? %s  "
					 "offset_p: %zu  offset_t1: %zu  offset_t2: %zu  "
					 "new_edit_distance: %u",
					 ref_res_2?"OK":"NO", offset_p_2, offset_t1_2, offset_t2_2, new_ed_2);
			 TRACE("Resulting alignment:");
			 TRACE("Factor of the EST sequence: %.*s%.*s%s", (int)(offset_p_2), efact+offset_p_1, (int)(offset_t2_2-offset_t1_2), DOT_STRING, efact+offset_p_1+offset_p_2);
			 TRACE("Factor of the genomic seq:  %.*s", (int)(allglen-offset_t2_1), allgfact+offset_t2_1);
			 if (!ref_res_2) {
				DEBUG("The border alignments of the second intron is not good enough. Small exon discarded.");
			 } else if ((int)offset_t2_2 - (int)offset_t1_2 < config->min_intron_length) {
				DEBUG("The second intron induced by the small exon is too short (%zunt). "
						"Small exon discarded.", offset_t2_2 - offset_t1_2);
			 } else if (offset_p_2 < _LB_SMALL_EXON_LENGTH_) {
				DEBUG("The small exon is too short (%zunt). "
						"Small exon discarded.", offset_p_2);
			 } else if (offset_p_2 > _UB_SMALL_EXON_LENGTH_) {
				DEBUG("The small exon is too long (%zunt). "
						"Small exon discarded.", offset_p_2);
			 } else {
// Compare the "Burset scores"
				const double prev_burset_freq=
				  getBursetFrequency_adaptor(genomic->EST_seq,
													  p1->GEN_end + 1,
													  p2->GEN_start);
				const double new_avg_burset_freq=
				  (getBursetFrequency_adaptor(allgfact, offset_t1_1, offset_t2_1) +
					getBursetFrequency_adaptor(allgfact+offset_t2_1, offset_t1_2, offset_t2_2)) / 2.0;
				DEBUG("Previous 'Burset score': %.1f", prev_burset_freq);
				DEBUG("New average 'Burset score': %.1f", new_avg_burset_freq);
				if (new_avg_burset_freq >= prev_burset_freq) {
				  INFO("A new small exon (%zunt) has been detected. "
						 "The new factorization does not decrease intron score, thus the small exon will be added "
						 "(old score: %.1f, new score %.1f).", offset_p_2, prev_burset_freq, new_avg_burset_freq);
				  INFO("Previous 'local' factorization on EST:     (%8d -%8d)   ...   ...   ...    (%8d -%8d).",
						 p1->EST_start, p1->EST_end,
						 p2->EST_start, p2->EST_end);
				  INFO("Previous 'local' factorization on genomic: (%8d -%8d)   ...   ...   ...    (%8d -%8d).",
						 p1->GEN_start, p1->GEN_end,
						 p2->GEN_start, p2->GEN_end);
				  pfactor pnew= factor_create();
				  pnew->EST_start= p1->EST_end + 1 - e1slen + offset_p_1;
				  pnew->EST_end  = p1->EST_end + 1 - e1slen + offset_p_1 + offset_p_2 - 1;
				  pnew->GEN_start= p1->GEN_end + 1 - g1slen + offset_t2_1;
				  pnew->GEN_end  = p1->GEN_end + 1 - g1slen + offset_t2_1 + offset_t1_2 - 1;
				  p2->EST_start= p1->EST_end + 1 - e1slen + offset_p_1 + offset_p_2;
				  p2->GEN_start= p1->GEN_end + 1 - g1slen + offset_t2_1 + offset_t2_2;
				  p1->EST_end  = p1->EST_end + 1 - e1slen + offset_p_1 - 1;
				  p1->GEN_end  = p1->GEN_end + 1 - g1slen + offset_t1_1 - 1;
				  INFO("Modified 'local' factorization on EST:     (%8d -%8d) (%8d -%8d) (%8d -%8d).",
						 p1->EST_start, p1->EST_end,
						 pnew->EST_start, pnew->EST_end,
						 p2->EST_start, p2->EST_end);
				  INFO("Modified 'local' factorization on genomic: (%8d -%8d) (%8d -%8d) (%8d -%8d).",
						 p1->GEN_start, p1->GEN_end,
						 pnew->GEN_start, pnew->GEN_end,
						 p2->GEN_start, p2->GEN_end);
				  list_add_before_iterator(pfactit, pfact, pnew);
				  *pp2= pnew;
				} else {
				  DEBUG("The new placement decreases intron score, thus the small exon will NOT be added.");
				}
			 }
		  }
		}
		pfree(efact);
		pfree(allgfact);
	 }

	 pfree(orige1sfact);
	 pfree(e2pfact);
	 pfree(origg1sfact);
	 pfree(g2pfact);
  }
}

static
void
search_for_new_small_exons(pEST_info genomic,
									pEST factorized_est,
									pconfiguration config) {
  DEBUG("Searching for possible small exons...");
  plistit pl_f_it= list_first(factorized_est->factorizations);
  while (listit_has_next(pl_f_it)) {
	 DEBUG("Analyzing a factorization...");
	 pfactorization pfact= listit_next(pl_f_it);

	 pfactor p1= NULL, p2= NULL;
	 plistit pl_factor_it= list_first(pfact);
	 plistit pl_factor_it_prev= list_first(pfact);
	 if (listit_has_next(pl_factor_it)) {
		p1= listit_next(pl_factor_it);
		if (p1->EST_start > _LB_SMALL_EXON_LENGTH_) {
		  DEBUG("A prefix of %dbp is not aligned. Check if a small exon is present...",
				  p1->EST_start);
		  search_small_exon_at_prefix(&p1, pl_factor_it, pfact,
												genomic, factorized_est, config);
		}
	 }
	 if (listit_has_next(pl_factor_it)) {
		p2= listit_next(pl_factor_it);
		listit_next(pl_factor_it_prev);
	 }
	 while (p2 != NULL) {
		search_small_exon(&p1, &p2, pl_factor_it, pfact,
								genomic, factorized_est, config);
		p1= p2;
		p2= NULL;
		if (listit_has_next(pl_factor_it)) {
		  p2= listit_next(pl_factor_it);
		  listit_next(pl_factor_it_prev);
		}
	 }
	 listit_destroy(pl_factor_it);
	 listit_destroy(pl_factor_it_prev);
  }
  listit_destroy(pl_f_it);
}



// ------------------------------------------------------------------------
//
// Removing false small exons
// Issue #3: Unnecessary exons
//

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
		(elen <= _UB_MED_EXON_LENGTH_)) {
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
  search_for_new_small_exons(genomic, factorized_est, config);
}

