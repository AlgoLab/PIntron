/**
 *
 *
 *                              PIntron
 *
 * A novel pipeline for computational gene-structure prediction based on
 * spliced alignment of expressed sequences (ESTs and mRNAs).
 *
 * Copyright (C) 2010  Yuri Pirola
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
** max-emb-graph.c
*/

#include "max-emb-graph.h"

#include "aug_suffix_tree.h"
#include "list.h"
#include "util.h"

#include "log.h"

inline static char
first_char_on_edge(const LST_Edge* const edge) {
  if (edge->range.start_index == edge->range.string->num_items-1)
	 return '\0';
  return ((char*)edge->range.string->data)[edge->range.start_index];
}

inline static size_t
get_lcp(const char * const s1, const size_t l1, const char * const s2) {
  size_t i= 0;
  while ((i<l1) && (s2[i]!='\0') && (s1[i]==s2[i]))
	 ++i;
  return i;
}


static void
find_deepest_common_node_rec(const char const* pattern,
									  LST_Node* const node,
									  const size_t already_matched,
									  const char avoid_prev_char,
									  LST_Edge** final, size_t* const matched_len) {
  NOT_NULL(node);
  if (pattern[0]=='\0') {
	 *final= node->up_edge;
	 *matched_len= lst_edge_get_length(node->up_edge);
  } else {
	 TRACE("The next character on the pattern is >%c<", pattern[0]);
	 LST_Edge *edge;
	 for (edge= node->kids.lh_first;
			(edge!=NULL);
			edge= edge->siblings.le_next) {
		TRACE("The first character of the edge is >%c<",
				first_char_on_edge(edge));
		if (first_char_on_edge(edge)==pattern[0]) {
		  if ((edge->dst_node->single_char == '\0') ||
				(edge->dst_node->single_char != avoid_prev_char)) {
			 TRACE("Edge found!");
		  } else {
			 TRACE("Edge ignored because the subtree contains "
					 "only occurrences preceeded by char %c that is the same "
					 "of %c.", edge->dst_node->single_char, avoid_prev_char);
			 edge= NULL;
		  }
		  break;
		}
	 }
	 if (edge==NULL) {
		*final= node->up_edge;
		*matched_len= lst_edge_get_length(node->up_edge);
	 } else {
		const size_t edge_length= lst_edge_get_length(edge);
		size_t tmp_lcp= 0;
		if (edge_length == 1) {
		  tmp_lcp= 1;
		} else if (already_matched > 0) {
		  if (already_matched >= edge_length) {
			 tmp_lcp= edge_length;
		  } else {
			 tmp_lcp= already_matched +
				get_lcp(((char*)edge->range.string->data)+edge->range.start_index+already_matched,
						  edge_length-already_matched,
						  pattern+already_matched);
		  }
		} else {
		  tmp_lcp= get_lcp(((char*)edge->range.string->data)+edge->range.start_index,
								 edge_length,
								 pattern);
		}
		const size_t lcp= tmp_lcp;
		TRACE("Edge length %zu-%zu", lcp, edge_length);
		if ((pattern[lcp] == '\0') || (lcp < edge_length)) {
		  *final= edge;
		  *matched_len= lcp;
		} else {
		  my_assert(lcp==edge_length);
		  size_t new_already_matched= 0;
		  if (already_matched > lcp)
			 new_already_matched= already_matched-lcp;
		  find_deepest_common_node_rec(pattern+edge_length,
												 edge->dst_node,
												 new_already_matched,
												 avoid_prev_char,
												 final,
												 matched_len);
		}
	 }
  }
}


static void
find_deepest_common_node(const char* const pattern,
								 LST_STree* const tree,
								 const char avoid_prev_char,
								 LST_Edge** final, size_t* const matched_len) {
  TRACE("Finding deepest common node.");
  find_deepest_common_node_rec(pattern, tree->root_node, 0, avoid_prev_char, final, matched_len);
}

static void
follow_suffix_link_and_fast_fwd(const char* pattern,
										  const LST_Edge* const prev_edge,
										  size_t matched_len,
										  const char avoid_prev_char,
										  LST_Edge** final,
										  size_t* const out_matched_len) {
  TRACE("Following the suffix link.");
  LST_Node* sl= NULL;
  if (lst_edge_get_length(prev_edge)==matched_len) {
	 sl= prev_edge->dst_node->suffix_link_node;
	 matched_len= 0;
  } else {
	 sl= prev_edge->src_node->suffix_link_node;
  }
  NOT_NULL(sl);
  TRACE("Matched len %zu and string-depth %zu", matched_len, sl->string_depth);
  find_deepest_common_node_rec(pattern + sl->string_depth,
										 sl, matched_len,
										 avoid_prev_char,
										 final, out_matched_len);
}



static void
fill_list_pairings(LST_Edge* const edge,
						 LST_Edge* const block_edge,
						 const size_t* const occ,
						 const size_t alph_size,
						 const size_t const symbol_k,
						 plist Vi,
						 const int p,
						 const size_t l) {
  size_t block_start, block_end;
  TRACE("Previous symbol key %zu", symbol_k);
  for (size_t k= 0; k<alph_size; ++k) {
	 if (k == symbol_k)
		continue;
	 TRACE("Considering symbol %zu.", k);
	 const size_t start= edge->dst_node->slices[2*k];
	 const size_t end= edge->dst_node->slices[2*k+1];
	 if ((block_edge == NULL) || (block_edge->dst_node->slices[2*k+1]==0)) {
		block_start= block_end= end;
	 } else {
		block_start= block_edge->dst_node->slices[2*k];
		block_end= block_edge->dst_node->slices[2*k+1];
	 }
	 TRACE("Considering slice %zu-%zu with block slice %zu-%zu.",
			 start, end, block_start, block_end);
	 for (size_t i= start;
			i<block_start; ++i) {
		const int t= occ[i];
		if ((t>0) || ((t==0) && ((k==0)||((k==1)&&(symbol_k==0))))) {
		  ppairing pairing= pairing_create();
		  pairing->p= p;
		  pairing->t= t;
		  pairing->l= l;
		  list_add_to_tail(Vi, pairing);
		  DEBUG("Found the new pairing (%d, %d, %zd).", p, t, l);
		}
	 }
	 for (size_t i= block_end;
			i<end; ++i) {
		int t= occ[i];
		ppairing pairing= pairing_create();
		pairing->p= p;
		pairing->t= t;
		pairing->l= l;
		list_add_to_tail(Vi, pairing);
		DEBUG("Found the new pairing (%d, %d, %zd).", p, t, l);
	 }
  }
}

pext_array
build_vertex_set(pEST_info pattern,
					  LST_STree* tree,
					  const ppreproc_gen const pg,
					  pconfiguration config) {
  INFO("Starting the build of vertex set of the MEG.");
  my_assert(pattern!=NULL);
  my_assert(pattern->EST_seq!=NULL);
  my_assert(pattern->EST_id!=NULL);
  my_assert(tree!=NULL);
  my_assert(config!=NULL);
  const size_t pattern_len= strlen(pattern->EST_seq);
  DEBUG("The pattern is %zd characters long.", pattern_len);
  pext_array V= EA_create();

// Creation of the source pairing
  plist Vi= list_create();
  plistit Vij= NULL;
  plistit Vii= NULL;
  ppairing pairing= pairing_create();
  pairing->p= SOURCE_PAIRING_START;
  pairing->t= SOURCE_PAIRING_START;
  pairing->l= SOURCE_PAIRING_LEN;
  list_add_to_tail(Vi, pairing);
  EA_insert(V, Vi);

  LST_Edge* prev_N= NULL;
  size_t prev_matched_len= 0;
  char prev_symbol= '\0';
  size_t prev_symbol_key= pg->alph_size;
  for (unsigned int i= 0; i<pattern_len; ++i) {
	 TRACE("Considering the %dth suffix of the pattern.", i);
	 LST_Edge* block_edge= NULL;
	 Vi= list_create();
	 EA_insert(V, Vi);
	 LST_Edge* N= NULL;
	 size_t matched_len;
	 if (prev_N==NULL || prev_N->src_node->suffix_link_node==NULL) {
		find_deepest_common_node(pattern->EST_seq+i, tree, prev_symbol, &N, &matched_len);
	 } else {
		follow_suffix_link_and_fast_fwd(pattern->EST_seq+i,
												  prev_N,
												  prev_matched_len,
												  prev_symbol,
												  &N, &matched_len);
	 }
	 if (N == NULL) {
		DEBUG("The suffix cannot be matched.");
		prev_N= NULL;
		prev_matched_len= 0;
	 } else {
		NOT_NULL(N);
		NOT_NULL(N->src_node);
		NOT_NULL(N->dst_node);
		TRACE("The deepest common node has string-depth %zd.",
				N->src_node->string_depth + matched_len);
#ifdef MUMMER_EMULATION
		size_t min_string_depth= config->min_factor_len;
#else
		size_t min_string_depth= MAX((N->src_node->string_depth + matched_len)*(config->min_string_depth_rate),
											  config->min_factor_len);
		TRACE("The minimum string-depth that will be considered is %zd.", min_string_depth);
#endif

// Salvo il nodo trovato per seguire poi il suffix link
		prev_N= N;
		prev_matched_len= matched_len;

		while (N->src_node->string_depth + matched_len
				 >= min_string_depth) {
		  TRACE("Analysing the common node at string-depth %zd.",
				  N->src_node->string_depth + matched_len);
		  fill_list_pairings(N,
									block_edge,
									tree->arr_occs,
									pg->alph_size,
									prev_symbol_key,
									Vi,
									i,
									N->src_node->string_depth + matched_len);

		  NOT_NULL(N->src_node->up_edge);
		  NOT_NULL(N->src_node->up_edge->src_node);

		  block_edge= N;
		  N= N->src_node->up_edge;
		  matched_len= lst_edge_get_length(N);
		}
		list_sort(Vi, (comparator)pairing_compare);
#ifndef MUMMER_EMULATION
		plist ltoremove= list_create();
		list_last_reuse(Vi, &Vij);
		while (listit_has_prev(Vij)) {
		  ppairing PJ= listit_prev(Vij);
		  listit_copy_reuse(Vij, &Vii);
		  bool rim= false;
		  while (!rim && listit_has_prev(Vii)) {
			 ppairing PI= listit_prev(Vii);
			 if ((PJ->t > PI->t) && (PJ->t+PJ->l <= PI->t+PI->l)) {
				DEBUG("Pairings (%d, %d, %d) and (%d, %d, %d) seem "
						"to be low-complexity repetitions.",
						PAIRING(PI), PAIRING(PJ));
				DEBUG("--> removing (%d, %d, %d)", PAIRING(PJ));
				list_add_to_head(ltoremove, PJ);
				rim= true;
			 }
			 if ((PJ->t == PI->t + 1) && (PJ->l == PI->l)) {
				DEBUG("Pairings (%d, %d, %d) and (%d, %d, %d) seem "
						"to be low-complexity repetitions.",
						PAIRING(PI), PAIRING(PJ));
				DEBUG("--> removing (%d, %d, %d)", PAIRING(PJ));
				list_add_to_head(ltoremove, PJ);
				rim= true;
			 }
		  }
		}
		list_first_reuse(Vi, &Vij);
		while (!list_is_empty(ltoremove)) {
		  ppairing PI= list_remove_from_head(ltoremove);
		  while (PI != listit_next(Vij)) ;
		  list_remove_at_iterator(Vij, (delete_function)pairing_destroy);
		}
		list_destroy(ltoremove, (delete_function)noop_free);
#endif
	 }
	 prev_symbol= pattern->EST_seq[i];
	 prev_symbol_key= get_key(pg, prev_symbol);
  }
  listit_destroy(Vij);
  listit_destroy(Vii);
  Vi= list_create();
  pairing= pairing_create();
  pairing->p= SINK_PAIRING_START;
  pairing->t= SINK_PAIRING_START;
  pairing->l= SINK_PAIRING_LEN;
  list_add_to_tail(Vi, pairing);
  EA_insert(V, Vi);

#ifndef MUMMER_EMULATION
  INFO("Cleaning low-complexity pairings starting at different positions on P.");
  {
	 const size_t n= EA_size(V);
	 plist Vi1= EA_get(V, n-2);
	 for (size_t i= n-3; i>0; --i) {
		TRACE("Analysis of the %zdth pairing list.", i);
		plist Vi= EA_get(V, i);
		plistit Vit1= list_first(Vi1);
		while (listit_has_next(Vit1)) {
		  ppairing I1= listit_next(Vit1);
		  plistit Vit= list_first(Vi);
		  bool rim= false;
		  while (!rim && listit_has_next(Vit)) {
			 ppairing I= listit_next(Vit);
			 TRACE("Pairings (%d, %d, %d) and (%d, %d, %d)", PAIRING(I), PAIRING(I1));
			 if ((I->t == I1->t) && (I->l >= I1->l)) {
				DEBUG("Pairing (%d, %d, %d) rimosso", PAIRING(I1));
				list_remove_at_iterator(Vit1, (delete_function)pairing_destroy);
				rim= true;
			 }
		  }
		  listit_destroy(Vit);
		}
		listit_destroy(Vit1);
		Vi1= Vi;
	 }
  }
#endif

  INFO("Build of vertex set of the MEG completed!");

  return V;
}

int compute_fl(const pconfiguration config) {
  my_assert(config!=NULL);
  return 2*(config->min_factor_len)+1;
}

int compute_gl(const pconfiguration config) {
  my_assert(config!=NULL);
  return 2*(config->min_factor_len)+3;
}


static bool
is_there_an_edge_strict(const ppairing I,
								const ppairing J,
								const int l,
								const int fl,
								pconfiguration config) {
  NOT_NULL( I );
  NOT_NULL( J );
  NOT_NULL( config );

  const double MAX_OVERLAP= 0.4;

//XXX: parametro arbitrario
  const bool I_is_long= (I->l >= 5*l);
  const bool J_is_long= (J->l >= 5*l);

  if (J->p <= I->p) {
	 TRACE("STRICT. "
			 "Pairing (%d, %d, %d) starts before (%d, %d, %d) on P. "
			 "Edge not possible.",
			 PAIRING(J), PAIRING(I));
	 return false;
  }
  if (J->t <= I->t) {
	 TRACE("STRICT. "
			 "Pairing (%d, %d, %d) starts before (%d, %d, %d) on T. "
			 "Edge not possible.",
			 PAIRING(J), PAIRING(I));
	 return false;
  }

  if (I->p + I->l <= J->p && J->p <= I->p + I->l + fl) {
// Simple-sequence on P
	 if ((I->t + I->l <= J->t) &&
		  ((config->max_intron_length==0) ||
			(J->t <= I->t + I->l + config->max_intron_length))) {
// and Simple-sequence on T
		return true;
	 } else if ((I->t + 2*l <= J->t + J->l) &&
					(J->t < I->t + I->l) &&
					(J->p + I->t - I->p - J->t <= fl)) {
// and Overlap on T
		if ((I_is_long) &&
			 (I->t + I->l - J->t > MAX_OVERLAP*I->l)) {
		  DEBUG("STRICT. "
				  "Pairing (%d, %d, %d) cuts too much of (%d, %d, %d) on T. "
				  "Edge not possible.",
				  PAIRING(J), PAIRING(I));
		  return false;
		}
		return true;
	 }
  } else if ((I->p + 2*l <= J->p + J->l) &&
				 (J->p < I->p + I->l)) {
// Overlap on P
	 if ((I->t + I->l <= J->t) &&
		  ((config->max_intron_length==0) ||
			(J->t <= I->t + I->l + config->max_intron_length))) {
// and Simple-sequence on T
		return true;
	 } else if ((I->t + 2*l <= J->t + J->l) &&
					(J->t < I->t + I->l) &&
					(J->p + I->t - I->p - J->t <= fl)) {
// and Overlap on T
		return true;
	 }
  }
  TRACE("Pairings (%d, %d, %d) and (%d, %d, %d) are "
		  "not connected.", PAIRING(I), PAIRING(J));
  return false;
}

static bool
is_there_an_edge(const ppairing I,
					  const ppairing J,
					  const int l,
					  const int fl,
					  pconfiguration config) {
  my_assert(I!=NULL);
  my_assert(J!=NULL);
  my_assert(config != NULL);

  if (I==J)
	 return false;

  //Modifica per evitare loop...
  if(J->p-I->p < 0 && (J->t-I->t > 0 && J->t-I->t < I->l))
	  return false;
  if(J->p-I->p <= 0 && J->t-I->t <= 0){
	  if((J->p-I->p < 0 || J->t-I->t < 0) || J->l < I->l)
		  return false;
  }

  if (I->p + I->l <= J->p && J->p <= I->p + I->l + fl) {
// Simple-sequence on P
	 if ((I->t + I->l <= J->t) &&
		  ((config->max_intron_length==0) ||
			(J->t <= I->t + I->l + config->max_intron_length))) {
// and Simple-sequence on T
		TRACE("Pairings (%d, %d, %d) and (%d, %d, %d) in "
				"Simple-Sequence on P and Simple-Sequence on T case.",
				PAIRING(I), PAIRING(J));
		return true;
	 } else if ((I->t + 2*l <= J->t + J->l) &&
					(J->t < I->t + I->l) &&
					(J->p + I->t - I->p - J->t <= fl)) {
// and Overlap on T
		TRACE("Pairings (%d, %d, %d) and (%d, %d, %d) in "
				"Simple-Sequence on P and Overlap on T case.",
				PAIRING(I), PAIRING(J));
		return true;
	 }
  } else if ((I->p + 2*l <= J->p + J->l) &&
				 (J->p < I->p + I->l)) {
// Overlap on P
	 if ((I->t + I->l <= J->t) &&
		  ((config->max_intron_length==0) ||
			(J->t <= I->t + I->l + config->max_intron_length))) {
// and Simple-sequence on T
		TRACE("Pairings (%d, %d, %d) and (%d, %d, %d) in "
				"Overlap on P and Simple-Sequence on T case.",
				PAIRING(I), PAIRING(J));
		return true;
	 } else if ((I->t + 2*l <= J->t + J->l) &&
					(J->t < I->t + I->l) &&
					(J->p + I->t - I->p - J->t <= fl)) {
// and Overlap on T
		TRACE("Pairings (%d, %d, %d) and (%d, %d, %d) in "
				"Overlap on P and Overlap on T case.",
				PAIRING(I), PAIRING(J));
		return true;
	 }
  }
  TRACE("Pairings (%d, %d, %d) and (%d, %d, %d) are "
		  "not connected.", PAIRING(I), PAIRING(J));
  return false;
}


static void
add_edges_from(const ppairing I, const pext_array V,
					const int l, const int fl, const int n,
					pconfiguration config) {
  TRACE("Adding edges from (%d, %d, %d)", PAIRING(I));
  TRACE("l= %d, fl= %d, n= %d", l, fl, n);
  const int ubound= MIN(I->p+I->l+fl+1, n-l);
  plistit Vjt= NULL;
  for (int j= 0; j<ubound; ++j) {
	 plist Vj= EA_get(V, j);
	 list_first_reuse(Vj, &Vjt);
	 while (listit_has_next(Vjt)) {
		ppairing J= listit_next(Vjt);
		if (is_there_an_edge_strict(I, J, l, fl, config)) {
			list_add_to_tail(I->adjs, J);
			list_add_to_tail(J->incs, I);
		}
	 }
  }
  listit_destroy(Vjt);
}

static void
add_edges_from_source(pext_array V,
							 pconfiguration config) {
  INFO("Adding the edges from the source pairing.");
  int p_len= EA_size(V)-2;
  const int L= config->min_factor_len;

  int max_p= (int)(((double)p_len) * config->max_prefix_discarded_rate);

// Get the source pairing
  plistit lit0= list_first(EA_get(V, 0));
  my_assert(listit_has_next(lit0));
  ppairing source= listit_next(lit0);
  listit_destroy(lit0);
// Iterate over the vertex set
  int i= 1;
  while (i<=max_p) {
	 plist Vi= EA_get(V, i);
	 plistit Vit= list_first(Vi);
	 while (listit_has_next(Vit)) {
		ppairing I= listit_next(Vit);
// Look for a non-overlapped incident
// If it has a non-overlapped incident, then it is not a source
		plistit incsit= list_first(I->incs);
		bool possible_source= true;
		while (possible_source && listit_has_next(incsit)) {
		  ppairing inc= listit_next(incsit);
		  possible_source= ! (((inc->p + inc->l <= I->p) || (I->p + I->l <= inc->p)) &&
									 ((inc->t + inc->l <= I->t) || (I->t + I->l <= inc->t)));
		  possible_source= possible_source &&
			 (((inc->p + L) > I->p) || ((inc->t + L) > I->t));
		}
		listit_destroy(incsit);
		if (possible_source) {
		  DEBUG("Pairing (%d, %d, %d) can be a source.", PAIRING(I));
		  list_add_to_tail(source->adjs, I);
		  list_add_to_tail(I->incs, source);
		} else {
		  TRACE("Pairing (%d, %d, %d) cannot be a source.", PAIRING(I));
		}
	 }
	 listit_destroy(Vit);
	 ++i;
  }
}

static void
add_edges_to_sink(pext_array V,
						pconfiguration config) {
  INFO("Adding the edges towards the sink pairing.");
  const int L= config->min_factor_len;
  int p_len= EA_size(V)-2;

  int min_p= (int)(((double)p_len) * (1.0-config->max_suffix_discarded_rate));

// Get the sink pairing
  plistit litN= list_first(EA_get(V, p_len+1));
  my_assert(listit_has_next(litN));
  ppairing sink= listit_next(litN);
  listit_destroy(litN);
// Iterate over the vertex set
  int i= 1;
  while (i<=p_len) {
	 plist Vi= EA_get(V, i);
	 plistit Vit= list_first(Vi);
	 while (listit_has_next(Vit)) {
		ppairing I= listit_next(Vit);
		if (I->p + I->l < min_p) continue;
// Look for a non-overlapped adjacent
// If it has a non-overlapped adjacent, then it is not a sink
		plistit adjsit= list_first(I->adjs);
		bool possible_sink= true;
		while (possible_sink && listit_has_next(adjsit)) {
		  ppairing adj= listit_next(adjsit);
		  possible_sink= ! (((adj->p + adj->l <= I->p) || (I->p + I->l <= adj->p)) &&
								  ((adj->t + adj->l <= I->t) || (I->t + I->l <= adj->t)));
		  possible_sink= possible_sink &&
			 (((I->p + I->l + L) > (adj->p + adj->l)) ||
			  ((I->t + I->l + L) > (adj->t + adj->l)));
		}
		listit_destroy(adjsit);
		if (possible_sink) {
		  DEBUG("Pairing (%d, %d, %d) can be a sink.", PAIRING(I));
		  list_add_to_tail(sink->incs, I);
		  list_add_to_tail(I->adjs, sink);
		} else {
		  TRACE("Pairing (%d, %d, %d) cannot be a sink.", PAIRING(I));
		}
	 }
	 listit_destroy(Vit);
	 ++i;
  }
}


void
build_edge_set(pext_array V,
					pconfiguration config
					) {
  INFO("Starting the build of the edge set of the MEG.");
  my_assert(V!=NULL);
  my_assert(config!=NULL);
  const size_t n= EA_size(V);
  const int fl= compute_fl(config);
  for (size_t i= 1; i<n-1; ++i) {
	 TRACE("Analysis of the %zdth pairing list.", i);
	 plist Vi= EA_get(V, i);
	 plistit Vit= list_first(Vi);
	 while (listit_has_next(Vit)) {
		ppairing I= listit_next(Vit);
		TRACE("Analysis of pairing (%d, %d, %d)", PAIRING(I));
		add_edges_from(I, V, config->min_factor_len, fl, n, config);
	 }
	 listit_destroy(Vit);
  }
  add_edges_from_source(V, config);
  add_edges_to_sink(V, config);
  INFO("Build of edge set of the MEG completed!");
}

//Etichetta come intronici gli archi piu' lunghi di questa quantita'
#define INTRONIC_EDGE 50
void
add_intronic_edges_to_file(FILE* f, pext_array V) {
  my_assert(V!=NULL);
  my_assert(f!=NULL);
  for (size_t i= 0; i< EA_size(V); ++i) {
	 plist Vi= EA_get(V, i);
	 plistit lit= list_first(Vi);
	 while (listit_has_next(lit)) {
		ppairing p= listit_next(lit);
		if ((p->p != SOURCE_PAIRING_START) && (p->p != SINK_PAIRING_START))  {
		  plistit ait= list_first(p->adjs);
		  while (listit_has_next(ait)) {
			 ppairing a= listit_next(ait);
			 if (a->p != SINK_PAIRING_START) {
				fprintf(f, "%d %d %d %d %d %d %d %d %d",
						  p->t + p->l, a->t,
						  p->p + p->l, a->p,
						  (a->t - p->t - p->l),
						  (a->p - p->p - p->l),
						  (a->t - p->t) - (a->p - p->p),
						  p->l, a->l);
				if ((a->t - p->t) - (a->p - p->p) >= INTRONIC_EDGE)
				  fprintf(f, " intronic");
				fprintf(f, "\n");
			 }
		  }
		  listit_destroy(ait);
		}
	 }
	 listit_destroy(lit);
  }
}
#undef INTRONIC_EDGE


#ifdef LOG_GRAPHS
// Funzione per la creazione del sorgente dot di un MEG.
// Abilitata solo in caso di debug.

// Colora (in giallo) i pairing di lunghezza almeno
// pari a quella specificata.
#define MIN_UNDERLINE_LEN 30
// Colora (in rosso) gli archi che corrispondono a gap
// su P di lunghezza superiore a quella specificata.
#define MAX_GAP_ON_P 4

void
save_meg_to_filename(pext_array V, const char* const filename) {
  NOT_NULL(V);
  NOT_NULL(filename);
  FILE* f= fopen(filename, "w");
  if (f==NULL) {
	 ERROR("Impossible to save the MEG to file '%s'.", filename);
  } else {
	 DEBUG("Saving the MEG to file '%s'.", filename);
	 print_meg(V, f);
	 fclose(f);
  }
}

void
print_meg(pext_array V, FILE* f) {
  my_assert(V!=NULL);
  my_assert(f!=NULL);
  fprintf(f, "digraph MEG {\n");
  for (size_t i= 0; i< EA_size(V); ++i) {
	 plist Vi= EA_get(V, i);
	 plistit lit= list_first(Vi);
	 while (listit_has_next(lit)) {
		ppairing p= listit_next(lit);
		if (p->p == SOURCE_PAIRING_START) {
		  fprintf(f, "n%d [label=\"source\"", p->id);
		} else if (p->p == SINK_PAIRING_START)  {
		  fprintf(f, "n%d [label=\"sink\"", p->id);
		} else {
		  fprintf(f, "n%d [label=\"%d (%d-%d, %d-%d)\"", p->id, p->id, p->p, p->p+p->l, p->t, p->t+p->l);
		}
// Evidenzia i nodi "abbastanza" lunghi
		if (p->l >= MIN_UNDERLINE_LEN) {
		  fprintf(f, ", style=filled, fillcolor=yellow");
		}
		fprintf(f,"];\n");
		plistit ait= list_first(p->adjs);
		while (listit_has_next(ait)) {
		  ppairing a= listit_next(ait);
		  fprintf(f, "\tn%d -> n%d[fontsize=12", p->id, a->id);
		  if ((p->p != SOURCE_PAIRING_START) && (a->p != SINK_PAIRING_START)) {
			 fprintf(f, ",label=\"P:%d\\nT:%d\\nD:%d\"",
						(a->p - p->p - p->l),
						(a->t - p->t - p->l),
						(a->t - p->t) - (a->p - p->p));
			 if ((a->t - p->t) - (a->p - p->p) < MAX_GAP_ON_P)
				fprintf(f, ",color=red");
			 else
				fprintf(f, ",color=blue");
		  }
		  fprintf(f, "];\n");
		}
		listit_destroy(ait);
	 }
	 listit_destroy(lit);
  }
  fprintf(f, "}\n");
}

#undef MIN_UNDERLINE_LEN
#undef MAX_GAP_ON_P

#else // if not defined LOG_GRAPHS

void
save_meg_to_filename(pext_array V, const char* const filename) {
  (void)V;
  (void)filename;
}

void
print_meg(pext_array V, FILE* f) {
  (void)V;
  (void)f;
}




#endif
