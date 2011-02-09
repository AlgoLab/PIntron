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
** aug_suffix_tree.c
*/

#include "aug_suffix_tree.h"

#include "util.h"

#include "log.h"

ppreproc_gen
PGen_create(void) {
  ppreproc_gen pg= PALLOC(struct _preproc_gen);
  pg->gen= NULL;
  pg->alph= NULL;
  pg->alph_occ= NULL;
  pg->keys= NULL;
  pg->alph_size= 0;
  pg->gen_len= 0;
  return pg;
}

void
PGen_destroy(ppreproc_gen pg) {
  if (pg == NULL)
	 return;
  if (pg->gen != NULL)
	 pfree(pg->gen);
  if (pg->alph != NULL)
	 pfree(pg->alph);
  if (pg->alph_occ != NULL)
	 pfree(pg->alph_occ);
  if (pg->keys != NULL)
	 pfree(pg->keys);
  pfree(pg);
}

void
preprocess_text(pEST_info gen, ppreproc_gen const pg) {
  NOT_NULL(gen);
  NOT_NULL(pg);
  const char* seq= gen->EST_seq;
  size_t i= 0;
  size_t n= 0;
  size_t* occ= NPALLOC(size_t, 256);
  for (i= 0; i<256; ++i)
	 occ[i]= 0;
  i= 0;
  while (seq[i] != '\0') {
	 if (occ[(unsigned int)seq[i]]==0)
		++n;
	 ++occ[(unsigned int)seq[i]];
	 ++i;
  }
  pg->gen= gen;
  pg->gen_len= i;
  char * alph= c_palloc(n+1);
  alph[n]= '\0';
  size_t* alph_occ= NPALLOC(size_t, n);
  size_t j= 0;
  for (i= 0; i<256; ++i) {
	 if (occ[i]>0) {
		alph[j]= (char)i;
		alph_occ[j]= occ[i]+1; // +1 is due to the fact that the first
									  // suffix is added for every previous character
		++j;
	 }
  }
  j= 0;
  for (i= 0; i<256; ++i) {
	 if (occ[i]>0) {
		occ[i]= j;
		++j;
	 } else {
		occ[i]= n;
	 }
  }
  pg->alph= alph;
  pg->alph_occ= alph_occ;
  pg->alph_size= n;
  pg->keys= occ;
#ifdef LOG_DEBUG_ENABLED
  DEBUG("Alphabet: |%s|", alph);
  size_t cumul= 0;
  for (i= 0; i<n; ++i) {
	 DEBUG("Occurrences of %c: %10zu  (%12zu)", alph[i], alph_occ[i], cumul);
	 cumul += alph_occ[i];
  }
#endif
}

static void
initialize_vector_list(LST_Node *node,
							  const size_t string_depth,
							  const size_t alph_size,
							  const size_t min_string_depth){
  node->string_depth= string_depth;
  node->single_char= '\0';
  if (string_depth>=min_string_depth)
	 node->slices= calloc(2*alph_size, sizeof(size_t));
  else
	 node->slices= NULL;
}

static void initialization(LST_Node *node,
									const size_t string_depth,
									const size_t alph_size,
									const size_t min_string_depth){
  LST_Edge *edge;
  initialize_vector_list(node, string_depth, alph_size, min_string_depth);
  for (edge= node->kids.lh_first;
		 edge!=NULL;
		 edge= edge->siblings.le_next) {
	 initialization(edge->dst_node,
						 string_depth+lst_edge_get_length(edge),
						 alph_size, min_string_depth);
  }
}


static void
fill_node_info(const char * const gen,
					const size_t alph_size,
					const size_t min_string_depth,
					size_t * const occ,
					const ppreproc_gen const pg,
					size_t * const next_el,
					LST_Node *node){
  if (node->kids.lh_first==NULL) {
	 node->single_char= '\0';
	 if (node->index>0) {
		const char prev_s= gen[node->index-1];
		node->single_char= prev_s;
		if (node->string_depth>=min_string_depth) {
		  const size_t k= get_key(pg, prev_s);
		  const size_t nek= next_el[k];
		  occ[nek]= node->index;
		  node->slices[2*k]= nek;
		  node->slices[2*k+1]= nek+1;
		  ++next_el[k];
		  TRACE("Found a leaf with index %d and previous "
				"character %c (key %zu) added at array position %zu.",
				  node->index, prev_s, k, nek);
		}
	 } else {
		for (size_t k= 0; k<pg->alph_size; ++k) {
		  const size_t nek= next_el[k];
		  occ[nek]= 0;
		  node->slices[2*k]= nek;
		  node->slices[2*k+1]= nek+1;
		  ++next_el[k];
		  TRACE("Found the leaf with index 0 added "
				  "at array position %zu.", nek);
		}
	 }
  } else {
	 LST_Edge *edge;
	 bool single= true;
	 char last= '\0';
	 for (edge= node->kids.lh_first; edge!=NULL;
			edge= edge->siblings.le_next) {
		LST_Node * const dnode= edge->dst_node;
		fill_node_info(gen, alph_size, min_string_depth, occ, pg, next_el, dnode);
		if (single) {
		  if (dnode->single_char == '\0') {
			 last= '\0';
			 single= false;
		  } else {
			 if (last == '\0')
				last= dnode->single_char;
			 else if (last != dnode->single_char) {
				single= false;
				last= '\0';
			 }
		  }
		}
		if (node->string_depth>=min_string_depth) {
		  for (size_t k= 0; k<alph_size; ++k) {
			 if (node->slices[2*k+1] == 0) {
				node->slices[2*k]= dnode->slices[2*k];
				node->slices[2*k+1]= dnode->slices[2*k+1];
			 } else if (node->slices[2*k+1] < dnode->slices[2*k+1]) {
				node->slices[2*k+1]= dnode->slices[2*k+1];
			 }
		  }
		}
	 }
	 node->single_char= last;
  }
}

void
stree_preprocess(LST_STree* stree, const ppreproc_gen const pg, const pconfiguration const conf) {
  NOT_NULL(stree);
  NOT_NULL(pg);
  my_assert(stree->num_strings == 1);
  INFO("Starting the suffix-tree preprocessing phase"
		 " with %d strings", stree->num_strings);
  stree->arr_occs= NPALLOC(size_t, pg->gen_len+pg->alph_size);
  size_t* next_el= NPALLOC(size_t, pg->alph_size);
  next_el[0]= 0;
  for (size_t i= 1; i<pg->alph_size; ++i)
	 next_el[i]= next_el[i-1] + pg->alph_occ[i-1];
  initialization(stree->root_node, 0, pg->alph_size, conf->min_factor_len);
  fill_node_info(pg->gen->EST_seq, pg->alph_size, conf->min_factor_len, stree->arr_occs,
					  pg, next_el, stree->root_node);
  pfree(next_el);
  INFO("Preprocessing complete!");
}

void
stree_info_destroy(LST_STree* stree) {
  INFO("Destroying the suffix-tree additional informations");
  pfree(stree->arr_occs);
  INFO("Destroying complete!");
}

static void
stree_print_node(LST_Node* node, FILE* f, int level) {
  print_repetitions(f, ' ', 3*level);
  fprintf(f, "<%p", (void*)node);
  if (node->suffix_link_node!=NULL){
	 fprintf(f, " suffix-link->%p", (void*)node->suffix_link_node);
  }
  if (node->up_edge!=NULL) {
	 fprintf(f, " label=\"%s\"",
				string_print_func_with_dollars(&node->up_edge->range));
  }
  if (node->kids.lh_first!=NULL) {
	 fprintf(f, ">\n");
	 for (LST_Edge* edge= node->kids.lh_first; edge!=NULL;
			edge= edge->siblings.le_next) {
		stree_print_node(edge->dst_node, f, level+1);
	 }
	 print_repetitions(f, ' ', 3*level);
	 fprintf(f, "</%p", (void*)node);
  } else {
	 fprintf(f, " occurrence=%d /", node->index);
  }
  fprintf(f, ">\n");
}


void
stree_print(LST_STree *stree, FILE* f) {
  stree_print_node(stree->root_node, f, 0);
}


char*
string_print_func_with_dollars(LST_StringIndex *index)
{
  static char s[4096];

  if (index->start_index == index->string->num_items - 1) {
	 sprintf(s, "$%d", index->string->id);
  } else {
	 TRACE("Stampa dell'intervallo %d-%d della stringa %d",
			 index->start_index, *(index->end_index),
			 index->string->id);
	 int tot= 1;
	 if (*(index->end_index)==index->string->num_items-1)
		tot= 0;
	 memcpy(s,
			  ((char *) index->string->data) + index->start_index,
			  *(index->end_index) - index->start_index + tot);
	 s[*(index->end_index) - index->start_index +tot] = '\0';
	 if (*(index->end_index)==index->string->num_items-1) {
		char myend[100];
		sprintf(myend, "$%d", index->string->id);
		strcat(s, myend);
	 }
  }
  return s;
}
