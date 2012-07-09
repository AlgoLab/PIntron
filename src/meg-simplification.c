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
** meg-simplification.c
**
** Made by (Yuri Pirola)
** Login   <yuri@dottxxii-12.dottorato.disco.unimib.it>
**
** Started on  Wed Dec  9 13:11:50 2009 Yuri Pirola
** Last update Sun May 12 01:17:25 2002 Speed Blue
*/

#include "meg-simplification.h"

#include "types.h"
#include "util.h"
#include "int_list.h"
#include "bit_vector.h"

#include "max-emb-graph.h"

#include "log.h"


void MEG_stats(pext_array V, size_t* tot_pairings, size_t* tot_edges) {
  *tot_pairings= 0;
  *tot_edges= 0;
  for (size_t i= 0; i< EA_size(V); ++i) {
	 plist Vi= EA_get(V, i);
	 plistit lit= list_first(Vi);
	 while (listit_has_next(lit)) {
		ppairing p= listit_next(lit);
		++(*tot_pairings);
		(*tot_edges)+= list_size(p->adjs);
	 }
	 listit_destroy(lit);
  }
}


bool
is_too_complex_for_compaction(pext_array V, pconfiguration config) {
  NOT_NULL(V);
  (void)config;

  size_t tot_pairings= 0;
  size_t tot_edges= 0;
  MEG_stats(V, &tot_pairings, &tot_edges);
  INFO("The MEG has %7zd vertices and %7zd edges.", tot_pairings, tot_edges);
//XXX: INSERITO PARAMETRO ARBITRARIO SUL NUMERO MASSIMO DI PAIRING IN UN MEG E DI ARCHI
  if (
		(tot_edges > 1000) ||
		(tot_pairings > 2000)
		) {
	 DEBUG("The MEG is over populated.");
	 return true;
  } else {
	 return false;
  }
}

bool
is_too_complex(pext_array V, pconfiguration config) {
  int min_len= 0;
  size_t freq_min_len= 0;
  size_t tot_pairings= 0;
  size_t tot_edges= 0;
  const size_t est_len= EA_size(V)-2;
  for (size_t i= 0; i< EA_size(V); ++i) {
	 plist Vi= EA_get(V, i);
	 plistit lit= list_first(Vi);
	 while (listit_has_next(lit)) {
		ppairing p= listit_next(lit);
		++tot_pairings;
		if (min_len==0 || p->l < min_len) {
		  min_len= p->l;
		  freq_min_len= 1;
		} else if (p->l==min_len) {
		  ++freq_min_len;
		}
		tot_edges+= list_size(p->adjs);
	 }
	 listit_destroy(lit);
  }
  INFO("The MEG has %7zd vertices and %7zd edges.", tot_pairings, tot_edges);
  DEBUG("The shortest pairings have length %4d, "
		  "they appear %8zd times over %8zd total pairings.",
		  min_len, freq_min_len, tot_pairings);
  if ((config->max_pairings_in_MEG!=0) &&
		(tot_pairings > config->max_pairings_in_MEG) &&
		(freq_min_len >
		 config->max_freq_shortest_pairing *
		 tot_pairings)) {
	 DEBUG("The MEG vertex set is over populated.");
	 return true;
//XXX: INSERITO PARAMETRO ARBITRARIO SUL NUMERO MASSIMO DI PAIRING IN UN MEG
  } else if (
				 (tot_edges > 5*tot_pairings) ||
				 (tot_pairings > ((2*est_len)/config->min_factor_len)) ||
				 (
				  (tot_pairings > ((est_len)/config->min_factor_len)) &&
				  (tot_pairings >= 50)
				 )
				) {
	 DEBUG("The MEG vertex set is over populated.");
	 return true;
  } else {
	 return false;
  }
}

void
remove_other_sources_and_sinks(pext_array V) {
  bool removed;
  int V_size= EA_size(V);
  plistit Vit= NULL;
  plistit adjs= NULL;
  plistit incs= NULL;
  do {
	 removed= false;
	 for (int i= 1; i<V_size-1; ++i) {
		plist Vi= EA_get(V, i);
		list_first_reuse(Vi, &Vit);
		while (listit_has_next(Vit)) {
		  ppairing I= listit_next(Vit);
		  if (list_is_empty(I->adjs) || list_is_empty(I->incs)) {
// Non ho adiacenti o incidenti -> rimuovo I
			 removed= true;
#ifdef LOG_DEBUG_ENABLED
			 if (list_is_empty(I->adjs)) {
				DEBUG("Remove (%d, %d, %d) because its adjacents are empty.", PAIRING(I));
			 } else {
				DEBUG("Remove (%d, %d, %d) because its incidents are empty.", PAIRING(I));
			 }
#endif
// Rimuovo I come incidente dei suoi adiacenti
			 { list_first_reuse(I->adjs, &adjs);
				while (listit_has_next(adjs)) {
				  ppairing adj= listit_next(adjs);
				  list_remove_element(adj->incs, I, (delete_function)noop_free);
				}
			 }
// Rimuovo I come adiacente dei suoi incidenti
			 {
				list_first_reuse(I->incs, &incs);
				while (listit_has_next(incs)) {
				  ppairing inc= listit_next(incs);
				  list_remove_element(inc->adjs, I, (delete_function)noop_free);
				}
			 }
// Rimuovo I dalla lista Vi
			 list_remove_at_iterator(Vit, (delete_function)pairing_destroy);
		  }
		}
	 }
  } while (removed);
  listit_destroy(adjs);
  listit_destroy(incs);
  listit_destroy(Vit);
}


void
remove_useless_edges(pext_array V, pconfiguration config) {
  INFO("Removing useless CMEG edges.");
  plistit lit= NULL;
  plistit ait= NULL;
  const int g= compute_gl(config);
  for (size_t i= 1; i< EA_size(V); ++i) {
	 plist Vi= EA_get(V, i);
	 list_first_reuse(Vi, &lit);
	 while (listit_has_next(lit)) {
		ppairing p= listit_next(lit);
		list_first_reuse(p->adjs, &ait);
		while (listit_has_next(ait)) {
		  ppairing a= listit_next(ait);
// Computing gap on genomic
		  if (a->t != SINK_PAIRING_START) {
			 const int gap= MAX(a->t - a->p - p->t + p->p, 0);
			 TRACE("Gap between pairing (%4d, %4d, %4d) and (%4d, %4d, %4d) "
					 "is %8d.", p->p, p->t, p->l, a->p, a->t, a->l, gap);
			 if (gap> g && gap<config->min_intron_length) {
				DEBUG("The edge between pairing (%4d, %4d, %4d) and "
						"(%4d, %4d, %4d) is useless because of the gap %8d.",
						p->p, p->t, p->l, a->p, a->t, a->l, gap);
				list_remove_at_iterator(ait, (delete_function)noop_free);
#ifndef NDEBUG
				bool ris=
#endif
				  list_remove_element(a->incs, p, (delete_function)noop_free);
#ifndef NDEBUG
				my_assert(ris);
#endif
			 }
		  }
		}
	 }
  }
  listit_destroy(ait);
  listit_destroy(lit);
}

static void
copy_adjacencies(ppairing new_v, ppairing old_v) {
  plistit old_adjs= list_first(old_v->adjs);
  while (listit_has_next(old_adjs)) {
	 ppairing a= (ppairing)listit_next(old_adjs);
	 DEBUG("Connecting (%4d, %4d, %4d) to (%4d, %4d, %4d).", PAIRING(new_v), PAIRING(a));
	 list_add_to_tail(new_v->adjs, a);
	 list_add_to_tail(a->incs, new_v);
  }
  listit_destroy(old_adjs);
}

static void
copy_incidencies(ppairing new_v, ppairing old_v) {
  plistit old_incs= list_first(old_v->incs);
  while (listit_has_next(old_incs)) {
	 ppairing i= (ppairing)listit_next(old_incs);
	 DEBUG("Connecting (%4d, %4d, %4d) to (%4d, %4d, %4d).", PAIRING(i), PAIRING(new_v));
	 list_add_to_tail(new_v->incs, i);
	 list_add_to_tail(i->adjs, new_v);
  }
  listit_destroy(old_incs);
}


void
compact_short_edges(pext_array V, pconfiguration config) {
  INFO("Compacting short edges.");
  plistit lit= NULL;
  plistit ait= NULL;
  bool removed= false;
  do {
	 removed= false;
	 for (size_t i= 1; i< EA_size(V); ++i) {
		plist Vi= EA_get(V, i);
		list_first_reuse(Vi, &lit);
		while (listit_has_next(lit)) {
		  ppairing p= listit_next(lit);
		  list_first_reuse(p->adjs, &ait);
		  while (listit_has_next(ait)) {
			 ppairing a= listit_next(ait);
// Computing gap on genomic
			 if (a->t != SINK_PAIRING_START) {
				bool compact= false;
				if (a->t + a->l - p->t == a->p + a->l - p->p) {
				  compact= (a->t >= p->t + p->l) && (a->t - p->t - p->l <= 3);
				}
				if (compact) {
				  DEBUG("Compacting edge between pairing (%4d, %4d, %4d) and "
						  "(%4d, %4d, %4d) because of the small gap.",
						  p->p, p->t, p->l, a->p, a->t, a->l);
				  removed= true;
				  list_remove_at_iterator(ait, (delete_function)noop_free);
// Removing edge
#ifndef NDEBUG
				  bool ris=
#endif
					 list_remove_element(a->incs, p, (delete_function)noop_free);
#ifndef NDEBUG
				  my_assert(ris);
#endif
// Replacing edge with a vertex
				  ppairing new_v= pairing_create();
				  new_v->p= p->p;
				  new_v->t= p->t;
				  new_v->l= a->p + a->l - p->p;
				  DEBUG("Created new vertex (%4d, %4d, %4d).", PAIRING(new_v));
				  copy_adjacencies(new_v, a);
				  copy_incidencies(new_v, p);
				  list_add_to_tail(EA_get(V, i), new_v);
				}
			 }
		  }
		}
	 }
	 remove_other_sources_and_sinks(V);
  } while (removed);
  listit_destroy(ait);
  listit_destroy(lit);
}

void
simplify_meg(pext_array V, pconfiguration config) {
  remove_useless_edges(V, config);
  remove_other_sources_and_sinks(V);
}


static void
add_vertex(pgraph graph, ppairing v) {
  static int next_id= 0;
  NOT_NULL(graph);
  NOT_NULL(v);

  v->id= next_id;
  ++next_id;

  EA_insert(graph, v);
}

pgraph
meg2graph(pext_array meg) {
  DEBUG("Transforming the MEG into a graph.");
  NOT_NULL(meg);

  pgraph graph= EA_create();

  unsigned int i;
  for (i= 0; i< EA_size(meg); ++i) {
	 plist Vi= EA_get(meg, i);
	 plistit lit= list_first(Vi);
	 while (listit_has_next(lit)) {
		ppairing p= listit_next(lit);
		add_vertex(graph, p);
	 }
	 listit_destroy(lit);
  }

  return graph;
}


void
graph_destroy(pgraph graph) {
  EA_destroy(graph, (delete_function)noop_free);
}

void
dfs_visit(pgraph G,
			 int** out_dtime,
			 int** out_ftime,
			 int** out_topological_ids,
			 bool* out_is_acyclic) {
  DEBUG("Starting the DFS visit of the graph.");
  NOT_NULL(G);
  NOT_NULL(out_is_acyclic);
  const size_t nv= EA_size(G);
  DEBUG("The graph has %zu vertices.", nv);
// Allocate the arrays if they are NULLs
  int* dtime= NULL;
  if (*out_dtime == NULL) {
	 dtime= NPALLOC(int, nv);
	 *out_dtime= dtime;
  } else {
	 dtime= *out_dtime;
  }
  int* ftime= NULL;
  if (*out_ftime == NULL) {
	 ftime= NPALLOC(int, nv);
	 *out_ftime= ftime;
  } else {
	 ftime= *out_ftime;
  }
  int* ids= NULL;
  if (*out_topological_ids == NULL) {
	 ids= NPALLOC(int, nv);
	 *out_topological_ids= ids;
  } else {
	 ids= *out_topological_ids;
  }
  bool is_acyclic= true;
// Ensure that the ids are correct and count the inbound edges
  int* color= NPALLOC(int, nv);  // 0 white, 1 grey, 2 black
  for (unsigned int i= 0; i<nv; ++i) {
	 color[i]= 0;
  }
  for (unsigned int i= 0; i<nv; ++i) {
	 ((ppairing)EA_get(G, i))->id= i;
  }
// Push source vertices
  pintlist S= intlist_create();
  for (unsigned int i= 0; i<nv; ++i) {
	 if (list_size(((ppairing)EA_get(G, i))->incs)==0)
		intlist_add_to_tail(S, i);
  }
  if (intlist_is_empty(S)) { // The graph is cyclic
	 is_acyclic= false;
  }
  unsigned int visited= 0;
  unsigned int time= 0;
  unsigned int progr_id= nv;
  plistit ait= NULL;
  do {
	 while (!intlist_is_empty(S)) {
		int id_v= intlist_remove_from_tail(S);
		TRACE("Extracted vertex %d.", id_v);
		ppairing v= EA_get(G, id_v);
		if (color[id_v]==0) {
		  TRACE("Discovered pairing (%4d, %5d, %4d)", PAIRING(v));
		  color[id_v]= 1;
		  dtime[id_v]= time;
		  ++time;
		  intlist_add_to_tail(S, id_v);
		  list_first_reuse(v->adjs, &ait);
		  while (listit_has_next(ait)) {
			 ppairing a= (ppairing)listit_next(ait);
			 my_assert(a->id < nv);
// Check the color of the adjacent vertex
			 if (color[a->id]==0) { // vertex white -> tree edge
				intlist_add_to_tail(S, a->id);
				TRACE("Added vertex %d.", a->id);
			 } else if (color[a->id]==1) { // vertex grey -> back edge
				is_acyclic= false;
			 } else { // vertex black or enqueued but not yet processed -> forward or cross edge
// do nothing
			 }
		  }
		} else if (color[id_v]==1) {
		  TRACE("Finished visiting pairing (%4d, %5d, %4d)", PAIRING(v));
		  color[id_v]= 2;
		  ftime[id_v]= time;
		  ++time;
		  --progr_id;
		  ids[id_v]= progr_id;
		  ++visited;
		} else {
// do nothing
		  my_assert(color[id_v]==2);
		  TRACE("The vertex has been already visited.");
		}
	 }
	 for (unsigned int i=0; i<nv && intlist_is_empty(S); ++i) {
		if (color[i]==0) {
		  is_acyclic= false;
		  intlist_add_to_tail(S, i);
		  TRACE("Added vertex %d.", i);
		}
	 }
  } while (!intlist_is_empty(S));

  my_assert(visited==nv);

  listit_destroy(ait);
  intlist_destroy(S);
  pfree(color);

  *out_is_acyclic= is_acyclic;
}

static int
compare_pairings_by_id(ppairing * p1, ppairing * p2) {
  return ((*p1)->id-(*p2)->id);
}

void
topological_sort(pgraph G, bool* out_is_acyclic) {
  DEBUG("Starting the topological sort of the graph.");
  NOT_NULL(G);
  NOT_NULL(out_is_acyclic);
  int* dtime= NULL;
  int* ftime= NULL;
  int* ids= NULL;
  dfs_visit(G, &dtime, &ftime, &ids, out_is_acyclic);
  if (*out_is_acyclic) {
	 const size_t nv= EA_size(G);
	 for (unsigned int i= 0; i<nv; ++i) {
		ppairing p= EA_get(G, i);
		p->id= ids[i];
	 }
// Sort the array
	 for (unsigned int i= 0; i<nv; ++i) {
		ppairing pa= EA_get(G, i);
		while (pa->id != i) {
		  TRACE("Swapping elements in position %4d and %4d.", pa->id, i);
		  ppairing pb= EA_get(G, pa->id);
		  EA_set(G, pa->id, pa);
		  EA_set(G, i, pb);
		  pa= pb;
		}
	 }
// Sort adjacencies and incidencies
	 for (unsigned int i= 0; i<nv; ++i) {
		ppairing pa= EA_get(G, i);
		list_sort(pa->adjs, (comparator)compare_pairings_by_id);
		list_sort(pa->incs, (comparator)compare_pairings_by_id);
	 }
  } else {
	 WARN("The graph was cyclic. Topological sort has not been actually performed.");
  }
  pfree(dtime);
  pfree(ftime);
  pfree(ids);
}


void
transitive_reduction(pgraph G) {
  NOT_NULL(G);
  INFO("Starting transitive reduction.");
  bool is_acyclic= false;
  topological_sort(G, &is_acyclic);
  if (!is_acyclic) {
	 FATAL("The graph is cyclic. Transitive reduction not possible! Terminating.");
	 fail();
  }
  const size_t nv= EA_size(G);
  plist* outs_star= NPALLOC(plist, nv);
  plist* outs_red= NPALLOC(plist, nv);
  plist* outs_red_inc= NPALLOC(plist, nv);
  for (size_t i= 0; i<nv; ++i) {
	 outs_star[i]= list_create();
	 outs_red[i]= list_create();
	 outs_red_inc[i]= list_create();
  }
  pbit_vect out_star_v= BV_create(nv);
  plistit ait= NULL;
  plistit wai= NULL;
  size_t removed_edges= 0;
  for (size_t i= nv; i>=1;) {
	 --i;
	 ppairing v= EA_get(G, i);
	 DEBUG("Analysis of pairing %4d (%4d,%7d,%4d).", v->id, PAIRING(v));
	 my_assert(v->id == i);
	 BV_clear(out_star_v);

	 BV_set(out_star_v, i, true);
	 list_add_to_tail(outs_star[i], v);

	 list_first_reuse(v->adjs, &ait);
	 while (listit_has_next(ait)) {
		ppairing w= listit_next(ait);
		TRACE("  adjacent pairing %4d (%4d,%7d,%4d)...", w->id, PAIRING(w));
		if ((! BV_get(out_star_v, w->id)) ||
			 ((w->p < v->p) || (w->t < v->t)) ||
			 ((w->p + w->l < v->p + v->l) || (w->t + w->l < v->t + v->l))) {
#if LOG_LEVEL_TRACE <= LOG_THRESHOLD
		  if (!BV_get(out_star_v, w->id)) {
			 TRACE("    ...NOT in the transitive closure.");
		  } else if ((w->p < v->p) || (w->t < v->t)) {
			 TRACE("    ...in the transitive closure BUT it starts earlier, so it will be kept.");
		  } else if ((w->p + w->l < v->p + v->l) || (w->t + w->l < v->t + v->l)) {
			 TRACE("    ...in the transitive closure BUT it ends earlier, so it will be kept.");
		  }
#elif LOG_LEVEL_INFO <= LOG_THRESHOLD
		  if ((w->p < v->p) || (w->t < v->t)) {
			 INFO("Edge (%4d,%7d,%4d) -> (%4d,%7d,%4d) "
					"is in the transitive closure but "
					"it will not be removed because "
					"(%4d,%7d,%4d) starts 'earlier' than "
					"(%4d,%7d,%4d).",
					PAIRING(v), PAIRING(w),
					PAIRING(w), PAIRING(v));
		  } else if ((w->p + w->l < v->p + v->l) || (w->t + w->l < v->t + v->l)) {
			 INFO("Edge (%4d,%7d,%4d) -> (%4d,%7d,%4d) "
					"is in the transitive closure but "
					"it will not be removed because "
					"(%4d,%7d,%4d) ends 'earlier' than "
					"(%4d,%7d,%4d).",
					PAIRING(v), PAIRING(w),
					PAIRING(w), PAIRING(v));
		  }
#endif
		  list_add_to_tail(outs_red[i], w);
		  list_add_to_tail(outs_red_inc[w->id], v);
		  if (!((w->p + w->l < v->p + v->l) || (w->t + w->l < v->t + v->l))) {
			 list_first_reuse(outs_star[w->id], &wai);
			 while (listit_has_next(wai)) {
				ppairing wa= listit_next(wai);
				if (!BV_get(out_star_v, wa->id)) {
				  if ((v->t <= wa->t) && (v->p <= wa->p) &&
						(v->t + v->l <= wa->t + wa->l) && (v->p + v->l <= wa->p + wa->l)) {
					 TRACE("    Adding %4d (%4d,%7d,%4d) to the "
							 "transitive closure of pairing %4zd.",
							 wa->id, PAIRING(wa), i);
					 BV_set(out_star_v, wa->id, true);
					 list_add_to_tail(outs_star[i], wa);
				  } else {
					 WARN("Edge (%4d,%7d,%4d) -> (%4d,%7d,%4d) "
							"not added to transitive closure because of "
							"pairing order.",
							PAIRING(v), PAIRING(wa));
				  }
				}
			 }
		  }
		} else {
		  TRACE("    ...in the transitive closure.");
		  DEBUG("Transitive reduction has removed edge between "
				  "(%4d,%7d,%4d) and (%4d,%7d,%4d).",
				  PAIRING(v), PAIRING(w));
		  ++removed_edges;
		}
	 }
  }
  BV_destroy(out_star_v);
  listit_destroy(ait);
  listit_destroy(wai);
  for (size_t i= 0; i<nv; ++i) {
	 list_destroy(outs_star[i], (delete_function)noop_free);
	 ppairing v= EA_get(G, i);
	 list_destroy(v->adjs, (delete_function)noop_free);
	 list_destroy(v->incs, (delete_function)noop_free);
	 v->adjs= outs_red[i];
	 v->incs= outs_red_inc[i];
  }
  pfree(outs_star);
  pfree(outs_red);
  pfree(outs_red_inc);
  INFO("Transitive reduction terminated. Removed %4zd edges.", removed_edges);
}
