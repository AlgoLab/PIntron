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
** meg-simplification.h
**
** Made by Yuri Pirola
** Login   <yuri@dottxxii-12.dottorato.disco.unimib.it>
**
** Started on  Wed Dec  9 13:12:10 2009 Yuri Pirola
** Last update Wed Dec  9 13:12:10 2009 Yuri Pirola
*/

#ifndef __MEG_SIMPLIFICATION_H__
#define __MEG_SIMPLIFICATION_H__

#include <stdbool.h>

#include "ext_array.h"
#include "configuration.h"

/**
 * This file contains procedures that simplifies the CMEG graph by removing unnecessary
 * edges (w.r.t. the spliced alignment problem).
 ***/

typedef pext_array pgraph;


/**
 *
 * Return the number of vertexes and edges of the CMEG @p V.
 *
 * @param V The CMEG to analyse.
 * @param tot_pairings (output) the number of vertexes
 * @param tot_edges (output) the number of edges
 *
 **/
void MEG_stats(pext_array V, size_t* tot_pairings, size_t* tot_edges);

/**
 *
 * Decide if the CMEG @p V is too complex to perform the compact short edges.
 *
 * @param V The CMEG to analyse.
 * @param config The configuration structure
 *
 **/
bool
is_too_complex_for_compaction(pext_array V, pconfiguration config);

/**
 *
 * Decide if the CMEG @p V is too complex.
 *
 * @param V The CMEG to analyse.
 * @param config The configuration structure
 *
 **/
bool
is_too_complex(pext_array V, pconfiguration config);

/**
 *
 * Remove edges and vertices from the CMEG @p V that cannot be reached from the
 * source or that cannot reach the sink.
 *
 * @param V The CMEG to simplify.
 *
 **/
void
remove_other_sources_and_sinks(pext_array V);

/**
 *
 * Remove edges from the CMEG @p V that cannot corresponds to mutations nor introns
 * because of their length.
 *
 * @param V The CMEG to simplify.
 * @param config The configuration structure
 *
 **/
void
remove_useless_edges(pext_array V, pconfiguration config);


void
compact_short_edges(pext_array V, pconfiguration config);

/**
 *
 * Apply simplification procedures of the CMEG @p V.
 *
 * @param V The CMEG to simplify.
 * @param config The configuration structure
 *
 **/
void
simplify_meg(pext_array V, pconfiguration config);

/**
 *
 * Transforms the (C)MEG @p meg into a "normal" graph intended as an extendible array
 * of ppairing.
 * It enforces arr[i]->id= i
 *
 * @param meg The (C)MEG to transform.
 * @return An extendible array containing the vertex set of the MEG graph.
 *
 **/
pgraph
meg2graph(pext_array meg);

/**
 *
 * Destroy the graph @p graph without freeing the memory of vertices, adjacencies, and
 * incidencies lists.
 * @param graph the graph to destroy.
 *
 **/
void
graph_destroy(pgraph graph);

/**
 *
 * Perform a DFS visit of the graph @p G.
 *
 * Remark: It may modify the field id of the pairings.
 * Note: Outbound array parameters are automatically allocated if NULLs.
 *
 * @param G The graph to visit.
 * @param out_dtime The discovering time of each vertex (out)
 * @param out_ftime The finishing time of each vertex (out)
 * @param out_topological_ids The topological order of the vertices (out) VALID ONLY IF
 * ACYCLIC.
 * @param out_is_acyclic The flag that indicates if the graph is acyclic (out)
 *
 **/
void
dfs_visit(pgraph G, int** out_dtime, int** out_ftime, int** out_topological_ids, bool* out_is_acyclic);

/**
 *
 * Perform a topological sort of the graph @p G.
 * The result is stored in the field id of the ppairing.
 *
 * @param G The graph to sort.
 * @param out_is_acyclic The flag that indicates if the graph is acyclic (out). The
 * sorting process is not valid if @p out_is_acyclic is false.
 *
 **/
void
topological_sort(pgraph G, bool* out_is_acyclic);

/**
 *
 * Perform a transitive reduction of the acyclic graph @p G.
 *
 * Remark: The digraph @p G must be acyclic, otherwise it fails.
 *
 * @param G The graph to sort.
 * @param out_is_acyclic The flag that indicates if the graph is acyclic (out). The
 * sorting process is not valid if @p out_is_acyclic is false.
 *
 **/
void
transitive_reduction(pgraph G);

#endif /* !__MEG_SIMPLIFICATION_H__ */
