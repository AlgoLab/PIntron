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
** max-emb-graph.h
*/

#ifndef _MAX_EMB_GRAPH_H_
#define _MAX_EMB_GRAPH_H_

#include <stdio.h>

#include "types.h"
#include "configuration.h"
#include "ext_array.h"
#include "aug_suffix_tree.h"

#define PAIRING( P ) P->p, P->t, P->l


int
compute_gl(pconfiguration config);

pext_array
build_vertex_set(pEST_info pattern,
					  LST_STree* tree,
					  const ppreproc_gen const pg,
					  pconfiguration config);

void
build_edge_set(pext_array V,
					pconfiguration config
					)
#ifndef __ICC
  __attribute__ ((nonnull))
#endif
  ;

void
add_intronic_edges_to_file(FILE* f, pext_array V);

// Funzione per la creazione del sorgente dot di un MEG.
// Abilitata solo se LOG_GRAPHS e' definito, altrimenti non compie nulla.

// Vedere file .c per parametri di formattazione

void
save_meg_to_filename(pext_array V, const char* const filename);

void
print_meg(pext_array V, FILE* f);


#endif /* _MAX-EMB-GRAPH_H_ */
