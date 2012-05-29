/**
 *
 *
 *                              PIntron
 *
 * A novel pipeline for computational gene-structure prediction based on
 * spliced alignment of expressed sequences (ESTs and mRNAs).
 *
 * Copyright (C) 2010  Yuri Pirola, Raffaella Rizzi
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
** configuration.h
*/

#ifndef _CONFIGURATION_H_
#define _CONFIGURATION_H_

#include <stdbool.h>

struct _configuration {

// Minimum lenght of a factor (l)
// Suggested value: at least 10
  unsigned int min_factor_len;

//The minimum length of an intron
//Suggested value:
  int min_intron_length;

//The maximum length of an intron
//Suggested value:
  int max_intron_length;

//When a deepest common node N was found, the algorithm
//processes common nodes whose string depth is greater
//than (min_string_depth_rate*string_depth[N])
//Rationale: if a long common factor was found, then
//considering a much
//shorter prefix does not seem useful.
//Suggested value: at least 0.2 (=20%)
//Valid values: (0, 1]
//The value 0 disables this feature.
  double min_string_depth_rate;

//The maximum rate of the EST lenght that can be discarded
//as prefix (or suffix, resp.) by a factorization.
//Suggested value: at most 0.5 (=50%)
//Valid values: (0, 1]
//The value 1.0 disables this feature. (sono inutili ma se li tolgo non funziona...)
  double max_prefix_discarded_rate;
  double max_suffix_discarded_rate;

//The maximum prefix/suffix that can be discarded from an EST (for
//linking to a source/sink node in the MEG))
//If set to -1, this check is not performed
// **NOTE**
// These options are used only during the constructions of
// the factorizations, NOT for the MEG.
  int  max_prefix_discarded;
  int  max_suffix_discarded;

  //The maximum difference (nt) for considering the same two splicing sites
  unsigned int max_site_difference;

  //The maximum number of factorizations in order to consider an
  //EST not an artifact
  //If set to -1, this check is not performed
  int max_number_of_factorizations;

  //The maximum coverage difference for accepting a factorization
  //wrt max coverage found
  double max_coverage_diff;

  //The maximum exon number difference for accepting a
  //factorization wrt min exon number found
  int max_exonNUM_diff;

  //The maximum gap length (on P) difference for accepting a
  //factorization wrt min gap length found
  int max_gapLength_diff;

  //If set to 1, external factors are retained in the output
  //factorizations, otherwise the first factor is deleted and the last factor is retained only if
  //the EST has a polyA chain
  char retain_externals;

  //The maximum number of pairings that can compose the vertex
  //set of a MEG
  // AND
  //The maximum frequency of the shortest pairings in a vertex
  //set of a MEG
  //Rationale: if more than max_freq*tot pairings have minimum
  //length and we have more than max pairings,
  //than it is plausible we are allowing too short pairings.
  unsigned int max_pairings_in_MEG;
  double max_freq_shortest_pairing;

  //Parameters for performin gap alignment in intron refinment
  int suffpref_length_on_est;
  int suffpref_length_for_intron;
  int suffpref_length_on_gen;

  bool trans_red;
  bool short_edge_comp;
};

typedef struct _configuration* pconfiguration;

pconfiguration config_create(int argc, char** argv);

pconfiguration config_clone(pconfiguration src);

void config_destroy(pconfiguration config);


#endif /* _CONFIGURATION_H_ */
