/**
 *
 *
 *                              PIntron
 *
 * A novel pipeline for computational gene-structure prediction based on
 * spliced alignment of expressed sequences (ESTs and mRNAs).
 *
 * Copyright (C) 2010  Raffaella Rizzi
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
/**
 *
 * @file agree-introns.h
 *
 * Procedures for agreeing an intron set.
 *
 **/

#ifndef _AGREE_INTRONS_H_
#define _AGREE_INTRONS_H_

#include "types.h"


/*
 * Try the agreement of a intron to some intron inside the list
 */

bool try_agreement_to_intron_list(char *, pintron, plist, int);
bool try_agreement_to_intron_list_on_single_site(char *, pintron, plist, plist);

bool try_agreement(char *, pintron, pgenomic_intron, int);

bool try_agreement_on_single_site(char *, pintron, pgenomic_intron, plist);

bool try_agreement_on_donor_site(char *, pintron, pgenomic_intron, plist);
bool try_agreement_on_acceptor_site(char *, pintron, pgenomic_intron, plist);

bool find_better_intron(char *, pintron, plist);

bool try_agreement_to_a_burset_frequency_list(char *, pintron, plist, plist, unsigned int);

/*
 * Set the two agree flags: (1) if the intron can be agreed and (2) what type of EST
 * supports it
 */
pintron set_agree_flags(pintron);

plist get_exon_composition_from_an_intron_composition(plist);

/*
 * Returns a pintron list from a pfactor list and update the gen_intron_list
 */
plist get_intron_composition_from_an_exon_composition(pEST_info, int, char *, plist, plist, bool);

/*
 * Adds a genomic intron to a list e returns the updated list.
 */
plist add_genomic_intron(char *, plist, int, int, pgenomic_intron *);

/*
 * Returns the error that the agreement (to an intron) would produce, but it does not perform any agreement
 */
unsigned int get_agreement_error(char *, pintron, pgenomic_intron);

/*
 * Returns the error that the agreement (to a start-end genomic region) would produce, but it does not perform any agreement
 */
unsigned int get_agreement_error_start_end(char *, pintron, int, int);

/*
 * Correct the EST alignment depending on the genomic intron
 */
void Correct_est_alignment(char *, pintron);

/*
 * Returns the Burset frequency of a genomic intron (0 it the pattern does not exist)
 */
int get_intron_Burset_frequency(char *, pgenomic_intron);

pgenomic_intron set_intron_Burset_frequency(pgenomic_intron);

int get_intron_Burset_frequency_start_end(char *, int, int);

#endif

