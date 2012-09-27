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
 * @file est-factorizations.h
 *
 * Funzioni per trovare i maximal embeddings delle EST.
 *
 **/

#ifndef _EST_FACTORIZATIONS_H_
#define _EST_FACTORIZATIONS_H_

#include "types.h"
#include "ext_array.h"
#include "int_list.h"
#include "configuration.h"
#include "my_time.h"

//Include

//Computa per data una EST tutte le fattorizzazioni ammissibili a partire
//dal grafo degli embedding massimali (pext_array)
//pEST get_EST_factorizations(pconfiguration, pEST_info, pext_array);
pEST get_EST_factorizations(pEST_info, pext_array, pconfiguration, pEST_info,
									 pmytime_timeout);

char* real_substring(int index, int length, const char* const string);

void print_split_string_on_stderr(const int, char*);

/*
 * Discards low complexity exons from a factorization and retains the best part
 */
plist clean_low_complexity_exons(plist, char *, char *);
plist clean_low_complexity_exons_2(plist, char *, char *);

plist clean_external_exons(plist, char *, char *);

plist clean_noisy_exons(plist, char *, char *, bool);

bool check_exon_start_end(plist);

bool check_small_exons(plist);

plist update_with_subfact_with_best_coverage(plist, pintlist);

plist add_if_not_exists(plist, plist, pconfiguration, bool*);

bool check_for_not_source_sink_factorization(plist, int);

plist handle_endpoints(plist, char *, char *);

bool check_est_coverage(plist, char *);

#endif
