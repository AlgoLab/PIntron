/**
 *
 *
 *                              PIntron
 *
 * A novel pipeline for computational gene-structure prediction based on
 * spliced alignment of expressed sequences (ESTs and mRNAs).
 *
 * Copyright (C) 2010,2012  Yuri Pirola, Raffaella Rizzi
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

#ifndef _COMPUTE_EST_FACT_H_
#define _COMPUTE_EST_FACT_H_


#include <stdio.h>

#include "types.h"
#include "my_time.h"
#include "configuration.h"

#include "aug_suffix_tree.h"


pEST
compute_est_fact(pEST_info gen,
					  pEST_info est,
					  LST_STree* tree,
					  ppreproc_gen pg,
					  FILE* floginfoext,
					  FILE* fmeg, FILE* fpmeg, FILE* ftmeg,
					  FILE* fintronic,
					  pmytime pt_alg, pmytime pt_comp, pmytime pt_io,
					  pconfiguration shared_config);


#endif
