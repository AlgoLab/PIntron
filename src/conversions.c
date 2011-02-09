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
#include "log.h"


int get_ABS_coord(int gen_abs_start, int gen_abs_end, int strand, int coord){

	if(strand == 1){
		return gen_abs_start+coord-1;
	}
	else{
		return gen_abs_end-coord+1;
	}
}

void get_ABS_region_start_end(int gen_abs_start, int gen_abs_end, int strand, int start, int end, int *abs_start, int *abs_end){

	if(strand == 1){
		*abs_start=get_ABS_coord(gen_abs_start, gen_abs_end, strand, start);
		*abs_end=get_ABS_coord(gen_abs_start, gen_abs_end, strand, end);
	}
	else{
		*abs_start=get_ABS_coord(gen_abs_start, gen_abs_end, strand, end);
		*abs_end=get_ABS_coord(gen_abs_start, gen_abs_end, strand, start);
	}
}
