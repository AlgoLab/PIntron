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
** refine.c
**
** Made by (Yuri)
** Login   <yuri@dottxxii-12.dottorato.disco.unimib.it>
**
** Started on  Wed Jul 22 14:15:30 2009 Yuri
** Last update Sun May 12 01:17:25 2002 Speed Blue
*/

#include "refine.h"

#include "util.h"
#include "log.h"
#include "refine-intron.h"

#define INDEL_COST 1
#define MISMATCH_COST 1
#define MATCH_COST 0

unsigned int*
edit_distance(const char* const s1, const size_t ls1,
				  const char* const s2, const size_t ls2) {
unsigned int* M= NPALLOC(unsigned int, (ls1+1)*(ls2+1));
  M[0]= 0;
  for (size_t i= 0; i<ls1; ++i) {
	 M[i+1]= M[i]+INDEL_COST;
  }
  for (size_t i= 0, j=0; i<ls2; ++i, j+=ls1+1) {
	 M[j+ls1+1]= M[j]+INDEL_COST;
  }
  size_t pos= ls1+2;
  for (size_t i2= 0; i2<ls2; ++i2) {
	 const char s2c= s2[i2];
	 for (size_t i1= 0; i1<ls1; ++i1) {
		const unsigned int min_cost_ul= M[pos-ls1-2]+
		  ((s2c==s1[i1])?MATCH_COST:MISMATCH_COST);
		const unsigned int min_cost_u= M[pos-ls1-1]+INDEL_COST;
		const unsigned int min_cost_l= M[pos-1]+INDEL_COST;
		M[pos]= MIN(min_cost_ul, MIN(min_cost_u, min_cost_l));
		++pos;
	 }
	 ++pos;
  }
// Print the matrix
#if defined (LOG_MSG) && (LOG_LEVEL_TRACE <= LOG_THRESHOLD)
  for (int i= 0; i<(ls1+1)*(ls2+1); ++i) {
	 printf("%3d ", M[i]);
	 if ((i+1)%(ls1+1)==0)
		printf("\n");
  }
#endif
  return M;
}

static char*
reverse(const char* const s, const size_t len) {
  char* const rs= c_palloc(len);
  for (size_t i= 0; i<len; ++i)
	 rs[len-1-i]= s[i];
  return rs;
}

bool
refine_borders(const char* const p,
					const size_t len_p,
					const char* const t,
					const size_t len_t,
					const unsigned int max_errs,
					size_t* out_offset_p,
					size_t* out_offset_t1,
					size_t* out_offset_t2,
					unsigned int* out_edit_distance) {
  return general_refine_borders(p, len_p, 0, len_p,
										  t, len_t,
										  max_errs,
										  out_offset_p,
										  out_offset_t1,
										  out_offset_t2,
										  out_edit_distance);
}


bool
general_refine_borders(const char* const p,
							  const size_t len_p,
							  const size_t min_p_cut,
							  const size_t max_p_cut,
							  const char* const t,
							  const size_t len_t,
							  const unsigned int max_errs,
							  size_t* out_offset_p,
							  size_t* out_offset_t1,
							  size_t* out_offset_t2,
							  unsigned int* out_edit_distance) {
  my_assert(min_p_cut <= max_p_cut);
  my_assert(max_p_cut <= len_p);
  const size_t t_win= MIN(len_p+max_errs, len_t);
  const unsigned int* Mp= edit_distance(t, t_win, p, len_p);
  const char* const rt= reverse(t, len_t);
  const char* const rp= reverse(p, len_p);
  const unsigned int* Ms= edit_distance(rt, t_win, rp, len_p);

  unsigned int* const min_pp= NPALLOC(unsigned int, len_p+1);
  unsigned int* const min_sp= NPALLOC(unsigned int, len_p+1);
  unsigned int* const min_pos_pp= NPALLOC(unsigned int, len_p+1);
  unsigned int* const min_pos_sp= NPALLOC(unsigned int, len_p+1);

  min_pp[0]= 0;
  min_pos_pp[0]= 0;
  size_t pos= t_win+1;
  for (size_t i= 0; i<len_p; ++i) {
	 min_pp[i+1]= Mp[pos];
	 min_pos_pp[i+1]= 0;
	 for (size_t j= 0; j<t_win; ++j) {
		++pos;
		if (min_pp[i+1]>Mp[pos]) {
		  min_pp[i+1]= Mp[pos];
		  min_pos_pp[i+1]= j+1;
		}
	 }
	 ++pos;
  }
  min_sp[0]= 0;
  min_pos_sp[0]= 0;
  pos= t_win+1;
  for (size_t i= 0; i<len_p; ++i) {
	 min_sp[i+1]= Ms[pos];
	 min_pos_sp[i+1]= 0;
	 for (size_t j= 0; j<t_win; ++j) {
		++pos;
		if (min_sp[i+1]>Ms[pos]) {
		  min_sp[i+1]= Ms[pos];
		  min_pos_sp[i+1]= j+1;
		}
	 }
	 ++pos;
  }

  size_t off_p= min_p_cut;
  size_t off_t1= min_pos_pp[min_p_cut];
  size_t off_t2= min_pos_sp[len_p-min_p_cut];
  unsigned int min= min_pp[min_p_cut]+min_sp[len_p-min_p_cut];
  int best_burset_freq= getBursetFrequency_adaptor(t, off_t1, len_t-off_t2);
  for (size_t i= min_p_cut+1; i<=max_p_cut; ++i) {
	 const int curr_burset_freq= getBursetFrequency_adaptor(t, min_pos_pp[i],
																			  len_t-min_pos_sp[len_p-i]);
	 const unsigned int curr= min_pp[i]+min_sp[len_p-i];
	 if ((min>curr) ||
		  ((min==curr) && (curr_burset_freq>best_burset_freq))) {
		min= curr;
		off_p= i;
		off_t1= min_pos_pp[i];
		off_t2= min_pos_sp[len_p-i];
		best_burset_freq= curr_burset_freq;
	 }
  }
  *out_offset_p= off_p;
  *out_offset_t1= off_t1;
  *out_offset_t2= len_t-off_t2;
  *out_edit_distance= min;
  pfree(Mp);
  pfree(rt);
  pfree(rp);
  pfree(Ms);
  pfree(min_pp);
  pfree(min_pos_pp);
  pfree(min_sp);
  pfree(min_pos_sp);
  return MAX(0, min)<=max_errs;
}

// Example main()
/*
int main() {
  const char* const s1= "AAAAABBBBBB";
  const char* const s2= "xAAAAxxxxxAxxxBBBBBBx";
  size_t p, t1, t2;
  refine_borders(s1, strlen(s1), s2, strlen(s2),
					  strlen(s1), &p, &t1, &t2);
  printf("P %s\nT %s\n", s1, s2);
  printf("out_p %3d  out_t1 %3d out_t2 %3d\n", p, t1, t2);
}
*/
