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

#include "compute-alignments.h"
#include "types.h"
#include <string.h>

#include "log.h"

//#define LOG_THRESHOLD LOG_LEVEL_TRACE

plist compute_alignment(char *EST_seq, char *genomic_seq, bool only_one_align){
	const size_t n=strlen(EST_seq);
	const size_t m=strlen(genomic_seq);

	//Per ora computa sempre un solo allineamento
	fail_if(!only_one_align);

	plist alignments=list_create();
// If the strings are equal, then compute the naive alignment
	if ((EST_seq==genomic_seq) ||
		 ((n==m) && (strcmp(EST_seq, genomic_seq)==0))) {
	  palignment alignment=alignment_create(n+1);
	  alignment->score= 0;
	  strncpy(alignment->EST_alignment, EST_seq, n+1);
	  strncpy(alignment->GEN_alignment, genomic_seq, n+1);
	  alignment->alignment_dim= n;

	  list_add_to_tail(alignments, alignment);
	  return alignments;
	}

	char *Mdir=NPALLOC(char, (n+1)*(m+1));

	memset(Mdir, 0u, (n+1)*(m+1));

//Per ora fare calcolare un solo allineamento
	if (only_one_align) {
		palignment alignment=alignment_create(n+m+1);
		alignment->score= ComputeAlignMatrix(EST_seq, n, genomic_seq, m, Mdir);
		TracebackAlignment(m, alignment, EST_seq, genomic_seq, Mdir, (int)n, (int) m);

		alignment->EST_alignment[alignment->alignment_dim]='\0';
		alignment->GEN_alignment[alignment->alignment_dim]='\0';

		list_add_to_tail(alignments, alignment);
	}
	else{
		//Computazione di piu' allineamenti ==> DA FARE
	  fail();
	}

	pfree(Mdir);

	return alignments;
}

unsigned int
ComputeAlignMatrix(const char * const EST_seq,
						 const size_t EST_len,
						 const char * const genomic_seq,
						 const size_t genomic_len,
						 char * const Mdir) {

	my_assert(EST_seq != NULL);
	my_assert(genomic_seq != NULL);
	my_assert(Mdir != NULL);

	const size_t n= EST_len;
	const size_t m= genomic_len;

	unsigned int *M1= NPALLOC(unsigned int, (m+1));
	unsigned int *M2= NPALLOC(unsigned int, (m+1));

	unsigned int i, j;

	for(j=0; j<m+1; j++) {
		M1[j]=j;
	}

	for(i=1; i<n+1; i++) {
	  const char est= EST_seq[i-1];
	  const bool is_est_a_N= ((est == 'n') || (est == 'N'));
	  M2[0]= i;
	  unsigned int left= i;
	  for(j=1; j<m+1; j++) {
		 unsigned int current= M1[j-1];
		 char current_dir= (char)0; //0 per allineamento caratteri in i-1 e j-1
		 if ((est == genomic_seq[j-1]) ||
			  is_est_a_N ||
			  (genomic_seq[j-1] == 'n') ||
			  (genomic_seq[j-1] == 'N')) {
			//Costo match a 0
		 } else {
			current += 1; //Costo mismatch a +1
		 }

		 if (current > M1[j]+1){
			current= M1[j]+1;  //Costo spazio a +1
			current_dir= (char)1;  //1 per cancellazione in genomic_seq; il carattere in i-1 matcha con -
		 }

		 if (current > left){ // E' gia' aumentato di 1
			current= left;  //Costo spazio a +1
			current_dir= (char)2;  //2 per cancellazione in EST_seq; il carattere in j-1 matcha con -
		 }
		 M2[j]= current;
		 Mdir[(i*m)+j]= current_dir;
		 left= current+1;
	  }
	  MY_SWAP(unsigned int*, M1, M2);
	}

	const unsigned int score= M1[m];

	pfree(M1);
	pfree(M2);

	return score;
}

void TracebackAlignment(const size_t m, palignment alignment, const char * const EST_seq, const char * const genomic_seq, const char * const Mdir, int i, int j){

	char direction=0;
	alignment->alignment_dim= 0;

// The alignment is computed backwards
	while (i>0 && j>0) {
		direction=Mdir[i*m+j];

		if (direction == 0){
			alignment->EST_alignment[alignment->alignment_dim]=EST_seq[i-1];
			alignment->GEN_alignment[alignment->alignment_dim]=genomic_seq[j-1];
			alignment->alignment_dim=alignment->alignment_dim+1;
			--i;
			--j;
		} else if (direction == 1){
		  alignment->EST_alignment[alignment->alignment_dim]=EST_seq[i-1];
		  alignment->GEN_alignment[alignment->alignment_dim]='-';
		  alignment->alignment_dim=alignment->alignment_dim+1;
		  --i;
		} else { // direction == 2
		  alignment->EST_alignment[alignment->alignment_dim]='-';
		  alignment->GEN_alignment[alignment->alignment_dim]=genomic_seq[j-1];
		  alignment->alignment_dim=alignment->alignment_dim+1;
		  --j;
		}
	}

	while (i > 0) {
	  alignment->EST_alignment[alignment->alignment_dim]=EST_seq[i-1];
	  alignment->GEN_alignment[alignment->alignment_dim]='-';
	  alignment->alignment_dim=alignment->alignment_dim+1;
	  --i;
	}
	while (j > 0) {
	  alignment->EST_alignment[alignment->alignment_dim]='-';
	  alignment->GEN_alignment[alignment->alignment_dim]=genomic_seq[j-1];
	  alignment->alignment_dim=alignment->alignment_dim+1;
	  --j;
	}
	alignment->EST_alignment[alignment->alignment_dim]='\0';
	alignment->GEN_alignment[alignment->alignment_dim]='\0';
// Reverse the alignment
	char *p1, *p2;
	for (p1= alignment->EST_alignment, p2= alignment->EST_alignment+alignment->alignment_dim-1;
		  p1 < p2;
		  ++p1, --p2) {
	  *p1 ^= *p2;
	  *p2 ^= *p1;
	  *p1 ^= *p2;
	}
	for (p1= alignment->GEN_alignment, p2= alignment->GEN_alignment+alignment->alignment_dim-1;
		  p1 < p2;
		  ++p1, --p2) {
	  *p1 ^= *p2;
	  *p2 ^= *p1;
	  *p1 ^= *p2;
	}
}


size_t*
edit_distance_matrix(const char* const s1, const size_t l1,
							const char* const s2, const size_t l2) {
  size_t* matrix= NPALLOC(size_t, (l1+1)*(l2+1));
  size_t i, j;
  for (i= 0; i<=l2; ++i) {
	 matrix[i]= i;
  }
  for (i= 0, j=0; i<=l1; ++i, j+= l2+1) {
	 matrix[j]= i;
  }
  size_t base= l2+1;
  for (i= 0; i<l1; ++i) {
	 const char cs1= s1[i];
	 for (j= 0; j<l2; ++j) {
		size_t val= matrix[base+j-l2-1];
		if (cs1 != s2[j]) {
		  val += 1;
		}
		val= MIN(matrix[base+j-l2]+1, val);
		val= MIN(matrix[base+j]+1, val);
		matrix[base+j+1]= val;
	 }
	 base+= (l2+1);
  }
  return matrix;
}


size_t
compute_edit_distance(const char* const s1, const size_t l1,
							 const char* const s2, const size_t l2) {
  if ((l1 == l2) && (strcmp(s1, s2)==0)) {
	 return 0;
  }
  size_t* matrix= edit_distance_matrix(s1, l1, s2, l2);
  const size_t dist= matrix[(l1+1)*(l2+1)-1];
  pfree(matrix);
  return dist;
}


/*
 * In edit is the real error if the procedures returns true
 */
bool K_band_edit_distance(char *seq1, char *seq2, unsigned int upper_bound, unsigned int *edit){
	my_assert(seq1 != NULL);
	my_assert(seq2 != NULL);
	my_assert(edit != NULL);

	size_t length1= strlen(seq1);
	size_t length2= strlen(seq2);

	/*If the two input sequences are the same sequence*/
	if ((length1==length2) &&
		 (strcmp(seq1,seq2)==0)) {
		*edit=0;
		return true;
	}

	if (upper_bound==0) {
	  if ((length1==length2) &&
			(strcmp(seq1,seq2)==0)) {
		 *edit= 0;
		 return true;
	  } else {
		 *edit= 1;
		 return false;
	  }
	}

	if (length1 < length2){
	  char *help=NULL;
	  help=seq1;
	  seq1=seq2;
	  seq2=help;
	  size_t tmp= length1;
	  length1= length2;
	  length2= tmp;
	}

	const size_t n= length1;
	const size_t m= length2;

	/*
	  If the two input sequences are not the same sequence and the edit
	  distance limit is 0, or the absolute length
	  difference is greater than the edit distance limit
	*/
	if ((n - m) > upper_bound) {
		*edit=n-m;
		return false;
	}

	const size_t k= upper_bound;

	if (2*k+1 >= n) {
	  *edit= compute_edit_distance(seq1, n, seq2, m);
	  return (*edit)<=upper_bound;
	}

	size_t * M1= NPALLOC(size_t, (2*k)+1);
	size_t * M2= NPALLOC(size_t, (2*k)+1);

	size_t r, c;
	for (c=0; c <= k; ++c)
	  M1[k+c]= c;

	for (c=0; c < 2*k+1; ++c)
	  M2[c]= k+1;

	size_t d;
	for (r= 1; r <= k; ++r) {
	  M2[k-r]= r;
	  for (c= 1; c < r+k; ++c) {
		 d= M1[k-r+c];
		 if (seq1[c-1] != seq2[r-1])
			d += 1;
		 d= MIN(d, M2[k-r+c-1]+1);
		 d= MIN(d, M1[k-r+c+1]+1);
		 M2[k-r+c]= d;
	  }
	  d= M1[2*k];
	  if (seq1[r+k-1] != seq2[r-1])
		 d += 1;
	  d= MIN(d, M2[2*k-1]+1);
	  M2[2*k]= d;
	  MY_SWAP(size_t*, M1, M2)
	}

	for (r= k+1; r<= n-k; ++r) {
	  M2[0]= M1[0];
	  if (seq1[r-k-1] != seq2[r-1])
		 M2[0] += 1;
	  M2[0]= MIN(M2[0], M1[1]+1);

	  for (c= r+1-k; c < r+k; ++c) {
		 d= M1[c+k-r];
		 if (seq1[c-1] != seq2[r-1])
			d += 1;
		 d= MIN(d, M2[c+k-r-1]+1);
		 d= MIN(d, M1[c+k-r+1]+1);
		 M2[c+k-r]= d;
	  }
	  d= M1[2*k];
	  if (seq1[r+k-1] != seq2[r-1])
		 d += 1;
	  d= MIN(d, M2[2*k-1]+1);
	  M2[2*k]= d;
	  MY_SWAP(size_t*, M1, M2)
	}

	for (r= n+1-k; r<= m; ++r) {
	  M2[0]= M1[0];
	  if (seq1[r-k-1] != seq2[r-1])
		 M2[0] += 1;
	  M2[0]= MIN(M2[0], M1[1]+1);

	  for (c= r+1-k; c <= n; ++c) {
		 d= M1[c+k-r];
		 if (seq1[c-1] != seq2[r-1])
			d += 1;
		 d= MIN(d, M2[c+k-r-1]+1);
		 d= MIN(d, M1[c+k-r+1]+1);
		 M2[c+k-r]= d;
	  }
	  MY_SWAP(size_t*, M1, M2)
	}

	const size_t result= M1[n+k-m];
	pfree(M1);
	pfree(M2);

	*edit=result;

	if (result > upper_bound)
		return false;
	else
		return true;
}

/*
 *	 MAIN TEST
 *
 *#include <stdio.h>
 *
 *int main(int argc, char** argv) {
 *  unsigned int ris;
 *  bool tmp;
 *  tmp= K_band_edit_distance(argv[1], argv[2], atoi(argv[3]), &ris);
 *  printf("%u %s |%s| |%s|\n", ris, tmp?"yes":"no", argv[1], argv[2]);
 *}
 *
 ***/
