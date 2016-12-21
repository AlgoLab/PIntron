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
#include "bit_vector.h"

#include "util.h"
#include <stdio.h>
#include <limits.h>
#include <string.h>

#define _ASSERT_VALID_BV( bv )						\
  my_assert(bv!=NULL);

#define _ASSERT_VALID_POS( bv, i )				\
  my_assert((i)<(bv)->n)
//  my_assert(0<=(i) && (i)<(bv)->n)

pbit_vect BV_create(const unsigned int n)
{
  my_assert(n>0u);
  const unsigned int ncells= (n/_LBTYPE)+1;
  pbit_vect bv= (pbit_vect)calloc(1, sizeof(struct _bit_vect)+sizeof(_BTYPE)*ncells);
  bv->n= n;
  bv->ncells= ncells;
  bv->arr= (_BTYPE*)(bv+1);
  return bv;
}

void BV_clear(pbit_vect bv) {
  _ASSERT_VALID_BV(bv);
  for (unsigned int i= 0; i<bv->ncells; ++i) {
	 bv->arr[i]= (_BTYPE)0;
  }
}

void BV_or(pbit_vect bv1, pbit_vect bv2) {
  _ASSERT_VALID_BV(bv1);
  _ASSERT_VALID_BV(bv2);
  my_assert(bv1->n == bv2->n);
  for (unsigned int i= 0; i<bv1->ncells; ++i) {
	 bv1->arr[i]= bv1->arr[i] | bv2->arr[i];
  }
}


void BV_destroy(pbit_vect bv)
{
  _ASSERT_VALID_BV(bv);
  bv->arr= NULL;
  pfree(bv);
}

#define _LINT unsigned long

static void
transform_coord(_LINT i, _LINT* cella, _BTYPE* mask)
{
  *cella= i/_LBTYPE;
  *mask= 1L<<(i%_LBTYPE);
}

void BV_set(pbit_vect bv, const unsigned int i, bool value)
{
  _ASSERT_VALID_BV(bv);
  _ASSERT_VALID_POS(bv, i);
  _LINT cella;
  _BTYPE mask;
  transform_coord(i, &cella, &mask);
  if (value) {
	 bv->arr[cella]= bv->arr[cella] | mask;
  } else {
	 bv->arr[cella]= bv->arr[cella] & (~mask);
  }
}

#define MODULO %


void BV_set_block(pbit_vect bv, const unsigned int i, _BTYPE block)
{
  _ASSERT_VALID_BV(bv);
  _ASSERT_VALID_POS(bv, i);
  my_assert((i MODULO (_LBTYPE))==0);
  _LINT cella;
  _BTYPE mask;
  transform_coord(i, &cella, &mask);
  bv->arr[cella]= block;
}

bool BV_get(pbit_vect bv, const unsigned int i)
{
  _ASSERT_VALID_BV(bv);
  _ASSERT_VALID_POS(bv, i);
  _LINT cella;
  _BTYPE mask;
  transform_coord(i, &cella, &mask);
  return (bv->arr[cella] & mask)!=0;
}

_BTYPE BV_get_block(pbit_vect bv, const unsigned int i)
{
  _ASSERT_VALID_BV(bv);
  _ASSERT_VALID_POS(bv, i);
  my_assert(i MODULO _LBTYPE==0);
  _LINT cella;
  _BTYPE mask;
  transform_coord(i, &cella, &mask);
  return bv->arr[cella];
}

#undef MODULO

_BTYPE
BV_get_unaligned_block(pbit_vect v, const unsigned int i, const unsigned int l)
{
  _ASSERT_VALID_BV(v);
  _ASSERT_VALID_POS(v, i);
  my_assert(l<=_LBTYPE);
  _LINT pb1= (i/_LBTYPE)*_LBTYPE;
  _LINT pb2= pb1+_LBTYPE;
  _LINT mod= i%_LBTYPE;
  _BTYPE br, tmp1, tmp2, mask;
  if (mod==0)
	 br= BV_get_block(v, i);
  else {
	 tmp1= BV_get_block(v, pb1);
	 if ((pb1+l<pb1)||(pb2>v->n))
		tmp2= 0L;
	 else
		tmp2= BV_get_block(v, pb2);
//	 print_block("tmp1 ", tmp1, _LBTYPE);
//	 print_block("tmp2 ", tmp2, _LBTYPE);
	 br= (tmp1>>mod) | (tmp2<<(_LBTYPE-mod));
  }
//  print_block("br   ", br, l);
  unsigned int j;
  mask= 0L;
  for (j= 0; j<l; ++j) {
	 mask= mask << 1;
	 mask= mask | 1L;
  }
//  printf("L= %d\n", l);
//  print_block("mask ", mask , _LBTYPE);
  br= br & mask;
//  print_block("br   ", br, _LBTYPE);
  return br;
}




void BV_print(pbit_vect bv)
{
  _ASSERT_VALID_BV(bv);
  unsigned int i;
  for (i= 0; i<bv->n; ++i) {
	 printf("%c", BV_get(bv, i)?'1':'0');
  }
  printf("\n");
}


pbit_vect BV_clone(pbit_vect bv)
{
  pbit_vect ris= PALLOC(struct _bit_vect);
  ris->n= bv->n;
  ris->ncells= bv->ncells;
  ris->arr= NPALLOC(_BTYPE, bv->ncells);
  memcpy(ris->arr, bv->arr, bv->ncells*sizeof(_BTYPE));
  return ris;
}

void BV_copy(pbit_vect ris, pbit_vect bv)
{
  memcpy(ris->arr, bv->arr, bv->ncells*sizeof(_BTYPE));
}

bool BV_contained(pbit_vect bv1, pbit_vect bv2)
{
  _BTYPE r=0L;
  unsigned int i=0;

  while(i<(bv1->n)){

	 _BTYPE bbv1=BV_get_block(bv1,i);
	 _BTYPE bbv2=BV_get_block(bv2,i);

	 r=~((~bbv1)|bbv2);
	 if(r!=0){ return false; }
	 i+=_LBTYPE;
  }
  return true;
}

bool BV_all_true(pbit_vect bv)
{
  for (unsigned int i=0; i<bv->n; ++i){
	 if (!BV_get(bv,i))
		return false;
  }
  return true;
}

#ifdef LOG_MSG
char* BV_to_string(pbit_vect bv)
{
  my_assert(bv->n<10001);
  static char string[10001];

  for (unsigned int i= 0; i<bv->n; ++i) {
	 if(BV_get(bv,i))
		string[i]='1';
	 else
		string[i]='0';
  }

  string[bv->n]='\0';
  return string;
}
#endif


#undef _ASSERT_VALID_BV
#undef _ASSERT_VALID_POS
