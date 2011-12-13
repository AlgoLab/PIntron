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
#include "min_factorization.h"
#include "color_matrix.h"
#include "list.h"
#include "log.h"
#include "bit_vector.h"
#include <assert.h>
#include "sempl_info.h"
#include <stdint.h>

bool valuate_list(pbit_vect comb, plist p)
{
  my_assert(p!=NULL);

  pbit_vect bv;
  plistit list_it_fact;

  list_it_fact=list_first(p);

  while(listit_has_next(list_it_fact)){
	 bv=listit_next(list_it_fact);
	 if(contained(bv,comb))
		{

		  listit_destroy(list_it_fact);
		  return true;
		}
  }
  listit_destroy(list_it_fact);
  return false;
}


bool valuate_combination(pbit_vect comb,plist color_matrix)
{

  my_assert(color_matrix!=NULL);

  bool factorized=true;
  pEST p;
  plistit list_it_est;

  list_it_est=list_first(color_matrix);

  while(listit_has_next(list_it_est)){
	 p=listit_next(list_it_est);

	 factorized=valuate_list(comb,p->bin_factorizations);

	 if(factorized==false){
		listit_destroy(list_it_est);
		return false;
	 }
  }
  listit_destroy(list_it_est);
  DEBUG("*** Combinazione trovata!! ***");
  return true;
}

int set_bin_value( pbit_vect test,psempl ps,int pos,bool value)
{
  int i=-1;
  do{
	 i=i+1;
	 if((BV_get(ps->factors_used,i)==false)&&(BV_get(ps->factors_not_used,i)==false)){
		pos=pos-1;
	 }
  }while(pos!=-1);
  BV_set(test,i,value);
  return i;
}


/*****************************************************
******************************************************/

bool create_combinations(int s, int k,pbit_vect comb,plist color_matrix)
{
  int card=(int)comb->n;
  bool factorized=false;
  if(k==1){
	 for(int cont=s;cont<card;cont++){

		BV_set(comb,cont,true);

		factorized=valuate_combination(comb,color_matrix);
		if(factorized==true){return true;}

		BV_set(comb,cont,false);
	 }
  }
  if(k>1){

	 for(int cont=s;cont<card-(k-1);cont++){

		BV_set(comb,cont,true);

		bool ris= create_combinations(cont+1,k-1,comb,color_matrix);

		if (ris) return true;
		BV_set(comb,cont,false);
	 }
  }
  return false;
}



bool all_false(pbit_vect bv)
{
  for(unsigned int i=0;i<bv->n;i++){
	 if (BV_get(bv,i)) return false;
  }
  return true;
}

void set_all_false(pbit_vect bv)
{
  for(unsigned int i=0;i<bv->n;i++){
	 BV_set(bv,i,false);
  }
}

int count_true(pbit_vect bv)
{
  int cont=0;
  for(unsigned int i=0;i<bv->n;i++){
	 if(BV_get(bv,i)==true){
		cont=cont+1;
	 }
  }
  return cont;
}

int min_number_of_factor(plist bin_factorizations)
{
  int min=0;
  int n_fact;
  plistit list_it;
  pbit_vect bv;

  list_it=list_first(bin_factorizations);

  while(listit_has_next(list_it)){

	 bv=listit_next(list_it);
	 n_fact=count_true(bv);

	 if((n_fact<min)||(min==0))min=n_fact;
  }

  listit_destroy(list_it);
  return min;
}



int max_of_min(plist color_matrix)
{
  int max=0;
  pEST est;
  plistit list_it;

  list_it=list_first(color_matrix);

  while(listit_has_next(list_it)){

	 est=listit_next(list_it);

	 if(min_number_of_factor(est->bin_factorizations)>max) max=min_number_of_factor(est->bin_factorizations);
  }
  listit_destroy(list_it);

  return max;
}


#if defined (LOG_MSG) && (LOG_LEVEL_DEBUG <= LOG_THRESHOLD)

void print_min_fact(pbit_vect bv, plist list_of_unique_fact)
{
  int cont=-1;
  plistit list_it;
  pfactor pf;

  DEBUG("Stampa del minimo numero di fattori necessari alla fattorizzazione di ogni EST");

  list_it=list_first(list_of_unique_fact);

  while(listit_has_next(list_it)){
	 cont=cont+1;
	 pf=listit_next(list_it);

	 if(BV_get(bv,cont)){
		 printf("%d)\t %d\t %d\t %d\t %d\n",cont+1,
				pf->EST_start, pf->EST_end,
				pf->GEN_start,pf->GEN_end);
	 }
  }
  listit_destroy(list_it);
}
#endif

//Ingloba il risultato nel vettore dei fattori certi

void inglobe(pbit_vect result,psempl p)
{
  int cont=0;
// int size=(int)p->factors_used->n;

  for(size_t i=0;i<p->factors_used->n;i++){

	 if((!BV_get(p->factors_used,i))&&(!BV_get(p->factors_not_used,i))){
		BV_set(p->factors_used,i,BV_get(result,cont));
		cont=cont+1;
	 }
	 }
}
void print_n_factorization(plist factorizations,int n)
{
  pfactor pf;
  pfactorization pfact;

  plistit plist_it_factor;
  plistit plist_it_factorizations;

  plist_it_factorizations=list_first(factorizations);
  int cont=0;

  while(listit_has_next(plist_it_factorizations)){

	 pfact=listit_next(plist_it_factorizations);
	 cont=cont+1;
	 if(cont==n){
		plist_it_factor=list_first(pfact);
		while(listit_has_next(plist_it_factor)) {
		  pf=listit_next(plist_it_factor);
		  printf("%d\t %d\t %d\t %d\n",
					 pf->EST_start, pf->EST_end,
					 pf->GEN_start,pf->GEN_end);
		}
		listit_destroy(plist_it_factor);
	 }
  }
  listit_destroy(plist_it_factorizations);
}

void print_n_factorization_complete(pEST est,int n)
{
  pfactor pf;
  pfactorization pfact;

  plist factorizations=est->factorizations;
  pboollist polya_signals=est->polyA_signals;
  pboollist polyadenil_signals=est->polyadenil_signals;

  plistit plist_it_factor;
  plistit plist_it_factorizations;
  pboollistit pboollist_it_polya;
  pboollistit pboollist_it_polyadenil;

  plist_it_factorizations=list_first(factorizations);
  pboollist_it_polya=boollist_first(polya_signals);
  pboollist_it_polyadenil=boollist_first(polyadenil_signals);

  int cont=0;

  while(listit_has_next(plist_it_factorizations)){
	 my_assert(boollistit_has_next(pboollist_it_polya));
	 my_assert(boollistit_has_next(pboollist_it_polyadenil));

	 pfact=listit_next(plist_it_factorizations);
	 bool polya=(bool)boollistit_next(pboollist_it_polya);
	 bool polyadenil=(bool)boollistit_next(pboollist_it_polyadenil);

	 cont=cont+1;
	 if(cont==n){
		printf("#polya=%d\n#polyad=%d\n", (polya == true)?(1):(0), (polyadenil == true)?(1):(0));
		plist_it_factor=list_first(pfact);
		while(listit_has_next(plist_it_factor)) {
		  pf=listit_next(plist_it_factor);
		  printf("%d\t %d\t %d\t %d\n",
					 pf->EST_start, pf->EST_end,
					 pf->GEN_start,pf->GEN_end);
		}
		listit_destroy(plist_it_factor);
	 }
  }
  listit_destroy(plist_it_factorizations);
}

// See issue #7
static void
compute_coverage_and_exons(pfactorization pfact,
									size_t * coverage, size_t * n_exons) {
  *coverage= 0;
  *n_exons= 0;
  pfactor pf;
  plistit plist_it_factor;
  plist_it_factor= list_first(pfact);
  while(listit_has_next(plist_it_factor)) {
	 pf= listit_next(plist_it_factor);
	 *coverage= (*coverage) + (pf->EST_end + 1 - pf->EST_start);
	 *n_exons= *n_exons + 1;
  }
  listit_destroy(plist_it_factor);
}


// See issue #7
void print_factorizations_result(pbit_vect min_factors, plist p,
											plist list_of_unique_fact, psempl psemp)
{
  pfactorization pfact;
  plistit list_it_est, list_it_bin, list_it_fact;
  pEST est;
  pbit_vect bv;
  int count_est;

  if(min_factors!=NULL)inglobe(min_factors,psemp);

  list_it_est=list_first(p);

  count_est=0;
  while(listit_has_next(list_it_est)){
	 est=listit_next(list_it_est);

	 size_t best_factorization= 0;
	 size_t best_coverage= 0;
	 size_t best_n_exons= SIZE_MAX;
	 size_t current_factorization= 0;

	 list_it_bin= list_first(est->bin_factorizations);
	 list_it_fact= list_first(est->factorizations);

	 while (listit_has_next(list_it_bin)){
		my_assert(listit_has_next(list_it_fact));

		current_factorization= current_factorization+1;

		bv= listit_next(list_it_bin);
		pfact= listit_next(list_it_fact);

		if (contained(bv, psemp->factors_used)){
		  size_t current_coverage= 0;
		  size_t current_n_exons= SIZE_MAX;
		  compute_coverage_and_exons(pfact, &current_coverage, &current_n_exons);
		  if ((best_coverage < current_coverage) ||
				((best_coverage == current_coverage) &&
				 (best_n_exons > current_n_exons))) {
			 DEBUG("Found a better factorization. Currently: coverage %zunt, no. of exons %zu.",
					 current_coverage, current_n_exons);
			 best_coverage= current_coverage;
			 best_n_exons= current_n_exons;
			 best_factorization= current_factorization;
		  }
		}
	 }
	 listit_destroy(list_it_bin);
	 listit_destroy(list_it_fact);

// Print the "best" factorization of the current EST
	 INFO("Saving factorization %zu (coverage: %zunt, no. of exons: %zu) for EST '%s'",
			best_factorization, best_coverage, best_n_exons, est->info->EST_id);
	 printf(">%s\n",est->info->EST_id);
	 print_n_factorization_complete(est, best_factorization);
	 count_est++;
  }
  listit_destroy(list_it_est);
}


bool all_true(pbit_vect bv)
{
  for(unsigned int i=0;i<bv->n;i++){
	 if (BV_get(bv,i)==0) return false;
  }
  return true;
}

/* Semplificazione ulteriore della matrice colorata*/
//*******************************************************************

void color_matrix_semplified_destroy(plist p)
{
  plistit list_it;
  pEST est;
  list_it=list_first(p);

  while(listit_has_next(list_it)){
		est=listit_next(list_it);
		list_destroy(est->bin_factorizations,(delete_function)BV_destroy);
  }
 pfree(est);
//  pfree(p);
  listit_destroy(list_it);
}



void addFactoriz(pbit_vect bv,plist bin_fact,psempl psemp)
{
  int dim=(int)bv->n;

  pbit_vect bin=BV_create(abs(dim-count_true(psemp->factors_used)-count_true(psemp->factors_not_used)));

  int pos=0;

  for(int i=0;i<dim;i++){

	 if((!BV_get(psemp->factors_used,i))&&(!BV_get(psemp->factors_not_used,i))){
		BV_set(bin,pos,BV_get(bv,i));
		pos=pos+1;
	 }
  }
  list_add_to_tail(bin_fact,bin);
}

void addEST(pEST p,plist col_mat_semp,psempl psemp)
{

  plistit list_it;
  pEST est=EST_create();

  est->info=p->info;


  pbit_vect bv;

  plist bin_fact=list_create();

  list_it=list_first(p->bin_factorizations);

  while(listit_has_next(list_it)){

	 bv=listit_next(list_it);
	 addFactoriz(bv,bin_fact,psemp);
  }

  est->bin_factorizations=bin_fact;

  list_add_to_tail(col_mat_semp,est);
  listit_destroy(list_it);
}

plist color_matrix_semplified_create(plist ests_factorizations,psempl psemp)
{

  my_assert(ests_factorizations!=NULL);

  pEST p;
  plistit plist_it_id;
  int cont=0;

  plist col_mat_semp=list_create();
  plist_it_id=list_first(ests_factorizations);

  while(listit_has_next(plist_it_id)){
	 p= listit_next(plist_it_id);

	 if(!BV_get(psemp->ests_ok,cont)){addEST(p,col_mat_semp,psemp);}
	 cont=cont+1;
  }
  listit_destroy(plist_it_id);
  return col_mat_semp;
}

/*****************************************************************************
*******************************************************************************/

pbit_vect min_fact(plist color_matrix)
{
  my_assert(color_matrix!=NULL);
  int result;

  int k=max_of_min(color_matrix);
  int start=0;

  bool factorized=false;

  pEST e=list_head(color_matrix);
  pbit_vect bv=list_head(e->bin_factorizations);

  pbit_vect test=BV_create(bv->n);
  start=k;

  int size=(int)bv->n;
  while(factorized==false){

	 INFO("Combinazione di %d fattori",start);
	 factorized=create_combinations(0,start,test,color_matrix);
	 start=start+1;
  }
  INFO("*** ricerca del numero minimo di fattori terminata ***");
  return test;
}








