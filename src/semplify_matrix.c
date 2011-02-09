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
#include "semplify_matrix.h"
#include "list.h"
#include "bit_vector.h"
#include "log.h"
#include "stdlib.h"

//Semplificazione delle colonne di tutti 1 all'interno
//della matrice colorata.

bool semplify_column(plist bin_fact,pbit_vect used)
{
  bool elimination=false;
  bool semplify=true;
  pbit_vect vect_fact;
  plistit list_it_fact;

  for(unsigned int i=0;i<used->n;i++){

	 semplify=true;
	 list_it_fact=list_first(bin_fact);

	 while((listit_has_next(list_it_fact))&&(semplify==true)&&(BV_get(used,i)==false)){

		vect_fact=listit_next(list_it_fact);

		if(BV_get(vect_fact,i)==false){
		  semplify=false;
		}
	 }
	 if((semplify==true)&&(BV_get(used,i)==false)) {
		//fprintf(stdout, "Used %d\n", i);
		BV_set(used,i,true);
		elimination=true;
	 }
	 listit_destroy(list_it_fact);
  }
  return elimination;
}

//semplificazione delle righe avente tutti 0
//nella matrice colorata.

bool semplify_row(plist bin_fact,pbit_vect ests_ok,pbit_vect factors_used,int number_est)
{
  bool elimination=false;
  bool all_zero=true;
  plistit list_it_fact;
  pbit_vect vect_fact;
  list_it_fact=list_first(bin_fact);

  while(listit_has_next(list_it_fact)){
	 vect_fact=listit_next(list_it_fact);

	 for(unsigned int i=0;i<vect_fact->n;i++){

		 //fprintf(stdout, "%d(%d,%d)-", i, BV_get(vect_fact,i), BV_get(factors_used,i));

		if((BV_get(vect_fact,i)==true)&&(BV_get(factors_used,i)!=true)) all_zero=false;
	 }

	 //fprintf(stdout, "All zero %d\n\n", all_zero);

	 if((all_zero==true)&&(BV_get(ests_ok,number_est)==false)){
		BV_set(ests_ok,number_est,true);
		elimination=true;
	 }
	 all_zero=true;
  }
  listit_destroy(list_it_fact);

  return elimination;
}


//semplificazione delle colonne di tutti 0
//nella matrice colorata

bool semplify_column_zero(plist bin_factorizations,int column)
{
  plistit list_it_fact;
  pbit_vect bv;


  list_it_fact=list_first(bin_factorizations);

  while(listit_has_next(list_it_fact)){
	 bv=listit_next(list_it_fact);

	 if(BV_get(bv,column)==true){
		listit_destroy(list_it_fact);
		return false;
	 }
  }
  listit_destroy(list_it_fact);

  return true;
}

//Funzione che iterativamente applica il processo di semplificazione della matrice
//colorata. Non avviene una vera e propria eliminazione in quanto i risultati del
//processo vengono memorizzati nella struttura _sempl_info.



psempl semplification(plist color_matrix, plist unique_factors)
{
  psempl p=psempl_create();

  p->factors_used=BV_create(list_size(unique_factors));
  p->factors_not_used=BV_create(list_size(unique_factors));
  p->ests_ok=BV_create(list_size(color_matrix));

  pEST est;
  bool el_column=false;
  bool el_row=false;
  bool el_col_zero=false;
  bool all_zero=true;

  plistit list_it_est;

  INFO("Inizio dell fase di semplificazione");
  do{

  list_it_est=list_first(color_matrix);

  el_column=false;

  while(listit_has_next(list_it_est)){
	 est=listit_next(list_it_est);

	 /*fprintf(stdout, "EST %s\n", est->info->EST_id);
	 plistit temp=list_first(est->bin_factorizations);
	 while(listit_has_next(temp)){
		 pbit_vect ltemp=(pbit_vect)listit_next(temp);
		 unsigned int i;
		 for(i=0; i<ltemp->n; i++)
			 fprintf(stdout, "%d", BV_get(ltemp,i));
		 fprintf(stdout, "\n");
	 }
	 listit_destroy(temp);*/


	 el_column= semplify_column(est->bin_factorizations,p->factors_used);
  }
  listit_destroy(list_it_est);


  /*unsigned int j;
  fprintf(stdout, "Used: ");
  for(j=0; j<p->factors_used->n; j++){
	  if(BV_get(p->factors_used,j))
		  fprintf(stdout, "%d-", j);
  }
  fprintf(stdout, "\n");
  exit(1);*/


  list_it_est=list_first(color_matrix);
  int number_est=0;

  el_row=false;

  while(listit_has_next(list_it_est)){
	 est=listit_next(list_it_est);

	 //fprintf(stdout, "EST %s\n", est->info->EST_id);

	 el_row=semplify_row(est->bin_factorizations,p->ests_ok,p->factors_used,number_est);
	 number_est=number_est+1;
  }

  listit_destroy(list_it_est);

  el_col_zero=false;

  for(size_t column=0;column<list_size(unique_factors);column++){

	 number_est=0;
	 list_it_est=list_first(color_matrix);
	 all_zero=true;

	 while(listit_has_next(list_it_est)&&(all_zero==true)){
		est=listit_next(list_it_est);

		if((BV_get(p->ests_ok,number_est)==false)&&(BV_get(p->factors_used,column)==false)){
		  all_zero= semplify_column_zero(est->bin_factorizations,column);
		}
		number_est=number_est+1;
		}

	 listit_destroy(list_it_est);

	 if((all_zero==true)&&(BV_get(p->factors_used,column)==false)&&(BV_get(p->factors_not_used,column)==false)){
		BV_set(p->factors_not_used,column,true);
		el_col_zero=true;
	 }
  }


 }while((el_column==true)||(el_row==true)||(el_col_zero==true));

  INFO("fine semplificazione");
  psempl_print(p);
  return p;
}
