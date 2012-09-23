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
#include "io-factorizations.h"
#include "list.h"
#include <string.h>
#include "util.h"
#include "log.h"
#include <ctype.h>

#define BUFFER 100000
#define mult 1

//La funzione print_factorization riceve una lista  rappresentante una
//fattorizzazione di un particolare est. La funzione stampa l'id dell'EST
//e le coordinate genomiche di ogni fattore viste come quadruple di numeri interi
// (est start,est end, gen start, gen end)
static void print_factorization(pfactorization pfact,char* EST_id,FILE* dest)
{
  my_assert(pfact!=NULL && EST_id!=NULL && dest!=NULL);
  pfactor pf;
  plistit plist_it_factor;
  plist_it_factor=list_first(pfact);
  fprintf(dest,">%s\n",EST_id);
  while(listit_has_next(plist_it_factor)) {
	 pf=listit_next(plist_it_factor);
	 fprintf(dest,"%d\t %d\t %d\t %d\n",
				pf->EST_start, pf->EST_end,
				pf->GEN_start,pf->GEN_end);
  }
  listit_destroy(plist_it_factor);
}

//La seguente funzione riceve in ingresso la lista degli ests. Ogni est ha
//associata una lista delle fattorizzazioni. Tale funzione si preoccupa di stampare
//su file L'id dell'est considerato e ogni sua fattorizzazione richiamando
//la funzione print_factorization
void  write_factorizations(plist ests_factorizations, FILE* dest)
{
  my_assert(ests_factorizations!=NULL);

  pEST p;
  plistit plist_it_id, plist_it_f;
  pfactorization pfact;
  plist_it_id=list_first(ests_factorizations);
  while(listit_has_next(plist_it_id)){
	 p= listit_next(plist_it_id);
	 plist_it_f=list_first(p->factorizations);
	 while(listit_has_next(plist_it_f)){
		pfact=listit_next(plist_it_f);
		print_factorization(pfact,p->info->EST_id,dest);
	 }
	 listit_destroy(plist_it_f);
  }
  listit_destroy(plist_it_id);
  DEBUG("Scrittura del file eseguita correttamente!!");
}



//funzione che aggiunge un fattore alla lista che descrive una
//singola fattorizzazione
static void addFactor(int EST_start, int EST_end,
							 int GEN_start, int GEN_end,
							 pfactorization pfact)
{
  pfactor fact=factor_create();
  fact->EST_start=EST_start;
  fact->EST_end=EST_end;
  fact->GEN_start=GEN_start;
  fact->GEN_end=GEN_end;
  list_add_to_tail(pfact,fact);
}

//aggiunge una fattorizzazione all'EST passato in ingresso alla funzione.
//Viene inizializzata la lista delle fattorizzazioni di quell'EST la prima volta
//che si aggiunge una fattorizzazione di quell'EST.
static void addFactorization(FILE *fp,pEST est,char* my_string)
{
  my_assert(my_string!=NULL && est!=NULL && fp!=NULL);

  int bytes_read;
  size_t n_bytes=BUFFER;
  int n;

  int EST_start,EST_end,GEN_start,GEN_end;
  int polya=0, polyadenil=0;

  if(est->factorizations==NULL){
	 est->factorizations=list_create();
  }
  if(est->polyA_signals==NULL){
	 est->polyA_signals=boollist_create();
  }
  if(est->polyadenil_signals==NULL){
	 est->polyadenil_signals=boollist_create();
  }

  pfactorization pfact=list_create();

  bytes_read= custom_getline(&my_string, &n_bytes, fp);
  n= sscanf(my_string,"%d %d %d %d", &EST_start, &EST_end,
				&GEN_start, &GEN_end);
   while((my_string[0]!='>')&&(!(feof(fp)))){
	 if(my_string[0]!= '#'){
		 n=sscanf(my_string,"%d %d %d %d",&EST_start,&EST_end,&GEN_start,&GEN_end);
		 if((isdigit(my_string[0]))&&(n==4)) {

			EST_start=(EST_start/mult)*mult;
			EST_end=(EST_end/mult)*mult;
			GEN_start=(GEN_start/mult)*mult;
			GEN_end=(GEN_end/mult)*mult;

			if(EST_start==0) EST_start=1;
			if(EST_end==0)   EST_end=1;
			if(EST_start==0) GEN_start=1;
			if(EST_start==0) GEN_end=1;

			addFactor(EST_start,EST_end,GEN_start,GEN_end,pfact);
		 }
	 }
	 else{
		 int polya_s, polyadenil_s;
		 n=sscanf(my_string,"#polya=%d",&polya_s);
		 if(n==1) {
			 polya=polya_s;
		 }
		 else{
			 n=sscanf(my_string,"#polyad=%d",&polyadenil_s);
			 if(n==1) {
				 polyadenil=polyadenil_s;
			 }
		 }
	 }
	 bytes_read= custom_getline(&my_string,&n_bytes,fp);
 }
  list_add_to_tail(est->factorizations,pfact);
  boollist_add_to_tail(est->polyA_signals, (BTYPE) (polya == 1)?(true):(false));
  boollist_add_to_tail(est->polyadenil_signals, (BTYPE) (polyadenil == 1)?(true):(false));
}

//funzione che aggiunge un nuovo EST alla lista. Tale funzione richiama la
//addFactorization in modo da aggiungere le fattorizzazioni dell'EST

void set_EST_id(FILE* fp,char* my_string,plist est_factorizations)
{
  my_assert(my_string!=NULL && est_factorizations!=NULL && fp!=NULL);

  pEST est;
  pEST_info iest;
  char* substr;

  my_string[strlen(my_string)-1]='\0';
  substr=substring(1,my_string);

  iest=EST_info_create();
  est=EST_create();
  iest->EST_id=substr;
  est->info=iest;
  list_add_to_tail(est_factorizations,est);
  addFactorization(fp,est,my_string);
 }

//la funzione ReadFile() legge il file delle fattorizzazioni e crea la struttura
//a liste concatenate che descrive tutte le possibili fattorizzazioni di ogni EST
//descritto all'interno del file letto.

plist read_factorizations(FILE* fp) {
  ssize_t bytes_read;
  size_t n_bytes=BUFFER;
  char* my_string;
  char* current_ID="/";
  plist est_factorizations;
  pEST p;
  char* substr;

  int numero_fattorizz;  //no commit!

  est_factorizations=list_create();
  my_string = c_palloc(n_bytes+1);
  bytes_read= custom_getline(&my_string,&n_bytes,fp);

  while(!(feof(fp))){
	 if((bytes_read>0)&&(my_string[0]=='>')){

		substr=substring(1,my_string);
		substr[strlen(substr)-1]='\0';

		if((strcmp(substr,current_ID)!=0)){

		  numero_fattorizz=1; //no commit!!

		  pfree(substr);
		  set_EST_id(fp,my_string,est_factorizations);
		  p= list_tail(est_factorizations);
		  current_ID=p->info->EST_id;
		}else{

		  numero_fattorizz=numero_fattorizz+1;  //no commit!
		  pfree(substr);
		  addFactorization(fp,list_tail(est_factorizations),my_string);  //no commit!
		}
	 } else {
		bytes_read= custom_getline(&my_string,&n_bytes,fp);
	 }
  }
  pfree(my_string);
  return est_factorizations;
}

#undef BUFFER







