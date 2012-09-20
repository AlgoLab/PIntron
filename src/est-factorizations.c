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
#include <stdlib.h>
#include <math.h>
#include "est-factorizations.h"
#include "list.h"
#include "types.h"
#include "double_list.h"
#include "int_list.h"
#include "bool_list.h"
#include "refine.h"
#include "refine-intron.h"
#include "exon-complexity.h"
#include "compute-alignments.h"
#include "detect-polya.h"

#include "log.h"
#include "max-emb-graph.h"
//#define LOG_THRESHOLD LOG_LEVEL_TRACE

// Define a long string composed only by spaces (1024, in particular) to be used to align the recursive calls of getsubtreeembeddings
#define _A_LOT_OF_SPACES_ "                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "

//Computa per un dato subtree
//di un grafo degli embedding (GEM) gli embeddings e restituisce una lista di ppairing
static plist get_subtree_embeddings(const int counter, ppairing root, pconfiguration config,
												pmytime_timeout ptt, const char* const GEN_seq);

//Prende in input una lista di embeddings (lista di liste di ppairing) e fornisce in output una
//lista di fattorizzazioni (eliminando eventualmente embedding non buoni)
static plist get_factorizations_from_embeddings(plist, pconfiguration, pEST_info, int);

//Aggiorna l'embedding con il nuovo nodo (solo se e' compatibile)
static plist update_embedding(plist embedding, ppairing node,
										const char* const GEN_seq, pconfiguration config);

static pfactor create_and_set_factor(int donor_EST_start, int donor_EST_end, int donor_GEN_start, int donor_GEN_end);

///Create a new subtree embeddings object and set its parameter
static psubtree_embeddings create_and_set_subtree_embeddings(ppairing, plist);
static plist get_computed_subtree_embeddings(ppairing);

//Funzione per copiare liste di pfactor
//static pfactor copy_pfactor(pfactor);

//Funzione per copiare liste di pfactor
static ppairing copy_ppairing(ppairing);

//Funzione di stampa di una lista di fattorizzazioni
static void print_factorizations(plist);
//Funzione di stampa di una fattorizzazione (lista di fattori pfactor)
static void print_factorization(plist);

//Funzione di stampa di una lista di embeddings
static void print_embeddings(plist);
//Funzione di stampa di una embedding (lista di fattori ppairing)
static void print_embedding(plist);

//L'argomento deve essere una lista di fattorizzazioni (lista di liste di fattori)
//static void factorization_list_destroy(plist);

//L'argomento deve essere una lista di embeddings (lista di liste di pairing)
static void embedding_list_destroy(plist);

//Funzione di stampa di una lista di fattorizzazioni
static void print_factorizations_on_log(const int log_level, plist, const char* const);
//Funzione di stampa di una fattorizzazione (lista di fattori pfactor)
static void print_factorization_on_log(const int log_level, plist, const char* const);

//Funzione di stampa di una lista di embeddings
//static void print_embeddings_for_info(plist);
//Funzione di stampa di una fattorizzazione (lista di fattori ppairing)
//static void print_embedding_for_info(plist);

//Controlla che la fattorizzazione passata come input non tagli un suffisso/prefisso eccessivo
//factorization deve essere una lista di pfactor
//static bool test_prefix_suffix(plist, int, pconfiguration, pEST_info);

//Comparator per confrontare due fattori (pfactor)
//static int factor_compare(const pfactor*, const pfactor*);

//Comparator per confrontare due pairing (ppairing)
//static int perfect_pairing_compare(const ppairing* pp1, const ppairing* pp2);

//Comparator per confrontare liste di fattori (pfactor)
//Due fattori sono uguali se hanno gli estremi (su EST e GEN) che differiscono di pochi nt)
//Il parametro type permette di specificare se il confronto e' su entrambi gli estremi (type=0),
//solo sui left end (-1) o solo sui right end (-2)
static int relaxed_factor_compare(const pfactor*, const pfactor*, int, int, plist);

//Calcola la copertura (rispetto alla lunghezza) sulla EST di una fattorizzazione (valore tra 0.0 e 1.0)
static double compute_coverage(plist, unsigned int);

//Calcola la lunghezza totale dei gap (su P) di una fattorizzazione (tagli prefisso/suffisso esclusi)
static int compute_gapLength(plist);

//Controlla la massimalita' tra due embedding
//Ritorna 2 se il primo e' massimale, 1 se entrambi sono massimali e 0 se il secondo e' massimale.
static char maximality_relation(plist, plist);

static bool check_gap_errors(plist, char*, char *, pconfiguration);

plist list_of_subtree_embeddings;

//Computa per una data EST tutte le fattorizzazioni ammissibili a partire
//dal grafo degli embedding massimali (pext_array)
pEST get_EST_factorizations(pEST_info pest_info, pext_array pext, pconfiguration config,
									 pEST_info gen_info, pmytime_timeout ptt)
{
  unsigned int i, pext_size;
  pEST est;
  plist pgem, factorization_list;
  plist subtree_fact_list;
  plist subtree_embedding_list;
  plistit pgem_iter;
  ppairing next_pairing;

  unsigned int counter;

  my_assert(pest_info != NULL);
  my_assert(pext != NULL);

  est = EST_create();
  est->info = pest_info;
  INFO("Getting factorizations of EST %s", est->info->EST_id);

  pext_size=EA_size(pext);
  unsigned int EST_length=pext_size-2;
  INFO("\tEA dimension: %d, EST length %d", pext_size, EST_length);

//Creazione della lista delle fattorizzazioni ammissibili vuota
  factorization_list=list_create();

//Settaggio dei campi visited e number_of_visits nei pairing del MEG
  for(i=0; i<pext_size; i++){
//Puntatore alla lista di pairings in posizione i
	 pgem =(plist)EA_get(pext, i);
	 pgem_iter=list_first(pgem);
	 while(listit_has_next(pgem_iter)){
		next_pairing=(ppairing) listit_next(pgem_iter);
		next_pairing->number_of_visits=0;
		next_pairing->visited=false;
	 }
	 listit_destroy(pgem_iter);
  }

  list_of_subtree_embeddings=list_create();

  for(i=0; i<pext_size; i++){
//Puntatore alla lista di pairings in posizione i
	 pgem =(plist)EA_get(pext, i);

	 if(!list_is_empty(pgem))
		DEBUG("\tPairings in position %d of EST", i);

//Iteratore su pgem
	 pgem_iter=list_first(pgem);
//Scansione dei pairings della lista
	 while(listit_has_next(pgem_iter)){
		next_pairing=(ppairing) listit_next(pgem_iter);

//if next_pairing has not been yet visited as a root o a factorization substree
		if(!next_pairing->visited){

		  counter=1;

		  DEBUG("\t\t%d) Path rooted in pairing (%d, %d, %d)",
				  counter, next_pairing->p, next_pairing->t, next_pairing->l);

		  subtree_embedding_list= get_subtree_embeddings(counter, next_pairing, config, ptt, gen_info->EST_seq);

		  if (subtree_embedding_list == NULL) {
			 return NULL;
		  }

		  DEBUG("\t\t...ALL THE EMBEDDINGS FOR THE PATH ARE OBTAINED!");
		  print_embeddings(subtree_embedding_list);

		  subtree_fact_list=get_factorizations_from_embeddings(subtree_embedding_list, config, est->info, pext_size-2);

// DISTRUGGERE subtree_embedding_list QUI

		  //printf("...ALL THE FACTORIZATIONS FOR THE PATH ARE OBTAINED %zu!\n", list_size(subtree_fact_list));

		  DEBUG("\t\t...ALL THE FACTORIZATIONS FOR THE PATH ARE OBTAINED!");
		  print_factorizations(subtree_fact_list);

		  plistit add_it=list_first(subtree_fact_list);
		  int counter_add=1;
		  while(listit_has_next(add_it)){
			 plist add_f=(plist)listit_next(add_it);

			 bool is_ok=check_for_not_source_sink_factorization(add_f, (int)EST_length);

			 if(is_ok)
				 is_ok=check_exon_start_end(add_f);

			 if(is_ok){
				 add_f=handle_endpoints(add_f, gen_info->EST_seq, est->info->EST_seq);
				 if(list_is_empty(add_f))
					 is_ok=false;
			 }

			 if(is_ok){
				 add_f=clean_external_exons(add_f, gen_info->EST_seq, est->info->EST_seq);
				 if(list_is_empty(add_f))
					 is_ok=false;
			 }

			 if(is_ok){
				 add_f=clean_low_complexity_exons_2(add_f, gen_info->EST_seq, est->info->EST_seq);
				 if(list_is_empty(add_f))
					 is_ok=false;
			 }

			 if(is_ok){
				 add_f=clean_noisy_exons(add_f, gen_info->EST_seq, est->info->EST_seq, false);
				 if(list_is_empty(add_f))
					 is_ok=false;
			 }

			 if(is_ok){
				 is_ok=check_est_coverage(add_f, est->info->EST_seq);
			 }

			 if(is_ok){
				 bool check_adding;
				 factorization_list=add_if_not_exists(add_f, factorization_list, config, &check_adding);
				 if(check_adding == false){
				 	 list_remove_at_iterator(add_it, (delete_function) factorization_destroy);
				 }
			 }
			 else
				 list_remove_at_iterator(add_it, (delete_function) factorization_destroy);

			 counter_add++;
		  }
		  listit_destroy(add_it);
		  list_destroy(subtree_fact_list, (delete_function)noop_free);
		}
		else{
		  DEBUG("\t...already visited! Cannot be a source of a path!");
		}
	 }
	 listit_destroy(pgem_iter);
  }

  plistit plist_add_factorization;

   /**Calcolo della coperture su P (non si tiene conto di gap su P ma solo del prefisso/suffisso tagliati)
 * e di quella massima**/
  int count_f=1;
  double max_coverage=0.0;
  pdoublelist coverage_list=doublelist_create();
  plistit it_for_cover=list_first(factorization_list);

  while(listit_has_next(it_for_cover)){
	 plist add_factorization=(plist)listit_next(it_for_cover);
	 DEBUG("\t\tThe factorization number %d", count_f);
	 bool isSourceSink=false;
	 if(list_size(add_factorization) == 1){
		pfactor head=(pfactor)list_head(add_factorization);
		if(head->EST_start < 0 || head->EST_start >= (int)EST_length){
		  DEBUG("\t\t\t...is a source-sink factorization and will be discarded!");
		  double coverage=-1.0;
		  doublelist_add_to_tail(coverage_list, coverage);
		  isSourceSink=true;
		}
	 }
	 if(isSourceSink == false){
		double coverage=compute_coverage(add_factorization, EST_length);
		doublelist_add_to_tail(coverage_list, coverage);
		if(max_coverage < coverage)
		  max_coverage=coverage;
		DEBUG("\t\t\t...has a coverage of %f!!", coverage);
	 }
	 count_f++;
  }
  listit_destroy(it_for_cover);

//FILTRO 1
//Prendo solo le fattorizzazioni che hanno una copertura che si avvicina a quella max
//entro una prefissata soglia. Vengono eliminate anche le source-sink
  count_f=1;
  plist_add_factorization=list_first(factorization_list);
  pdoublelistit d_it=doublelist_first(coverage_list);
  while(listit_has_next(plist_add_factorization)){
	 plist add_factorization=(plist)listit_next(plist_add_factorization);
	 double coverage=(double)doublelistit_next(d_it);
	 DEBUG("\t\tThe factorization number %d", count_f);
	 print_factorization(add_factorization);

	 if(coverage == -1.0 || max_coverage-coverage > config->max_coverage_diff){
		list_remove_at_iterator(plist_add_factorization, (delete_function)factorization_destroy);
		DEBUG("\t\t\t...has been discarded!");
	 }
	 else{
		 //XXX
		 if((max_coverage-coverage)*((int)strlen(est->info->EST_seq)) > 100){
			list_remove_at_iterator(plist_add_factorization, (delete_function)factorization_destroy);
			DEBUG("\t\t\t...has been discarded! The abs coverage difference is too much!");
		 }
	 }

	 count_f++;
 }
 doublelistit_destroy(d_it);
 doublelist_destroy(coverage_list);
 listit_destroy(plist_add_factorization);

/**Calcolo del numero minimo di esoni sulle fattorizzazioni rimaste dopo FILTRO 1**/
  /*count_f=1;
  int min_exonNUM=0;
  pintlist exonNUM_list=intlist_create();
  plistit it_for_exonNUM=list_first(factorization_list);
  while(listit_has_next(it_for_exonNUM)){
	 plist add_factorization=(plist)listit_next(it_for_exonNUM);
	 DEBUG("\t\tThe factorization number %d", count_f);
	 int exonNUM=list_size(add_factorization);
	 intlist_add_to_tail(exonNUM_list, exonNUM);
	 if(min_exonNUM == 0 || min_exonNUM > exonNUM)
		min_exonNUM=exonNUM;
	 DEBUG("\t\t\t...has %d exon(s)!", exonNUM);
	 count_f++;
  }
  listit_destroy(it_for_exonNUM);

//FILTRO 2
//Prendo solo le fattorizzazioni che hanno un numero di esoni che si avvicina a quello minimo
//entro una soglia prestabilita
  count_f=1;
  plist_add_factorization=list_first(factorization_list);
  pintlistit i_it=intlist_first(exonNUM_list);
  while(listit_has_next(plist_add_factorization)){
	 plist add_factorization=(plist)listit_next(plist_add_factorization);
	 int exonNUM=(int)intlistit_next(i_it);
	 DEBUG("\t\tThe factorization number %d", count_f);
	 print_factorization(add_factorization);

	 if ((config->max_exonNUM_diff != -1) && (exonNUM-min_exonNUM > config->max_exonNUM_diff)){
		list_remove_at_iterator(plist_add_factorization, (delete_function)factorization_destroy);
		print_factorizations(factorization_list);
		DEBUG("\t\t\t...has been discarded!");
	 }
	 count_f++;
  }
  intlistit_destroy(i_it);
  intlist_destroy(exonNUM_list);
  listit_destroy(plist_add_factorization);*/

 /**Calcolo della minima lunghezza totale di gap su P (esclusi i tagli di prefisso/suffisso) per le
 * fattorizzazioni rimaste dopo FILTRO 2. ATTENZIONE: i gap su P rimasti sono solo quelli in corrispondenza
 * di un introne su T**/
  count_f=1;
  int min_gapLength=-1;
  pintlist gapLength_list=intlist_create();
  plistit it_for_gapLength=list_first(factorization_list);

  while(listit_has_next(it_for_gapLength)){
	 plist add_factorization=(plist)listit_next(it_for_gapLength);
	 DEBUG("\t\tThe factorization number %d", count_f);
	 int gapLength=compute_gapLength(add_factorization);
	 intlist_add_to_tail(gapLength_list, gapLength);
	 DEBUG("\t\t\t...has a total gap length on P of %d nt (%zd exons)!", gapLength, list_size(add_factorization));
	 if(min_gapLength == -1 || min_gapLength > gapLength)
		min_gapLength=gapLength;
	 count_f++;
  }
  listit_destroy(it_for_gapLength);

  //FILTRO 3
  //Prendo solo le fattorizzazioni che hanno una lunghezza totale di gap su P che si avvicina a quella minima
  //entro una soglia prestabilita
   count_f=1;
   plist_add_factorization=list_first(factorization_list);
   pintlistit i2_it=intlist_first(gapLength_list);
   while(listit_has_next(plist_add_factorization)){
 	 plist add_factorization=(plist)listit_next(plist_add_factorization);
 	 int gapLength=(int)intlistit_next(i2_it);
 	 DEBUG("\t\tThe factorization number %d", count_f);
 	 print_factorization(add_factorization);

 	 if ((config->max_gapLength_diff != -1) && (gapLength-min_gapLength > config->max_gapLength_diff)) {
 		list_remove_at_iterator(plist_add_factorization, (delete_function)factorization_destroy);
		print_factorizations(factorization_list);
 		DEBUG("\t\t\t...has been discarded!");
 	 }
 	 count_f++;
   }
   intlistit_destroy(i2_it);
   intlist_destroy(gapLength_list);
   listit_destroy(plist_add_factorization);

  //FILTRO 4
  //I gap vengono controllati e le fattorizzazioni, con anche solo un gap che supera una certa soglia e che
  //non puo' essere riempito (oppure quelle che hanno un livello totale di errore sui gap troppo elevato),
  //vengono eliminate
   count_f=1;
   plist_add_factorization=list_first(factorization_list);
     while(listit_has_next(plist_add_factorization)){
  	 plist add_factorization=(plist)listit_next(plist_add_factorization);
  	 DEBUG("\t\tThe factorization number %d", count_f);
  	 print_factorization(add_factorization);
  	 if(!check_gap_errors(add_factorization, est->info->EST_seq, gen_info->EST_seq, config)){
  		list_remove_at_iterator(plist_add_factorization, (delete_function)factorization_destroy);
 		print_factorizations(factorization_list);
 		DEBUG("\t\t\t...has been discarded!");
  	 }
  	 count_f++;
    }
    listit_destroy(plist_add_factorization);

  DEBUG("Preliminary factorizations of EST %s => %zd", est->info->EST_id, list_size(factorization_list));
  print_factorizations_on_log(LOG_LEVEL_DEBUG, factorization_list, gen_info->EST_seq);

  if(config->max_number_of_factorizations != 0 && ((int) list_size(factorization_list)) > config->max_number_of_factorizations){
	 INFO("\tbut it is an artifact! There are too many factorizations!");
	 factorization_list_destroy(factorization_list);
	 factorization_list=list_create();
  }

  //Refine factorizations
  DEBUG("Refine factorizations...");
  plistit plist_fact_to_be_refined=list_first(factorization_list);
  while(listit_has_next(plist_fact_to_be_refined)){
	 plist fact_to_be_refined=(plist)listit_next(plist_fact_to_be_refined);
	 if(!list_is_empty(fact_to_be_refined)){
		 DEBUG("New factorization to be refined:");
		 print_factorization_on_log(LOG_LEVEL_DEBUG, fact_to_be_refined, gen_info->EST_seq);
		 plistit plist_f_t_r=list_first(fact_to_be_refined);
		 pfactor donor=(pfactor)listit_next(plist_f_t_r);
		 bool first_intron=true;
		 while(listit_has_next(plist_f_t_r)){
			 pfactor acceptor=(pfactor)listit_next(plist_f_t_r);

			 DEBUG(" - Refining intron:  gen. coord.=[%7d-%7d],  pattern=[ %.2s - %.2s],  EST cut=%6d",
					 donor->GEN_end+1, acceptor->GEN_start-1,
					 gen_info->EST_seq+donor->GEN_end+1, gen_info->EST_seq+acceptor->GEN_start-2,
					 acceptor->EST_start);

			 refine_intron(config, gen_info, pest_info, donor, acceptor, first_intron);

			 first_intron=false;

			 DEBUG("   Resulting intron: gen. coord.=[%7d-%7d],  pattern=[ %.2s - %.2s],  EST cut=%6d",
					 donor->GEN_end+1, acceptor->GEN_start-1,
					 gen_info->EST_seq+donor->GEN_end+1, gen_info->EST_seq+acceptor->GEN_start-2,
					 acceptor->EST_start);

			 donor=acceptor;
		 }
		 listit_destroy(plist_f_t_r);

		 plist_f_t_r=list_first(fact_to_be_refined);
		 pfactor first_exon=(pfactor)listit_next(plist_f_t_r);
		 if(listit_has_next(plist_f_t_r)){
			 pfactor second_exon=(pfactor)listit_next(plist_f_t_r);
			 if(first_exon->EST_start == second_exon->EST_start){
				 DEBUG("\tThe first and the second exons in the list are the same");
				 DEBUG("\t...and the first one will be removed");
				 list_remove_from_head(fact_to_be_refined);
			 }
		 }
		 listit_destroy(plist_f_t_r);
		 DEBUG("Refined factorization:");
		 print_factorization_on_log(LOG_LEVEL_DEBUG, fact_to_be_refined, gen_info->EST_seq);
	 }
  }

  listit_destroy(plist_fact_to_be_refined);

  plist final_factorization_list=list_create();
  //plist final_factorization_list=factorization_list;

  plistit plist_fact_to_be_corrected=list_first(factorization_list);
  while(listit_has_next(plist_fact_to_be_corrected)){
	  plist fact_to_be_corrected=(plist)listit_next(plist_fact_to_be_corrected);
	  fact_to_be_corrected=clean_noisy_exons(fact_to_be_corrected, gen_info->EST_seq, est->info->EST_seq, false);
	  fact_to_be_corrected=clean_external_exons(fact_to_be_corrected, gen_info->EST_seq, est->info->EST_seq);
	  /*bool is_ok=false;
	  if(!list_is_empty(fact_to_be_corrected))
		  is_ok=check_exon_start_end(fact_to_be_corrected);*/
	  //if(is_ok == false || (list_is_empty(fact_to_be_corrected) || !check_small_exons(fact_to_be_corrected)))
	  if(list_is_empty(fact_to_be_corrected))// || !check_small_exons(fact_to_be_corrected))
			list_remove_at_iterator(plist_fact_to_be_corrected, (delete_function)factorization_destroy);
	  else{
			bool check_adding;
			final_factorization_list=add_if_not_exists(fact_to_be_corrected, final_factorization_list, config, &check_adding);
			if(check_adding == false){
				list_remove_at_iterator(plist_fact_to_be_corrected, (delete_function) factorization_destroy);
			}
	  }
  }

  listit_destroy(plist_fact_to_be_corrected);
  list_destroy(factorization_list, (delete_function)noop_free);

  /*STAMPA QUALITA' FATTORIZZAZIONI*****************************************/
    /*plistit print_it=list_first(final_factorization_list);
    while(listit_has_next(print_it)){
  		 plist print_fact=(plist)listit_next(print_it);
  		 plistit print_it2=list_first(print_fact);
  		 printf("Factorization:\n");
  		 unsigned int tot_edit=0;
  		 int count_f=1;
  		 while(listit_has_next(print_it2)){
  			 pfactor print_f=(pfactor)listit_next(print_it2);
  			 printf("%d %d %d %d ", print_f->EST_start+1, print_f->EST_end+1, print_f->GEN_start+1, print_f->GEN_end+1);
  			 if(count_f == 1)
  				 printf(". ");
  			 else
  				 printf("%c%c ", gen_info->EST_seq[print_f->GEN_start-2], gen_info->EST_seq[print_f->GEN_start-1]);

  			 if(count_f == (int)list_size(print_fact))
  				 printf(". ");
  			 else
  				 printf("%c%c ", gen_info->EST_seq[print_f->GEN_end+1], gen_info->EST_seq[print_f->GEN_end+2]);

  			 int exon_length=print_f->GEN_end-print_f->GEN_start+1;
  				//XXX
  			int max_allowed_error=(int)((double)exon_length*(3.0/100.0));

  			char *gen_exon_seq=real_substring(print_f->GEN_start, print_f->GEN_end-print_f->GEN_start+1, gen_info->EST_seq);
  			char *est_exon_seq=real_substring(print_f->EST_start, print_f->EST_end-print_f->EST_start+1, est->info->EST_seq);

  			unsigned int edit;
  			bool ok=false;
  			ok=K_band_edit_distance(gen_exon_seq, est_exon_seq, max_allowed_error, &edit);
  			if(!ok)
  				printf("NOT-OK\n");
  			else
  				printf("%d\n", edit);

  			tot_edit+=edit;

  			pfree(gen_exon_seq);
  			pfree(est_exon_seq);

  			 count_f++;
  		 }
  		 listit_destroy(print_it2);

  		 printf("*******Error %d\n", tot_edit);
    }
    listit_destroy(print_it);*/
    /*STAMPA QUALITA' FATTORIZZAZIONI*****************************************/

  //Detect polyA signal
  plist_fact_to_be_corrected=list_first(final_factorization_list);
  est->polyA_signals=boollist_create();
  est->polyadenil_signals=boollist_create();
  while(listit_has_next(plist_fact_to_be_corrected)){
		  plist fact_to_be_corrected=(plist)listit_next(plist_fact_to_be_corrected);
		  DEBUG("Correcting tail exon");
		  fact_to_be_corrected=correct_composition_tail(fact_to_be_corrected, gen_info->EST_seq, est->info->original_EST_seq);
		  bool polyadenil=false;
		  bool polyA=detect_polyA_signal(fact_to_be_corrected, gen_info->EST_seq, est->info->original_EST_seq, &polyadenil);
		  boollist_add_to_tail(est->polyA_signals, polyA);
		  boollist_add_to_tail(est->polyadenil_signals, polyadenil);
 }
  listit_destroy(plist_fact_to_be_corrected);

  est->factorizations=final_factorization_list;

  INFO("Refined factorizations of EST %s => %zd",
		 est->info->EST_id, list_size(final_factorization_list));
  print_factorizations_on_log(LOG_LEVEL_INFO, final_factorization_list, gen_info->EST_seq);

  return est;
}

//Computa per un dato subtree tutti gli embedding e restituisce una lista di liste di ppairing
static plist get_subtree_embeddings(const int counter, ppairing root, pconfiguration config,
												pmytime_timeout ptt, const char* const GEN_seq)
{
  plist embedding_list; 	//Lista degli embedding
  plist updated_embedding_list;
  plist subtree_embedding_list; 	//Lista degli embedding relativi ad un subtree
  plist embedding; 	//Embedding: lista di ppairing
  plist next_embedding;
  plist adj_list;
  plistit adj_list_iter, subtree_embedding_iter;
  ppairing pair;
  ppairing next_adj_pairing;

  my_assert(root != NULL);

  DEBUG("\t\t%.*s->%d) Pairing node (%d, %d, %d)", counter, _A_LOT_OF_SPACES_, counter,
		  root->p, root->t, root->l);

  plist computed_sub_e=get_computed_subtree_embeddings(root);
  if(computed_sub_e != NULL){
	 DEBUG("\t\tThe subtree embeddings are retrieved!");

	 print_embeddings(computed_sub_e);

	 return computed_sub_e;
  }
  DEBUG("\t\t  %.*sThe subtree embeddings are to be computed!",
		  counter, _A_LOT_OF_SPACES_);

// Check timeout
  if (MYTIME_timeout_expired(ptt)) {
	 return NULL;
  }

//Recupero la lista di adiacenza del nodo in input
  adj_list=root->adjs;

//Creazione della lista delle fattorizzazioni vuota
  embedding_list=embedding_create();

  root->visited=true;
  root->number_of_visits++;

//Se il nodo root e' una foglia
  if(list_is_empty(adj_list)){
	 DEBUG("\t\t\t%.*s...IS A LEAF!", counter, _A_LOT_OF_SPACES_);

//Creazione dell'embedding (lista di ppairing vuota)
	 embedding=list_create();

//Creazione del pairing
	 pair=pairing_create();

	 pair->p=root->p;
	 pair->t=root->t;
	 pair->l=root->l;

//Aggiunta del pairing all'embedding
	 list_add_to_head(embedding, pair);

//Aggiunta dell'embedding alla lista degli embedding
	 list_add_to_head(embedding_list, embedding);
  }
  else{
	 adj_list_iter=list_first(adj_list);

//Scansione dei pairings della lista
	 while(listit_has_next(adj_list_iter)){
		next_adj_pairing=(ppairing) listit_next(adj_list_iter);

		subtree_embedding_list= get_subtree_embeddings(counter+1, next_adj_pairing, config, ptt, GEN_seq);
		if (subtree_embedding_list == NULL) {
		  return NULL;
		}

		DEBUG("\t\t\t%.*sEmbeddings of subtree rooted in (%d, %d, %d) obtained!",
				counter, _A_LOT_OF_SPACES_,
				next_adj_pairing->p, next_adj_pairing->t, next_adj_pairing->l);

		print_embeddings(subtree_embedding_list);

		//Aggiungo (se compatibile) il root node ad ogni embedding in subtree_embedding_list
		subtree_embedding_iter=list_first(subtree_embedding_list);

		DEBUG("\t\t\t%.*sAdding the node (%d, %d, %d) to the embeddings above...",
				counter, _A_LOT_OF_SPACES_, root->p, root->t, root->l);

		int count_f=1;

		unsigned int time_limit_check= 0;
		while(listit_has_next(subtree_embedding_iter)){

		  next_embedding=(plist)listit_next(subtree_embedding_iter);

		  DEBUG("\t\t\t\t%.*s...adding the node to the embedding %d...",
				  counter, _A_LOT_OF_SPACES_, count_f);
		  print_embedding(next_embedding);

//Adds the root node; ritorna una lista di embedding (attualmente uno solo)
		  updated_embedding_list=update_embedding(next_embedding, root, GEN_seq, config);

		  DEBUG("\t\t\t\t%.*s...node added!", counter, _A_LOT_OF_SPACES_);
		  print_embeddings(updated_embedding_list);

		  if(!list_is_empty(updated_embedding_list)){
			 plistit add_it=list_first(updated_embedding_list);
			 while(listit_has_next(add_it)){

// Check timeout (but not often)
				if (!time_limit_check && MYTIME_timeout_expired(ptt)) {
				  return NULL;
				}
				time_limit_check += 1;
				time_limit_check &= 1023u;

				plist add_emb=(plist)listit_next(add_it);

				plistit cmp_it=list_first(embedding_list);

//Valori di ritorno della funzione maximality_relation(add_emb, cmp_emb):
//2 ==> add_emb e' massimale
//1 ==> entrambe sono massimali
//0 ==> cmp_emb e' massimale
				char is_maximal=2;

				while(listit_has_next(cmp_it) && is_maximal >= 1){
				  plist cmp_emb=(plist)listit_next(cmp_it);
				  is_maximal=maximality_relation(add_emb, cmp_emb);

				  if(is_maximal == 2){
					 list_remove_at_iterator(cmp_it, (delete_function)pairing_destroy_2);
				  }
				}
				listit_destroy(cmp_it);

				if(is_maximal >= 1)
				  list_add_to_tail(embedding_list, add_emb);
				else{
				  list_remove_at_iterator(add_it, (delete_function)pairing_destroy_2);
				}
			 }

			 if(list_is_empty(updated_embedding_list)){
				embedding_list_destroy(updated_embedding_list);
			 }

			 listit_destroy(add_it);
		  }
		  else{
			 embedding_list_destroy(updated_embedding_list);
		  }
		  count_f++;
		}

		listit_destroy(subtree_embedding_iter);
	 }
	 listit_destroy(adj_list_iter);

	 DEBUG("\t\t\t%.*s...node (%d, %d, %d) added to all the embeddings of all its adjacent nodes!",
			 counter, _A_LOT_OF_SPACES_, root->p, root->t, root->l);
	 print_embeddings(embedding_list);
  }

  list_add_to_tail(list_of_subtree_embeddings, create_and_set_subtree_embeddings(root, embedding_list));

  return embedding_list;
}

//Aggiorna l'embedding con il nuovo nodo (solo se e' compatibile)
static plist update_embedding(plist embedding, ppairing node,
										const char* const GEN_seq, pconfiguration config){
  plist copy_embedding;
  plist sink_embedding;
  int node_copy_l;
  int head_copy_l, head_copy_p, head_copy_t;

  my_assert(!list_is_empty(embedding));

  ppairing head=(ppairing)list_head(embedding);
  plist return_embedding_list=list_create();

  if(head->p == SINK_PAIRING_START){
	 if(node->p >= 0){
		sink_embedding=list_create();
		ppairing node_copy=pairing_simple_copy(node);
		list_add_to_head(sink_embedding, node_copy);
		list_add_to_head(return_embedding_list, sink_embedding);
	 }
	 return return_embedding_list;
  }

  if(node->p < 0){
	 copy_embedding=list_copy(embedding, (copy_item) copy_ppairing);
	 list_add_to_head(return_embedding_list, copy_embedding);

	 return return_embedding_list;
  }

  TRACE("Adding node (%d, %d, %d) to the embeddings starting with head (%d, %d, %d)",
		  PAIRING(node), PAIRING(head));
//Verifico la compatibilita' di node con embedding, cioe' node e' il primo pairing (head) di embedding
//devono soddisfare gli stessi criteri di linking (node->head) con cui il MEG e' stato costruito
  const int small_delta=(head->p+head->l)-node->p;
  const int big_delta=(head->t+head->l)-node->t;
  TRACE("Deltas: small=%d, big=%d", small_delta, big_delta);

//L'embedding dei due nodi deve dare fattori di una certa lunghezza minima
  const int min_fl= (int)config->min_factor_len;
  const int fl=2*min_fl;
  if(small_delta >= fl && big_delta >= fl){
//L'eventuale gap su P tra node e head deve essere limitato
	 if(small_delta-(node->l+head->l) <= fl){
//L'embedding dei due nodi non deve aprire un gap su P maggiore di una certa soglia (su T e' libero)
//perche' ammettiamo gli introni e l'introne su T deve essere >= ad una certa soglia se questa e' specificata in config
		if(small_delta-big_delta <= fl){

//Esiste un gap su P e su T
		  if(small_delta >= (node->l+head->l) && big_delta >= (node->l+head->l)){
			 head_copy_p=head->p;
			 head_copy_t=head->t;
			 head_copy_l=head->l;
			 node_copy_l=node->l;
		  }
//Esiste un overlap o su P o su T o su entrambi
		  else{
			 TRACE("There is an overlap on P and/or T.");
			 const int ref_delta= MIN(small_delta, big_delta);

			 int temp_length_node=ref_delta/2;
			 int temp_length_head=ref_delta-temp_length_node;
			 if(temp_length_node > node->l){
				temp_length_node=node->l;
				temp_length_head=ref_delta-temp_length_node;
			 }
			 else{
				if(temp_length_head > head->l){
				  temp_length_head=head->l;
				  temp_length_node=ref_delta-temp_length_head;
				}
			 }
			 head_copy_l=temp_length_head;
			 head_copy_p=head->p+head->l-head_copy_l;
			 head_copy_t=head->t+head->l-head_copy_l;
			 node_copy_l=temp_length_node;
		  }

		  const bool is_overlap_on_p= ( small_delta < (node->l + head->l) );

		  const int gap_length_on_p= head_copy_p - node->p - node_copy_l - 1;
		  const int gap_length_on_t= head_copy_t - node->t - node_copy_l - 1;
		  const int possible_intron_length= gap_length_on_t - MAX(0, gap_length_on_p);
		  const bool is_intron_on_t= ( (possible_intron_length >= 0) &&
												 (config->min_intron_length == 0 ||
												  possible_intron_length >= config->min_intron_length) );

// Ho un introne e avevo un overlap su P, cerco il miglior introne secondo Burset
		  if (is_overlap_on_p && is_intron_on_t) {
			 DEBUG("There is an overlap on P and an intron on T. "
					 "Searching the best cut according Burset...");
			 int best_burset_freq= -1;
			 int best_P_cut= 0;
			 const int min_P_cut= MAX(node->p + min_fl, head->p);
			 const int max_P_cut= MIN(head->p + head->l - min_fl, node->p + node->l);
			 TRACE("Trying to cut between %d and %d (inclusive)...", min_P_cut, max_P_cut);
			 for (int current_P_cut= min_P_cut;
					current_P_cut <= max_P_cut;
					++current_P_cut) {
				TRACE("Trying P cut %d...", current_P_cut);
				const int current_burset_freq=
				  getBursetFrequency_adaptor(GEN_seq,
													  current_P_cut - node->p + node->t,
													  current_P_cut - head->p + head->t);
				if (current_burset_freq >= best_burset_freq) {
				  best_burset_freq= current_burset_freq;
				  best_P_cut= current_P_cut;
				}
			 }
			 DEBUG("Best intron placement (according to Burset frequency): "
					 "P cut=%d, Burset frequency=%d.",
					 best_P_cut, best_burset_freq);
			 DEBUG("Original pairings on P: [%9d --%10d)   -->   [%9d --%10d)",
					 node->p, node->p + node->l,
					 head->p, head->p + head->l);
			 DEBUG("Original pairings on T: [%9d --%10d)   -->   [%9d --%10d)",
					 node->t, node->t + node->l,
					 head->t, head->t + head->l);
			 const int tmpdeltaH= best_P_cut - head->p;
			 head_copy_l= head->l - tmpdeltaH;
			 head_copy_p= head->p + tmpdeltaH;
			 head_copy_t= head->t + tmpdeltaH;
			 const int tmpdeltaN= node->p + node->l - best_P_cut;
			 node_copy_l= node->l - tmpdeltaN;
			 DEBUG("Refined pairings on P:  [%9d --%10d)   -->   [%9d --%10d)",
					 node->p, node->p + node_copy_l,
					 head_copy_p, head_copy_p + head_copy_l);
			 DEBUG("Refined pairings on T:  [%9d --%10d) %.2s...%.2s [%9d --%10d)",
					 node->t, node->t + node_copy_l,
					 GEN_seq+node->t + node_copy_l, GEN_seq+head_copy_t-2,
					 head_copy_t, head_copy_t + head_copy_l);
		  }

// Se ( ho un piccolo gap "netto" su T  OR ho un introne di lunghezza >= al minimo )
//  --> aggiungo
		  if((gap_length_on_t <= fl) || is_intron_on_t ) {
			 copy_embedding=list_copy(embedding, (copy_item) copy_ppairing);
			 ppairing head_copy=(ppairing)list_head(copy_embedding);
			 head_copy->p=head_copy_p;
			 head_copy->t=head_copy_t;
			 head_copy->l=head_copy_l;

			 ppairing node_copy=pairing_simple_copy(node);
			 node_copy->l=node_copy_l;

			 list_add_to_head(copy_embedding, node_copy);
			 list_add_to_head(return_embedding_list, copy_embedding);
		  }
		}
	 }
  }

  return return_embedding_list;
}

//Create a new factor and set its parameter
static pfactor create_and_set_factor(int donor_EST_start, int donor_EST_end, int donor_GEN_start, int donor_GEN_end)
{
  pfactor new_factor;

  new_factor=factor_create();

//Set the new factor with input parameters
  new_factor->EST_start=donor_EST_start;
  new_factor->EST_end=donor_EST_end;
  new_factor->GEN_start=donor_GEN_start;
  new_factor->GEN_end=donor_GEN_end;

  return new_factor;
}


//Create a new subtree embeddings object and set its parameter
static psubtree_embeddings create_and_set_subtree_embeddings(ppairing root, plist embeddings)
{
  psubtree_embeddings new_sub_e;

  new_sub_e=subtree_embeddings_create();

  new_sub_e->root=root;
  new_sub_e->embeddings=embeddings;

  return new_sub_e;
}

static plist get_computed_subtree_embeddings(ppairing root)
{
  plistit list_iter=list_first(list_of_subtree_embeddings);
  plist ret_e;

  bool found=false;
  while(listit_has_next(list_iter) && !found){
	 psubtree_embeddings sub_e=(psubtree_embeddings)listit_next(list_iter);

	 if(root == sub_e->root){
		found=true;
		ret_e=sub_e->embeddings;
	 }
  }

  listit_destroy(list_iter);

  if(found==true)
	 return ret_e;
  else
	 return NULL;
}

//UNUSED
/*
//Funzione per copiare liste di pfactor
static pfactor copy_pfactor(pfactor arg)
{
  pfactor copy = PALLOC(struct _factor);

  *copy = *arg;

  return copy;

}
*/

//Funzione per copiare liste di ppairing
static ppairing copy_ppairing(ppairing arg)
{
  ppairing copy = PALLOC(struct _pairing);

  *copy = *arg;

  return copy;

}

static void print_factorizations(plist factorization_list){
  plistit plist_it_id;

  TRACE("\t\tFactorizations->");

  plist_it_id=list_first(factorization_list);
  while(listit_has_next(plist_it_id)){
	 plist factorization=(plist)listit_next(plist_it_id);
	 print_factorization(factorization);
  }
  listit_destroy(plist_it_id);
}

static void print_factorization(plist factorization){
  plistit plist_it_f;

  TRACE("\t\tFactorization->");

  plist_it_f=list_first(factorization);
  while(listit_has_next(plist_it_f)){
	 pfactor pfact=listit_next(plist_it_f);
	 TRACE("\t\t\tE_START=%d E_END=%d G_START=%d G_END=%d", pfact->EST_start, pfact->EST_end, pfact->GEN_start, pfact->GEN_end);
  }
  listit_destroy(plist_it_f);
}

static void print_embeddings(plist embedding_list){
  plistit plist_it_id;

  TRACE("\t\tEmbeddings->");

  plist_it_id=list_first(embedding_list);
  while(listit_has_next(plist_it_id)){
	 plist embedding=(plist)listit_next(plist_it_id);
	 print_embedding(embedding);
  }
  listit_destroy(plist_it_id);
}

static void print_embedding(plist embedding){
  plistit plist_it_e;

  TRACE("\t\tEmbedding->");

  plist_it_e=list_first(embedding);
  while(listit_has_next(plist_it_e)){
	 ppairing pair=(ppairing)listit_next(plist_it_e);
	 TRACE("\t\t\tp=%d t=%d l=%d", pair->p, pair->t, pair->l);
  }
  listit_destroy(plist_it_e);
}

static void
print_factorization_on_log(const int log_level,
									plist factorization, const char* const gen_seq){
  plistit plist_it_f;

  plist_it_f=list_first(factorization);
  size_t exon_n= 0;
  while(listit_has_next(plist_it_f)){
	 pfactor pfact=(pfactor)listit_next(plist_it_f);
	 ++exon_n;
	 LOG(log_level, "   %3zu) E_START=%8d  E_END=%8d  --  G_START=%10d  G_END=%10d  "
		  "%.2s...%.2s",
		  exon_n, pfact->EST_start, pfact->EST_end, pfact->GEN_start, pfact->GEN_end,
		  (exon_n==1)?"  ":gen_seq+pfact->GEN_start-2,
			!listit_has_next(plist_it_f)?"  ":gen_seq+pfact->GEN_end+1);
  }
  listit_destroy(plist_it_f);
}

static void
print_factorizations_on_log(const int log_level,
									 plist factorization_list, const char* const gen_seq){
  plistit plist_it_id= list_first(factorization_list);
  size_t factorization_no= 1;
  while(listit_has_next(plist_it_id)){
	 plist factorization=(plist)listit_next(plist_it_id);
	 LOG(log_level, "  Factorization %-3zu [%2zu exon(s)]",
		  factorization_no, list_size(factorization));
	 print_factorization_on_log(log_level, factorization, gen_seq);
	 ++factorization_no;
  }
  listit_destroy(plist_it_id);
}


/*
//UNUSED
static void print_embeddings_for_info(plist embedding_list){
  plistit plist_it_id;

  INFO("\t\tEmbeddings>");

  plist_it_id=list_first(embedding_list);
  while(listit_has_next(plist_it_id)){
	 plist embedding=(plist)listit_next(plist_it_id);
	 print_embedding_for_info(embedding);
  }
  listit_destroy(plist_it_id);
}

static void print_embedding_for_info(plist embedding){
  plistit plist_it_e;

  INFO("\t\tEmbedding->");

  plist_it_e=list_first(embedding);
  while(listit_has_next(plist_it_e)){
	 ppairing pair=(ppairing)listit_next(plist_it_e);
	 INFO("\t\t\tp=%d t=%d l=%d", pair->p, pair->t, pair->l);
  }
  listit_destroy(plist_it_e);
}
*/

//L'argomento deve essere una lista di fattorizzazioni (lista di liste di fattori)
/*static void factorization_list_destroy(plist factorization_list){
  if(factorization_list != NULL)
	 list_destroy(factorization_list,(delete_function)factorization_destroy);
}*/

//L'argomento deve essere una lista di embeddings (lista di liste di pairing)
static void embedding_list_destroy(plist embedding_list){
  if(embedding_list != NULL)
	 list_destroy(embedding_list,(delete_function)embedding_destroy);
}

//UNUSED
/*
//Controlla che la fattorizzazione passata come input non tagli un suffisso/prefisso eccessivo
//factorization deve essere una lista di pfactor
static bool test_prefix_suffix(plist factorization, int est_length, pconfiguration config, pEST_info info){
  pfactor head_factor;
  pfactor tail_factor;

  my_assert(info->pref_polyA_length == -1 || info->pref_polyT_length == -1);
  my_assert(info->suff_polyA_length == -1 || info->suff_polyT_length == -1);

  head_factor=(pfactor) list_head(factorization);
  int actual_cut_prefix=head_factor->EST_start;
  int pref_poly_red=(info->pref_polyA_length != -1)?(info->pref_polyA_length):((info->pref_polyT_length != -1)?(info->pref_polyT_length):(0));
  actual_cut_prefix=actual_cut_prefix-pref_poly_red;
  if(config->max_prefix_discarded != -1 && actual_cut_prefix > config->max_prefix_discarded){
	 DEBUG("\tThe prefix cut is too long");
	 return false;
  }

  tail_factor=(pfactor) list_tail(factorization);
  int actual_cut_suffix=est_length-tail_factor->EST_end-1;
  int suff_poly_red=(info->suff_polyA_length != -1)?(info->suff_polyA_length):((info->suff_polyT_length != -1)?(info->suff_polyT_length):(0));
  actual_cut_suffix=actual_cut_suffix-suff_poly_red;
  if(config->max_suffix_discarded != -1 && actual_cut_suffix > config->max_suffix_discarded){
	 DEBUG("\tThe suffix cut is too long");
	 return false;
  }

  return true;
}
*/

//UNUSED
/*
//Comparator per confrontare due fattori (pfactor)
static int factor_compare(const pfactor* pf1, const pfactor* pf2) {
  my_assert(pf1!=NULL);
  my_assert(pf2!=NULL);
  pfactor p1= *pf1;
  pfactor p2= *pf2;
  my_assert(p1!=NULL);
  my_assert(p2!=NULL);

  if (p1->EST_start == p2->EST_start && p1->EST_end == p2->EST_end) {
	 if (p1->GEN_start == p2->GEN_start && p1->GEN_end == p2->GEN_end) {
		return 0;
	 }
  }

  return 1;
}

//Comparator per confrontare due pairing (ppairing)
static int perfect_pairing_compare(const ppairing* pp1, const ppairing* pp2) {
  my_assert(pp1!=NULL);
  my_assert(pp2!=NULL);
  ppairing p1= *pp1;
  ppairing p2= *pp2;
  my_assert(p1!=NULL);
  my_assert(p2!=NULL);

  if (p1->p == p2->p && p1->t == p2->t) {
	 if (p1->l == p2->l) {
		return 0;
	 }
  }

  return 1;
}
*/

//Comparator per confrontare liste di fattori (pfactor)
//Due fattori sono uguali se la differenza tra gli ss (solo sulla genomica) e' contenuta in allowed_diff
//cfr_type = 0 ==> confronto su entrambi gli estremi
//cfr_type = -2 ==> confronto sul right end
//cfr_type = -1 ==> confronto sul left end
//cfr_type = 1 ==> confronto sul left (il right di pf2 puo' essere < del right di pf1, ma non puo'
//superare il right di pf1 di 5 volte allowed_diff)
//cfr_type = 2 ==> confronto sul right (il left di pf2 puo' essere > del left di pf1, ma non puo'
//arretrare rispetto al left di pf1 di 5 volte allowed_diff)
//DA FARE: l1 e' la lista relativa a pf1 per fare in modo che con cfr_type=1 (=2) esista anche un controllo sul
//superamento (arretramento) che (nel caso in cui verifichi le condizioni sopra descritte) non deve essere
//dello stesso ordine della lunghezza complessiva dei fattori di genomica che precedono (seguono) pf1.
//NOTA: e' il caso di limitare il controllo sulla sola genomica??
static int relaxed_factor_compare(const pfactor* pf1, const pfactor* pf2, int cfr_type, int allowed_diff, plist l1) {
  my_assert(pf1!=NULL);
  my_assert(pf2!=NULL);
  pfactor p1= *pf1;
  pfactor p2= *pf2;
  my_assert(p1!=NULL);
  my_assert(p2!=NULL);
  my_assert(l1!=NULL);

  my_assert(cfr_type == 0 || abs(cfr_type) == 1 || abs(cfr_type) == 2);

  //Se i due fattori non si sovrappongono e' inutile il controllo.
  if(p1->GEN_start < p2->GEN_start && p1->GEN_end < p2->GEN_start){
	  return 1;
  }
  else{
	  if(p2->GEN_start < p1->GEN_start && p2->GEN_end < p1->GEN_start)
		  return 1;
  }

//printf("\tCFR %d %d-%d %d-%d\n", cfr_type, p1->GEN_start+1, p1->GEN_end+1, p2->GEN_start+1, p2->GEN_end+1);

  //int max_unconf_diff=5*allowed_diff;
  int max_unconf_diff=20;

  if(max_unconf_diff <= 0)
	 max_unconf_diff=5;

  if(cfr_type == 0){
//Confronto solo sulla genomica
	 if (abs(p1->GEN_end-p2->GEN_end)<=allowed_diff){
		if(abs(p1->GEN_start-p2->GEN_start)<=allowed_diff){
		  return 0;
		}
	 }
  }

  if(abs(cfr_type) == 2){
	 if(abs(p1->GEN_end-p2->GEN_end)<=allowed_diff){
		if(cfr_type == 2){
		  if(p1->GEN_start-p2->GEN_start > max_unconf_diff){
//printf("\tCFR %d %d-%d %d-%d\n", cfr_type, p1->GEN_start+1, p1->GEN_end+1, p2->GEN_start+1, p2->GEN_end+1);
			 return 1;
		  }
		  else{
			 if(p1->GEN_start-p2->GEN_start > 0){
				plistit it1=list_first(l1);
				int tot_l=0;
				bool stop=false;
				while(listit_has_next(it1) && !stop){
				  pfactor f=(pfactor)listit_next(it1);
				  if(p1->GEN_start == f->GEN_start){
					 stop=true;
				  }
				  else{
					 tot_l+=(f->GEN_end-f->GEN_start+1);
				  }
				}
				my_assert(stop==true);
				if(abs(p1->GEN_start-p2->GEN_start-tot_l) < 10){
//printf("\tCFR %d %d-%d %d-%d\n", cfr_type, p1->GEN_start+1, p1->GEN_end+1, p2->GEN_start+1, p2->GEN_end+1);
				  return 1;
				}
			 }
		  }
		}
		return 0;
	 }
  }

  if(abs(cfr_type) == 1){
	 if (abs(p1->GEN_start-p2->GEN_start)<=allowed_diff){
		if(cfr_type == 1){
		  if(p2->GEN_end-p1->GEN_end > max_unconf_diff){
//printf("\tCFR %d %d-%d %d-%d\n", cfr_type, p1->GEN_start+1, p1->GEN_end+1, p2->GEN_start+1, p2->GEN_end+1);
			 return 1;
		  }
		  else{
			 if(p2->GEN_end-p1->GEN_end > 0){
				plistit it1=list_last(l1);
				int tot_l=0;
				bool stop=false;
				while(listit_has_prev(it1) && !stop){
				  pfactor f=(pfactor)listit_prev(it1);
				  if(p1->GEN_start == f->GEN_start){
					 stop=true;
				  }
				  else{
					 tot_l+=(f->GEN_end-f->GEN_start+1);
				  }
				}
				my_assert(stop==true);
				if(abs(p2->GEN_end-p1->GEN_end-tot_l) < 20){
//printf("\tCFR %d %d-%d %d-%d\n", cfr_type, p1->GEN_start+1, p1->GEN_end+1, p2->GEN_start+1, p2->GEN_end+1);
				  return 1;
				}
			 }
		  }
		}
		return 0;
	 }
  }

//printf("\tCFR %d %d-%d %d-%d\n", cfr_type, p1->GEN_start+1, p1->GEN_end+1, p2->GEN_start+1, p2->GEN_end+1);
  return 1;
}

//Calcola la copertura (rispetto alla lunghezza) su P di una fattorizzazione (valore tra 0.0 e 1.0)
//Eventuali gap interni su P sono ignorati
static double compute_coverage(plist factorization, unsigned int length){
  pfactor head=list_head(factorization);
  pfactor tail=list_tail(factorization);

  int cut_prefix_length=head->EST_start;
  int cut_suffix_length=length-tail->EST_end-1;

  int cover=length-(cut_prefix_length+cut_suffix_length);

  return ((double)cover)/(double)length;
}

//Calcola la lunghezza totale dei gap (su P) di una fattorizzazione (tagli prefisso/suffisso esclusi)
static int compute_gapLength(plist factorization){

  my_assert(list_size(factorization) >= 1);

  if(list_size(factorization) == 1)
	 return 0;

  int gapLength=0;

  plistit plist_it_id=list_first(factorization);
  pfactor donor_factor=(pfactor)listit_next(plist_it_id);
  while(listit_has_next(plist_it_id)){
	 pfactor accept_factor=(pfactor)listit_next(plist_it_id);
	 gapLength+=accept_factor->EST_start-donor_factor->EST_end-1;
	 donor_factor=accept_factor;
  }
  listit_destroy(plist_it_id);

  return gapLength;
}

static plist get_factorizations_from_embeddings(plist embedding_list, pconfiguration config, pEST_info info, int est_length){
  plist return_factorization_list=list_create();
  int fl=2*(config->min_factor_len);
  unsigned int count=1;
  plist factorization;

  plistit plist_it_id=list_first(embedding_list);
  while(listit_has_next(plist_it_id)){
	 plist embedding=(plist)listit_next(plist_it_id);
	 DEBUG("The embedding %d", count);
	 print_embedding(embedding);
	 ppairing head=(ppairing)list_head(embedding);
	 int actual_cut_prefix=head->p;
	 int pref_poly_red=(info->pref_polyA_length != -1)?(info->pref_polyA_length):((info->pref_polyT_length != -1)?(info->pref_polyT_length):(0));
	 actual_cut_prefix=actual_cut_prefix-pref_poly_red;

	 ppairing tail=(ppairing)list_tail(embedding);
	 int actual_cut_suffix=est_length-(tail->p+tail->l);
	 int suff_poly_red=(info->suff_polyA_length != -1)?(info->suff_polyA_length):((info->suff_polyT_length != -1)?(info->suff_polyT_length):(0));
	 actual_cut_suffix=actual_cut_suffix-suff_poly_red;
	 factorization=list_create();
	 plistit plist_it_e=list_first(embedding);
	 pfactor last_factor= NULL;
	 bool stop=false;
	 while(listit_has_next(plist_it_e) && stop == false){
		ppairing pair=(ppairing)listit_next(plist_it_e);
		if(list_is_empty(factorization)){
		  last_factor=create_and_set_factor(pair->p, pair->p+pair->l-1, pair->t, pair->t+pair->l-1);
		  list_add_to_tail(factorization, last_factor);
		  DEBUG("\t...first factor %d-%d (%d-%d) added!", last_factor->EST_start, last_factor->EST_end, last_factor->GEN_start, last_factor->GEN_end);
		} else {
//Si genera un introne su T
		  if ((pair->t-last_factor->GEN_end-1) > fl) {
			 pfactor tail_factor=(pfactor)list_tail(factorization);
			 last_factor= create_and_set_factor(pair->p, pair->p+pair->l-1, pair->t, pair->t+pair->l-1);
//L'eventuale gap su P viene lasciato e poi risolto in fase di post-processing
			 if (0 && (last_factor->EST_start-tail_factor->EST_end-1) > 0){
				tail_factor->EST_end=tail_factor->EST_end+((last_factor->EST_start-tail_factor->EST_end-1)/2);
				last_factor->EST_start=tail_factor->EST_end+1;
				DEBUG("\t...right end of the last factor, updated on P to %d!", tail_factor->EST_end);
			 }
			 list_add_to_tail(factorization, last_factor);
			 DEBUG("\t...new factor %d-%d (%d-%d) added!", last_factor->EST_start, last_factor->EST_end, last_factor->GEN_start, last_factor->GEN_end);
		  } else {
//I due esoni vengono fusi (e quindi l'eventuale gap su P scompare)
			 last_factor=(pfactor)list_tail(factorization);
			 last_factor->EST_end=pair->p+pair->l-1;
			 last_factor->GEN_end=pair->t+pair->l-1;

			 DEBUG("\t...updating last factor to %d-%d (%d-%d)!",
					 last_factor->EST_start, last_factor->EST_end,
					 last_factor->GEN_start, last_factor->GEN_end);
		  }
		}
	 }
	 if(stop == false){
		list_add_to_tail(return_factorization_list, factorization);
	 }
	 count++;
	 listit_destroy(plist_it_e);
  }
  listit_destroy(plist_it_id);

  return return_factorization_list;
}

//Controlla la massimalita' tra due embedding
//Ritorna 2 se il primo e' massimale, 1 se entrambi sono massimali e 0 se il secondo e' massimale.
//Il confronto parte dal primo pairing per entrambi gli embedding in quanto e' una procedura applicata
//agli embedding relativi ad un sottoalbero del MEG radicato in un determinato nodo.
static char maximality_relation(plist add_emb, plist cmp_emb){

  if(list_size(add_emb) > list_size(cmp_emb)){
	 plistit add_pair_it=list_first(add_emb);
	 plistit cmp_pair_it=list_first(cmp_emb);

	 bool check=true;
	 while(listit_has_next(add_pair_it) && listit_has_next(cmp_pair_it) && check){
		ppairing add_pair=(ppairing)listit_next(add_pair_it);
		ppairing cmp_pair=(ppairing)listit_next(cmp_pair_it);

		if(cmp_pair->p < add_pair->p || (cmp_pair->p+cmp_pair->l > add_pair->p+add_pair->l)){
		  check=false;
		}
		else{
		  if(cmp_pair->t < add_pair->t || (cmp_pair->t+cmp_pair->l > add_pair->t+add_pair->l)){
			 check=false;
		  }
		}
	 }

	 listit_destroy(add_pair_it);
	 listit_destroy(cmp_pair_it);

	 if(check)
		return 2;
	 else
		return 1;
  }

  if(list_size(add_emb) < list_size(cmp_emb)){
	 plistit add_pair_it=list_first(add_emb);
	 plistit cmp_pair_it=list_first(cmp_emb);

	 bool check=true;
	 while(listit_has_next(add_pair_it) && listit_has_next(cmp_pair_it) && check){
		ppairing add_pair=(ppairing)listit_next(add_pair_it);
		ppairing cmp_pair=(ppairing)listit_next(cmp_pair_it);

		if(add_pair->p < cmp_pair->p || (add_pair->p+add_pair->l > cmp_pair->p+cmp_pair->l)){
		  check=false;
		}
		else{
		  if(add_pair->t < cmp_pair->t || (add_pair->t+add_pair->l > cmp_pair->t+cmp_pair->l)){
			 check=false;
		  }
		}
	 }

	 listit_destroy(add_pair_it);
	 listit_destroy(cmp_pair_it);

	 if(check)
		return 0;
	 else
		return 1;
  }

  plistit add_pair_it_2=list_first(add_emb);
  plistit cmp_pair_it_2=list_first(cmp_emb);

  bool check_2=true;
  while(listit_has_next(add_pair_it_2) && listit_has_next(cmp_pair_it_2) && check_2){
	 ppairing add_pair=(ppairing)listit_next(add_pair_it_2);
	 ppairing cmp_pair=(ppairing)listit_next(cmp_pair_it_2);

	 if(add_pair->p < cmp_pair->p || (add_pair->p+add_pair->l > cmp_pair->p+cmp_pair->l)){
		check_2=false;
	 }
	 else{
		if(add_pair->t < cmp_pair->t || (add_pair->t+add_pair->l > cmp_pair->t+cmp_pair->l)){
		  check_2=false;
		}
	 }
  }
  listit_destroy(add_pair_it_2);
  listit_destroy(cmp_pair_it_2);

  if(check_2){
	 return 0;
  }

  add_pair_it_2=list_first(add_emb);
  cmp_pair_it_2=list_first(cmp_emb);

  check_2=true;
  while(listit_has_next(add_pair_it_2) && listit_has_next(cmp_pair_it_2) && check_2){
	 ppairing add_pair=(ppairing)listit_next(add_pair_it_2);
	 ppairing cmp_pair=(ppairing)listit_next(cmp_pair_it_2);

	 if(cmp_pair->p < add_pair->p || (cmp_pair->p+cmp_pair->l > add_pair->p+add_pair->l)){
		check_2=false;
	 }
	 else{
		if(cmp_pair->t < add_pair->t || (cmp_pair->t+cmp_pair->l > add_pair->t+add_pair->l)){
		  check_2=false;
		}
	 }
  }

  listit_destroy(add_pair_it_2);
  listit_destroy(cmp_pair_it_2);

  if(check_2)
	 return 2;

  return 1;
}

static bool check_gap_errors(plist factorization, char *est_seq, char *gen_seq, pconfiguration config){

	 //Parametro da mettere in config
	 unsigned int threshold_ed=20;

	 size_t out_offset_p, out_offset_t1, out_offset_t2;
	 unsigned int out_edit_distance, tot_out_edit_distance;
	 bool ok=true;

	 DEBUG("... checking of the gap errors");

	 tot_out_edit_distance=0;

	  plistit plist_it_id=list_first(factorization);
	  pfactor donor_factor=(pfactor)listit_next(plist_it_id);
	  while(listit_has_next(plist_it_id) && ok == true){
		 pfactor accept_factor=(pfactor)listit_next(plist_it_id);

		 size_t gapLengthOnP=accept_factor->EST_start-donor_factor->EST_end-1;

		 if(gapLengthOnP > 0){
			 size_t gapLengthOnT=accept_factor->GEN_start-donor_factor->GEN_end-1;

			 if(gapLengthOnP > gapLengthOnT)
				 FATAL("...the gap on P cannot be greater than the gap on T!");

			 int max_errs=gapLengthOnP;
			 char *p=real_substring(donor_factor->EST_end+1, gapLengthOnP, est_seq);
			 char *t=real_substring(donor_factor->GEN_end+1, gapLengthOnT, gen_seq);
			 DEBUG("...the gap from %d %d on P", donor_factor->EST_end+1, accept_factor->EST_start-1);
			 ok= refine_borders(p, gapLengthOnP,
									 t, gapLengthOnT,
									 max_errs,
									 &out_offset_p,
									 &out_offset_t1, &out_offset_t2,
									 &out_edit_distance);
			 pfree(p);
			 pfree(t);
			 if(ok == true){
				 DEBUG("\t...is ok!");
				 tot_out_edit_distance+=out_edit_distance;
				 donor_factor->EST_end+=out_offset_p;
				 accept_factor->EST_start=donor_factor->EST_end+1;
				 donor_factor->GEN_end+=out_offset_t1;
				 accept_factor->GEN_start-=gapLengthOnT-out_offset_t2;
			 }
			 else{
				 DEBUG("\t...is not ok and this factorization will be deleted!");
			 }
		 }

		 donor_factor=accept_factor;
	  }
	  listit_destroy(plist_it_id);

	  if(ok == true && tot_out_edit_distance > threshold_ed){
		  DEBUG("...The factorization will be deleted, since there are too many errors!");
		  ok=false;
	  }

	  if(ok == true){
		  plist_it_id=list_first(factorization);
		  pfactor d_factor=(pfactor)listit_next(plist_it_id);
		  while(listit_has_next(plist_it_id)){
			 pfactor a_factor=(pfactor)listit_next(plist_it_id);

			 //XXX
			 if(a_factor->GEN_start-d_factor->GEN_end-1 <= 3){
				 d_factor->EST_end=a_factor->EST_end;
				 d_factor->GEN_end=a_factor->GEN_end;
				 list_remove_at_iterator(plist_it_id, (delete_function)pairing_destroy_2);
			 }
			 else{
				 if(a_factor->GEN_start-d_factor->GEN_end-1 < config->min_intron_length){
					  DEBUG("...The intron %d-%d is too small!\n", d_factor->GEN_end+1, a_factor->GEN_start-1);
				 }
				 d_factor=a_factor;
			 }
		  }
		  listit_destroy(plist_it_id);
	  }

	  return ok;
}

char* real_substring(const int index, const int length, char* string){

  my_assert(index>=0);
  my_assert(string != NULL);
  my_assert(length >= 0);

  //TRACE("index value is: %d", index);
  //TRACE("string value is: |%s|", string);
  //TRACE("length value is: %d", length);

  char * const ris= c_palloc(length+1);
  memcpy(ris, string+index, length*sizeof(char));
  ris[length]= '\0';
  return ris;
}//end real_substring


void print_split_string_on_stderr(const int token_dim, char* string){

  my_assert(token_dim > 0);
  my_assert(string != NULL);

  TRACE("token dim is: %d", token_dim);
  TRACE("string value is: |%s|", string);

  size_t length=strlen(string);
  unsigned int i=0;
  while(i < length){
	  TRACE("token from %d", i);
	  char *token=real_substring((int)i, token_dim, string);
	  fprintf(stderr, "%s\n", token);
	  pfree(token);
	  i+=token_dim;
  }
}//end real_substring

//inutile...
plist clean_low_complexity_exons(plist factorization, char *genomic_sequence, char *est_sequence){
	my_assert(factorization != NULL);
	my_assert(!list_is_empty(factorization));

	DEBUG("Check complexity:");
	print_factorization(factorization);

	plist help_list=list_create();

	int best_left_index=-1, best_right_index=-1;
	int best_est_cover=-1;

	int left_index=1;
	int right_index=1;
	plistit plist_f_t_r=list_first(factorization);
	while(listit_has_next(plist_f_t_r)){
		pfactor exon=(pfactor)listit_next(plist_f_t_r);

		double gendscore=0.0f;
		double estdscore=0.0f;
		//Provvisorio solo per evitare problemi (da sistemare prima)
		if(exon->GEN_start <= exon->GEN_end){
			gendscore=dustScoreByLeftAndRight(genomic_sequence, exon->GEN_start, exon->GEN_end);
			estdscore=dustScoreByLeftAndRight(est_sequence, exon->EST_start, exon->EST_end);
		}

		if(gendscore > 1.0 || estdscore > 1.0){
			TRACE("\t exon %d-%d (%d-%d) has a low complexity", exon->GEN_start, exon->GEN_end, exon->EST_start, exon->EST_end);
			if(left_index < right_index){
				my_assert(!list_is_empty(help_list));

				pfactor first_exon=list_head(help_list);
				pfactor last_exon=list_tail(help_list);
				int est_cover=last_exon->EST_end-first_exon->EST_start+1;
				if(est_cover > best_est_cover){
					best_left_index=left_index;
					best_right_index=right_index-1;
					best_est_cover=est_cover;
				}

				list_destroy(help_list, noop_free);
				help_list=list_create();
			}
			left_index=right_index+1;
		}
		else{
			list_add_to_tail(help_list, exon);
		}
		right_index++;
	}
	listit_destroy(plist_f_t_r);

	my_assert(right_index-1 == (int)list_size(factorization));

	if(left_index < right_index){
		if(left_index != 1){
			pfactor first_exon=list_head(help_list);
			pfactor last_exon=list_tail(help_list);
			int est_cover=last_exon->EST_end-first_exon->EST_start+1;
			if(est_cover > best_est_cover){
				best_left_index=left_index;
				best_right_index=right_index-1;
				best_est_cover=est_cover;
			}
		}
		else{
			best_left_index=left_index;
			best_right_index=right_index-1;
		}
	}

	list_destroy(help_list, noop_free);

	int size=(int)list_size(factorization);

	if(best_left_index == -1 || best_right_index == -1){
		my_assert(best_left_index == -1 && best_right_index == -1);
		int remove=size;
		while(remove > 0){
			list_remove_from_tail(factorization);
			remove--;
		}
	}
	else{
		int remove_head=best_left_index-1;
		while(remove_head > 0){
			list_remove_from_head(factorization);
			remove_head--;
		}

		int remove_tail=best_right_index+1;
		while(remove_tail <= size){
			list_remove_from_tail(factorization);
			remove_tail++;
		}
	}

	return factorization;
}

plist clean_low_complexity_exons_2(plist factorization, char *genomic_sequence, char *est_sequence){
	my_assert(genomic_sequence != NULL);
	my_assert(est_sequence != NULL);
	my_assert(factorization != NULL);
	my_assert(!list_is_empty(factorization));

	DEBUG("Check complexity:");
	print_factorization(factorization);

	pintlist split_list=intlist_create();
	plistit plist_f_t_r=list_first(factorization);

	int index=1;
	while(listit_has_next(plist_f_t_r)){
		pfactor exon=(pfactor)listit_next(plist_f_t_r);

		double gendscore=0.0f;
		double estdscore=0.0f;
		//Provvisorio solo per evitare problemi (da sistemare prima)
		if(exon->GEN_start <= exon->GEN_end){
			gendscore=dustScoreByLeftAndRight(genomic_sequence, exon->GEN_start, exon->GEN_end);
			estdscore=dustScoreByLeftAndRight(est_sequence, exon->EST_start, exon->EST_end);
		}

		//XXX
		if(gendscore > 1.0 || estdscore > 1.0){
			TRACE("\t exon %d-%d (%d-%d) has a low complexity", exon->GEN_start, exon->GEN_end, exon->EST_start, exon->EST_end);
			intlist_add_to_tail(split_list, index);
		}
		index++;
	}
	listit_destroy(plist_f_t_r);

	factorization=update_with_subfact_with_best_coverage(factorization, split_list);

	intlist_destroy(split_list);

	return factorization;
}

plist clean_external_exons(plist factorization, char *genomic_sequence, char *est_sequence){
	my_assert(genomic_sequence != NULL);
	my_assert(est_sequence != NULL);
	my_assert(factorization != NULL);

	if(list_is_empty(factorization)){
		return factorization;
	}

	pfactor head=list_remove_from_head(factorization);

	int head_length=head->GEN_end-head->GEN_start+1;
	bool head_is_ok=true;

	//XXX
	if(head_length < 10){
		head_is_ok=false;
	}

	//XXX
	if(head_is_ok && head_length < 20){
		if(genomic_sequence[head->GEN_end+1] != 'G' && genomic_sequence[head->GEN_end+1] != 'g')
			head_is_ok=false;
		else{
			if((genomic_sequence[head->GEN_end+2] != 'T' && genomic_sequence[head->GEN_end+2] != 't') && (genomic_sequence[head->GEN_end+2] != 'C' && genomic_sequence[head->GEN_end+2] != 'c'))
				head_is_ok=false;
			else{
				if((int)list_size(factorization) >= 1){
					pfactor next_exon=list_head(factorization);
					if(genomic_sequence[next_exon->GEN_start-2] != 'A' && genomic_sequence[next_exon->GEN_start-2] != 'a')
						head_is_ok=false;
					else{
						if(genomic_sequence[next_exon->GEN_start-1] != 'G' && genomic_sequence[next_exon->GEN_start-1] != 'g')
							head_is_ok=false;
					}
				}
				else{
					head_is_ok=false;
				}
			}
		}
		if(head_is_ok){
			char *gen_exon_seq=real_substring(head->GEN_start, head->GEN_end-head->GEN_start+1, genomic_sequence);
			char *est_exon_seq=real_substring(head->EST_start, head->EST_end-head->EST_start+1, est_sequence);
			size_t l1=strlen(gen_exon_seq);
			size_t l2=strlen(est_exon_seq);
			unsigned int* M=edit_distance(gen_exon_seq, l1, est_exon_seq, l2);
			int error=M[(l1+1)*(l2+1)-1];
			pfree(M);
			pfree(gen_exon_seq);
			pfree(est_exon_seq);
			if(error > 0)
				head_is_ok=false;
		}
	}

	if(head_is_ok == true)
		list_add_to_head(factorization, head);
	else{
		DEBUG("First exon %d-%d removed!", head->GEN_start, head->GEN_end);
	}

	if(list_is_empty(factorization))
		return factorization;

	pfactor tail=list_remove_from_tail(factorization);

	int tail_length=tail->GEN_end-tail->GEN_start+1;
	bool tail_is_ok=true;

	//XXX
	if(tail_length < 10){
		tail_is_ok=false;
	}

	//XXX
	if(tail_is_ok && tail_length < 20){
		if(genomic_sequence[tail->GEN_start-2] != 'A' && genomic_sequence[tail->GEN_start-2] != 'a')
			tail_is_ok=false;
		else{
			if(genomic_sequence[tail->GEN_start-1] != 'G' && genomic_sequence[tail->GEN_start-1] != 'g')
				tail_is_ok=false;
			else{
				if((int)list_size(factorization) >= 1){
					pfactor prev_exon=list_tail(factorization);
					if(genomic_sequence[prev_exon->GEN_end+1] != 'G' && genomic_sequence[prev_exon->GEN_end+1] != 'g')
						tail_is_ok=false;
					else{
						if((genomic_sequence[prev_exon->GEN_end+2] != 'T' && genomic_sequence[prev_exon->GEN_end+2] != 't') && (genomic_sequence[prev_exon->GEN_end+2] != 'C' && genomic_sequence[prev_exon->GEN_end+2] != 'c'))
							tail_is_ok=false;
					}
				}
				else{
					tail_is_ok=false;
				}
			}
		}
		if(tail_is_ok){
			char *gen_exon_seq=real_substring(tail->GEN_start, tail->GEN_end-tail->GEN_start+1, genomic_sequence);
			char *est_exon_seq=real_substring(tail->EST_start, tail->EST_end-tail->EST_start+1, est_sequence);
			size_t l1=strlen(gen_exon_seq);
			size_t l2=strlen(est_exon_seq);
			unsigned int* M=edit_distance(gen_exon_seq, l1, est_exon_seq, l2);
			int error=M[(l1+1)*(l2+1)-1];
			pfree(M);
			pfree(gen_exon_seq);
			pfree(est_exon_seq);
			if(error > 0)
				tail_is_ok=false;
		}
	}

	if(tail_is_ok == true)
		list_add_to_tail(factorization, tail);
	else{
		DEBUG("Last exon %d-%d removed!", tail->GEN_start, tail->GEN_end);
	}

	return factorization;
}

static unsigned int
compute_maximum_edit_distance_for_exons(const size_t exon_length) {
  double max_error_rate= 0.0;
  if (exon_length > 100)
	 max_error_rate= 0.030;
  else if (exon_length > 50)
	 max_error_rate= 0.035;
  else
	 max_error_rate= 0.040;
  const unsigned int max_error= (unsigned int)MAX(1.0, ceil(exon_length*max_error_rate));
  TRACE("Exon length: %5zu.   Maximum edit distance: %ud", exon_length, max_error);
  return max_error;
}


plist clean_noisy_exons(plist factorization, char *genomic_sequence, char *est_sequence, bool only_internals){
	my_assert(genomic_sequence != NULL);
	my_assert(est_sequence != NULL);
	my_assert(factorization != NULL);
	my_assert(!list_is_empty(factorization));

	DEBUG("Check noise:");
	print_factorization(factorization);

	pintlist split_list=intlist_create();
	plistit plist_f_t_r=list_first(factorization);

	size_t size=list_size(factorization);

	int index=(only_internals == true)?(2):(1);
	int last_index=(only_internals == true)?((int)(size-1)):((int)size);
	if(only_internals == true)
		listit_next(plist_f_t_r);
	while(listit_has_next(plist_f_t_r) && index <= last_index){
		pfactor exon=(pfactor)listit_next(plist_f_t_r);

		int exon_length=exon->GEN_end-exon->GEN_start+1;

// See issue #5
		unsigned int max_allowed_error= compute_maximum_edit_distance_for_exons(exon_length);

		bool ok=false;

		//Provvisorio solo per evitare problemi (da sistemare prima)
		if(exon->GEN_start <= exon->GEN_end){
			char *gen_exon_seq=real_substring(exon->GEN_start, exon->GEN_end-exon->GEN_start+1, genomic_sequence);
			char *est_exon_seq=real_substring(exon->EST_start, exon->EST_end-exon->EST_start+1, est_sequence);

			unsigned int edit;
			ok=K_band_edit_distance(gen_exon_seq, est_exon_seq, max_allowed_error, &edit);

			pfree(gen_exon_seq);
			pfree(est_exon_seq);
		}

		if(!ok){
			TRACE("\t exon %d-%d (%d-%d) has a high error and the maximum allowed is %u",
					exon->GEN_start, exon->GEN_end,
					exon->EST_start, exon->EST_end,
					max_allowed_error);
			intlist_add_to_tail(split_list, index);
		}
		index++;
	}
	listit_destroy(plist_f_t_r);

	factorization=update_with_subfact_with_best_coverage(factorization, split_list);

	intlist_destroy(split_list);

	return factorization;
}

plist update_with_subfact_with_best_coverage(plist factorization, pintlist split_list){
	my_assert(factorization != NULL);
	my_assert(!list_is_empty(factorization));

	if(intlist_is_empty(split_list)){
		return factorization;
	}

	int best_left_index=-1, best_right_index=-1;
	int best_est_cover=-1;

	plistit f_it=list_first(factorization);
	pintlistit int_it=intlist_first(split_list);

	int left_index=1;
	while(intlistit_has_next(int_it)){
		int right_index=(int)intlistit_next(int_it);
		my_assert(listit_has_next(f_it));
		pfactor left_exon=(pfactor)listit_next(f_it);
		pfactor right_exon=left_exon;
		if(left_index < right_index){
			int times=right_index-left_index-1;
			while(times > 0){
				my_assert(listit_has_next(f_it));
				right_exon=(pfactor)listit_next(f_it);
				times--;
			}
			int est_cover=right_exon->EST_end-left_exon->EST_start+1;
			if(est_cover > best_est_cover){
				best_left_index=left_index;
				best_right_index=right_index-1;
				best_est_cover=est_cover;
			}
			my_assert(listit_has_next(f_it));
			listit_next(f_it);
		}
		left_index=right_index+1;
	}

	intlistit_destroy(int_it);

	int size=(int)list_size(factorization);

	if(left_index <= (int)size){
		my_assert(left_index != 1);
		my_assert(listit_has_next(f_it));
		pfactor left_exon=(pfactor)listit_next(f_it);
		pfactor right_exon=left_exon;
		int times=list_size(factorization)-left_index;
		while(times > 0){
			my_assert(listit_has_next(f_it));
			right_exon=(pfactor)listit_next(f_it);
			times--;
		}
		int est_cover=right_exon->EST_end-left_exon->EST_start+1;
		if(est_cover > best_est_cover){
			best_left_index=left_index;
			best_right_index=list_size(factorization);
			best_est_cover=est_cover;
		}
	}

	listit_destroy(f_it);

	if(best_left_index == -1 || best_right_index == -1){
		my_assert(best_left_index == -1 && best_right_index == -1);
		int remove=size;
		while(remove > 0){
			list_remove_from_tail(factorization);
			remove--;
		}
	}
	else{
		int remove_head=best_left_index-1;
		while(remove_head > 0){
			list_remove_from_head(factorization);
			remove_head--;
		}

		int remove_tail=best_right_index+1;
		while(remove_tail <= size){
			list_remove_from_tail(factorization);
			remove_tail++;
		}
	}

	return factorization;
}

bool check_exon_start_end(plist factorization){
	my_assert(factorization != NULL);
	my_assert(!list_is_empty(factorization));

	plistit f_it=list_first(factorization);

	bool exon_ok=true;
	int prev_EST_end=-1;
	int prev_GEN_end=-1;

	while(listit_has_next(f_it) && exon_ok == true){
		pfactor exon_lrcheck=(pfactor)listit_next(f_it);
		if(exon_lrcheck->EST_start > exon_lrcheck->EST_end || exon_lrcheck->GEN_start > exon_lrcheck->GEN_end){
			DEBUG("\t\t\t...exon left > exon right on the EST or on the genomic sequence!");
			exon_ok=false;
		}
		else{
			if(exon_lrcheck->EST_start < prev_EST_end || exon_lrcheck->GEN_start < prev_GEN_end){
				DEBUG("\t\t\t...exon left < previous exon right on the EST or on the genomic sequence!");
				exon_ok=false;
			}
			else{
				prev_EST_end=exon_lrcheck->EST_end;
				prev_GEN_end=exon_lrcheck->GEN_end;
			}
		}
	}
	listit_destroy(f_it);

	return exon_ok;
}

bool check_small_exons(plist factorization){
	my_assert(factorization != NULL);
	my_assert(!list_is_empty(factorization));

	plistit f_it=list_first(factorization);

	bool exon_ok=true;

	while(listit_has_next(f_it) && exon_ok == true){
		pfactor exon_lrcheck=(pfactor)listit_next(f_it);
		if(exon_lrcheck->EST_end-exon_lrcheck->EST_start+1 < 6){
			DEBUG("\t\t\t...exon too small!");
			exon_ok=false;
		}
	}
	listit_destroy(f_it);

	return exon_ok;
}

plist add_if_not_exists(plist factorization_to_be_added, plist factorization_list, pconfiguration config, bool *check_adding){
	my_assert(factorization_to_be_added != NULL);
	my_assert(factorization_list != NULL);
	my_assert(!list_is_empty(factorization_to_be_added));

	plistit cmp_it=list_first(factorization_list);
	bool found=false;
	while(listit_has_next(cmp_it) && found == false){
		plist cmp_f=(plist)listit_next(cmp_it);

		print_factorization(cmp_f);

		int cont_result=0;
		//Gestione del caso in cui entrambe le fattorizzazioni siano composte da un solo esone
		if(list_size(cmp_f) == list_size(factorization_to_be_added) && list_size(cmp_f) == 1){
			pfactor head1=(pfactor)list_head(factorization_to_be_added);
			pfactor head2=(pfactor)list_head(cmp_f);
			if(head1->GEN_start == head2->GEN_start && head1->GEN_end == head2->GEN_end)
				cont_result=-2;
			else{
				if(head1->GEN_start >= head2->GEN_start && head1->GEN_end <= head2->GEN_end)
					cont_result=-1;
				else{
					if(head1->GEN_start <= head2->GEN_start && head1->GEN_end >= head2->GEN_end)
						cont_result=1;
				}
			}
		}
		else{
			//Il confronto e' perfetto sui siti confermati
			cont_result=relaxed_list_contained(factorization_to_be_added, cmp_f, (relaxed_comparator)relaxed_factor_compare, config->max_site_difference);
		}

		//factorization_to_be_added e' uguale (-2) o contenuta (-1) in cmp_f
		if(cont_result < 0){
			if(cont_result == -2){
				pfactor head1=(pfactor)list_head(factorization_to_be_added);
				pfactor head2=(pfactor)list_head(cmp_f);
				if(head1->EST_start < head2->EST_start){
					head2->EST_start=head1->EST_start;
					head2->GEN_start=head1->GEN_start;
				}
				pfactor tail1=(pfactor)list_tail(factorization_to_be_added);
				pfactor tail2=(pfactor)list_tail(cmp_f);
				if(tail1->EST_end > tail2->EST_end){
					tail2->EST_end=tail1->EST_end;
					tail2->GEN_end=tail1->GEN_end;
				}
			}
			found=true;
		}
		else{
			//cmp_f e' contenuta in factorization_to_be_added e ha meno esoni.
			if(cont_result == 1){
				DEBUG("\t\t\t...removing the previous one!");
				list_remove_at_iterator(cmp_it, noop_free);
				factorization_destroy(cmp_f);
			}
		}
	}
	listit_destroy(cmp_it);

	if(found == false)
		list_add_to_tail(factorization_list, factorization_to_be_added);

	*check_adding=!found;

	return factorization_list;
}

bool check_for_not_source_sink_factorization(plist factorization, int EST_length){
	my_assert(factorization != NULL);
	my_assert(!list_is_empty(factorization));

	if((int)list_size(factorization) > 1){
		return true;
	}

	pfactor head=list_head(factorization);

	if(head->EST_start < 0 || head->EST_start >= (int)EST_length)
		return false;
	else
		return true;
}

plist handle_endpoints(plist factorization, char *genomic_sequence, char *est_sequence){
	my_assert(genomic_sequence != NULL);
	my_assert(est_sequence != NULL);
	my_assert(factorization != NULL);
	my_assert(!list_is_empty(factorization));

	DEBUG("Handle end points:");
	print_factorization(factorization);

	size_t gen_length=strlen(genomic_sequence);
	size_t est_length=strlen(est_sequence);

	pfactor head=list_head(factorization);
	my_assert(head->GEN_start >= 0 && head->GEN_end < (int)gen_length);
	my_assert(head->EST_start >= 0 && head->EST_end < (int)est_length);

	char *gen_exon_seq=real_substring(head->GEN_start, head->GEN_end-head->GEN_start+1, genomic_sequence);
	char *est_exon_seq=real_substring(head->EST_start, head->EST_end-head->EST_start+1, est_sequence);
	plist alignments=compute_alignment(est_exon_seq, gen_exon_seq, true);
	palignment alignment=list_head(alignments);

	int j=0;
	int matches=0;
	int cut_factor=head->EST_start;
	int cut_exon=head->GEN_start;
	bool stop=false;
	while(j < alignment->alignment_dim && !stop){
		//Il primo esone deve iniziare con 5 match almeno.
		//XXX
		if(matches > 5){
			stop=true;
		}
		else{
			if(alignment->EST_alignment[j] == alignment->GEN_alignment[j]){
				cut_factor++;
				cut_exon++;
				matches++;
			}
			else{
				if(alignment->EST_alignment[j] != '-')
					cut_factor++;

				if(alignment->GEN_alignment[j] != '-')
					cut_exon++;

				matches=0;
			}
			j++;
		}
	}

	if(stop == false){
		DEBUG("The first exon does not have at least 5 matches!");
		list_remove_from_head(factorization);
	}
	else{
		if(head->GEN_start < cut_exon-matches){
			DEBUG("Left end of the first exon reduced!");
		}
		head->EST_start=cut_factor-matches;
		head->GEN_start=cut_exon-matches;
	}

	alignments_destroy(alignments);

	pfree(gen_exon_seq);
	pfree(est_exon_seq);

	if(list_is_empty(factorization))
		return factorization;

	pfactor tail=list_tail(factorization);
	my_assert(tail->GEN_start >= 0 && tail->GEN_end < (int)gen_length);
	my_assert(tail->EST_start >= 0 && tail->EST_end < (int)est_length);

	gen_exon_seq=real_substring(tail->GEN_start, tail->GEN_end-tail->GEN_start+1, genomic_sequence);
	est_exon_seq=real_substring(tail->EST_start, tail->EST_end-tail->EST_start+1, est_sequence);
	alignments=compute_alignment(est_exon_seq, gen_exon_seq, true);
	alignment=list_head(alignments);

	j=alignment->alignment_dim-1;
	matches=0;
	cut_factor=tail->EST_end;
	cut_exon=tail->GEN_end;
	stop=false;
	while(j >= 0 && stop == false){
		//XXX
		if(matches > 10){
			stop=true;
		}
		else{
			if(alignment->EST_alignment[j] == alignment->GEN_alignment[j]){
				cut_factor--;
				cut_exon--;
				matches++;
			}
			else{
				if(alignment->EST_alignment[j] != '-')
					cut_factor--;

				if(alignment->GEN_alignment[j] != '-')
					cut_exon--;

				matches=0;
			}

			j--;
		}
	}

	int est_cleavage=cut_factor+matches;
	int gen_cleavage=cut_exon+matches;

	//Correzione del cleavage
	int cursor_cut=j+matches+1;
	stop=false;
	while(((alignment->EST_alignment[cursor_cut] == '-' || alignment->GEN_alignment[cursor_cut] == '-') && cursor_cut < alignment->alignment_dim-1) && stop == false){
		if(alignment->EST_alignment[cursor_cut] == '-'){
			int try=cursor_cut+1;
			while(alignment->EST_alignment[try] == '-'){
				try++;
			}
			if(try < alignment->alignment_dim){
				if(alignment->EST_alignment[try] == alignment->GEN_alignment[cursor_cut]){
					alignment->EST_alignment[cursor_cut]=alignment->EST_alignment[try];
					alignment->EST_alignment[try]='-';
					est_cleavage++;
					gen_cleavage++;
				}
				else
					stop=true;
			}
			else
				stop=true;
		}
		else{
			int try=cursor_cut+1;
			while(alignment->GEN_alignment[try] == '-'){
				try++;
			}
			if(try < alignment->alignment_dim){
				if(alignment->GEN_alignment[try] == alignment->EST_alignment[cursor_cut]){
					alignment->GEN_alignment[cursor_cut]=alignment->GEN_alignment[try];
					alignment->GEN_alignment[try]='-';
					est_cleavage++;
					gen_cleavage++;
				}
				else
					stop=true;
			}
			else
				stop=true;
		}
		cursor_cut++;
	}

	if(gen_cleavage >= tail->GEN_start){
		if(gen_cleavage < tail->GEN_end){
			DEBUG("Cleavage reduced!");
		}
		tail->EST_end=est_cleavage;
		tail->GEN_end=gen_cleavage;
	}
	else{
		DEBUG("The last exon does not have at least 10 matches!");
		list_remove_from_tail(factorization);
	}

	alignments_destroy(alignments);

	pfree(gen_exon_seq);
	pfree(est_exon_seq);

	return factorization;
}

bool check_est_coverage(plist factorization, char *est_seq){
	my_assert(factorization != NULL);
	my_assert(est_seq != NULL);
	my_assert(!list_is_empty(factorization));

	size_t est_length=strlen(est_seq);

	pfactor head=(pfactor)list_head(factorization);
	pfactor tail=(pfactor)list_tail(factorization);

	int coverage_length=tail->EST_end-head->EST_start+1;

	double coverage=((double)coverage_length/(double)est_length);

	if(coverage >= 0.35f)
		return true;
	else
		return false;
}
