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
#include "color_matrix.h"
#include "list.h"
#include "types.h"
#include "bit_vector.h"
#include "log.h"

//Stampa la lista dei fattori unici

#if defined (LOG_MSG) && (LOG_LEVEL_DEBUG <= LOG_THRESHOLD)
/*
 * Questa procedura considera la lista passata come argomento come lista di fattori unici (sulla
 * genomica).
 * E quindi e' inutile stampare anche EST_start e EST_end, o no?
 */
void print_factors_list(plist list_of_factors, bool is_not_window) {

  my_assert(list_of_factors!=NULL);

  if (is_not_window) {
	 DEBUG("List of the unique factors:");
  } else {
	 DEBUG("List of the unique genomic windows:");
  }

  int number_of_factor= 0;

  pfactor pf;
  plistit list_it;

  list_it=list_first(list_of_factors);

  INFO("Factors (%zd)", list_size(list_of_factors));

  while(listit_has_next(list_it)){
	 pf=listit_next(list_it);
	 number_of_factor=number_of_factor+1;
	 INFO("%d\t%d\t%d", number_of_factor, pf->GEN_start,pf->GEN_end);
  }
  listit_destroy(list_it);
}
#endif

/*
 * Procedura per confrontare sulla genomica un fattore (cfactor) con un altro (factor).
 * Sono ammessi due tipi di confronto: (1) di uguaglianza se is_not_window e' true,
 * (2) di contenimento se is_not_window e' false
 */
bool equal_factor(pfactor cfactor, pfactor factor, bool is_not_window)
{

	//fprintf(stdout, "CFR %d-%d\n", cfactor->GEN_start, cfactor->GEN_end);

  if(is_not_window){
	  if((cfactor->GEN_start==factor->GEN_start)&&
			  (cfactor->GEN_end==factor->GEN_end)){
		  return true;
	  }else{
		  return false;
	  }
  }
  else{
	  if((cfactor->GEN_start <= factor->GEN_start)&&
			  (cfactor->GEN_end >= factor->GEN_end)){
		  return true;
	  }else{
		  return false;
	  }
  }
}

//Verifica se il fattore in ingresso è presente nella lista di fattori unici list_of_factors

bool verify_factor(plist list_of_factors,pfactor factor)
{
  plistit list_it;
  pfactor current_factor;
  list_it=list_first(list_of_factors);
  bool already_exist=false;

  while(listit_has_next(list_it)){
	 current_factor=listit_next(list_it);
	 if(equal_factor(current_factor,factor,true)){
		already_exist= true;
	 }
  }
  listit_destroy(list_it);

  return already_exist;
}

/*
 * Aggiorna le finestre sulla genomica di list_of_factors con il fattore factor,
 * costruendo una lista ordinata per left crescente
 */
plist update_windows(plist list_of_factors, pfactor factor)
{
	pfactor copy;

	if(list_is_empty(list_of_factors)){
		copy=factor_create();
		copy->GEN_start=factor->GEN_start;
		copy->GEN_end=factor->GEN_end;
		copy->EST_start=-1;
		copy->EST_end=-1;
		list_add_to_head(list_of_factors, copy);

		return list_of_factors;
	}

	plistit list_it_for_start, list_it_for_end;
	pfactor start_current_window, end_current_window;
	bool stop_for_start, stop_for_end;
	bool start_inside, end_inside;

	list_it_for_start=list_first(list_of_factors);
	stop_for_start=false;
	start_inside=false;
	while(stop_for_start == false && listit_has_next(list_it_for_start)){
		 start_current_window=listit_next(list_it_for_start);
		 if(factor->GEN_start <= start_current_window->GEN_end){
			 stop_for_start=true;
			 //lo start di factor cade all'interno di current_factor (sulla genomica)
			 if(factor->GEN_start >= start_current_window->GEN_start)
				 start_inside=true;
			 //altrimenti lo start di factor cade prima dello start di current_factor (sulla genomica)
		 }
	}
	//All'uscita da questo ciclo:
		 //se stop_for_start e' true
			 //list_it_for_start punta alla finestra in cui cade lo start di factor se start_inside e' true,
			 //altrimenti lo start di factor cade prima dello start del fattore puntato da list_it_for_start
			 //e sicuramente dopo l'end della eventuale finestra precedente
		 //se stop_for_start e' false
			 //significa che lo start di factor cade dopo l'end dell'ultima finestra di list_of_factors


	//Se lo start di factor cade dopo l'end dell'ultima finestra
	if(stop_for_start == false){
		copy=factor_create();
		copy->GEN_start=factor->GEN_start;
		copy->GEN_end=factor->GEN_end;
		copy->EST_start=-1;
		copy->EST_end=-1;
		list_add_to_tail(list_of_factors, copy);
		listit_destroy(list_it_for_start);

		return list_of_factors;
	}

	//fprintf(stdout, "%d %d (%d-%d)\n", stop_for_start, start_inside, start_current_window->GEN_start, start_current_window->GEN_end);

	list_it_for_end=list_first(list_of_factors);
	stop_for_end=false;
	end_inside=false;
	int count_window=0;
	while(stop_for_end == false && listit_has_next(list_it_for_end)){
		 end_current_window=listit_next(list_it_for_end);
		 count_window++;
		 if(factor->GEN_end <= end_current_window->GEN_end){
			 stop_for_end=true;
			 //l'end di factor cade all'interno di current_factor (sulla genomica)
			 if(factor->GEN_end >= end_current_window->GEN_start){
				 end_inside=true;
			 }
			 //altrimenti l'end di factor cade prima dello start di current_factor (sulla genomica)
		 }
	}
	//All'uscita da questo ciclo:
		 //se stop_for_end e' true
			 //list_it_for_end punta alla finestra in cui cade l'end di factor se end_inside e' true,
			 //altrimenti l'end factor cade prima dello start del fattore puntato da list_it_for_end
			 //e sicuramente dopo l'end della eventuale finestra precedente
		 //se stop_for_end e' false
			 //significa che l'end di factor cade dopo l'end dell'ultima finestra di list_of_factors

	//Se l'end di factor cade prima dello start della prima finestra
	if(stop_for_end == true &&  end_inside == false && count_window == 1){
		copy=factor_create();
		copy->GEN_start=factor->GEN_start;
		copy->GEN_end=factor->GEN_end;
		copy->EST_start=-1;
		copy->EST_end=-1;
		list_add_to_head(list_of_factors, copy);
		listit_destroy(list_it_for_start);
		listit_destroy(list_it_for_end);

		return list_of_factors;
	}

	//fprintf(stdout, "%d %d (%d-%d)\n", stop_for_start, start_inside, start_current_window->GEN_start, start_current_window->GEN_end);

	if(start_inside == false && end_inside == false){
		//factor sta in mezzo a due finestre senza sovrapporsi ad alcuna di esse (va quindi aggiunto a list_of_factors)
		if(stop_for_end == true && start_current_window->GEN_start == end_current_window->GEN_start){
			//Aggiungo factor a list_of_factors appena prima di list_it_for_start
			copy=factor_create();
			copy->GEN_start=factor->GEN_start;
			copy->GEN_end=factor->GEN_end;
			copy->EST_start=-1;
			copy->EST_end=-1;
			list_add_before_iterator(list_it_for_start, list_of_factors, copy);
		}
		//factor ricopre una o piu' finestre consecutive estendendole a sx e a dx
		else{
			//Aggiorno la finestra in list_it_for_start con lo start e l'end di factor
			start_current_window->GEN_start=factor->GEN_start;
			start_current_window->GEN_end=factor->GEN_end;
			//Rimuovo le finestre iterando dalla successiva a list_it_for_start
			//a quella prima di list_it_for_end
			int end_curr_win_start=end_current_window->GEN_start;
			bool stop=false;
			while(listit_has_next(list_it_for_start) && stop == false){
				pfactor current_window=listit_next(list_it_for_start);
				if(current_window->GEN_start >= end_curr_win_start)
					stop=true;
				else
					list_remove_at_iterator(list_it_for_start, (delete_function)factor_destroy);
			}
		}
	}
	else{
		if(start_inside == true){
			if(end_inside == true){
				//factor ha lo start in una finestra e l'end in un'altra successiva
				if(start_current_window->GEN_start != end_current_window->GEN_start){
					//Aggiorno l'end della finestra in list_it_for_start con l'end della finestra in list_it_for_end
					start_current_window->GEN_end=end_current_window->GEN_end;

					//Rimuovo le finestre iterando dalla successiva a list_it_for_start
					//fino alla list_it_for_end (compresa)
					int end_curr_win_start=end_current_window->GEN_start;
					bool stop=false;
					while(stop == false && listit_has_next(list_it_for_start)){
						pfactor current_window=listit_next(list_it_for_start);
						/*list_remove_at_iterator(list_of_factors, list_it_for_start, (delete_function)factor_destroy);
						if(current_window->GEN_start == end_current_window->GEN_start)
							stop=true;*/
						if(current_window->GEN_start > end_curr_win_start)
							stop=true;
						else {
							list_remove_at_iterator(list_it_for_start, (delete_function)factor_destroy);
						}
					}
				}
			}
			//factor ha lo start in una finestra e l'end in mezzo a due finestre o dopo l'end
			//dell'ultima
			else{
				//Aggiorno l'end della finestra in list_it_for_start con l'end di factor
				start_current_window->GEN_end=factor->GEN_end;
				//Rimuovo le finestre iterando dalla successiva a list_it_for_start
				//a quella prima di list_it_for_end (rimuovo anche quella di list_it_for_end
				//se stop_for_end e' false (cioe' se l'end di factor e' oltre l'end
				//dell'ultima finestra)
				int end_curr_win_start=end_current_window->GEN_start;
				bool stop=false;
				while(listit_has_next(list_it_for_start) && stop == false){
					pfactor current_window=listit_next(list_it_for_start);
					if(stop_for_end == true && current_window->GEN_start >= end_curr_win_start)
						stop=true;
					else
						list_remove_at_iterator(list_it_for_start, (delete_function)factor_destroy);
				}
			}
		}
		//factor ha lo start in mezzo a due finestre (o prima dello start della prima)
		//e l'end in una finestra
		else{
			//Aggiorno lo start della finestra in list_it_for_start con lo start di factor
			start_current_window->GEN_start=factor->GEN_start;

			//Aggiorno l'end della finestra in list_it_for_start con l'end di list_it_for_end
			start_current_window->GEN_end=end_current_window->GEN_end;
			//Rimuovo le finestre iterando dalla successiva a list_it_for_start
			//fino alla list_it_for_end (compresa)
			int end_curr_win_start=end_current_window->GEN_start;
			bool stop=false;
			while(listit_has_next(list_it_for_start) && stop == false){
				pfactor current_window=listit_next(list_it_for_start);
				/*list_remove_at_iterator(list_of_factors, list_it_for_start, (delete_function)factor_destroy);
				if(current_window->GEN_start == end_current_window->GEN_start)
					stop=true;*/
				if(current_window->GEN_start > end_curr_win_start)
						stop=true;
				else
					list_remove_at_iterator(list_it_for_start, (delete_function)factor_destroy);
			}
		}
	}

	listit_destroy(list_it_for_start);
	listit_destroy(list_it_for_end);

	return list_of_factors;
}

/*
 * Funzione che restituisce la posizione del fattore all'interno della lista di fattori unici.
 * Se is_not_window e' true la ricerca e' precisa, altrimenti la ricerca e' di contenimento (da usare
 * nel caso in cui i fattori unici siano finestre di genomica (vedi procedura windows_list_create())
 */
int get_factor_position(pfactor factor,plist factors, bool is_not_window)
{

  my_assert((factor!=NULL)&&(factors!=NULL));


  pfactor pf;
  plistit plist_it_factor;
  int pos=-1;

  //fprintf(stdout, "\t\tSearch for %d-%d\n", factor->GEN_start, factor->GEN_end);

  plist_it_factor=list_first(factors);

  while(listit_has_next(plist_it_factor)) {

	 pos=pos+1;
	 pf=listit_next(plist_it_factor);

	 //fprintf(stdout, "\t\t\t...%d-%d\n", pf->GEN_start, pf->GEN_end);

	 if(equal_factor(pf,factor,is_not_window)){
		break;
	 }
  }
  listit_destroy(plist_it_factor);

  //fprintf(stdout, "\t\t\tPos=%d\n", pos);

  return pos;
}


//Creazione della lista di fattori unici

plist factors_list_create(plist ests_factorizations){

  my_assert(ests_factorizations!=NULL);

  DEBUG("Creating the list of unique factors...");


  pEST p;
  plistit plist_it_id, plist_it_f;
  pfactorization pfact;
  pfactor pf;
  plistit plist_it_factor;

  plist list_of_factors=list_create();
  plist_it_id=list_first(ests_factorizations);

  while(listit_has_next(plist_it_id)){
	 p= listit_next(plist_it_id);
	 plist_it_f=list_first(p->factorizations);

	 while(listit_has_next(plist_it_f)){
		pfact=listit_next(plist_it_f);
		plist_it_factor=list_first(pfact);

		while(listit_has_next(plist_it_factor)) {
		  pf=listit_next(plist_it_factor);
		  if(!verify_factor(list_of_factors,pf))list_add_to_tail(list_of_factors,pf);
		}
		listit_destroy(plist_it_factor);
	 }
	 listit_destroy(plist_it_f);
  }
  listit_destroy(plist_it_id);

  return list_of_factors;
}

/*
 *
 */
plist windows_list_create(plist ests_factorizations){

  my_assert(ests_factorizations!=NULL);

  DEBUG("Creating the list of genomic windows that will be considered as unique factors...");

  pEST p;
  plistit plist_it_id, plist_it_f;
  pfactorization pfact;
  pfactor pf;
  plistit plist_it_factor;

  plist list_of_factors=list_create();

  plist_it_id=list_first(ests_factorizations);

  while(listit_has_next(plist_it_id)){
	 p= listit_next(plist_it_id);
	 plist_it_f=list_first(p->factorizations);

	 while(listit_has_next(plist_it_f)){
		pfact=listit_next(plist_it_f);
		plist_it_factor=list_first(pfact);

		//fprintf(stdout, "EST %s\n", p->info->EST_id);

		while(listit_has_next(plist_it_factor)) {
		  pf=listit_next(plist_it_factor);

		  list_of_factors=update_windows(list_of_factors,pf);

		 /*fprintf(stdout, "add %d-%d\n", pf->GEN_start, pf->GEN_end);
		  fprintf(stdout, "Windows: ");
		  plistit t_it=list_first(list_of_factors);
		  while(listit_has_next(t_it)){
			  pfactor win=(pfactor)listit_next(t_it);
			  fprintf(stdout, "%d-%d ", win->GEN_start, win->GEN_end);
		  }
		  fprintf(stdout, "\n");
		  listit_destroy(t_it);*/

		}
		listit_destroy(plist_it_factor);
	 }
	 listit_destroy(plist_it_f);
  }

  listit_destroy(plist_it_id);

  return list_of_factors;
}

//Funzione che legge la fattorizzazione dell'est nel campo factorizations ed esprime tale
//fattorizzazione come vettore binario per l'inserimento di tale vettore nella lista bin_factorizations di
//quell'est.
/*
 * is_not_window deve essere true, se il confronto e' di uguaglianza (cioe' factors e' stata costruita con
 * factors_list_create(), e deve essere false se factors e' stata costruita con
 * windows_list_create()
 */
void add_factoriz(pEST est,plist factorization,plist factors, bool is_not_window)
{

  my_assert((est!=NULL)&&(factorization!=NULL)&&(factors!=NULL));

  pfactor pf;
  plistit plist_it_factor;

  int factor_pos;
  pbit_vect bv=BV_create(list_size(factors));

  if(est->bin_factorizations==NULL){
	 est->bin_factorizations=list_create();
  }

  plist_it_factor=list_first(factorization);

	 while(listit_has_next(plist_it_factor)) {

		pf=listit_next(plist_it_factor);
		  factor_pos=get_factor_position(pf,factors,is_not_window);

		  //fprintf(stdout, "\t\tFactor position %d\n", factor_pos);

		  BV_set(bv,factor_pos,true);
	 }
	 list_add_to_tail(est->bin_factorizations,bv);
	 listit_destroy(plist_it_factor);
}


/*
 * is_not_window deve essere true, se il confronto e' di uguaglianza (cioe' factors e' stata costruita con
 * factors_list_create(), e deve essere false se factors e' stata costruita con
 * windows_list_create()
 */
static void
add_EST(pEST p,plist factors, bool is_not_window)
{
  my_assert((p!=NULL)&&(factors!=NULL));

  plistit list_it;
  plist factorization;

  //fprintf(stdout, "EST %s\n", p->info->EST_id);

  list_it=list_first(p->factorizations);

  while(listit_has_next(list_it)){

	 factorization=listit_next(list_it);

	 /*plistit temp=list_first(factorization);
	  fprintf(stdout, "\tFactorization\n");
	  while(listit_has_next(temp)){
		  pfactor f_t=(pfactor)listit_next(temp);
		  fprintf(stdout, "\t\t%d-%d\n", f_t->GEN_start, f_t->GEN_end);
	  }
	  listit_destroy(temp);*/

	 add_factoriz(p,factorization,factors,is_not_window);
  }
  listit_destroy(list_it);
}

//Stampa la matrice colorata

#if defined (LOG_MSG) && (LOG_LEVEL_DEBUG <= LOG_THRESHOLD)
void color_matrix_print(plist color_matrix)
{
  my_assert(color_matrix!=NULL);

  DEBUG("Colored matrix. Gives the ID of each sequence and the releated binary factorization.");

  plistit list_it_est, list_it_factoriz;
  pEST est;

  pbit_vect bv;
  list_it_est=list_first(color_matrix);

  while(listit_has_next(list_it_est)){
	 est=listit_next(list_it_est);

	 list_it_factoriz=list_first(est->bin_factorizations);
	 printf("id:%s\n", est->info->EST_id);

	 while(listit_has_next(list_it_factoriz)){
		bv=listit_next(list_it_factoriz);
		BV_print(bv);
	 }
	  listit_destroy(list_it_factoriz);
  }
  listit_destroy(list_it_est);
}
#endif

//creazione della matrice colorata. La funzione prende in ingresso la lista di liste
//rappresentante le fattorizzazioni di tutti gli est e aggiunge alla lista bin_factorizations
//di ogni est le fattorizzazioni viste come vettori binari
/*
 * is_not_window deve essere true, se il confronto e' di uguaglianza (cioe' factors e' da costruire con
 * factors_list_create(), e deve essere false se factors e' da costruire con
 * windows_list_create()
 */
plist color_matrix_create(plist ests_factorizations, bool is_not_window)
{

  my_assert(ests_factorizations!=NULL);

  pEST p;
  plistit plist_it_id;

  plist factors;
  if (is_not_window) {
	 factors=factors_list_create(ests_factorizations);
  } else {
//Da sostituire
	 factors=windows_list_create(ests_factorizations);
  }

  INFO("Total factors (before simplification): %zu", list_size(factors));
  DEBUG("Creation of the list of unique factors successful!");
  plist_it_id=list_first(ests_factorizations);

  while(listit_has_next(plist_it_id)){
	 p= listit_next(plist_it_id);
	 add_EST(p,factors,is_not_window);
  }
  listit_destroy(plist_it_id);

  return factors;

}





