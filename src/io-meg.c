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
#include<stdlib.h>
#include "io-meg.h"
#include "util.h"
#include <string.h>

#include "log.h"
#define LEN_BUFFER 1001

static ppairing
get_pairing(char* pairing_as_string, int id_pairing) {
  my_assert(pairing_as_string != NULL);

  ppairing pairing = pairing_create();
  int p=0, t=0, l=0;
  int n_par_read= sscanf(pairing_as_string, "(%d,%d,%d)",
								 &p, &t, &l);
  if (n_par_read!=3) {
	 FATAL("Read a pairing with %d parameter instead of 3",
			 n_par_read);
	 fail();
  }
  DEBUG("Read pairing (%d,%d,%d)", p, t, l);

  pairing->p = p;
  pairing->t = t;
  pairing->l = l;
  pairing->id = id_pairing;
  return pairing;
}

pext_array meg_read(FILE * source){

  my_assert(source != NULL);
  int flag = 0;
  int dim = 0;
  int i;
  int id_pairing;
  id_pairing = 0;
  size_t size = LEN_BUFFER;
  char* buffer = c_palloc(size);
  pext_array array = EA_create();
  pext_array graph;
  DEBUG("Starting the read of the vertex-set");
  while( flag != 1 ){

	 do {
		my_getline(&buffer,&size,source);
		TRACE("Read line >%s<", buffer);
	 } while (buffer[0]=='\0');

	 if((strcmp(buffer,"#adj#")) == 0){
		DEBUG("Vertex-set read!");
		flag = 1;
	 } else {
		TRACE("Got a new pairing");
		ppairing pairing = get_pairing(buffer, id_pairing);
		EA_insert(array,pairing);
		id_pairing++;
		if (pairing->p>dim) {
		  dim= pairing->p;
		}
	 }
  }

  TRACE("Creating a vector called 'graph'of %d element. Each element is an empty list", dim+1);
  graph = EA_create();
  for(i = 0; i <= dim; i++) {
	 EA_insert(graph,list_create());
  }
  TRACE("Filling graph's lists");
  const int max_id_pairing= EA_size(array);
  for(i = 0; i < max_id_pairing; i++){
	 ppairing pairing = EA_get(array,i);
	 list_add_to_tail( (plist) EA_get(graph, pairing->p),
							 pairing );
  }

  DEBUG("Starting the read of the edge-set");

  while (!feof(source)) {
	 do {
		my_getline(&buffer,&size,source);
		TRACE("Read line >%s<",buffer);
	 } while (buffer[0]=='\0' && !feof(source));

     if(strcmp(buffer,"#\\#") == 0) break;

	 if (!feof(source)) {

		int vs, vt;
		int num_p= sscanf(buffer, "%d-%d", &vs, &vt);
		if (num_p!=2) {
		  fail();
		}
		if (vs<0 || vs>=max_id_pairing) {
		  FATAL("The id %d of the source vertex is not valid.", vs);
		  fail();
		}
		if (vt<0 || vt>=max_id_pairing) {
		  FATAL("The id %d of the target vertex is not valid.", vt);
		  fail();
		}
		ppairing ps= (ppairing)EA_get(array, vs);
		ppairing pt= (ppairing)EA_get(array, vt);
		list_add_to_tail(ps->adjs, pt);
		TRACE("Added edge (%d, %d, %d)->(%d, %d, %d)",
				ps->p, ps->t, ps->l,
				pt->p, pt->t, pt->l);
	 }
  }
  pfree(buffer);
  EA_destroy(array,noop_free);
  return graph;
}//Endread_meg


void meg_write(FILE* dest, pext_array graph ) {
  my_assert( dest != NULL);
  my_assert( graph != NULL);
  pext_arrayit array_it = EA_begin(graph);
  plist list;
  plistit list_it, adjs_it;
  ppairing pairing;
  int index = 0;
  DEBUG("Writing the MEG");
  DEBUG("Vertex-set");
  while (ext_arrayit_has_next(array_it)){
	 list = ext_arrayit_next(array_it);
	 list_it = list_first(list);
	 while(listit_has_next(list_it)){
		pairing= (ppairing)listit_next(list_it);
		fprintf(dest, "(%d,%d,%d)\n",
				  pairing->p,
				  pairing->t,
				  pairing->l);
		pairing->id= index;
		++index;
	 }
	 listit_destroy(list_it);
  }
  ext_arrayit_destroy(array_it);
  DEBUG("Edge-set");
  fprintf(dest,"#adj#\n");
  array_it = EA_begin(graph);
  while(ext_arrayit_has_next(array_it)){
	 list = ext_arrayit_next(array_it);
	 list_it = list_first(list);
	 while(listit_has_next(list_it)){
		pairing= (ppairing)listit_next(list_it);
		adjs_it= list_first(pairing->adjs);
		while(listit_has_next(adjs_it)){
		  ppairing aux = (ppairing)listit_next(adjs_it);
		  fprintf(dest,"%d-%d\n", pairing->id, aux->id);
		}
		listit_destroy(adjs_it);
	 }
	 listit_destroy(list_it);
  }
  ext_arrayit_destroy(array_it);
}//end meg_write


#undef LEN_BUFFER

