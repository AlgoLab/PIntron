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
#include "io-gen-ests.h"
#include"io-meg.h"
#include"io-multifasta.h"
#include"log.h"
#include <string.h>

#define LEN_BUFFER 100

static void write_EST_MEG(FILE* source, pEST_MEG EST_MEG){

       my_assert(source != NULL);
       my_assert(EST_MEG != NULL);

       fprintf(source,"#EST#\n");
       DEBUG("The pointer to est info is: %p",(void*)EST_MEG->est);
       write_single_EST_info(source,EST_MEG->est);
       fprintf(source,"#\\#\n");
       fprintf(source,"#MEG#\n");
       meg_write(source,EST_MEG->meg);
       fprintf(source,"#\\#\n");

       }

pGEN_ESTS GEN_ESTS_read(FILE* source){

 my_assert(source != NULL);

 pGEN_ESTS GEN_ESTS = GEN_ESTS_create();
 size_t size = LEN_BUFFER;
 char* buffer = c_palloc(size);
 my_getline(&buffer,&size,source);
 DEBUG("Read the line >%s<",buffer);
      if(strcmp(buffer,"#GENOMICA#")==0){
                  pEST_info tmp = read_single_EST_info(source);
                  GEN_ESTS->gen = tmp;
                  }
      else{
            DEBUG("The input file must start with genomic's data");
            }

 pEST_MEG EST_MEG;

 while(!feof(source)){

                      do {
                          my_getline(&buffer,&size,source);
                          DEBUG("Read the line >%s<",buffer);
                      } while(buffer[0]=='\0' && !feof(source));
                      if (!feof(source)) {
                      if(strcmp(buffer,"#EST#")==0){
                                                 EST_MEG = EST_MEG_create();
                                                 DEBUG("Creating an est");
                                                 EST_MEG->est = read_single_EST_info(source);
                      } else if(strcmp(buffer,"#MEG#")==0){
                                               DEBUG("Creating the MEG");
                                               EST_MEG->meg = meg_read(source);
                                               list_add_to_tail(GEN_ESTS->ests,EST_MEG);
                                               DEBUG("MEG created");
                                               }
                      }
                      }

 DEBUG("GEN_ESTS structure has been created");
 pfree(buffer);
 return GEN_ESTS;

}

void GEN_ESTS_write(FILE* source, pGEN_ESTS data){

     my_assert(source != NULL);
     my_assert(data != NULL);
     my_assert(data->gen != NULL);
     my_assert(data->ests != NULL);
     fprintf(source,"#GENOMICA#\n");
     write_single_EST_info(source,data->gen);
     fprintf(source,"#\\#\n");
     plistit it = list_first(data->ests);
     while(listit_has_next(it)){
                                DEBUG("Writing an EST_MEG");
                                write_EST_MEG(source,listit_next(it));
                                }
     listit_destroy(it);
     }

#undef LEN_BUFFER
