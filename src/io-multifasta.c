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
//This file provides an IO for multi-fasta format

#include "io-multifasta.h"
#include "log.h"
#include "util.h"
#include <string.h>
#define LEN_BUFFER 10000000
#define LEN_STRAND_ARRAY 10
#define LEN_ABSCOORD_ARRAY 100

static char get_Complement(char car){
  if(car == 'A')
	 return 'T';
  if(car == 'a')
	 return 't';
  if(car == 'T')
	 return 'A';
  if(car == 't')
	 return 'a';
  if(car == 'C')
	 return 'G';
  if(car == 'c')
	 return 'g';
  if(car == 'G')
	 return 'C';
  if(car == 'g')
	 return 'c';
  if(car == 'A')
	 return 'T';

  if(car == 'R')
	 return 'Y';
  if(car == 'r')
	 return 'y';
  if(car == 'Y')
	 return 'R';
  if(car == 'y')
	 return 'r';
  if(car == 'M')
	 return 'K';
  if(car == 'm')
	 return 'k';
  if(car == 'K')
	 return 'M';
  if(car == 'k')
	 return 'm';
  if(car == 'B')
	 return 'V';
  if(car == 'b')
	 return 'v';
  if(car == 'V')
	 return 'B';
  if(car == 'v')
	 return 'b';
  if(car == 'D')
	 return 'H';
  if(car == 'd')
	 return 'h';
  if(car == 'H')
	 return 'D';
  if(car == 'h')
	 return 'd';

  return car;
}

static char* getData(char *buffer, size_t n, FILE *f_in, int* lenght){
  my_assert(buffer != NULL);
  my_assert( n > 0);
  my_assert(f_in != NULL);
  my_assert(lenght != NULL);

  char* res;
  plist l= list_create();
  int len= 0;
  int leng= *lenght;

  while(buffer[0] != '>' && leng != EOF && strcmp(buffer,"#\\#")){
	 if (leng>0) {
		char* readline = alloc_and_copy(buffer);
		FINETRACE("the read line is |%s|", readline);
		FINETRACE("the length of the read line is %d", leng);
		list_add_to_tail(l, readline);
		len = len + leng;
	 }
	 leng = my_getline(&buffer,&n,f_in);
	 *lenght = leng;
  }//end-While

  plistit li = list_first(l);
  res= c_palloc(len+1);

  int ind = 0;
  while(listit_has_next(li)){
	 char* tmp = (char*) listit_next(li);
	 memcpy(res + ind, tmp, strlen(tmp));
	 ind = ind + strlen(tmp);
  }
  res[ind] = '\0';
  list_destroy(l, (delete_function)pfree_function);
  listit_destroy(li);
  FINETRACE("the concatenated string is |%s|", res);
  return res;
}//end-getData


plist read_multifasta(FILE* inputFile){
  my_assert(inputFile != NULL);
  plist data_list = list_create();
  size_t size;
  int length;
  size= LEN_BUFFER;
  char* string = c_palloc(size);
  char* tmp;
  length= my_getline(&string, &size, inputFile);
  while(length>0 && string[length]<32) {
	 string[length]= '\0';
	 length--;
  }
  while(length != EOF){
	 if(string[0] == '>'){
		pEST_info obj = EST_info_create();
		tmp = substring(1, string);
		obj->EST_id = tmp;
		length= my_getline(&string, &size, inputFile);
		tmp= getData(string, size, inputFile, &length);
		obj->EST_seq = tmp;
		char *original_seq=substring(0, obj->EST_seq);
		obj->original_EST_seq = original_seq;
		list_add_to_tail(data_list,obj);
	 } else {
		length = my_getline(&string, &size, inputFile);
		while(length>0 && string[length]<32) {
		  string[length]= '\0';
		  length--;
		}
	 }
  }//end-while
  pfree(string);
  return data_list;
}


void write_multifasta(plist inputList, FILE* output_file){
  my_assert(inputList != NULL);
  my_assert(output_file != NULL);
  plistit li = list_first(inputList);
  my_assert(li != NULL);

  while(listit_has_next(li)){
	 pEST_info tmp = (pEST_info) listit_next(li);
	 fprintf(output_file,">%s\n",tmp->EST_id );
	 fprintf(output_file,"%s\n",tmp->original_EST_seq);
  }//End-While

  listit_destroy(li);

}//End-print_list_in_file


void write_multifasta_output(pEST_info gen, pEST est, FILE* output_file, char retain_externals){
  my_assert(est != NULL);
  my_assert(output_file != NULL);

  if(!(est->factorizations == NULL || list_is_empty(est->factorizations))){
	 plistit f_it=list_first(est->factorizations);
	 pboollistit polya_it=boollist_first(est->polyA_signals);
	 pboollistit polyadenil_it=boollist_first(est->polyadenil_signals);

	 while(listit_has_next(f_it)){
		my_assert(boollistit_has_next(polya_it));
		my_assert(boollistit_has_next(polyadenil_it));

		plist factorization=(plist)listit_next(f_it);
		bool polya=(bool) boollistit_next(polya_it);
		bool polyadenil=(bool) boollistit_next(polyadenil_it);

		if(retain_externals || (list_size(factorization) > 2 || (list_size(factorization) == 2 && est->info->suff_polyA_length != -1))){
		  fprintf(output_file,">%s\n",est->info->EST_id);

		  if(!retain_externals){
			  polya=0;
			  polyadenil=0;
		  }
		  fprintf(output_file,"#polya=%d\n#polyad=%d\n", polya, polyadenil);

		  plistit factor_it;
		  factor_it=list_first(factorization);

		  unsigned int counter=1;
		  unsigned int l_index=(retain_externals == 0)?(1):(0);
		  unsigned int r_index=(retain_externals == 0)?((est->info->suff_polyA_length == -1)?(list_size(factorization)):(list_size(factorization)+1)):(list_size(factorization)+1);

		  while(listit_has_next(factor_it)){
                    pfactor factor=(pfactor)listit_next(factor_it);
                    if(counter > l_index && counter < r_index){
                      fprintf(output_file, "%d %d %d %d %.*s %.*s\n",
                              factor->EST_start + 1, factor->EST_end + 1,
                              gen->pref_N_length + factor->GEN_start + 1, gen->pref_N_length + factor->GEN_end + 1,
                              factor->EST_end + 1 - factor->EST_start,
                              est->info->original_EST_seq + factor->EST_start,
                              factor->GEN_end + 1 - factor->GEN_start,
                              gen->original_EST_seq + gen->pref_N_length + factor->GEN_start);
                    }
                    counter++;
		  }
		  listit_destroy(factor_it);
		}
	 }

	 listit_destroy(f_it);
	 boollistit_destroy(polya_it);
	 boollistit_destroy(polyadenil_it);
 }
}

pEST_info read_single_EST_info(FILE* source){

  my_assert(source != NULL);
  size_t size = LEN_BUFFER;
  int length;
  char* buffer = c_palloc(size);
  char* tmp;
  my_getline(&buffer,&size,source);

  pEST_info obj = EST_info_create();

  while( (!(feof(source))) && (strcmp(buffer,"#\\#") != 0)){
	 tmp = substring(1, buffer);
	 obj->EST_id = tmp;
	 length = my_getline(&buffer, &size, source);
	 tmp = getData(buffer, size, source, &length);
	 obj->EST_seq = tmp;
	 char *original_seq=substring(0, obj->EST_seq);
	 obj->original_EST_seq = original_seq;
  }
  DEBUG("ID is %s",obj->EST_id);
  DEBUG("The sequence created is >%s<",obj->original_EST_seq);
  DEBUG("The pointer's value is: %p",(void*) obj);
  pfree(buffer);
  return obj;
}

void write_single_EST_info(FILE* source, pEST_info EST_info){

  my_assert(source != NULL);
  my_assert(EST_info != NULL);

  fprintf(source,">%s\n",EST_info->EST_id);
  fprintf(source,"%s\n",EST_info->original_EST_seq);
}

void set_EST_GB_identification(pEST_info EST_info){
  my_assert(EST_info != NULL);
  my_assert(EST_info->EST_id != NULL);

  char* p=strstr(EST_info->EST_id, "/gb=");
  if(p == NULL){
	 p=strstr(EST_info->EST_id, "/GB=");
  }
  if(p != NULL){
	 p+=4;
// Calculate GB len
	 size_t len= 0;
	 while ((*(p+len) != ' ') &&
			  (*(p+len) != '/') &&
			  (*(p+len) != '\0')) {
		++len;
	 }
	 EST_info->EST_gb=c_palloc(len+1);
	 strncpy(EST_info->EST_gb, p, len);
	 EST_info->EST_gb[len]='\0';
	 DEBUG("GB ==> %s", EST_info->EST_gb);
  }
  else{
	 DEBUG("GB ==> unknown");
  }
}

static bool
parse_genomic_header_standard(pEST_info gen) {
// Parse the genomic header in the following format
// >chrNN:NNNNN:NNNN:+-1
// return true if successful
  char* header= strdup(gen->EST_id);
  char* bakheader= header;
  if (header == NULL) return false;

  char* chr= strsep(&header, ":");
  if (chr == NULL) {
	 free(bakheader);
	 return false;
  }
  chr= strdup(chr);

  char* start= strsep(&header, ":");
  if (start == NULL) {
	 free(bakheader);
	 free(chr);
	 return false;
  }

  char* end= strsep(&header, ":");
  if (end == NULL) {
	 free(bakheader);
	 free(chr);
	 return false;
  }

  char* strand= strsep(&header, ":");
  if (strand == NULL) {
	 free(bakheader);
	 free(chr);
	 return false;
  }
  strand= strdup(strand);

  char* last= strsep(&header, ":");
  if (last != NULL) {
	 free(bakheader);
	 free(chr);
	 free(strand);
	 return false;
  }

// Check that everything is correct
  int abs_start= atoi(start);
  int abs_end= atoi(end);
  int EST_strand= atoi(strand);
  if ((abs_start < 1) || (abs_end < 1) ||
		((EST_strand != -1) && (EST_strand != +1))) {
	 free(bakheader);
	 free(chr);
	 free(strand);
	 return false;
  }

  INFO("Genomic chromosome: %s", chr);
  INFO("Genomic abs start:  %d", abs_start);
  INFO("Genomic abs end:    %d", abs_end);
  INFO("Genomic strand:     %s", strand);

  gen->EST_chr= chr;
  gen->abs_start= abs_start;
  gen->abs_end= abs_end;
  gen->EST_strand= EST_strand;
  gen->EST_strand_as_read= strand;

  free(bakheader);
  return true;
}

static bool
genomic_header_defaults(pEST_info gen) {
// Parse the genomic header in the following format
// >chrNN:NNNNN:NNNN:+-1
// return true if successful
  char* chr= strdup("unknown");
  char* strand= strdup("+1");

  int abs_start= 1;
  int abs_end= strlen(gen->EST_seq);
  int EST_strand= atoi(strand);

  ERROR("The header of the genomic file is not in the correct standard!");
  ERROR("This may lead to prediction errors.");
  ERROR("Please correct the header if possible.");
  ERROR("Guessed values:");
  INFO("Genomic chromosome: %s", chr);
  INFO("Genomic abs start:  %d", abs_start);
  INFO("Genomic abs end:    %d", abs_end);
  INFO("Genomic strand:     %s", strand);

  gen->EST_chr= chr;
  gen->abs_start= abs_start;
  gen->abs_end= abs_end;
  gen->EST_strand= EST_strand;
  gen->EST_strand_as_read= strand;

  return true;
}


void parse_genomic_header(pEST_info gen) {
  my_assert(gen != NULL);
  my_assert(gen->EST_id != NULL);
  bool ris= parse_genomic_header_standard(gen);
  if (!ris) {
// Compute resonable defaults
	 ris= genomic_header_defaults(gen);
  }
  my_assert(ris);
  if (!ris) {
	 ERROR("Something went wrong when parsing the genomic header.");
	 fail();
  }
}


void set_EST_Strand_and_RC(pEST_info EST_info, pEST_info gen){
  my_assert(EST_info != NULL);
  my_assert(EST_info->EST_id != NULL);
  my_assert(gen->EST_strand_as_read != NULL);

  EST_info->EST_strand_as_read= c_palloc(LEN_STRAND_ARRAY+1);

  if (EST_info->EST_gb != NULL && EST_info->EST_gb[0] == 'N' && EST_info->EST_gb[1] == 'M') {
    strcpy(EST_info->EST_strand_as_read, "1");
    EST_info->EST_strand= 1;
    DEBUG("STRAND RefSeq==> %s", EST_info->EST_strand_as_read);
  } else {
    char* p=strstr(EST_info->EST_id, "/clone_end=");
    if (p == NULL) {
      p=strstr(EST_info->EST_id, "/CLONE_END=");
    }
    if (p != NULL) {
      p+=11;
      int i=0;
      while ((i<LEN_STRAND_ARRAY) &&
            (*p != '\0') &&
            (*p != '\'')){
        EST_info->EST_strand_as_read[i]= *p;
        p++;
        i++;
      }
      EST_info->EST_strand_as_read[i]='\0';
      DEBUG("STRAND as read ==> %s", EST_info->EST_strand_as_read);

      if (!strcmp(EST_info->EST_strand_as_read, "3")) {
        EST_info->EST_strand= 1;
      } else if (!strcmp(EST_info->EST_strand_as_read, "5")) {
        EST_info->EST_strand= -1;
      } else {
        WARN("Failed to interpret the strand. Setting to '1' by default. Read: '%s'", EST_info->EST_strand_as_read);
        EST_info->EST_strand= 1;
      }
    } else {
      WARN("Failed to find a strand. Setting to '1' by default.");
      EST_info->EST_strand= 1;
      EST_info->EST_strand_as_read[0]= '\0';
    }
  }

  DEBUG("Interpreted strand ==> %d", EST_info->EST_strand);

  if (EST_info->EST_strand == -1) {
    reverse_and_complement(EST_info);
  }
}

void reverse_and_complement(pEST_info EST_info) {
//Faccio il reverse and complement della sequenza
  size_t right= MAX(strlen(EST_info->EST_seq), 1)-1;
  DEBUG("Reverse and complement of EST sequence (%zubp)", strlen(EST_info->EST_seq));
  size_t left= 0;

  while (left<=right) {
	 char new_right=get_Complement(EST_info->EST_seq[left]);
	 char new_left=get_Complement(EST_info->EST_seq[right]);
	 EST_info->EST_seq[right]= new_right;
	 EST_info->EST_seq[left]= new_left;
	 EST_info->original_EST_seq[right]= new_right;
	 EST_info->original_EST_seq[left]= new_left;
	 ++left;
	 --right;
  }
  FINETRACE("==>%s", EST_info->EST_seq);
}


static void old_polyAT_substitution(pEST_info EST_info) {
  my_assert(EST_info != NULL);
  my_assert(EST_info->EST_seq != NULL);

  char c;
  size_t i;
  size_t est_len= strlen(EST_info->EST_seq);

  my_assert(est_len>0);

  EST_info->pref_polyA_length=-1;
  EST_info->suff_polyA_length=-1;
  EST_info->pref_polyT_length=-1;
  EST_info->suff_polyT_length=-1;

// Check the beginning
// Look for a sequence of at least
// POLYA_MIN_LEN characters, otherwise it does not consider
// it as a polyAT
  size_t count_A= 0;
  size_t count_T= 0;
  size_t last_A= 0;
  size_t last_T= 0;
  size_t last_A_count= 0;
  size_t last_T_count= 0;
  for (i= 0;
		 i<_POLYA_MIN_LEN &&
			i<est_len;
		 ++i) {
	 if (EST_info->EST_seq[i]=='A') {
		++count_A;
		last_A= i;
		last_A_count= count_A;
	 }
	 if (EST_info->EST_seq[i]=='T') {
		++count_T;
		last_T= i;
		last_T_count= count_T;
	 }
  }
  DEBUG("Found %zu A and %zu T", count_A, count_T);
  while ((i<est_len) &&
			((count_A>=(_POLYA_MIN_FRACTION*i)) ||
			 (count_T>=(_POLYA_MIN_FRACTION*i)))) {
	 if (EST_info->EST_seq[i]=='A') {
		++count_A;
		last_A= i;
		last_A_count= count_A;
	 }
	 if (EST_info->EST_seq[i]=='T') {
		++count_T;
		last_T= i;
		last_T_count= count_T;
	 }
	 ++i;
  }
  last_A= (last_A<_POLYA_MIN_LEN-1)?_POLYA_MIN_LEN-1:last_A;
  last_T= (last_T<_POLYA_MIN_LEN-1)?_POLYA_MIN_LEN-1:last_T;
  if ((last_A_count>=(_POLYA_MIN_FRACTION*(last_A+1))) ||
		(last_T_count>=(_POLYA_MIN_FRACTION*(last_T+1)))) {
	 if ((((double)last_A_count)/(last_A+1))>=(((double)last_T_count)/(last_T+1)))
		c= 'A';
	 else
		c= 'T';
// A polyA/T has been found
	 char sc= (c=='A')?_POLYA_CHR:_POLYT_CHR;
	 size_t mlen= ((c=='A')?last_A+1:last_T+1);
	 INFO("Found a %zubp long initial poly%c.", mlen, c);
	 for (i= 0; i<mlen; ++i)
		EST_info->EST_seq[i]= sc;
	 if(c=='A'){
		EST_info->pref_polyA_length=mlen;
	 } else {
		EST_info->pref_polyT_length=mlen;
	 }
  }

// Check the end
  count_A= 0;
  count_T= 0;
  last_A= 0;
  last_T= 0;
  last_A_count= 0;
  last_T_count= 0;
  for (i= 0;
		 i<_POLYA_MIN_LEN &&
			i<est_len;
		 ++i) {
	 if (EST_info->EST_seq[est_len-i-1]=='A') {
		++count_A;
		last_A= i;
		last_A_count= count_A;
	 }
	 if (EST_info->EST_seq[est_len-i-1]=='T') {
		++count_T;
		last_T= i;
		last_T_count= count_T;
	 }
  }
  DEBUG("Found %zu A and %zu T", count_A, count_T);
  while ((i<est_len) &&
			((count_A>=(_POLYA_MIN_FRACTION*i)) ||
			 (count_T>=(_POLYA_MIN_FRACTION*i)))) {
	 if (EST_info->EST_seq[est_len-i-1]=='A') {
		++count_A;
		last_A= i;
		last_A_count= count_A;
	 }
	 if (EST_info->EST_seq[est_len-i-1]=='T') {
		++count_T;
		last_T= i;
		last_T_count= count_T;
	 }
	 ++i;
  }
  last_A= (last_A<_POLYA_MIN_LEN-1)?_POLYA_MIN_LEN-1:last_A;
  last_T= (last_T<_POLYA_MIN_LEN-1)?_POLYA_MIN_LEN-1:last_T;
  if ((last_A_count>=(_POLYA_MIN_FRACTION*(last_A+1))) ||
		(last_T_count>=(_POLYA_MIN_FRACTION*(last_T+1)))) {
	 if ((((double)last_A_count)/(last_A+1))>=(((double)last_T_count)/(last_T+1)))
		c= 'A';
	 else
		c= 'T';
// A polyA/T has been found
	 char sc= (c=='A')?_POLYA_CHR:_POLYT_CHR;
	 size_t mlen= ((c=='A')?last_A+1:last_T+1);
	 INFO("Found a %zubp long final poly%c.", mlen, c);
	 for (i= 0; i<mlen; ++i)
		EST_info->EST_seq[est_len-i-1]= sc;
	 if(c=='A'){
		EST_info->suff_polyA_length=mlen;
	 } else {
		EST_info->suff_polyT_length=mlen;
	 }
  }
}

void polyAT_substitution(pEST_info EST_info) {
  my_assert(EST_info != NULL);
  my_assert(EST_info->EST_seq != NULL);

  char c;
  size_t i;
  const size_t est_len= strlen(EST_info->EST_seq);

  my_assert(est_len>0);

  EST_info->pref_polyA_length=-1;
  EST_info->suff_polyA_length=-1;
  EST_info->pref_polyT_length=-1;
  EST_info->suff_polyT_length=-1;

  if (est_len<_POLYA_MIN_LEN) {
	 return;
  }

// Check the beginning
// Look for a sequence of at least
// POLYA_MIN_LEN characters, otherwise it does not consider
// it as a polyAT
  size_t count_A= 0;
  size_t count_T= 0;
  size_t running_count_A= 0;
  size_t running_count_T= 0;
  size_t last_A= 0;
  size_t last_T= 0;
  size_t last_A_count= 0;
  size_t last_T_count= 0;
  for (i= 0;
		 i<_POLYA_MIN_LEN &&
			i<est_len;
		 ++i) {
	 if (EST_info->EST_seq[i]=='A') {
		++count_A;
		last_A= i;
		last_A_count= count_A;
	 }
	 if (EST_info->EST_seq[i]=='T') {
		++count_T;
		last_T= i;
		last_T_count= count_T;
	 }
  }
  DEBUG("Found %zu A and %zu T at the beginning of the sequence.", count_A, count_T);
  running_count_A= count_A;
  running_count_T= count_T;
  while ((i<est_len) &&
			((running_count_A>=(_POLYA_MIN_FRACTION*_POLYA_MIN_LEN)) ||
			 (running_count_T>=(_POLYA_MIN_FRACTION*_POLYA_MIN_LEN)))) {
	 if (EST_info->EST_seq[i-_POLYA_MIN_LEN]=='A') {
		--running_count_A;
	 }
	 if (EST_info->EST_seq[i-_POLYA_MIN_LEN]=='T') {
		--running_count_T;
	 }
	 if (EST_info->EST_seq[i]=='A') {
		++count_A;
		++running_count_A;
		last_A= i;
		last_A_count= count_A;
	 }
	 if (EST_info->EST_seq[i]=='T') {
		++count_T;
		++running_count_T;
		last_T= i;
		last_T_count= count_T;
	 }
	 ++i;
  }
  last_A= (last_A<_POLYA_MIN_LEN-1)?_POLYA_MIN_LEN-1:last_A;
  last_T= (last_T<_POLYA_MIN_LEN-1)?_POLYA_MIN_LEN-1:last_T;
  if ((last_A_count>=(_POLYA_MIN_FRACTION*(last_A+1))) ||
		(last_T_count>=(_POLYA_MIN_FRACTION*(last_T+1)))) {
	 if ((((double)last_A_count)/(last_A+1))>=(((double)last_T_count)/(last_T+1)))
		c= 'A';
	 else
		c= 'T';
// A polyA/T has been found
	 char sc= (c=='A')?_POLYA_CHR:_POLYT_CHR;
	 size_t mlen= ((c=='A')?last_A+1:last_T+1);
	 INFO("Found a %zubp long initial poly%c.", mlen, c);
	 for (i= 0; i<mlen; ++i)
		EST_info->EST_seq[i]= sc;
	 if(c=='A'){
		EST_info->pref_polyA_length=mlen;
	 } else {
		EST_info->pref_polyT_length=mlen;
	 }
  } else {
	 DEBUG("No polyA/T has been found at the beginning of the sequence.");
  }

// Check the end
  count_A= 0;
  count_T= 0;
  last_A= 0;
  last_T= 0;
  last_A_count= 0;
  last_T_count= 0;
  for (i= 0;
		 i<_POLYA_MIN_LEN &&
			i<est_len;
		 ++i) {
	 if (EST_info->EST_seq[est_len-i-1]=='A') {
		++count_A;
		last_A= i;
		last_A_count= count_A;
	 }
	 if (EST_info->EST_seq[est_len-i-1]=='T') {
		++count_T;
		last_T= i;
		last_T_count= count_T;
	 }
  }
  DEBUG("Found %zu A and %zu T at the end of the sequence.", count_A, count_T);
  running_count_A= count_A;
  running_count_T= count_T;
  while ((i<est_len) &&
			((running_count_A>=(_POLYA_MIN_FRACTION*_POLYA_MIN_LEN)) ||
			 (running_count_T>=(_POLYA_MIN_FRACTION*_POLYA_MIN_LEN)))) {
	 if (EST_info->EST_seq[est_len-i-1+_POLYA_MIN_LEN]=='A') {
		--running_count_A;
	 }
	 if (EST_info->EST_seq[est_len-i-1+_POLYA_MIN_LEN]=='T') {
		--running_count_T;
	 }
	 if (EST_info->EST_seq[est_len-i-1]=='A') {
		++count_A;
		++running_count_A;
		last_A= i;
		last_A_count= count_A;
	 }
	 if (EST_info->EST_seq[est_len-i-1]=='T') {
		++count_T;
		++running_count_T;
		last_T= i;
		last_T_count= count_T;
	 }
	 ++i;
  }
  last_A= (last_A<_POLYA_MIN_LEN-1)?_POLYA_MIN_LEN-1:last_A;
  last_T= (last_T<_POLYA_MIN_LEN-1)?_POLYA_MIN_LEN-1:last_T;
  if ((last_A_count>=(_POLYA_MIN_FRACTION*(last_A+1))) ||
		(last_T_count>=(_POLYA_MIN_FRACTION*(last_T+1)))) {
	 if ((((double)last_A_count)/(last_A+1))>=(((double)last_T_count)/(last_T+1)))
		c= 'A';
	 else
		c= 'T';
// A polyA/T has been found
	 char sc= (c=='A')?_POLYA_CHR:_POLYT_CHR;
	 size_t mlen= ((c=='A')?last_A+1:last_T+1);
	 INFO("Found a %zubp long final poly%c.", mlen, c);
	 for (i= 0; i<mlen; ++i)
		EST_info->EST_seq[est_len-i-1]= sc;
	 if(c=='A'){
		EST_info->suff_polyA_length=mlen;
	 } else {
		EST_info->suff_polyT_length=mlen;
	 }
  } else {
	 DEBUG("No polyA/T has been found at the end of the sequence.");
  }
}

void Ntails_removal(pEST_info EST_info) {
  int pref= 0;
  int suff= 0;
  int est_len= strlen(EST_info->EST_seq);

  my_assert(est_len>0);
  EST_info->pref_N_length= 0;
  if (EST_info->EST_seq[0]=='N') {
	 for (pref= 0;
			EST_info->EST_seq[pref]=='N';
			++pref) {
// Do nothing
	 }
	 int i= pref;
	 int j= 0;
	 while (i<=est_len) {
		EST_info->EST_seq[j]=
		  EST_info->EST_seq[i];
		++i;
		++j;
	 }
	 EST_info->pref_N_length= pref;
	 INFO("Removed a prefix of %d N.", pref);
  }
  est_len= strlen(EST_info->EST_seq);
  while ((suff<est_len) &&
			(EST_info->EST_seq[est_len-1-suff] == 'N')) {
	 ++suff;
  }
  if (suff == est_len) {
	 FATAL("The sequence is only composed by Ns.");
	 fail();
  }
  EST_info->EST_seq[est_len-suff]= '\0';
  EST_info->suff_N_length= suff;
  if (suff>0) {
	 INFO("Removed a suffix of %d N.", suff);
  }
}

#undef LEN_BUFFER
