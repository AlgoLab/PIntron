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
  list_destroy(l, (delete_function)pfree);
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
				fprintf(output_file, "%d %d %d %d ", factor->EST_start+1, factor->EST_end+1, factor->GEN_start+1, factor->GEN_end+1);
				int i;
				for(i=factor->EST_start; i<=factor->EST_end; i++){
				  fprintf(output_file, "%c", est->info->original_EST_seq[i]);
				}
				fprintf(output_file, " ");
				for(i=factor->GEN_start; i<=factor->GEN_end; i++){
				  fprintf(output_file, "%c", gen->original_EST_seq[i]);
				}
				fprintf(output_file, "\n");
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

void set_Chromosome(pEST_info EST_info){
  my_assert(EST_info != NULL);
  my_assert(EST_info->EST_id != NULL);

  char* p=strstr(EST_info->EST_id, "chr");
  if(p == NULL){
	 p=strstr(EST_info->EST_id, "CHR");
  }
  if(p != NULL){
	 p+=3;
	 EST_info->EST_chr=c_palloc(10);
	 int i=0;
	 //fare una versione safe!!
	 while(p!=NULL && *p != ':' && i < 9){
		EST_info->EST_chr[i]=*p;
		p++;
		i++;
	 }
	 EST_info->EST_chr[i]='\0';
	 DEBUG("CHR ==> %s", EST_info->EST_chr);
 }
  else{
	 DEBUG("CHR ==> unknown");
  }
}

void set_Chromosome_coordinates(pEST_info EST_info){
	  my_assert(EST_info != NULL);
	  my_assert(EST_info->EST_id != NULL);

	  char *p=strchr(EST_info->EST_id, ':');
	  char *l=strrchr(EST_info->EST_id, ':');

	  if(p != NULL && p != l){
		 my_assert(p < l);

		 char *inner_str=c_palloc(2*LEN_ABSCOORD_ARRAY+2);
		 strcpy(inner_str,"");
		 p++;
		 int index=0;
		 while(p != l){
			 inner_str[index]=*p;
			 p++;
			 index++;
		 }
		 inner_str[index]='\0';

		 if(strcmp(inner_str, "")){

			 char *c1=strchr(inner_str, ':');
			 char *c2=strrchr(inner_str, ':');

			 if(c1 == NULL || c1 != c2){
				 INFO("Genomic strand ==> unknown");
				 fail();
			 }
			 c2++;
			 *c1='\0';
			 EST_info->abs_start=atoi(inner_str);
			 EST_info->abs_end=atoi(c2);
		 }
	  }
	  else{
		 INFO("Genomic strand ==> unknown");
		 fail();
	  }
}

void set_Genomic_Strand(pEST_info EST_info){
  my_assert(EST_info != NULL);
  my_assert(EST_info->EST_id != NULL);

  char* p=strrchr(EST_info->EST_id, ':');

  if(p != NULL){
	 EST_info->EST_strand_as_read=c_palloc(LEN_STRAND_ARRAY+1);
	 int i=0;
	 do{
		p++;
		EST_info->EST_strand_as_read[i]=*p;
		i++;
	 }while((*p != '1') && (i<LEN_STRAND_ARRAY));
	 EST_info->EST_strand_as_read[i]='\0';
	 INFO("Genomic strand %s", EST_info->EST_strand_as_read);
  }
  else{
	 INFO("Genomic strand ==> unknown");
	 fail();
  }
}

void set_EST_Strand_and_RC(pEST_info EST_info, pEST_info gen){
  my_assert(EST_info != NULL);
  my_assert(EST_info->EST_id != NULL);
  my_assert(gen->EST_strand_as_read != NULL);

  EST_info->EST_strand_as_read= c_palloc(LEN_STRAND_ARRAY+1);

  if(EST_info->EST_gb != NULL && EST_info->EST_gb[0] == 'N' && EST_info->EST_gb[0] == 'M'){
	 strcpy(EST_info->EST_strand_as_read, "1");
	 EST_info->EST_strand= 1;
	 DEBUG("STRAND RefSeq==> %s", EST_info->EST_strand_as_read);
  } else {
	 char* p=strstr(EST_info->EST_id, "/clone_end=");
	 if(p == NULL){
		p=strstr(EST_info->EST_id, "/CLONE_END=");
	 }
	 if(p != NULL){
		p+=11;
		int i=0;
		while((i<LEN_STRAND_ARRAY) &&
				(*p != '\0') &&
				(*p != '\'')){
		  EST_info->EST_strand_as_read[i]= *p;
		  p++;
		  i++;
		}
		EST_info->EST_strand_as_read[i]='\0';
		DEBUG("STRAND ==> %s", EST_info->EST_strand_as_read);

		if(!strcmp(gen->EST_strand_as_read, "1") || !strcmp(gen->EST_strand_as_read, "+1")){
		  if(!strcmp(EST_info->EST_strand_as_read, "3"))
			 EST_info->EST_strand= 1;
		  else
			 EST_info->EST_strand= -1;
		}
		else{
		  if(!strcmp(EST_info->EST_strand_as_read, "5"))
			 EST_info->EST_strand= 1;
		  else
			 EST_info->EST_strand= -1;
		}
	 }
	 else{
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


void polyAT_substitution(pEST_info EST_info) {
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
