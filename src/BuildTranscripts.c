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
//Gestione del flag -4
//Riduzione: riduzione degli esterni (flag 2 e 1) con vincolo che solo l'esone precedente o successivo coincida
//Elimina le forme alternative in cui almeno un esone e' stato ridotto
//Output esteso
//Catena polya

#include <stdio.h>
#include <stdlib.h>

#define M 100
#define INCLUDE_EXTERNALS 1
#define READ_ABS_COORD

#define PRINT_POLYA

//19dic06
//Se questa macro e' definita, i trascritti associati ad un refseq (NM_*) non vengono estesi
#define DONT_EXTEND_REFSEQ


/**
	ex-file list.c
**/

//file contenente le procedure per la gestione delle strutture dati
//#define MESSAGES 1 //macro per la visualizzazione di messaggi in caso di errore
#include <string.h>
#include "my_time.h"
#include "log.h"
#include "util.h"
#include "log-build-info.h"

typedef struct structSimpleList{//lista semplice, cioe' una lista che contiene solo una stringa
  char *value;
  struct structSimpleList*prev;//puntatore all'elemento precedente
  struct structSimpleList*next;//puntatore all'elemento successivo
} structSimpleList;

typedef struct structExonList{//lista degli esoni:contiene 3 stringhe,start end e sequence e una lista semplice di ESTs
  char *start;//start dell'esone
  char *end;//end dell'esone
  char *sequence;//stringa dell'esone
  structSimpleList *est_ids;//insieme di est

//RAFFA 30mag04
//Lista delle prime coord dei trascritti in cui l'esone compare nella struttura transcripts (prima di ogni riduzione)
  structSimpleList *first_tr_coord;//first_tr_coord=i se il trascritto inizia con l'esone i
//Lista delle seconde coord dei trascritti in cui l'esone compare nella struttura transcripts (prima di ogni riduzione)
  structSimpleList *second_tr_coord;//second_tr_coord=j se il trascritto è il (j+1)-esimo tra quelli che iniziano con l'esone i
//Lista delle terze coord dei trascritti in cui l'esone compare nella struttura transcripts (prima di ogni riduzione)
  structSimpleList *third_tr_coord;//third_tr_coord=p se l'esone è il (p+1)-esimo del trascritto

  int reduce;//l'esone e' una riduzione di un altro?

//RAFFA 4giu04
  int external;	//l'esone è un esterno?

  int polya;	//Compare in almeno una EST con catena polya?

//22mar07
  int used;	//-2 l'esone non compare in nessun trascritto, -1 compare in qualche trascritto e >= 0 compare in qualche tr e il valore fornisce
//un indice d'ordine esclusi gli esoni che non vengono usati

  struct structExonList *prev;//puntatore all'elemento precedente
  struct structExonList *next;//puntatore all'elemento successivo

} structExonList;

structExonList *initializeExonList(){
//inizializza la lista degli esoni
  structExonList *temp=(structExonList*)calloc(1,sizeof(structExonList));
  temp->start=NULL;
  temp->end=NULL;
  temp->sequence=NULL;
  temp->est_ids=NULL;

//RAFFA 30mag04
  temp->first_tr_coord=NULL;
  temp->second_tr_coord=NULL;
  temp->third_tr_coord=NULL;

  temp->prev=NULL;
  temp->next=NULL;

//RAFFA 30mag04
//temp->reduce=0;
  temp->reduce=1;

//RAFFA 4giu04
  temp->external=1;

  temp->polya=0;

  return temp;
}

structSimpleList *initializeSimpleList(){
//inizializza una lista semplice
  structSimpleList *temp=(structSimpleList*)calloc(1,sizeof(structSimpleList));
  temp->value=NULL;
  temp->prev=NULL;
  temp->next=NULL;

  return temp;
}

/* static void printElementExonList(structExonList *s){ */
/* //stampa un elemento della lista degli esoni */
/*   if(s==NULL || s->start==NULL || s->sequence==NULL || s->end==NULL) { */
/* #ifdef MESSAGES */
/* 	 printf("Elemento vuoto\n"); */
/* #endif */
/* 	 return; */
/*   } */
/*   printf("start    %s\n",s->start); */
/*   printf("end      %s\n",s->end); */
/*   printf("sequence %s\n",s->sequence); */
/*   printf("\n"); */
/* } */

/* static void printElementSimpleList(structSimpleList *s){ */
/* //stampa un elemento della lista semplice */
/*   if(s==NULL || s->value==NULL) { */
/* #ifdef MESSAGES */
/* 	 printf("Elemento vuoto\n"); */
/* #endif */
/* 	 return; */
/*   } */
/*   printf("value %s\n",s->value); */
/* } */

static int lengthSimpleList(structSimpleList *s){
//ritorna la lunghezza della lista semplice
  int temp=0;
  if(s==NULL || s->value==NULL) return 0;
  while(s!=NULL){
	 temp++;
	 s=s->next;
  }
  return temp;
}

/* static void printSimpleList(structSimpleList *s){ */
/* //stampa l'intera lista semplice */
/*   if(lengthSimpleList(s)==0){ */
/* #ifdef MESSAGES */
/* 	 printf("printSimleList, lista vuota!\n"); */
/* #endif */
/* 	 return ; */
/*   } */
/*   while(s!=NULL){ */
/* 	 printElementSimpleList(s); */
/* 	 s=s->next; */
/*   } */
/* } */

/* static void printExonList(structExonList *s){ */
/* //stampa l'intera lista degli esoni */
/*   if (s==NULL) { */
/* #ifdef MESSAGES */
/* 	 printf("Elemento vuoto \n"); */
/* #endif */
/* 	 return; */
/*   } */
/*   while(s!=NULL){ */
/* 	 printElementExonList(s); */
/* 	 s=s->next; */
/*   } */
/* } */

static char *itoa(int n){
//restituisce una stringa contenente l'intero n
  char *temp=(char*)calloc(20,sizeof(char));
  int l=sprintf(temp, "%d", n);
  char *result=(char*)calloc(l+1,sizeof(char));
  strcpy(result,temp);
  free(temp);
  return result;
}

static structSimpleList *pushSimpleList(structSimpleList *s,const char *value){
//aggiunge un elemento alla lista semplice
  structSimpleList *temp=s;
  if(s==NULL){
	 temp=initializeSimpleList();
	 temp->value=(char*)calloc(strlen(value)+1,sizeof(char));
	 strcpy(temp->value,value);
	 temp->next=NULL;
	 temp->prev=NULL;
	 return temp;
  }
  while(s->next!=NULL)
	 s=s->next;
  s->next=(structSimpleList*)calloc(1,sizeof(structSimpleList));
  s->next->prev=s;
  s->next->value=(char*)calloc(strlen(value)+1,sizeof(char));
  strcpy(s->next->value,value);
  s->next->next=NULL;
  return temp;
}

//RAFFA 30mag04
//structExonList *pushExonList(structExonList *s,const char *start,char *end,char *sequence,structSimpleList *est_ids,int reduce){
//RAFFA 4giu04
//structExonList *pushExonList(structExonList *s,const char *start,char *end,char *sequence,structSimpleList *est_ids,structSimpleList *first_tr_coord,structSimpleList *second_tr_coord,structSimpleList *third_tr_coord,int reduce){
static structExonList *pushExonList(structExonList *s,const char *start,char *end,char *sequence,structSimpleList *est_ids,structSimpleList *first_tr_coord,structSimpleList *second_tr_coord,structSimpleList *third_tr_coord,int reduce, int external, int polya){
//aggiunge un elemento alla lista degli esoni
//start,end e sequence vengono copiati, est_ids invece no!
  structExonList *temp=s;
  if(s==NULL){
	 temp=initializeExonList();
	 temp->start=(char*)calloc(strlen(start)+1,sizeof(char));
	 strcpy(temp->start,start);
	 temp->end=(char*)calloc(strlen(end)+1,sizeof(char));
	 strcpy(temp->end,end);
	 temp->sequence=(char*)calloc(strlen(sequence)+1,sizeof(char));
	 strcpy(temp->sequence,sequence);
	 temp->est_ids=est_ids;

//RAFFA 30mag04
	 temp->first_tr_coord=first_tr_coord;
	 temp->second_tr_coord=second_tr_coord;
	 temp->third_tr_coord=third_tr_coord;

	 temp->next=NULL;
	 temp->prev=NULL;
	 temp->reduce=reduce;

//RAFFA 4giu04
	 temp->external=external;

	 temp->polya=polya;

	 return temp;
  }
  while(s->next!=NULL)
	 s=s->next;
  s->next=(structExonList*)calloc(1,sizeof(structExonList));
  s->next->prev=s;
  s->next->start=(char*)calloc(strlen(start)+1,sizeof(char));
  strcpy(s->next->start,start);
  s->next->end=(char*)calloc(strlen(end)+1,sizeof(char));
  strcpy(s->next->end,end);
  s->next->sequence=(char*)calloc(strlen(sequence)+1,sizeof(char));
  strcpy(s->next->sequence,sequence);
  s->next->est_ids=est_ids;

//RAFFA 30mag04
  s->next->first_tr_coord=first_tr_coord;
  s->next->second_tr_coord=second_tr_coord;
  s->next->third_tr_coord=third_tr_coord;

  s->next->next=NULL;

//RAFFA 30mag04
//s->next->reduce=0;
  s->next->reduce=1;

//RAFFA 4giu04
  s->next->external=1;

  s->next->polya=0;

  return temp;
}

static int lengthExonList(structExonList *s){
//ritorna la lunghezza di una exonlist
  int temp=0;
  if(s==NULL) return 0;
  while(s!=NULL){
	 temp++;
	 s=s->next;
  }
  return temp;
}

static structSimpleList *elementAtSimpleList(structSimpleList*s,int index){
//ritorna l'ELEMENTO index-esimo si una lista semplice
  int i;
  int l=lengthSimpleList(s);
  if(l==0) {
#ifdef MESSAGES
	 printf("impossibile considerare l'elemento %d, lista vuota\n",index);
#endif
	 return NULL;
  }
  if(index>l-1) {
#ifdef MESSAGES
	 printf("elementatsimplelist :Out of index\n");
#endif
	 return NULL;
  }
  for(i=0;i<index;i++)
	 s=s->next;
  return s;
}

static structExonList *elementAtExonList(structExonList *s,int index){
//ritorna l'ELEMENTO index-esimo di una exon list
  int i;
  int l=lengthExonList(s);
  if(index>l-1) {
#ifdef MESSAGES
	 printf("Out of index %d exonlist\n",index);
#endif
	 return NULL;
  }
  for(i=0;i<index;i++)
	 s=s->next;
  return s;
}

static int isEmptySimpleList(structSimpleList *s){
//ritorna 1 se la lista semplice è vuota
  if(lengthSimpleList(s)==0) return 1;else return 0;
}

/* static int isEmptyExonList(structExonList *s){ */
/* //ritorna 1 se la exonlist è vuota */
/*   if(lengthExonList(s)==0) return 1;else return 0; */
/* } */

static int findSimpleList(structSimpleList *s,const char *value){
//ritorna il primo indice nella lista il cui elemento = value altrimenti -1
  int i=0;
  if(s==NULL){
#ifdef MESSAGES
	 printf("Impossibile cercare elementi,simpleList vuota\n");
#endif
	 return -1;
  }
  if(value==NULL) {
#ifdef MESSAGES
	 printf("impossibile cercare parametro nullo!");
#endif
	 return -1;
  }
  while(s!=NULL){
	 if(s->value!=NULL && !strcmp(value,s->value)) return i;
	 i++;
	 s=s->next;
  }
  return -1;
}

/* static int equalsSimpleList(structSimpleList *s,structSimpleList *t){ */
/* //ritorna 1 se i value di s e t sono uguali(solo per elementi non per l'intera lista!) */
/*   if(s->value==NULL || s==NULL || t->value==NULL || t==NULL) return 0; */
/*   else return !strcmp(s->value,t->value); */
/* } */

/* static char *getValueSimpleList(structSimpleList *s){ */
/* //ritorna il value dell'elemento puntato da s in una lista semplice */
/*   return s->value; */
/* } */

static structSimpleList *difference(structSimpleList *s,structSimpleList *t){
//ritorna una nuova lista che contiene s\t
  int i;
  char *value;
  int ls=lengthSimpleList(s);
  structSimpleList *temp=NULL;
  for(i=0;i<ls;i++){
	 value=elementAtSimpleList(s,i)->value;
	 if(findSimpleList(t,value)==-1)
		temp=pushSimpleList(temp,value);
  }
  return temp;
}

static structSimpleList *intersection(structSimpleList *s,structSimpleList *t){
//ritorna una nuova lista che contiene s intersecato t
  int i;
  char *value;
  int ls=lengthSimpleList(s);
  structSimpleList *temp=NULL;
  for(i=0;i<ls;i++){
	 value=elementAtSimpleList(s,i)->value;
	 if(findSimpleList(t,value)!=-1)
		temp=pushSimpleList(temp,value);
  }
  return temp;
}

static structSimpleList *deleteSimpleList(structSimpleList *s,const char *value){
//elimina l'elemento che ha come valore value

  int index;
  int l=lengthSimpleList(s);

  if (l==0) {
#ifdef MESSAGES
	 printf("deletesimplelist! lista vuota\n");
#endif
	 return s;
  }

  index=findSimpleList(s,value);
  if(index>=l) {
#ifdef MESSAGES
	 printf("deletesimplelist! impossibile cancellare! out of index\n");
#endif
	 return s;
  }

  if(index==-1) {
#ifdef MESSAGES
	 printf("deletesimplelist! elemento non trovato!\n");
#endif
	 return s;
  }
  if(index==0){//sto cancellando il primo elemento
	 if(l==1) {
//la lista e' formata da un solo elemento che viene cancellato,ritorno una lista nulla
//			free(s);
		return NULL;
	 }
	 else{
//la lista ha almeno due elementi, quindi esiste il next
		structSimpleList *temp=s->next;
//			free(s);
		temp->prev=NULL;
		return temp;
	 }
  }
  if(index<l-1){//sto cancellando un elemento tra la testa e la coda
	 structSimpleList *prev=elementAtSimpleList(s,index)->prev;
	 structSimpleList *next=elementAtSimpleList(s,index)->next;
	 prev->next=next;
	 next->prev=prev;
//		free(elementAtSimpleList(s,index));
	 return s;
  }
  if(index==l-1){//sto cancellando l'ultimo elemento
	 structSimpleList *last=elementAtSimpleList(s,index)->prev;
	 last->next=NULL;
//		free(last->next);
	 return s;
  }
// by Yuri: Aggiunte le prossime istruzioni perche' il programma
// originale non prevedeva questo flusso di esecuzione
  fprintf(stderr, "ERROR - An impossible thing is happened!\n");
  exit(EXIT_FAILURE);
  return 0; // to make the compiler happy
}

static structSimpleList *copySimpleList(structSimpleList *s){//ritorna una copia della lista
  structSimpleList *temp=NULL;
  int l=lengthSimpleList(s);
  int i;
  for(i=0;i<l;i++)
	 temp=pushSimpleList(temp,elementAtSimpleList(s,i)->value);
  return temp;
}

static structSimpleList *unionSimpleList(structSimpleList *s,structSimpleList *t){
//ritorna l'unione di due liste semplici
//i duplicati non sono aggiunti
  int ls=lengthSimpleList(s);
  int lt=lengthSimpleList(t);
  int i;
  structSimpleList *temp=NULL;
  if(ls==0){
	 temp=copySimpleList(t);
	 return temp;
  }
  if(lt==0){
	 temp=copySimpleList(s);
	 return temp;
  }
  temp=copySimpleList(s);
  for(i=0;i<lt;i++){
	 char *value=elementAtSimpleList(t,i)->value;
	 if(findSimpleList(s,value)==-1) temp=pushSimpleList(temp,value);
  }
  return temp;
}

/**
	end ex-file list.c
**/

typedef struct structure2{
  structSimpleList *transcript;
  structSimpleList *coord1;
  structSimpleList *coord2;
} structure2;

typedef struct structure{
  int v1;
  int v2;

//RAFFA 30mag04
//Esone a cui si riduce
  int v3;

  struct structure *next;
} structure;

typedef struct structHash{
  char *id;
  int flag;
  int polya;
  int conf_ex_num;
  struct structHash* next;
} structHash ;

typedef struct struct_matx{
  int actual;
  int flag;
  structSimpleList *est_ids;
} struct_matx;

typedef struct struct_matx_pref_suff{
  int flag;
  int actual;
} struct_matx_pref_suff;

typedef struct listOfList{//una lista in cui ogni elemento e' una lista
  structSimpleList *local;
  struct listOfList*next;
} listOfList;

#ifdef READ_ABS_COORD
long int gen_start, gen_end;
int strand;
long int boundary;
#endif

static int lengthStructure(structure *s){//ritorna la lunghezza di una struttura
  int temp=0;
  if(s==NULL)
	 return 0; 	while(s!=NULL){
	 temp++;
	 s=s->next;
  }
  return temp;
}

/* static void printElementHash(structHash *s){//stampa la tabella hash */
/*   if(s==NULL ){ */
/* 	 printf("Elemento vuoto\n"); */
/* 	 return; */
/*   } */
/*   while(s!=NULL){ */
/* 	 printf("id          %s\n",s->id); */
/* 	 printf("flag        %d\n",s->flag); */
/* 	 printf("conf_ex_num %d\n",s->conf_ex_num); */
/* 	 printf("local       %d\n",s); */
/* 	 printf("\n"); */
/* 	 s=s->next; */
/*   } */
/*   return; */
/* } */

/* static structure *copyStructure(structure *s){//ritorna una copia di una struttura */
/*   if(s==NULL) { */
/* 	 printf("copyStructure:stringa da copiare nulla\n"); */
/* 	 return NULL; */
/*   } */
/*   structure *temp=NULL; */
/*   temp=(structure*)calloc(1,sizeof(structure)); */
/*   temp->v1=s->v1; */
/*   temp->v2=s->v2; */

/* //RAFFA 30mag04 */
/*   temp->v3=s->v3; */

/*   temp->next=NULL; */
/*   return temp; */
/* } */

//RAFFA 30mag04
static structure *pushStructure(structure *s,int v1,int v2,int v3){//aggiunge in coda una struttura
//structure *pushStructure(structure *s,int v1,int v2){//aggiunge in coda una struttura
  structure *temp=s;
  if(temp==NULL){
	 temp=(structure*)calloc(1,sizeof(structure));
	 temp->v1=v1;
	 temp->v2=v2;

//RAFFA 30mag04
	 temp->v3=v3;

	 temp->next=NULL;
	 return temp;
  }
  while(s->next!=NULL)
	 s=s->next;
  s->next=(structure*)calloc(1,sizeof(structure));
  s->next->v1=v1;
  s->next->v2=v2;

//RAFFA 30mag04
  s->next->v3=v3;

  s->next->next=NULL;
  return temp;
}

static structure *elementAtStructure(structure *s,int index){//ritorna un elemento della struttra
  if(index>=lengthStructure(s))
	 return NULL;
  int i;
  for(i=0;i<index;i++)
	 s=s->next; 	return s;
}

/* static structure *getStructure(structure *s,int index){//ritorna l'elemento index-esimo di una struttura */
/*   structure *temp=NULL; */
/*   if(index>=lengthStructure(s)){ */
/* 	 printf("getStructure: out of index!\n"); */
/* 	 return NULL; */
/*   } */
/*   temp=elementAtStructure(s,index); */
/*   return temp; */
/* } */

//RAFFA 30mag04
static void setStructure(structure *s,int index,int v1,int v2,int v3){//modifica i dati di una struttura
//void setStructure(structure *s,int index,int v1,int v2){//modifica i dati di una struttura
  if(index>=lengthStructure(s))
	 return;
  elementAtStructure(s,index)->v1=v1;
  elementAtStructure(s,index)->v2=v2;

//RAFFA 30mag04
  elementAtStructure(s,index)->v3=v3;
}

//RAFFA 30mag04
static int equalsStructure(structure *s,int index,int v1,int v2,int v3){//la struttura ha i valori v1 e v2?
//int equalsStructure(structure *s,int index,int v1,int v2){//la struttura ha i valori v1 e v2?
  if(index>=lengthStructure(s))
	 return 0;
  structure *temp=elementAtStructure(s,index);

//RAFFA 30mag04
  if (temp->v1==v1 && temp->v2==v2 && temp->v3==v3)
//if (temp->v1==v1 && temp->v2==v2)
	 return 1;

  else return 0;
}

/* static void printStructure(structure *s){//stampa una struttura */
/*   int i; */
/*   for(i=0;i<lengthStructure(s);i++) */
/* //RAFFA 30mag04 */
/* 	 printf("(%d,%d,%d)\n",elementAtStructure(s,i)->v1,elementAtStructure(s,i)->v2,elementAtStructure(s,i)->v3); */
/* //printf("(%d,%d)\n",elementAtStructure(s,i)->v1,elementAtStructure(s,i)->v2); */
/* } */

static structHash*initializeHash(){//inizializza hash
  structHash *temp=(structHash*)calloc(1,sizeof(structHash));
  temp->id=NULL;
  temp->flag=-100;
  temp->conf_ex_num=-100;
  temp->next=NULL;
  return temp;
}

/* static int emptyHash(structHash *h){//ritorna 1 se h e' NULL */
/*   if (h==NULL) return 1;else return 0; */
/* } */

static char *read_single_char(FILE *source){//legge da file un parametro
//05dic05
  char *temp=(char*)calloc(300000,sizeof(char));//variabile temporanea
  fscanf(source,"%s",temp);
  char *result=(char*)calloc(strlen(temp)+1,sizeof(char));
  strcpy(result,temp);
  free(temp);
  return result;
}

static long int Hashing_word(char *w){//funzione hash
  long int h;
  for(h=0;*w!='\0';w++)
	 h=(4*h+*w) % M; /*M global variable already filled*/ 	return h;
}

static int findHash(structHash *s,char *id){//ritorna la posizione che contiene id nell'hash altrimenti -1
  int temp=0;
  if(s==NULL){
	 printf("hash vuoto!\n");
	 return -1;
  }
  while(s!=NULL){
	 if(!strcmp(s->id,id)) return temp;
	 temp++;
	 s=s->next;
  }
  return -1; }

static int lengthHash(structHash *s){//ritorna la lunghezza dell'hash
  int temp=0;
  while(s!=NULL){
	 temp++;
	 s=s->next;
  }     return temp;
}

static structHash *elementAtHash(structHash *s,int n){//ritorna l'n-esimo elemento dell'hash
  int i;
  if(n>=lengthHash(s)) {
	 printf("Hash : out of index\n");
	 return NULL;
  }
  for(i=0;i<n;i++)
	 s=s->next;
  return s;
}
//variabili globali
structHash *est_hash[M];
struct_matx* matx;
struct_matx_pref_suff *matx_pref_suff;
structExonList *exon_list;
structExonList *exon_list_copy;

/* static structSimpleList* findExon(structExonList *exon_list,const char *start,const char *end){ */
/*   int i; */
/*   structSimpleList *temp=NULL; */
/*   for(i=0;i<lengthExonList(exon_list);i++){ */
/* 	 if(!strcmp(elementAtExonList(exon_list,i)->start,start) && !strcmp(elementAtExonList(exon_list,i)->end,end)) */
/* 		temp=pushSimpleList(temp,itoa(i)); */
/*   } */
/*   return temp; */
/* } */
/* static structExonList *deleteExonList(structExonList *exon_list,int pos){ */
/* //cancella un elemento dalla exon_list */
/* //cancella il pos-esimo elemento */
/*   structExonList *temp=exon_list; */
/*   int l=lengthExonList(exon_list); */
/*   if(l==1) */
/* 	 return NULL; */
/*   if (pos>=l){ */
/* 	 printf("deleteExonList: out out index\n"); */
/* 	 return temp; */
/*   } */
/*   if(pos==0) { */
/* //sto cancellando la testa */
/* 	 exon_list=exon_list->next; */
/* 	 exon_list->prev=NULL; */
/* 	 return exon_list; */
/*   } */
/*   if(pos==l-1){ */
/* //sto cancellando l'ultimo elemento */
/* 	 elementAtExonList(exon_list,pos-1)->next=NULL; */
/* 	 return temp; */
/*   } */
/* //sto cancellando un elemento tra la testa e la coda */
/*   structExonList *prev=elementAtExonList(exon_list,pos)->prev; */
/*   structExonList *next=elementAtExonList(exon_list,pos)->next; */
/*   elementAtExonList(exon_list,pos-1)->next=elementAtExonList(exon_list,pos)->next; */
/*   elementAtExonList(exon_list,pos+1)->prev=elementAtExonList(exon_list,pos-1); */
/*   return temp; */
/* } */

/* static structExonList *processing(structExonList *exon_list){ //elimina dalla lista gli esoni che hanno lo stesso start e lo stesso end */
/*   structSimpleList *found=NULL; */
/*   structExonList *result=NULL; */
/*   int i; */
/*   int j; */
/*   structSimpleList *est_ids=NULL; */

/* //RAFFA 30mag04 */
/*   structSimpleList *first_tr_coord=NULL; */
/*   structSimpleList *second_tr_coord=NULL; */
/*   structSimpleList *third_tr_coord=NULL; */


/*   while(lengthExonList(exon_list)>0){ */
/* 	 structExonList *exon=elementAtExonList(exon_list,0); */
/* 	 char *start=exon->start; */
/* 	 char *end=exon->end; */

/* 	 char *sequence=exon->sequence; */
/* 	 structSimpleList *found=findExon(exon_list,start,end); */
/* 	 int l=lengthSimpleList(found); */
/* 	 int pos[l]; */
/* 	 for(i=0;i<l;i++) */
/* 		pos[i]=atoi(elementAtSimpleList(found,i)->value); */
/* 	 for(i=0;i<l;i++){ */
/* 		est_ids=unionSimpleList(est_ids,elementAtExonList(exon_list,pos[i])->est_ids); */
/* 	 } */

/* //RAFFA 30mag04 (posso commentare?) */
/* 	 for(i=0;i<l;i++){ */
/* 		first_tr_coord=unionSimpleList(first_tr_coord,elementAtExonList(exon_list,pos[i])->first_tr_coord); */
/* 	 } */
/* 	 for(i=0;i<l;i++){ */
/* 		second_tr_coord=unionSimpleList(second_tr_coord,elementAtExonList(exon_list,pos[i])->second_tr_coord); */
/* 	 } */
/* 	 for(i=0;i<l;i++){ */
/* 		third_tr_coord=unionSimpleList(third_tr_coord,elementAtExonList(exon_list,pos[i])->third_tr_coord); */
/* 	 } */

/* 	 for(i=0;i<l;i++){ */
/* 		exon_list=deleteExonList(exon_list,pos[i]-i); */
/* 	 } */

/* //RAFFA 30mag04 */
/* //result=pushExonList(result,start,end,sequence,copySimpleList(est_ids),0); */
/* //RAFFA 4giu04 */
/* //result=pushExonList(result,start,end,sequence,copySimpleList(est_ids),copySimpleList(first_tr_coord),copySimpleList(second_tr_coord),copySimpleList(third_tr_coord),1); */
/* 	 result=pushExonList(result,start,end,sequence,copySimpleList(est_ids),copySimpleList(first_tr_coord),copySimpleList(second_tr_coord),copySimpleList(third_tr_coord),1,1,0); */

/* 	 est_ids=NULL; */

/* //RAFFA 30mag04 */
/* 	 first_tr_coord=NULL; */
/* 	 second_tr_coord=NULL; */
/* 	 third_tr_coord=NULL; */
/*   } */
/*   return result; */
/* } */

static structExonList *build_structures(char *filename){//costuisce l'hash e la exon_list
  int k;
  int p;
  char *temp;
  //char *temp2;
  char *start;
  char *end;
  char *sequence;
  structExonList *result=NULL;
  structSimpleList *l1=NULL;
  structSimpleList *l2=NULL;
  structSimpleList *l3=NULL;
  structSimpleList *est_ids=NULL;

//RAFFA 30mag04
  structSimpleList *first_tr_coord=NULL;
  structSimpleList *second_tr_coord=NULL;
  structSimpleList *third_tr_coord=NULL;

  //structExonList *exon_list=NULL;
  //int tt=0;
  //int h=0;
  int key;
  int stop=0;
  FILE *source=fopen(filename,"r");
  for(k=0;k<M;k++){
	 est_hash[k]=NULL;
  }

#ifdef READ_ABS_COORD
  fscanf(source,"%ld",&gen_start);
  fscanf(source,"%ld",&gen_end);
  fscanf(source,"%d",&strand);
  fscanf(source,"%ld",&boundary);
#endif

  while (!feof(source) && !stop){
	 fscanf(source,"%d",&p);
	 if(p==0){
		stop=-1;
	 }
	 else{
		l1=NULL;
		l2=NULL;
		l3=NULL;
		for(k=0;k<p;k++){
		  temp=read_single_char(source);
//fprintf(stdout, "ID %s ", temp);
		  l1=pushSimpleList(l1,temp);
		  free(temp);
		  temp=read_single_char(source);
//fprintf(stderr, "flag %s ", temp);
		  l2=pushSimpleList(l2,temp);
		  free(temp);
		  temp=read_single_char(source);
//fprintf(stderr, "polya %s ", temp);
		  l3=pushSimpleList(l3,temp);
		  free(temp);
		}

		start=read_single_char(source);
		end=read_single_char(source);

//fprintf(stdout, "start %s end %s\n", start, end);

		sequence=read_single_char(source);

//fprintf(stdout, "%s\n", sequence);

		est_ids=NULL;
		for(k=0;k<p;k++){
		  int flag;
		  int polya;
		  char *id;
		  flag=atoi(elementAtSimpleList(l2,k)->value);//contiene i flags
		  id=elementAtSimpleList(l1,k)->value;//contiene gli id degli est
		  polya=atoi(elementAtSimpleList(l3,k)->value);//contiene i flag polya
//printf("ID %s flag %d polya %d\n", id, flag, polya);
		  key=Hashing_word(id);	//RAFFA
//printf("ID %s key %d\n", id, key);
		  int kk;
#ifdef INCLUDE_EXTERNALS
		  est_ids=pushSimpleList(est_ids,id);
#endif
#ifndef INCLUDE_EXTERNALS
		  if (flag==0){
			 est_ids=pushSimpleList(est_ids,id);
#endif
			 if(est_hash[key]==NULL){//l'elemento e' vuoto, lo inserisco
				est_hash[key]=initializeHash();
				est_hash[key]->id=(char*)calloc(strlen(id)+1,sizeof(char));
				strcpy(est_hash[key]->id,id);
				est_hash[key]->flag=flag;
				est_hash[key]->polya=polya;
				est_hash[key]->conf_ex_num=1;
			 }
			 else{//non e' vuoto, guardo se esiste gia' l'est
				kk=findHash(est_hash[key],id);
				if(kk==-1){//non e' presenmte l'est , lo aggiungo come se fosse una lista
//scorro tutto l'hash
				  structHash *temp=est_hash[key];
				  while(temp->next!=NULL)
					 temp=temp->next;
//e accodo
				  temp->next=initializeHash();
				  temp->next->id=(char*)calloc(strlen(id)+1,sizeof(char));
				  strcpy(temp->next->id,id);
				  temp->next->flag=flag;
				  temp->next->polya=polya;
				  temp->next->conf_ex_num=1;
				}
				else{//kk != -1, l'est esiste gia' incremento solo il conf_Ex_num
#ifdef INCLUDE_EXTERNALS
				  int value=elementAtHash(est_hash[key],kk)->flag;
				  if(value==-1)
					 if(flag==-2)
						elementAtHash(est_hash[key],kk)->flag=-3;
					 else
						elementAtHash(est_hash[key],kk)->flag=-1;
				  else
					 elementAtHash(est_hash[key],kk)->flag=flag;
#endif

				  elementAtHash(est_hash[key],kk)->conf_ex_num++;
				}
			 }
		  }
#ifndef INCLUDE_EXTERNALS
		  else{//il flag e' diverso da 0
			 if(est_hash[key]==NULL){//l'elemento e' vuoto, lo inserisco
				est_hash[key]=initializeHash();
				est_hash[key]->id=(char*)calloc(strlen(id)+1,sizeof(char));
				strcpy(est_hash[key]->id,id);
				est_hash[key]->flag=flag;
				est_hash[key]->polya=polya;
				est_hash[key]->conf_ex_num=0;
			 }
			 else{//non e' vuoto, guardo se esiste gia' l'est
				int kk;
				kk=findHash(est_hash[key],id);
				structHash *temp=est_hash[key];
				if(kk==-1){//non e' presenmte l'est , lo aggiungo come se fosse una lista
//scorro tutto l'hash
				  while(temp->next!=NULL)
					 temp=temp->next;
//e accodo l'est
				  temp->next=initializeHash();
				  temp->next->id=(char*)calloc(strlen(id)+1,sizeof(char));
				  strcpy(temp->next->id,id);
				  temp->next->flag=flag;
				  temp->next->polya=polya;
				  temp->next->conf_ex_num=0;
				}
				else{//l'est esiste gia'
//est_hash{list1(k)}->flag=(est_hash{list1(k)}->flag=0)?(list2(k)):(-3);
				  if (elementAtHash(est_hash[key],kk)->flag==0)elementAtHash(est_hash[key],kk)->flag=flag;
				  else elementAtHash(est_hash[key],kk)->flag=-3;
				}
			 }
		  }
		}
#endif
		if(!isEmptySimpleList(est_ids)){
//RAFFA 30mag04
//result=pushExonList(result,start,end,sequence,est_ids,0);
//RAFFA 4giu04
//result=pushExonList(result,start,end,sequence,est_ids,first_tr_coord, second_tr_coord,third_tr_coord,1);
		  result=pushExonList(result,start,end,sequence,est_ids,first_tr_coord, second_tr_coord,third_tr_coord,1,1,0);
		}
	 }
  }

//05feb07
//structExonList *result2=processing(result);
//int exon_number=lengthExonList(result2);
  int exon_number=lengthExonList(result);

  int index;

//05feb07
/*for(index=0;index<exon_number;index++){
  structSimpleList *temp_ids=elementAtExonList(result2,index)->est_ids;
//printf("Start %s End %s poly %d\n", elementAtExonList(result2,index)->start, elementAtExonList(result2,index)->end, elementAtExonList(result2,index)->polya);
char stop=0;
int z=0;
while(z<lengthSimpleList(temp_ids) && !stop){
char *ident=elementAtSimpleList(temp_ids,z)->value;
int key=Hashing_word(ident);
//printf("ID %s key %d\n", ident, key);
structHash *temp_hash=est_hash[key];
//Cerco la EST
char stop2=0;
while(temp_hash != NULL && ! stop2){
//printf("	ID %s polya %d\n", temp_hash->id, temp_hash->polya);
if(!strcmp(temp_hash->id, ident))
stop2=1;
else
temp_hash=temp_hash->next;
}
//printf("		poly %d\n", temp_hash->polya);
if(temp_hash->polya == 1){
elementAtExonList(result2,index)->polya=1;
stop=1;
}
else
z++;
}
}

int i;

//	printExonList(exon_list);
//	exit(0);
return result2;*/
  for(index=0;index<exon_number;index++){
	 structSimpleList *temp_ids=elementAtExonList(result,index)->est_ids;
//printf("Start %s End %s poly %d\n", elementAtExonList(result,index)->start, elementAtExonList(result,index)->end, elementAtExonList(result,index)->polya);
	 char stop=0;
	 int z=0;
	 while(z<lengthSimpleList(temp_ids) && !stop){
		char *ident=elementAtSimpleList(temp_ids,z)->value;
		int key=Hashing_word(ident);
//printf("ID %s key %d\n", ident, key);
		structHash *temp_hash=est_hash[key];
//Cerco la EST
		char stop2=0;
		while(temp_hash != NULL && ! stop2){
//printf("	ID %s polya %d\n", temp_hash->id, temp_hash->polya);
		  if(!strcmp(temp_hash->id, ident))
			 stop2=1;
		  else
			 temp_hash=temp_hash->next;
		}
//printf("		poly %d\n", temp_hash->polya);
		if(temp_hash->polya == 1){
		  elementAtExonList(result,index)->polya=1;
		  stop=1;
		}
		else
		  z++;
	 }
  }

 // int i;

//	printExonList(exon_list);
//	exit(0);
  return result;
}

static struct_matx_pref_suff *getMatxPrefSuff(struct_matx_pref_suff *a,int i,int j,int n){//ritorna l'elemento [i,j] di matx_preff_suff
  a+=i*n+j;
  return a;
}

static struct_matx *getMatx(struct_matx *a,int i,int j,int n){//ritorna l'elemento [i,j] di matx
  a+=i*n+j;
  return a;
}

static listOfList*pushListOfList(listOfList *s,structSimpleList *t){//aggiunge un elemento (lista) ad una listOfList
  listOfList *temp=s;
  if (s==NULL){
	 temp=(listOfList*)calloc(1,sizeof(listOfList));
	 temp->local=t;
	 temp->next=NULL;
	 return temp;
  }
  while(s->next!=NULL)
	 s=s->next;
  s->next=(listOfList*)calloc(1,sizeof(listOfList));
  (s->next)->local=t;
  (s->next)->next=NULL;
  return temp;
}

static int lengthListOfList(listOfList *s){//ritorna la lunghezza si una listOfList(cioe' il numero di liste che contiene)
  int temp=0;
  if(s==NULL) return 0;
  while(s!=NULL){
	 temp++;
	 s=s->next;
  }
  return temp;
}

static listOfList*elementAtListOfList(listOfList*s,int n){//ritorna l'elemento n-esimo (una lista)
  int i;
  if(n>=lengthListOfList(s)) {
	 printf("list of list! Out of index\n");
	 return NULL;
  }
  for (i=0;i<n;i++)
	 s=s->next;
  return s;
}

/* static void printElementListOfList(listOfList*s,int index){//stampa un elemento (lista) di una listOfList */
/*   printSimpleList(elementAtListOfList(s,index)->local); */
/* } */

static int findListOfList(listOfList*s,char *value){	//ritorna l'indice della prima lista che ha al suo interno value
  int i;
  int k;
  int l=lengthListOfList(s);
  if (l==0) {
//printf("List of list vuota\n");
	 return -1;
  }
  for(i=0;i<l;i++){
	 structSimpleList *temp=elementAtListOfList(s,i)->local;
	 k=findSimpleList(temp,value);
	 if (k!=-1)
		return i;
  }
  return -1;
}

static listOfList *deleteListOfList(listOfList *s,int index){//cancella un elemento (lista) da una listOfList
  int l=lengthListOfList(s);
  if(index>=l){
	 printf("deleteListOfList: out of index!\n");
	 return s;
  }
  if(index==0){
//stocancellando la testa
	 if(l==1)
		return NULL;
	 else {
		return s->next;
	 }
  }
  if(index==l-1){//sto cancellando l'ultimo elemento
//imposto il next dell'elemento precedente a null
	 elementAtListOfList(s,index-1)->next=NULL;
	 return s;
  }
//se entro qui significa che 0<index<l-1
  elementAtListOfList(s,index-1)->next=elementAtListOfList(s,index)->next;
  return s;
}

/* static char *getLastTranscript(char *s){//ritorna l'ultimo esone di un trascritto */
/*   int l=strlen(s); */
/*   s+=l-1; */
/*   *s--; */
/*   while(*s!='.') */
/* 	 s--; */
/*   s++; */
/*   int size=strlen(s)+1; */
/*   char *temp=(char*)calloc(size,sizeof(char)); */
/*   strcpy(temp,s); */
/*   return temp; */
/* } */
/* static int findFlagFromHash(char *id,int flag){//ritorna true se esiste un est_hash[id]->flag==flag */
/*   int n=Hashing_word(id); */
/*   structHash *temp=est_hash[n]; */
/*   while(temp!=NULL){ */
/* 	 if(temp->flag==flag && !strcmp(temp->id,id)) */
/* 		return 1; */
/* 	 temp=temp->next; */
/*   } */
/*   return 0; */
/* } */

static int verify_competing53(listOfList *end_tr_est_ids[],int i,int j,int n){
  int matx_pref_suff_ij=getMatxPrefSuff(matx_pref_suff,i,j,n)->flag;

  if(matx_pref_suff_ij==-1 || matx_pref_suff_ij==-3){
	 if(elementAtExonList(exon_list_copy, j)->polya != 0 || elementAtExonList(exon_list_copy, j)->external == -1 || elementAtExonList(exon_list_copy, j)->external == 0){
		return 1;
	 }
	 else{
		//int exist=0;
#ifndef INCLUDE_EXTERNALS
		structSimpleList *end_tr_est_id=end_tr_est_ids[j]->local;
		int ii;
		for(ii=0;ii<lengthSimpleList(end_tr_est_id);ii++){
		  char *id=elementAtSimpleList(end_tr_est_id,ii)->value;
		  if( findFlagFromHash(id,-2) ||findFlagFromHash(id,-3)){
			 return 1;
		  }
		}
#endif
		return 0;
	 }
  }

  if(matx_pref_suff_ij==1 || matx_pref_suff_ij==3){
//if(lengthSimpleList(end_transcripts[i]->transcript)==0){
	 if(elementAtExonList(exon_list_copy, i)->polya != 0 || elementAtExonList(exon_list_copy, i)->external == -1 || elementAtExonList(exon_list_copy, i)->external == 0){
		return 1;
	 }
	 else{
#ifndef INCLUDE_EXTERNALS
		structSimpleList *end_tr_est_id=end_tr_est_ids[i]->local;
		int ii;
		for(ii=0;ii<lengthSimpleList(end_tr_est_id);ii++){
		  char *id=elementAtSimpleList(end_tr_est_id,ii)->value;
		  if( findFlagFromHash(id,-2) || findFlagFromHash(id,-3)){
			 return 1;
		  }
		}
#endif
		return 0;
	 }
  }

  if (matx_pref_suff_ij==-4){
//if (lengthSimpleList(end_transcripts[i]->transcript)==0){
	 if(elementAtExonList(exon_list_copy, i)->polya != 0 || elementAtExonList(exon_list_copy, i)->external == -1 || elementAtExonList(exon_list_copy, i)->external == 0){
		return 1;
	 }
	 else{
//RICONTROLLARE
#ifndef INCLUDE_EXTERNALS
		structSimpleList *end_tr_est_id=end_tr_est_ids[i];
		int ii;
		for(ii=0;ii<lengthSimpleList(end_tr_est_id);ii++){
		  char *id=elementAtSimpleList(end_tr_est_id,ii)->value;
		  if( findFlagFromHash(id,-2)||findFlagFromHash(id,-3)){
			 return 1;
		  }
		}
#endif
		return 0;
	 }
  }

  return 0;
}
static int verify_competing35(listOfList *tr_est_ids[],int i,int j,int n){
  int matx_pref_suff_ij=getMatxPrefSuff(matx_pref_suff,i,j,n)->flag;

  if(matx_pref_suff_ij==-2 || matx_pref_suff_ij==-3){
	 if(elementAtExonList(exon_list_copy, j)->external == -2 || elementAtExonList(exon_list_copy, j)->external == 0){
		return 1;
	 }
	 else{
#ifndef INCLUDE_EXTERNALS
		structSimpleList *tr_est_id=tr_est_ids[j]->local;
		int ii;
		for(ii=0;ii<lengthSimpleList(tr_est_id);ii++){
		  char *id=elementAtSimpleList(tr_est_id,ii)->value;
		  if(findFlagFromHash(id,-1) ||findFlagFromHash(id,-3)){
			 return 1;
		  }
		}
#endif
		return 0;
	 }
  }

  if(matx_pref_suff_ij==2 || matx_pref_suff_ij==3){
	 if(elementAtExonList(exon_list_copy, i)->external == -2 || elementAtExonList(exon_list_copy, i)->external == 0){
		return 1;
	 }
	 else{
#ifndef INCLUDE_EXTERNALS
		structSimpleList *tr_est_id=tr_est_ids[i]->local;
		int ii;
		int exist=0;
		for(ii=0;ii<lengthSimpleList(tr_est_id);ii++){
		  char *id=elementAtSimpleList(tr_est_id,ii)->value;
		  if( findFlagFromHash(id,-1)||findFlagFromHash(id,-3)){
			 return 1;
		  }
		}
#endif
		return 0;
	 }
  }
  if (matx_pref_suff_ij==-4){
//if (lengthSimpleList(transcripts[j])==0){
	 if(elementAtExonList(exon_list_copy, j)->external == -2 || elementAtExonList(exon_list_copy, j)->external == 0){
		return 1;
	 }
	 else{
//RICONTROLLARE
#ifndef INCLUDE_EXTERNALS
		structSimpleList *end_tr_est_id=end_tr_est_ids[i];
		int ii;
		for(ii=0;ii<lengthSimpleList(tr_est_id);ii++){
		  char *id=elementAtSimpleList(tr_est_id,ii)->value;
		  if( findFlagFromHash(id,-2)||findFlagFromHash(id,-3)){
			 return 1;
		  }
		}
#endif
		return 0;
	 }
  }
  return 0;
}

static int lengthTranscript( char *s){//ritorna il numero di esoni presenti nel trascritto
  int i;
  int temp=1;
  int l=strlen(s);
  for(i=0;i<l;i++){
	 if (*s=='.')
		temp++;
	 s++;
  }
  return temp;
}

static char *elementAtTranscript(char *s,int index){//ritorna l'esone index-esimo
  char *t=(char*)calloc(20,sizeof(char));
   int l=lengthTranscript(s);
  char *temp2=s;
  if(index>=l )
	 return NULL;
  int temp=0;
  int pos[l];
  int i;
  for(i=0;i<l;i++)
	 pos[i]=-1;
  for(i=0;i<(int)strlen(temp2);i++){
	 if(*s=='.')
		pos[temp++]=i;
	 s++;
  }
  s=temp2;
  if(index==0) {
	 if(pos[index]>0)
		strncpy(t,s,pos[index]);
	 else
		strcpy(t,s);
	 return t;
  }
  if(index==lengthTranscript(s)-1){
	 s+=pos[index-1]+1;
	 strcpy(t,s);
	 return t;
  }
  if(index>=1){
	 s+=pos[index-1]+1;
	 strncpy(t,s,pos[index]-pos[index-1]-1);
  }
  return t;
}

static char *firstExon(structSimpleList *s,int index){//ritorna il primo esone
  char *temp=elementAtSimpleList(s,index)->value;
  return elementAtTranscript(temp,0);
}

static char *lastExon(structSimpleList *s,int index){//ritorna l'ultimo esone
  char *temp=elementAtSimpleList(s,index)->value;
  int l=lengthTranscript(temp);
  return elementAtTranscript(temp,l-1);
}

/* static char *getInternTranscript(char *s){//ritorna un trascritto che non contiene il primo e l'ultimo esone */
/*   int i; */
/*   int l=lengthTranscript(s); */
/*   if(l<=2) { */
/* 	 return NULL; */
/*   } 	char *temp=(char*)calloc(strlen(s),sizeof(char)); */
/*   strcpy(temp,""); */
/*   for(i=1;i<l-1;i++){ */
/* 	 strcat(temp,elementAtTranscript(s,i)); */
/* 	 if(i!=l-2) strcat(temp,"."); */
/*   } */
/*   return temp; */
/* } */
/* static void printAlternativesForms(FILE *out,structSimpleList *tr_list_i,structSimpleList *tr_list_j,int i,int j,int flag){ */
/*   fprintf(out,"%d,%d,%d\n",i,j,flag); */
/*   int ii; */
/*   int length=lengthSimpleList(tr_list_i); */
/*   if(length==0) */
/* 	 fprintf(out,"%s\n","-1"); */
/*   else{ */
/* 	 for(ii=0;ii<length-1;ii++) */
/* //fprintf(out,"%s,",elementAtSimpleList(tr_list_i,ii)->value); */
/* 		fprintf(out,"%d,",atoi(elementAtSimpleList(tr_list_i,ii)->value)+1);//giorgio2luglio */

/* //fprintf(out,"%s\n",elementAtSimpleList(tr_list_i,ii)->value); */
/* 	 fprintf(out,"%d\n",atoi(elementAtSimpleList(tr_list_i,ii)->value)+1);//giorgio2luglio */
/*   } */

/*   length=lengthSimpleList(tr_list_j); */
/*   if(length==0) */
/* 	 fprintf(out,"%s\n","-1"); */
/*   else{ */
/* 	 for(ii=0;ii<length-1;ii++) */
/* //fprintf(out,"%s,",elementAtSimpleList(tr_list_j,ii)->value); */
/* 		fprintf(out,"%d,",atoi(elementAtSimpleList(tr_list_j,ii)->value)+1);//giorgio2luglio */

/* //fprintf(out,"%s\n",elementAtSimpleList(tr_list_j,ii)->value); */
/* 	 fprintf(out,"%d\n",atoi(elementAtSimpleList(tr_list_j,ii)->value)+1);//giorgio2luglio */
/*   } */

/* } */


static void transcripts_finder(char *inputfile,char *outputfile){
  int n;
  int i,j;
  int t;
  char cmp;

  exon_list=build_structures(inputfile);

  exon_list_copy=build_structures(inputfile);

  n=lengthExonList(exon_list);

/*for(i=0;i<n;i++){
  printf("%d) %s:%s\n",i, elementAtExonList(exon_list_copy,i)->start,elementAtExonList(exon_list_copy,i)->end);
  printSimpleList(elementAtExonList(exon_list_copy,i)->est_ids);
  printf("%s\n", elementAtExonList(exon_list_copy,i)->sequence);
  }*/
//exit(0);

  matx=(struct_matx*)calloc(n*n,sizeof(struct_matx));
  matx_pref_suff=(struct_matx_pref_suff *)calloc(n*n,sizeof(struct_matx_pref_suff));
  for(i=0;i<n;i++){
	 for(j=0;j<n;j++){
		getMatx(matx,i,j,n)->actual = -1;
		getMatx(matx,i,j,n)->flag=0;
		getMatx(matx,i,j,n)->est_ids=NULL;
	 }
  }

  structSimpleList *inters=NULL;
  structSimpleList *sub=NULL;
  for(i=0;i<n;i++){
	 for(j=i+1;j<n;j++){
		inters=intersection(elementAtExonList(exon_list,i)->est_ids,elementAtExonList(exon_list,j)->est_ids);

/*if(!isEmptySimpleList(inters)){
  fprintf(stdout, "INT %d and %d\n	(%s-%s %s-%s)\n", i, j, elementAtExonList(exon_list_copy,i)->start,elementAtExonList(exon_list_copy,i)->end, elementAtExonList(exon_list_copy,j)->start,elementAtExonList(exon_list_copy,j)->end);
  printSimpleList(inters);
  }*/

/*if(i == 59 && j == 138){
  fprintf(stdout, "INTERS*******\n");
  printSimpleList(inters);
  fprintf(stdout, "*******\n");
  fprintf(stdout, "EST LIST 1*******\n");
  printSimpleList(elementAtExonList(exon_list,i)->est_ids);
  fprintf(stdout, "*******\n");
  fprintf(stdout, "EST LIST 2*******\n");
  printSimpleList(elementAtExonList(exon_list,j)->est_ids);
  fprintf(stdout, "*******\n");
  }*/

		if(!isEmptySimpleList(inters)){
		  getMatx(matx,i,j,n)->flag=1;
		  getMatx(matx,i,j,n)->est_ids=inters;
		  sub=difference(elementAtExonList(exon_list,j)->est_ids,inters);
		  elementAtExonList(exon_list,j)->est_ids=sub;
		}
	 }
  }

//riempimento di matx_preff_suff
  structExonList *ei;
  structExonList *ej;
  for(i=0;i<n;i++){
	 for(j=0;j<n;j++){
		getMatxPrefSuff(matx_pref_suff,i,j,n)->flag=0;
		getMatxPrefSuff(matx_pref_suff,i,j,n)->actual=1;
	 }
  }

  i=0;
  while( i < n-1){
	 j=i+1;
	 int j_start,j_end,i_start,i_end;
	 ei=elementAtExonList(exon_list,i);
	 ej=elementAtExonList(exon_list,j);
	 j_start=atoi(ej->start);
	 j_end=atoi(ej->end);
	 i_start=atoi(ei->start);
	 i_end=atoi(ei->end);
	 while (j_start <= i_end && j<n){
		if (getMatx(matx,i,j,n)->flag == 0){
		  if (j_start == i_start){
			 if (j_end < i_end)
				getMatxPrefSuff(matx_pref_suff,i,j,n)->flag=-1;
			 else
				getMatxPrefSuff(matx_pref_suff,i,j,n)->flag=1;
		  }
		  if (j_end == i_end){
			 if (j_start > i_start)
				getMatxPrefSuff(matx_pref_suff,i,j,n)->flag=-2;
			 else
				getMatxPrefSuff(matx_pref_suff,i,j,n)->flag=2;
		  }
		  if (j_end < i_end && j_start > i_start)
			 getMatxPrefSuff(matx_pref_suff,i,j,n)->flag=-3;
		  if (j_end > i_end && j_start < i_start)
			 getMatxPrefSuff(matx_pref_suff,i,j,n)->flag=3;
		  if (j_end > i_end && j_start > i_start)
			 getMatxPrefSuff(matx_pref_suff,i,j,n)->flag=-4;
		}
		j++;
		if(j<n){
		  ej=elementAtExonList(exon_list,j);
		  j_start=atoi(ej->start);
		  j_end=atoi(ej->end);
		  i_start=atoi(ei->start);
		  i_end=atoi(ei->end);
		}
	 }
	 i++;
  }

/*	printf("\n");
	for(i=0;i<n;i++){
	for(j=0;j<n;j++)
	printf("%d ",getMatx(matx,i,j,n)->flag);
	printf("\t ");
	for(j=0;j<n;j++){
	printf("%d ",getMatxPrefSuff(matx_pref_suff,i,j,n)->flag);
	printf("\t ");
	}
	printf("\n");
	for(j=0;j<n;j++)
	printf("%d ",getMatxPrefSuff(matx_pref_suff,i,j,n)->actual);

	printf("\n"); 	}*/


/*for(i=0; i<n; i++)
  fprintf(stdout, " Exon %d left %s right %s\n", i, (elementAtExonList(exon_list,i)->start), (elementAtExonList(exon_list,i)->end));
  exit(0);*/

/*for(i=0; i<n; i++)
  for(j=i+1; j<n; j++)
//fprintf(stdout, "HERE (%d,%d) = %d\n", i,j,getMatxPrefSuff(matx_pref_suff,i,j,n)->flag);
fprintf(stdout, "HERE (%d,%d) = %d\n", i,j,getMatx(matx,i,j,n)->flag);
*/

/*Ricostruzione dei trascritti effettivi*/
  structSimpleList *transcripts[n];
  listOfList *tr_est_ids[n];

//19dic06
  structSimpleList *isRefSeqTranscript[n];

  char *tr_string=(char*)calloc(10000,sizeof(char));
  char *id=(char*)calloc(1000,sizeof(char));
  int k;
  int p;
  structSimpleList *tempList=NULL;
  for(i=0;i<n;i++){
	 strcpy(tr_string,"");     		//Ricostruzione di tutti i trascritti che iniziano dall'i-esimo esone e non sono composti dal solo esone i
	 transcripts[i]=NULL;
	 tr_est_ids[i]=NULL;

//19dic06
	 isRefSeqTranscript[i]=NULL;

	 for(j=i+1;j<n;j++){
		if(getMatx(matx,i,j,n)->flag==1){
		  inters = intersection(elementAtExonList(exon_list,i)->est_ids,getMatx(matx,i,j,n)->est_ids); //deve essere non vuota
		  for (t=0;t<lengthSimpleList(inters);t++){
			 strcpy(id,elementAtSimpleList(inters,t)->value);
			 k=findListOfList(tr_est_ids[i],id);
			 if(k==-1){
				strcpy(tr_string,itoa(i));
				strcat(tr_string,".");
				strcat(tr_string,itoa(j));
				if (findSimpleList(transcripts[i],tr_string)==-1){
				  transcripts[i]=pushSimpleList(transcripts[i],tr_string);
				  tempList=NULL;
				  tempList=pushSimpleList(tempList,id);
				  tr_est_ids[i]=pushListOfList(tr_est_ids[i],tempList);
				}
				else{
				  p=findSimpleList(transcripts[i],tr_string);
				  elementAtListOfList(tr_est_ids[i],p)->local=pushSimpleList(elementAtListOfList(tr_est_ids[i],p)->local,id);
				}
			 }
			 else{
				strcpy(tr_string,elementAtSimpleList(transcripts[i],k)->value);
				strcat(tr_string,".");
				strcat(tr_string,itoa(j));
				if (findSimpleList(transcripts[i],tr_string)==-1){
				  transcripts[i]=pushSimpleList(transcripts[i],tr_string);
				  tempList=NULL;
				  tempList=pushSimpleList(tempList,id);
				  tr_est_ids[i]=pushListOfList(tr_est_ids[i],tempList);
				}
				else{
				  p=findSimpleList(transcripts[i],tr_string);
				  elementAtListOfList(tr_est_ids[i],p)->local=pushSimpleList(elementAtListOfList(tr_est_ids[i],p)->local,id);
				}
				elementAtListOfList(tr_est_ids[i],k)->local=deleteSimpleList((elementAtListOfList(tr_est_ids[i],k)->local),id);
				if(lengthSimpleList(elementAtListOfList(tr_est_ids[i],k)->local)==0){
				  char *value=elementAtSimpleList(transcripts[i],k)->value;
				  transcripts[i]=deleteSimpleList(transcripts[i],value);
				  tr_est_ids[i]=deleteListOfList(tr_est_ids[i],k);
				}
			 }
		  }
		} 		}
	 structSimpleList *one_ex_ids=NULL;
	 structSimpleList *unione=NULL;
	 int ii;
	 for(ii=0;ii<lengthListOfList(tr_est_ids[i]);ii++)
		unione=unionSimpleList(unione,elementAtListOfList(tr_est_ids[i],ii)->local);
	 one_ex_ids=difference(elementAtExonList(exon_list,i)->est_ids,unione);
	 if (lengthSimpleList(one_ex_ids)!=0){
		strcpy(tr_string,itoa(i));
		transcripts[i]=pushSimpleList(transcripts[i],tr_string);
		tr_est_ids[i]=pushListOfList(tr_est_ids[i],one_ex_ids);
	 }
  }
  structure2 *end_transcripts[n];
  listOfList *end_tr_est_ids[n];

//19dic06
  structSimpleList *end_isRefSeqTranscript[n];

  char *last;  	char *transcript;
  for(i=0;i<n;i++){
	 end_transcripts[i]=(structure2*)calloc(1,sizeof(structure2));
	 end_transcripts[i]->transcript=NULL;
	 end_transcripts[i]->coord1=NULL;
	 end_transcripts[i]->coord2=NULL;
	 end_tr_est_ids[i]=NULL;

//19dic06
	 end_isRefSeqTranscript[i]=NULL;
  }

  for(i=0;i<n;i++){
	 int l=lengthSimpleList(transcripts[i]);
	 for(j=0;j<l;j++){
		transcript=elementAtSimpleList(transcripts[i],j)->value;

		int lunghezza=lengthTranscript(transcript)-1;

		if(lunghezza == 0 && elementAtExonList(exon_list_copy, i)->external == 1)
		  elementAtExonList(exon_list_copy, i)->external=-3;

//RAFFA 4giu04
		if(elementAtExonList(exon_list_copy, i)->external != 0 && lunghezza != 0){
		  if(elementAtExonList(exon_list_copy, i)->external == -2)
			 elementAtExonList(exon_list_copy, i)->external=0;
		  else
			 elementAtExonList(exon_list_copy, i)->external=-1;
		}

//RAFFA 30mag04
		int p;
		int exon_index;
		for(p=0; p<=lunghezza; p++){
		  exon_index=atoi(elementAtTranscript(transcript,p));
		  elementAtExonList(exon_list_copy, exon_index)->first_tr_coord=pushSimpleList(elementAtExonList(exon_list_copy, exon_index)->first_tr_coord, itoa(i));
		  elementAtExonList(exon_list_copy, exon_index)->second_tr_coord=pushSimpleList(elementAtExonList(exon_list_copy, exon_index)->second_tr_coord, itoa(j));
		  elementAtExonList(exon_list_copy, exon_index)->third_tr_coord=pushSimpleList(elementAtExonList(exon_list_copy, exon_index)->third_tr_coord, itoa(p));

		  if(p != 0 && p != lunghezza){
			 elementAtExonList(exon_list_copy, exon_index)->external=0;
		  }
		}

		last=elementAtTranscript(transcript,lunghezza);

		if(elementAtExonList(exon_list_copy, atoi(last))->polya == 1)
		  elementAtExonList(exon_list_copy, atoi(last))->polya=2;

//RAFFA 4giu04
		if(elementAtExonList(exon_list_copy, atoi(last))->external != 0 && lunghezza != 0){
		  if(elementAtExonList(exon_list_copy, atoi(last))->external == -1)
			 elementAtExonList(exon_list_copy, atoi(last))->external=0;
		  else
			 elementAtExonList(exon_list_copy, atoi(last))->external=-2;
		}

		end_transcripts[atoi(last)]->transcript=pushSimpleList(end_transcripts[atoi(last)]->transcript,transcript);
		end_transcripts[atoi(last)]->coord1=pushSimpleList(end_transcripts[atoi(last)]->coord1,itoa(i));
		end_transcripts[atoi(last)]->coord2=pushSimpleList(end_transcripts[atoi(last)]->coord2,itoa(j));
		if(lengthListOfList(tr_est_ids[i])!=0)
		  end_tr_est_ids[atoi(last)]=pushListOfList(end_tr_est_ids[atoi(last)],tr_est_ids[i]->local);


//19dic06
		structSimpleList *s_EST_ids=tr_est_ids[i]->local;
		char isRefSeq=0;
		while(s_EST_ids != NULL && !isRefSeq){
		  if(s_EST_ids->value != NULL) {
		  	/* UPDATE for noncoding RefSeq
		  	*/
			/*if(s_EST_ids->value[0] == 'N' && s_EST_ids->value[1] == 'M'){
				isRefSeq=1;
			 }*/
		  	if(s_EST_ids->value[0] == 'N' && s_EST_ids->value[2] == '_'){
		  		if(s_EST_ids->value[1] == 'M' || s_EST_ids->value[1] == 'R'){
					isRefSeq=1;
		  		}
		  	}	  
		  }
		  s_EST_ids=s_EST_ids->next;
		}
		if(isRefSeq)
		  end_isRefSeqTranscript[atoi(last)]=pushSimpleList(end_isRefSeqTranscript[atoi(last)], "1");
		else
		  end_isRefSeqTranscript[atoi(last)]=pushSimpleList(end_isRefSeqTranscript[atoi(last)], "0");
	 }
  }

/*for(i=0;i<n;i++){
  if(elementAtExonList(exon_list_copy,i)->polya == 2)
  printf("%d) %s:%s\n",i, elementAtExonList(exon_list_copy,i)->start,elementAtExonList(exon_list_copy,i)->end);
  }*/

//Verifica delle forme di competing
  for(i=0;i<n;i++)
	 for(j=i+1;j<n;j++) {
		if (getMatxPrefSuff(matx_pref_suff,i,j,n)->flag!=0){
		  if (abs(getMatxPrefSuff(matx_pref_suff,i,j,n)->flag)==1)
			 if(!verify_competing53(end_tr_est_ids,i,j,n))
				getMatxPrefSuff(matx_pref_suff,i,j,n)->actual=0;

		  if(abs(getMatxPrefSuff(matx_pref_suff,i,j,n)->flag)==2){
			 if(!verify_competing35(tr_est_ids,i,j,n))
				getMatxPrefSuff(matx_pref_suff,i,j,n)->actual=0;
		  }


		  if(abs(getMatxPrefSuff(matx_pref_suff,i,j,n)->flag)==3){
			 if(!verify_competing53(end_tr_est_ids,i,j,n))
				getMatxPrefSuff(matx_pref_suff,i,j,n)->actual=2;

			 if(verify_competing35(tr_est_ids,i,j,n)){
				if(getMatxPrefSuff(matx_pref_suff,i,j,n)->actual==1)
				  getMatxPrefSuff(matx_pref_suff,i,j,n)->actual=3;
			 }
			 else
				if(getMatxPrefSuff(matx_pref_suff,i,j,n)->actual==2)
				  getMatxPrefSuff(matx_pref_suff,i,j,n)->actual=0;

		  }

		  if (getMatxPrefSuff(matx_pref_suff,i,j,n)->flag == -4) {
			 if (!verify_competing53(tr_est_ids,i,j,n))
				getMatxPrefSuff(matx_pref_suff,i,j,n)->actual=2;
			 if (verify_competing35(tr_est_ids,i,j,n)){
				if (getMatxPrefSuff(matx_pref_suff,i,j,n)->actual==1)
				  getMatxPrefSuff(matx_pref_suff,i,j,n)->actual=3;
			 }
			 else{
				if (getMatxPrefSuff(matx_pref_suff,i,j,n)->actual==2)
				  getMatxPrefSuff(matx_pref_suff,i,j,n)->actual=0;
			 }

		  }
		}
	 }

/*for(i=0; i<n; i++)
  for(j=i+1; j<n; j++)
  fprintf(stderr, "HERE (%d,%d) = (%d,%d)\n", i,j,getMatxPrefSuff(matx_pref_suff,i,j,n)->flag, getMatxPrefSuff(matx_pref_suff,i,j,n)->actual);
*/

//19dic06
  for(i=0;i<n;i++){
	 for(k=0;k<lengthSimpleList(transcripts[i]);k++){
		structSimpleList *s_EST_ids=elementAtListOfList(tr_est_ids[i],k)->local;
		char isRefSeq=0;
		while(s_EST_ids != NULL && !isRefSeq){
		  if(s_EST_ids->value != NULL) {
			/* UPDATE for noncoding RefSeq
		  	*/
		  	/*if(s_EST_ids->value[0] == 'N' && s_EST_ids->value[1] == 'M'){
				isRefSeq=1;
			 }*/
		  	if(s_EST_ids->value[0] == 'N' && s_EST_ids->value[2] == '_'){
		  		if(s_EST_ids->value[1] == 'M' || s_EST_ids->value[1] == 'R'){
					isRefSeq=1;
		  		}
		  	}	  
		  }
		  s_EST_ids=s_EST_ids->next;
		}
		if(isRefSeq)
		  isRefSeqTranscript[i]=pushSimpleList(isRefSeqTranscript[i], "1");
		else
		  isRefSeqTranscript[i]=pushSimpleList(isRefSeqTranscript[i], "0");
	 }
  }


/*STAMPA DEI TRASCRITTI*/
/*	for(i=0;i<n;i++){
	printf("Trascritti che iniziano con %d\n", i);
	for(k=0;k<lengthSimpleList(transcripts[i]);k++){
	char *transcript=elementAtSimpleList(transcripts[i],k)->value;
	printf("%s\n",transcript);
	for(j=0;j<lengthTranscript(transcript);j++){
	int index=atoi(elementAtTranscript(transcript,j));
	char *sequence=elementAtExonList(exon_list,index)->sequence;
//			printf("%s \n",sequence);
}
printSimpleList(elementAtListOfList(tr_est_ids[i],k)->local);

//19dic06
char *isMrna=elementAtSimpleList(isRefSeqTranscript[i],k)->value;
printf("IS refSeq %s\n", isMrna);

//19dic06
//structSimpleList *s_EST_ids=elementAtListOfList(tr_est_ids[i],k)->local;
//while(s_EST_ids != NULL){
//	if(s_EST_ids->value != NULL) {
//		printf("value raffa %s\n",s_EST_ids->value);
//	}
//	s_EST_ids=s_EST_ids->next;
//}
}
}
exit(0);*/
/*	i=111;
	for(k=0;k<lengthSimpleList(transcripts[i]);k++){
	char *transcript=elementAtSimpleList(transcripts[i],k)->value;
	printf("%s\n",transcript);
	for(j=0;j<lengthTranscript(transcript);j++){
	int index=atoi(elementAtTranscript(transcript,j));
	char *sequence=elementAtExonList(exon_list,index)->sequence;
//			printf("%s \n",sequence);
}
//		printSimpleList(elementAtListOfList(tr_est_ids[i],k)->local);
}
//exit(0);*/

  structure *left_keep_transcript[n];
  structure *right_keep_transcript[n];

/*nuovocodice2*/

  for(i=0;i<n;i++){
	 left_keep_transcript[i]=NULL;
	 right_keep_transcript[i]=NULL;
  }
  for(i=0;i<n;i++){
	 for(j=0;j<lengthSimpleList(transcripts[i]);j++){
//RAFFA 30mag04
		left_keep_transcript[i]=pushStructure(left_keep_transcript[i],-1,-1,-1);
		right_keep_transcript[i]=pushStructure(right_keep_transcript[i],-1,-1,-1);
//left_keep_transcript[i]=pushStructure(left_keep_transcript[i],-1,-1);
//right_keep_transcript[i]=pushStructure(right_keep_transcript[i],-1,-1);
	 }

  }

//RAFFA 30mag04 PARTE COMMENTATA
/*	for(i=0;i<n;i++)
	for(j=i+1;j<n;j++){
	int index;
	int matx_pref_suff_ij=getMatxPrefSuff(matx_pref_suff,i,j,n)->flag;
	int actual=getMatxPrefSuff(matx_pref_suff,i,j,n)->actual;
	if(abs(matx_pref_suff_ij)==2){
	if(actual==0){
	if(matx_pref_suff_ij==2)
	index=i;
	else
	index=j;
	for(k=0;k<lengthSimpleList(transcripts[i]);k++){
	if (equalsStructure(left_keep_transcript[i],k,-1,-1)){
	int p;
	for(p=0;p<lengthSimpleList(transcripts[j]);p++){
	char *transcript_ik=elementAtSimpleList(transcripts[i],k)->value;
	char *transcript_jp=elementAtSimpleList(transcripts[j],p)->value;
	char *intern_ik=getInternTranscript(transcript_ik);
	char *intern_jp=getInternTranscript(transcript_jp);
	cmp=0;
	if(intern_ik == NULL && intern_jp == NULL)
	cmp=1;
	else
	if(intern_ik != NULL & intern_jp != NULL)
	cmp=(!strcmp(intern_ik,intern_jp))?(1):(0);
	if(cmp){
// porre reduce a true dell'indeex-esimo esone in ex list cpy
elementAtExonList(exon_list_copy,index)->reduce=1;

//if(intern_ik!=NULL && intern_jp!=NULL && !strcmp(intern_ik,intern_jp)){
if(index==i)
setStructure(left_keep_transcript[index],k,j,p);
else
setStructure(left_keep_transcript[index],p,i,k);
if(!strcmp(lastExon(transcripts[i],k),lastExon(transcripts[j],p))){
if(index==i){
structure *temp3=getStructure(left_keep_transcript[index],k);
int v1=temp3->v1;
int v2=temp3->v2;
setStructure(right_keep_transcript[index],k,v1,v2);
}
else{
structure *temp3=getStructure(left_keep_transcript[index],p);
int v1=temp3->v1;
int v2=temp3->v2;
setStructure(right_keep_transcript[index],p,v1,v2);
}
}
}
}
}
}
}
}
}


//for(i=0;i<n;i++)
//for(j=0;j<lengthStructure(left_keep_transcript[i]);j++)
for(i=0;i<n;i++)
for(j=i+1;j<n;j++){
int index;
int matx_pref_suff_ij=getMatxPrefSuff(matx_pref_suff,i,j,n)->flag;
int actual=getMatxPrefSuff(matx_pref_suff,i,j,n)->actual;
if(abs(matx_pref_suff_ij)==1){
if(actual==0){
if(matx_pref_suff_ij==1)
index=i;
else
index=j;
for(k=0;k<lengthSimpleList(end_transcripts[i]->transcript);k++){
int w=atoi(elementAtSimpleList(end_transcripts[i]->coord1,k)->value);
int z=atoi(elementAtSimpleList(end_transcripts[i]->coord2,k)->value);
if (equalsStructure(right_keep_transcript[w],z,-1,-1)){
int p;
for(p=0;p<lengthSimpleList(end_transcripts[j]->transcript);p++){
int x=atoi(elementAtSimpleList(end_transcripts[j]->coord1,p)->value);
int y=atoi(elementAtSimpleList(end_transcripts[j]->coord2,p)->value);
char *end_transcript_ik=elementAtSimpleList(end_transcripts[i]->transcript,k)->value;
char *end_transcript_jp=elementAtSimpleList(end_transcripts[j]->transcript,p)->value;
char *intern_ik=getInternTranscript(end_transcript_ik);
char *intern_jp=getInternTranscript(end_transcript_jp);
cmp=0;
if(intern_ik == NULL && intern_jp == NULL)
cmp=1;
else
if(intern_ik != NULL & intern_jp != NULL)
cmp=(!strcmp(intern_ik,intern_jp))?(1):(0);


if(cmp){
//if(intern_ik!=NULL && intern_jp!=NULL && !strcmp(intern_ik,intern_jp)){
// porre riduce a true dell'indeex-esimo esone in ex list cpy
elementAtExonList(exon_list_copy,index)->reduce=1;
if(index==i){
setStructure(right_keep_transcript[w],z,x,y);
}
else
setStructure(right_keep_transcript[x],y,w,z);

if(!strcmp(firstExon(end_transcripts[i]->transcript,k),firstExon(end_transcripts[j]->transcript,p))){ 											if(index==i){
structure *temp3=getStructure(right_keep_transcript[w],z);
int v1=temp3->v1;
int v2=temp3->v2;
setStructure(left_keep_transcript[w],z,v1,v2);
}
else{
structure *temp3=getStructure(right_keep_transcript[x],y);
int v1=temp3->v1;
int v2=temp3->v2; 												setStructure(left_keep_transcript[x],y,v1,v2); 											 											}
}
}
}
}
}
}
}
} */
//RAFFA 30mag04 PARTE COMMENTATA

//RAFFA 30mag04 MODIFICA DELLA PARTE COMMENTATA
  for(i=0;i<n;i++)
	 for(j=i+1;j<n;j++){

		int index_to_be_red;	//Indice dell'esone da ridurre
		int index_of_red;		//Indice dell'esone al quale ridurre

		int matx_pref_suff_ij=getMatxPrefSuff(matx_pref_suff,i,j,n)->flag;
		int actual=getMatxPrefSuff(matx_pref_suff,i,j,n)->actual;
		if(abs(matx_pref_suff_ij)==2){
		  if(actual==0){
			 if(matx_pref_suff_ij==2){
				index_to_be_red=i;
				index_of_red=j;
			 }
			 else{
				index_to_be_red=j;
				index_of_red=i;
			 }

//Considero tutti i trascritti che iniziano con index_to_be_red
			 for(k=0;k<lengthSimpleList(transcripts[index_to_be_red]);k++){

//19dic06
				char extends_first_exon=1;
#ifdef DONT_EXTEND_REFSEQ
				char *isMrna=elementAtSimpleList(isRefSeqTranscript[index_to_be_red],k)->value;
				if(!strcmp(isMrna, "1"))
				  extends_first_exon=0;
#endif

//19dic06
//if (equalsStructure(left_keep_transcript[index_to_be_red],k,-1,-1,-1)){
				if(extends_first_exon && equalsStructure(left_keep_transcript[index_to_be_red],k,-1,-1,-1)){
				  char *transcript_ik=elementAtSimpleList(transcripts[index_to_be_red],k)->value;
				  int p;
				  int exon_to_be_compared;
				  if(lengthTranscript(transcript_ik) > 1){
					 exon_to_be_compared=atoi(elementAtTranscript(transcript_ik,1));

//Considero tutti i trascritti che contengono index_of_red
					 structSimpleList *first_tr_coord=elementAtExonList(exon_list_copy, index_of_red)->first_tr_coord;
					 structSimpleList *second_tr_coord=elementAtExonList(exon_list_copy, index_of_red)->second_tr_coord;
					 structSimpleList *third_tr_coord=elementAtExonList(exon_list_copy, index_of_red)->third_tr_coord;

					 for(p=0;p<lengthSimpleList(first_tr_coord);p++){
						int first=atoi(elementAtSimpleList(first_tr_coord,p)->value);
						int second=atoi(elementAtSimpleList(second_tr_coord,p)->value);
						int third=atoi(elementAtSimpleList(third_tr_coord,p)->value);
						char *transcript_jp=elementAtSimpleList(transcripts[first],second)->value;

//21feb07
						char extends_first_exon2=1;
#ifdef DONT_EXTEND_REFSEQ
						char *isMrna2=elementAtSimpleList(isRefSeqTranscript[first],second)->value;
						if(!strcmp(isMrna2, "1"))
						  extends_first_exon2=0;
#endif

//21feb07
						if(extends_first_exon2){

						  cmp=0;

						  int lgth=lengthTranscript(transcript_jp)-1;
//Se l'esone index_of_red non e' l'ultimo del trascritto
						  if(third<lgth){
							 int third_following=third+1;
							 int exon_following=atoi(elementAtTranscript(transcript_jp,third_following));
							 if(exon_following==exon_to_be_compared)
								cmp=1;
							 else{
								int first_index;
								int second_index;
								if(exon_following < exon_to_be_compared){
								  first_index=exon_following;
								  second_index=exon_to_be_compared;
								}
								else{
								  first_index=exon_to_be_compared;
								  second_index=exon_following;
								}
								int matx_pref_suff_ij_1=getMatxPrefSuff(matx_pref_suff,first_index,second_index,n)->flag;
								int actual_1=getMatxPrefSuff(matx_pref_suff,first_index,second_index,n)->actual;
								if(abs(matx_pref_suff_ij_1)==1){
								  if(actual_1==0){
									 cmp=1;
									 if(actual_1 == 1){
										if(first_index == exon_following)
										  setStructure(right_keep_transcript[first],second,index_to_be_red,k,1);
										else
										  setStructure(right_keep_transcript[index_to_be_red],k,first,second,third_following);
									 }
									 else{
										if(second_index == exon_following){
										  setStructure(right_keep_transcript[first],second,index_to_be_red,k,1);
										}
										else
										  setStructure(right_keep_transcript[index_to_be_red],k,first,second,third_following);
									 }
								  }
								}
							 }
						  }

						  if(cmp){
							 setStructure(left_keep_transcript[index_to_be_red],k,first,second,third);

							 if(!strcmp(lastExon(transcripts[index_to_be_red],k),lastExon(transcripts[first],second)))
								setStructure(right_keep_transcript[index_to_be_red],k,first,second,lgth);
						  }
//21feb07
						}
					 }
				  }
				}
			 }
		  }
		}
	 }

  for(i=0;i<n;i++)
	 for(j=i+1;j<n;j++){

		int index_to_be_red;	//Indice dell'esone da ridurre
		int index_of_red;		//Indice dell'esone al quale ridurre

		int matx_pref_suff_ij=getMatxPrefSuff(matx_pref_suff,i,j,n)->flag;
		int actual=getMatxPrefSuff(matx_pref_suff,i,j,n)->actual;
		if(abs(matx_pref_suff_ij)==1){
		  if(actual==0){
			 if(matx_pref_suff_ij==1){
				index_to_be_red=i;
				index_of_red=j;
			 }
			 else{
				index_to_be_red=j;
				index_of_red=i;
			 }

			 for(k=0;k<lengthSimpleList(end_transcripts[index_to_be_red]->transcript);k++){
				int w=atoi(elementAtSimpleList(end_transcripts[index_to_be_red]->coord1,k)->value);
				int z=atoi(elementAtSimpleList(end_transcripts[index_to_be_red]->coord2,k)->value);

//19dic06
				char extends_last_exon=1;
#ifdef DONT_EXTEND_REFSEQ
				char *isMrna=elementAtSimpleList(end_isRefSeqTranscript[index_to_be_red],k)->value;
				if(!strcmp(isMrna, "1"))
				  extends_last_exon=0;
#endif

//19dic06
//if (equalsStructure(right_keep_transcript[w],z,-1,-1,-1)){
				if(extends_last_exon && equalsStructure(right_keep_transcript[w],z,-1,-1,-1)){
				  char *end_transcript_ik=elementAtSimpleList(end_transcripts[index_to_be_red]->transcript,k)->value;
				  int p;
				  int exon_to_be_compared;
				  if(lengthTranscript(end_transcript_ik) > 1){
					 exon_to_be_compared=atoi(elementAtTranscript(end_transcript_ik,lengthTranscript(end_transcript_ik)-2));

//Considero tutti i trascritti che contengono index_of_red
					 structSimpleList *first_tr_coord=elementAtExonList(exon_list_copy, index_of_red)->first_tr_coord;
					 structSimpleList *second_tr_coord=elementAtExonList(exon_list_copy, index_of_red)->second_tr_coord;
					 structSimpleList *third_tr_coord=elementAtExonList(exon_list_copy, index_of_red)->third_tr_coord;

					 for(p=0;p<lengthSimpleList(first_tr_coord);p++){
						int first=atoi(elementAtSimpleList(first_tr_coord,p)->value);
						int second=atoi(elementAtSimpleList(second_tr_coord,p)->value);
						int third=atoi(elementAtSimpleList(third_tr_coord,p)->value);
						char *transcript_jp=elementAtSimpleList(transcripts[first],second)->value;

//21feb07
						char extends_last_exon2=1;
#ifdef DONT_EXTEND_REFSEQ
						char *isMrna2=elementAtSimpleList(isRefSeqTranscript[first],second)->value;
						if(!strcmp(isMrna2, "1"))
						  extends_last_exon2=0;
#endif

//21feb07
						if(extends_last_exon2){

//6ott06 (Per disattivare la verifica sugli esoni precedenti)
//cmp=0;
						  cmp=1;

						  //int lgth=lengthTranscript(transcript_jp)-1;
//Se l'esone index_of_red non e' il primo del trascritto
						  if(third>0){
							 int third_prev=third-1;
							 int exon_prev=atoi(elementAtTranscript(transcript_jp,third_prev));
							 if(exon_prev==exon_to_be_compared)
								cmp=1;
						  }

						  if(cmp){
							 setStructure(right_keep_transcript[w],z,first,second,third);

							 if(!strcmp(firstExon(transcripts[w],z),firstExon(transcripts[first],second)))
								setStructure(left_keep_transcript[w],z,first,second,0);
						  }
//21feb07
						}
					 }
				  }
				}
			 }
		  }
		}
	 }
//RAFFA 30mag04 MODIFICA DELLA PARTE COMMENTATA

/*	i=57;
//for(i=0;i<n;i++)
for(k=0;k<lengthSimpleList(transcripts[i]);k++){
fprintf(stderr, "HERE2 %d-%d = %d,%d\n", i,k,elementAtStructure(left_keep_transcript[i],k)->v1, elementAtStructure(left_keep_transcript[i],k)->v2);
fprintf(stderr, "HERE3 %d-%d = %d,%d\n", i,k,elementAtStructure(right_keep_transcript[i],k)->v1, elementAtStructure(right_keep_transcript[i],k)->v2);
fprintf(stderr, "HERE4 %d-%d = %d,%d\n", i,k,elementAtStructure(left_keep_transcript[i],k)->v3, elementAtStructure(right_keep_transcript[i],k)->v3);
}
exit(0);*/

  listOfList *global_tr_est_ids=NULL;
  structSimpleList *global_transcripts=NULL;

//19dic06
#ifdef DONT_EXTEND_REFSEQ
  structSimpleList *isRefSeqTranscriptForGlobal=NULL;
#endif

  structSimpleList *temp_global_transcripts=NULL;
  listOfList *temp_global_tr_est_ids=NULL; 	structSimpleList *duplicate=NULL;

//RAFFA 30mag04 PARTE COMMENTATA
/*	for(i=0;i<n;i++){
	for(k=0;k<lengthSimpleList(transcripts[i]);k++){
	tr_string=elementAtSimpleList(transcripts[i],k)->value;
	int coord1=elementAtStructure(left_keep_transcript[i],k)->v1;
	int coord2=elementAtStructure(left_keep_transcript[i],k)->v2;
	int coord1_1=i;
	int coord2_1=k;
	while(coord1!=-1 && coord2 !=-1){
	coord1_1=coord1;
	coord2_1=coord2;
	int tempCoord1;
	tempCoord1=elementAtStructure(left_keep_transcript[coord1],coord2)->v1;
	coord2=elementAtStructure(left_keep_transcript[coord1],coord2)->v2;
	coord1=tempCoord1;
	}

	char *gen_tr_string=elementAtSimpleList(transcripts[coord1_1],coord2_1)->value;
	char *tempString=(char*)calloc(1000,sizeof(char));
	strcpy(tempString,"");
	int l=lengthTranscript(tr_string);
	int ii;
	char *value=elementAtTranscript(gen_tr_string,0);
	strcat(tempString,value);
	if (l>1) strcat(tempString,".");
	for(ii=1;ii<l;ii++){
	strcat(tempString,elementAtTranscript(tr_string,ii));
	if (ii<l-1) strcat(tempString,".");
	}
	strcpy(tr_string,tempString);

	coord1=elementAtStructure(right_keep_transcript[i],k)->v1;
	coord2=elementAtStructure(right_keep_transcript[i],k)->v2;
	coord1_1=i;
	coord2_1=k;

	while(coord1!=-1 && coord2 !=-1){
	coord1_1=coord1;
	coord2_1=coord2;
	int tempCoord1;
	tempCoord1=elementAtStructure(right_keep_transcript[coord1],coord2)->v1;
	coord2=elementAtStructure(right_keep_transcript[coord1],coord2)->v2;
	coord1=tempCoord1;
	}

	gen_tr_string=elementAtSimpleList(transcripts[coord1_1],coord2_1)->value;
	l=lengthTranscript(tr_string);
	strcpy(tempString,"");
	value=elementAtTranscript(gen_tr_string,lengthTranscript(gen_tr_string)-1);
	for(ii=0;ii<l-1;ii++){
	strcat(tempString,elementAtTranscript(tr_string,ii));
	strcat(tempString,".");
	}
	strcat(tempString,value);
	strcpy(tr_string,tempString);

	elementAtSimpleList(transcripts[i],k)->value=tr_string;
	temp_global_transcripts=pushSimpleList(temp_global_transcripts,tr_string);
	temp_global_tr_est_ids=pushListOfList(temp_global_tr_est_ids,elementAtListOfList(tr_est_ids[i],k)->local);
	duplicate=pushSimpleList(duplicate,"0");
	}
	}*/
//RAFFA 30mag04 PARTE COMMENTATA

//RAFFA 30mag04 MODIFICA DELLA PARTE COMMENTATA
  for(i=0;i<n;i++){
	 for(k=0;k<lengthSimpleList(transcripts[i]);k++){
		tr_string=elementAtSimpleList(transcripts[i],k)->value;

		int coord1=elementAtStructure(left_keep_transcript[i],k)->v1;
		int coord2=elementAtStructure(left_keep_transcript[i],k)->v2;
		int coord3=elementAtStructure(left_keep_transcript[i],k)->v3;

		int coord1_1=(coord1 == -1)?(i):(coord1);
		int coord2_1=(coord2 == -1)?(k):(coord2);
		int coord3_1=(coord3 == -1)?(0):(coord3);

		while(coord1!=-1 && coord2 !=-1 && coord3 == 0){
		  coord1_1=coord1;
		  coord2_1=coord2;
		  coord3_1=coord3;
		  int tempCoord1;
		  int tempCoord2;
		  tempCoord1=elementAtStructure(left_keep_transcript[coord1],coord2)->v1;
		  tempCoord2=elementAtStructure(left_keep_transcript[coord1],coord2)->v2;
		  coord3=elementAtStructure(left_keep_transcript[coord1],coord2)->v3;
		  coord1=tempCoord1;
		  coord2=tempCoord2;
		  if(coord3 > 0){
			 coord1_1=coord1;
			 coord2_1=coord2;
			 coord3_1=coord3;
		  }
		}

		char *gen_tr_string=elementAtSimpleList(transcripts[coord1_1],coord2_1)->value;
		char *tempString=(char*)calloc(1000,sizeof(char));
		strcpy(tempString,"");
		int l=lengthTranscript(tr_string);
		int ii;
		char *value=elementAtTranscript(gen_tr_string,coord3_1);
		strcat(tempString,value);
		if (l>1)
		  strcat(tempString,".");
		for(ii=1;ii<l;ii++){
		  strcat(tempString,elementAtTranscript(tr_string,ii));
		  if (ii<l-1) strcat(tempString,".");
		}
		strcpy(tr_string,tempString);

		coord1=elementAtStructure(right_keep_transcript[i],k)->v1;
		coord2=elementAtStructure(right_keep_transcript[i],k)->v2;
		coord3=elementAtStructure(right_keep_transcript[i],k)->v3;

		coord1_1=(coord1 == -1)?(i):(coord1);
		coord2_1=(coord2 == -1)?(k):(coord2);
		coord3_1=(coord3 == -1)?(lengthTranscript(elementAtSimpleList(transcripts[i],k)->value)-1):(coord3);

		while(coord1!=-1 && coord2 !=-1 && coord3 == lengthTranscript(elementAtSimpleList(transcripts[coord1_1],coord2_1)->value)-1){
		  coord1_1=coord1;
		  coord2_1=coord2;
		  coord3_1=coord3;
		  int tempCoord1;
		  int tempCoord2;
		  tempCoord1=elementAtStructure(right_keep_transcript[coord1],coord2)->v1;
		  tempCoord2=elementAtStructure(right_keep_transcript[coord1],coord2)->v2;
		  coord3=elementAtStructure(right_keep_transcript[coord1],coord2)->v3;
		  coord1=tempCoord1;
		  coord2=tempCoord2;
		  if((coord1 != -1 && coord2 != -1) && coord3 < lengthTranscript(elementAtSimpleList(transcripts[coord1],coord2)->value)-1 && coord3 != -1){
			 coord1_1=coord1;
			 coord2_1=coord2;
			 coord3_1=coord3;
		  }
		}

		gen_tr_string=elementAtSimpleList(transcripts[coord1_1],coord2_1)->value;

		l=lengthTranscript(tr_string);
		strcpy(tempString,"");
		value=elementAtTranscript(gen_tr_string,coord3_1);
		for(ii=0;ii<l-1;ii++){
		  strcat(tempString,elementAtTranscript(tr_string,ii));
		  strcat(tempString,".");
		}
		strcat(tempString,value);
		strcpy(tr_string,tempString);

		elementAtSimpleList(transcripts[i],k)->value=tr_string;
		temp_global_transcripts=pushSimpleList(temp_global_transcripts,tr_string);
		temp_global_tr_est_ids=pushListOfList(temp_global_tr_est_ids,elementAtListOfList(tr_est_ids[i],k)->local);
		duplicate=pushSimpleList(duplicate,"0");
	 }
  }
//RAFFA 30mag04 MODIFICA DELLA PARTE COMMENTATA

  for(i=0;i<lengthSimpleList(temp_global_transcripts);i++){
	 if(!strcmp(elementAtSimpleList(duplicate,i)->value,"0")){
		structSimpleList *est_ids=NULL;
		char *value=elementAtSimpleList(temp_global_transcripts,i)->value;
		global_transcripts=pushSimpleList(global_transcripts,value);

		structSimpleList *tempList=elementAtListOfList(temp_global_tr_est_ids,i)->local;
		for(j=0;j<lengthSimpleList(tempList);j++)
		  est_ids=pushSimpleList(est_ids,elementAtSimpleList(tempList,j)->value);
		for(k=i+1;k<lengthSimpleList(temp_global_transcripts);k++){
		  char *t_g_t_k=elementAtSimpleList(temp_global_transcripts,k)->value;
		  char *t_g_t_i=elementAtSimpleList(temp_global_transcripts,i)->value;
		  if(!strcmp(elementAtSimpleList(duplicate,k)->value,"0") && !strcmp(t_g_t_k,t_g_t_i)){
			 tempList=elementAtListOfList(temp_global_tr_est_ids,k)->local;
			 for(j=0;j<lengthSimpleList(tempList);j++)
				est_ids=pushSimpleList(est_ids,elementAtSimpleList(tempList,j)->value);
			 strcpy(elementAtSimpleList(duplicate,k)->value,"1");
		  }
		}
		global_tr_est_ids=pushListOfList(global_tr_est_ids,est_ids);

//19dic06
#ifdef DONT_EXTEND_REFSEQ
		structSimpleList *s_EST_ids=est_ids;
		char isRefSeq=0;
		while(s_EST_ids != NULL && !isRefSeq){
		  if(s_EST_ids->value != NULL) {
			/* UPDATE for noncoding RefSeq
		  	*/
		  	/*if(s_EST_ids->value[0] == 'N' && s_EST_ids->value[1] == 'M'){
				isRefSeq=1;
			 }*/
		  	if(s_EST_ids->value[0] == 'N' && s_EST_ids->value[2] == '_'){
		  		if(s_EST_ids->value[1] == 'M' || s_EST_ids->value[1] == 'R'){
					isRefSeq=1;
		  		}
		  	}	  
		  }
		  s_EST_ids=s_EST_ids->next;
		}
		if(isRefSeq)
		  isRefSeqTranscriptForGlobal=pushSimpleList(isRefSeqTranscriptForGlobal, "1");
		else
		  isRefSeqTranscriptForGlobal=pushSimpleList(isRefSeqTranscriptForGlobal, "0");
#endif
	 }
  }

//output
//per francesca
  FILE *out=fopen(outputfile,"w");

#ifdef READ_ABS_COORD
  fprintf(out,"%ld\n",gen_start);
  fprintf(out,"%ld\n",gen_end);
  fprintf(out,"%d\n",strand);
  fprintf(out,"%ld\n",boundary);
#endif

//22mar07
  for(i=0;i<n;i++){
	 elementAtExonList(exon_list_copy,i)->used=-2;
  }
  for(i=0;i<lengthSimpleList(global_transcripts);i++){
	 char *tr=elementAtSimpleList(global_transcripts,i)->value;//stringa del trascritto

	 for(j=0;j<lengthTranscript(tr);j++){
		int exon_number=atoi(elementAtTranscript(tr,j));
		elementAtExonList(exon_list_copy,exon_number)->used=-1;
	 }
  }
  int mapping=0;
  for(i=0;i<n;i++){
	 if(elementAtExonList(exon_list_copy,i)->used == -1){
		elementAtExonList(exon_list_copy,i)->used=mapping;
		mapping++;
	 }
  }

  fprintf(out,"%d\n",lengthSimpleList(global_transcripts));//numero dei trascritti

//22mar07
//fprintf(out,"%d\n",lengthExonList(exon_list_copy));//numero esoni
  fprintf(out,"%d\n", mapping);//numero esoni effettivamente usati

  if(n != 0)
	 fprintf(out,"%s\n",elementAtExonList(exon_list_copy,lengthExonList(exon_list_copy)-1)->end);//lunghezza genomica
  else
	 fprintf(out,"0\n");//lunghezza genomica

  for(i=0;i<n;i++){
//22mar07
	 if(elementAtExonList(exon_list_copy,i)->used >= 0){
		fprintf(out,"%s:%s",elementAtExonList(exon_list_copy,i)->start,elementAtExonList(exon_list_copy,i)->end);//start:end
#ifdef PRINT_POLYA
		fprintf(out,":%d",(elementAtExonList(exon_list_copy,i)->polya == 2)?(1):(0));//start:end
#endif
		fprintf(out,"\n");//start:end
//22mar07
	 }
  }

  for(i=0;i<lengthSimpleList(global_transcripts);i++){
	 char *tr=elementAtSimpleList(global_transcripts,i)->value;//stringa del trascritto

//19dic06
//fprintf(out,".%d\n",lengthSimpleList(elementAtListOfList(global_tr_est_ids,i)->local));
	 fprintf(out,".%d",lengthSimpleList(elementAtListOfList(global_tr_est_ids,i)->local));
#ifdef DONT_EXTEND_REFSEQ
	 char *isRefSeq=elementAtSimpleList(isRefSeqTranscriptForGlobal,i)->value;
	 if(!strcmp(isRefSeq, "1")){
		structSimpleList *s_EST_ids=elementAtListOfList(global_tr_est_ids,i)->local;
		while(s_EST_ids != NULL){
		  if(s_EST_ids->value != NULL) {
			/* UPDATE for noncoding RefSeq
		  	*/
			 /*if(s_EST_ids->value[0] == 'N' && s_EST_ids->value[1] == 'M'){
				fprintf(out,".%s", s_EST_ids->value);
			 }*/
		  	if(s_EST_ids->value[0] == 'N' && s_EST_ids->value[2] == '_'){
		  		if(s_EST_ids->value[1] == 'M' || s_EST_ids->value[1] == 'R'){
					fprintf(out,".%s", s_EST_ids->value);
		  		}
		  	}	  
		  }
		  s_EST_ids=s_EST_ids->next;
		}
	 }
#endif
	 fprintf(out,"\n");

//22mar07
//fprintf(out,"%s\n",tr);
	 for(j=0;j<lengthTranscript(tr);j++){
		int exon_number=atoi(elementAtTranscript(tr,j));
		fprintf(out,"%d", elementAtExonList(exon_list_copy,exon_number)->used);
		if(j < lengthTranscript(tr)-1)
		  fprintf(out,".");
		else
		  fprintf(out,"\n");
	 }

	 for(j=0;j<lengthTranscript(tr);j++){
		int exon_number=atoi(elementAtTranscript(tr,j));

//RAFFA 30mag04 setto reduce = 0 per gli esoni che non sono stati ridotti (che compaiono ancora)
		elementAtExonList(exon_list_copy,exon_number)->reduce=0;

		fprintf(out,"%s\n",elementAtExonList(exon_list_copy,exon_number)->sequence);//sequenza dell'esone
	 }
  }

  fprintf(out,"#\n");

//22mar07
  fprintf(out,"*\n");
 // exit(0);


/*   structSimpleList *tr_list[n]; */
/*   for(i=0;i<n;i++) */
/* 	 tr_list[i]=NULL; */

/*   for(i=0;i<lengthSimpleList(global_transcripts);i++){ */
/* 	 char *tr=elementAtSimpleList(global_transcripts,i)->value; */
/* 	 for(j=0;j<lengthTranscript(tr);j++){ */
/* 		int exon_number=atoi(elementAtTranscript(tr,j)); */
/* 		tr_list[exon_number]=pushSimpleList(tr_list[exon_number],itoa(i)); */
/* 	 } */
/*   } */
/* /\*fprintf(stderr,"HERE %s:%s\n",elementAtExonList(exon_list_copy,309)->start,elementAtExonList(exon_list_copy,309)->end);//start:end */
/*   fprintf(stderr,"HERE %s:%s\n",elementAtExonList(exon_list_copy,406)->start,elementAtExonList(exon_list_copy,406)->end);//start:end */
/*   fprintf(stderr,"HERE %d\n",getMatxPrefSuff(matx_pref_suff,309,406,n)->flag);//start:end */
/*   fprintf(stderr,"HERE %d\n",getMatxPrefSuff(matx_pref_suff,309,406,n)->actual);//start:end */
/*   exit(0);*\/ */

/* /\*or(i=0;i<lengthSimpleList(global_transcripts);i++){ */
/*   char *transcript=elementAtSimpleList(global_transcripts,i)->value; */
/*   printf("%s\n",transcript); */
/*   for(k=0;k<lengthTranscript(transcript);k++){ */
/*   char *ch=elementAtTranscript(transcript,k); */
/* //			printf("%s\n",elementAtExonList(exon_list,atoi(ch))->sequence); */
/* printf("Start %s End %s\n", elementAtExonList(exon_list,atoi(ch))->start, elementAtExonList(exon_list,atoi(ch))->end); */
/* } */
/* //		printSimpleList(elementAtListOfList(global_tr_est_ids,i)->local); */
/* } */
/* *\/ */
/* //	printf("forme di competing\n"); */
/* /\* */
/*   flag=abs(1) 3'5' */
/*   flag=abs(2) 5'3' */
/*   flag=abs(3) -> actual=1 5'3' */
/*   actual=2 3'5' */
/*   actual=3 5'3' 3'5'*\/ */
/*   int number=0; */
/* //ATTENZIONE CHE per la mod 22mar07 gli indici degli esoni che vengono stampati da qui in poi non sono piu' gli stessi di prima (TOGLIERE!!!) */
/* //calcolo per il conteggio del numero di forme alternative */
/*   for(i=0;i<n;i++){ */
/* 	 for(j=i+1;j<n;j++){ */
/* 		int actual=getMatxPrefSuff(matx_pref_suff,i,j,n)->actual; */
/* 		int flag=getMatxPrefSuff(matx_pref_suff,i,j,n)->flag; */
/* 		if (flag!=0){ */
/* 		  int i_reduce=elementAtExonList(exon_list_copy,i)->reduce; */
/* 		  int j_reduce=elementAtExonList(exon_list_copy,j)->reduce; */

/* //RAFFA 4giu04 */
/* 		  int i_external=elementAtExonList(exon_list_copy,i)->external; */
/* 		  int j_external=elementAtExonList(exon_list_copy,j)->external; */

/* /\*if(abs(flag)==1 && actual==1 && i_reduce==0 && j_reduce==0) number++; */
/*   if(abs(flag)==2 && actual==1 && i_reduce==0 && j_reduce==0) number++; */
/*   if(abs(flag)==3 && actual==1 && i_reduce==0 && j_reduce==0) number++; */
/*   if(abs(flag)==3 && actual==2 && i_reduce==0 && j_reduce==0) number++; */
/*   if(abs(flag)==3 && actual==3 && i_reduce==0 && j_reduce==0) number++; */
/*   if(flag==-4 && actual==1 && i_reduce==0 && j_reduce==0) number++; */
/*   if(flag==-4 && actual==2 && i_reduce==0 && j_reduce==0) number++; */
/*   if(flag==-4 && actual==3 && i_reduce==0 && j_reduce==0) number++;*\/ */

/* //RAFFA 4giu04 */
/* 		  if(abs(flag)==1 && actual==1 && i_reduce==0 && j_reduce==0 && i_external==0 && j_external==0) number++; */
/* 		  if(abs(flag)==2 && actual==1 && i_reduce==0 && j_reduce==0 && i_external==0 && j_external==0) number++; */
/* 		  if(abs(flag)==3 && actual==1 && i_reduce==0 && j_reduce==0 && i_external==0 && j_external==0) number++; */
/* 		  if(abs(flag)==3 && actual==2 && i_reduce==0 && j_reduce==0 && i_external==0 && j_external==0) number++; */
/* 		  if(abs(flag)==3 && actual==3 && i_reduce==0 && j_reduce==0 && i_external==0 && j_external==0) number++; */
/* 		  if(flag==-4 && actual==1 && i_reduce==0 && j_reduce==0 && i_external==0 && j_external==0) number++; */
/* 		  if(flag==-4 && actual==2 && i_reduce==0 && j_reduce==0 && i_external==0 && j_external==0) number++; */
/* 		  if(flag==-4 && actual==3 && i_reduce==0 && j_reduce==0 && i_external==0 && j_external==0) number++; */
/* 		} */
/* 	 } */
/*   } */
/*   fprintf(out,"%d\n",number); */
/*   for(i=0;i<n;i++){ */
/* 	 for(j=i+1;j<n;j++){ */
/* 		int actual=getMatxPrefSuff(matx_pref_suff,i,j,n)->actual; */
/* 		int flag=getMatxPrefSuff(matx_pref_suff,i,j,n)->flag; */
/* 		if (flag!=0){ */

/* 		  int i_reduce=elementAtExonList(exon_list_copy,i)->reduce; */
/* 		  int j_reduce=elementAtExonList(exon_list_copy,j)->reduce; */

/* //RAFFA 4giu04 */
/* 		  int i_external=elementAtExonList(exon_list_copy,i)->external; */
/* 		  int j_external=elementAtExonList(exon_list_copy,j)->external; */

/* //RAFFA 4giu04 */
/* //if(abs(flag)==1 && actual==1 && i_reduce==0 && j_reduce==0) //5'3' */
/* 		  if(abs(flag)==1 && actual==1 && i_reduce==0 && j_reduce==0 && i_external==0 && j_external==0) //5'3' */
/* 			 printAlternativesForms(out,tr_list[i],tr_list[j],i,j,-2); */

/* //RAFFA 4giu04 */
/* //if(abs(flag)==2 && actual==1 && i_reduce==0 && j_reduce==0) */
/* 		  if(abs(flag)==2 && actual==1 && i_reduce==0 && j_reduce==0 && i_external==0 && j_external==0) */
/* 			 printAlternativesForms(out,tr_list[i],tr_list[j],i,j,-1);//3'5' */

/* //RAFFA 4giu04 */
/* //if(abs(flag)==3 && actual==1 && i_reduce==0 && j_reduce==0) */
/* 		  if(abs(flag)==3 && actual==1 && i_reduce==0 && j_reduce==0 && i_external==0 && j_external==0) */
/* 			 printAlternativesForms(out,tr_list[i],tr_list[j],i,j,-2);//5'3' */

/* //RAFFA 4giu04 */
/* //if(abs(flag)==3 && actual==2 && i_reduce==0 && j_reduce==0) */
/* 		  if(abs(flag)==3 && actual==2 && i_reduce==0 && j_reduce==0 && i_external==0 && j_external==0) */
/* 			 printAlternativesForms(out,tr_list[i],tr_list[j],i,j,-1);//3'5' */

/* //RAFFA 4giu04 */
/* //if(abs(flag)==3 && actual==3 && i_reduce==0 && j_reduce==0) */
/* 		  if(abs(flag)==3 && actual==3 && i_reduce==0 && j_reduce==0 && i_external==0 && j_external==0) */
/* 			 printAlternativesForms(out,tr_list[i],tr_list[j],i,j,0);//3'5' & 5'3' */

/* //RAFFA 4giu04 */
/* //if(flag==-4 && actual==1 && i_reduce==0 && j_reduce==0) */
/* 		  if(flag==-4 && actual==1 && i_reduce==0 && j_reduce==0 && i_external==0 && j_external==0) */
/* 			 printAlternativesForms(out,tr_list[i],tr_list[j],i,j,-2);//5'3' */

/* //RAFFA 4giu04 */
/* //if(flag==-4 && actual==2 && i_reduce==0 && j_reduce==0) */
/* 		  if(flag==-4 && actual==2 && i_reduce==0 && j_reduce==0 && i_external==0 && j_external==0) */
/* 			 printAlternativesForms(out,tr_list[i],tr_list[j],i,j,-1);//3'5' */

/* //RAFFA 4giu04 */
/* //if(flag==-4 && actual==3 && i_reduce==0 && j_reduce==0) */
/* 		  if(flag==-4 && actual==3 && i_reduce==0 && j_reduce==0  && i_external==0 && j_external==0) */
/* 			 printAlternativesForms(out,tr_list[i],tr_list[j],i,j,0);//3'5' & 5'3' */

/* 		} */
/* 	 } */
/*   } */
/*   fprintf(out,"*\n"); */
/* //peril debug */
/* //crea un file con tutti i trascritti numerati */
/* /\*FILE *filetemp=fopen("temp.txt","w"); */
/*   for(i=0;i<lengthSimpleList(global_transcripts);i++) */
/*   fprintf(filetemp,"%d %s\n",i,elementAtSimpleList(global_transcripts,i)->value); */
/*   fclose(filetemp);*\/ */
}

int main(int argc,char **argv){
  if (argc<3){
	 printf("Usage: %s <inputfile> <outputfile>\n",argv[0]);
	 exit(0);
  }

  INFO("BUILD-TRANSCRIPTS");
  PRINT_LICENSE_INFORMATION;
  PRINT_SYSTEM_INFORMATION;

  pmytime pt_tot= MYTIME_create_with_name("Total");

  MYTIME_start(pt_tot);

  transcripts_finder(argv[1],argv[2]);

  MYTIME_stop(pt_tot);
  MYTIME_LOG(INFO, pt_tot);

  MYTIME_destroy(pt_tot);

  INFO("End");
  resource_usage_log();

   return 0;
}
