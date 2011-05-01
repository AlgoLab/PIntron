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
//Versione attualmente sul sito (i criteri di matching degli esoni sono quelli di Pesole)
//I trascritti vengono compattati per polyA e vengono evidenziate le diverse terminazioni
//Esiste un filtraggio per introni finale

//ATTENZIONE!!!!!: alcuni percorsi di composizione vengono persi quando viene scandita la coda

/*ATTENZIONE: non funziona se vengono anche letti trascritti di un solo esone*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "my_time.h"
#include "log.h"
#include "util.h"
#include "log-build-info.h"

//#define HALT_EXIT_MODE //If defined, exits with halt(), otherwise with exit(EXIT_FAILURE)

//31ott07 Modifica per sperimentazione trascritti
#define STRONG_FIRST_LAST_MATCH

#define MAX_EXONS 500           //Massimo numero di esoni gestibili per un trascritto

//07ott05
#define MAX_POLYA_END 20                //Massimo numero di terminazioni polyA gestibili per un trascritto

#define MAX_TRANSCRIPTS 5000    //Massimo numero di trascritti gestibili (sia per l'input r che per l'output)
#define MAX_NLD 300000          //Massimo numero di nucleotidi gestibili in una sequenza

#define EXT_EDIT 2              //Distanza di edit per la sostituzione degli esterni

#define MAX_EXON_IN_TRANS 100   //Dimensione della lista degli indici dei trascritti in cui compare un esone

//21set05
#define MAX_DIFF_FOR_REDUCING 20        //Differenza massima tra i left (right) end tra un esterno sinistro (destro)
//ed un interno per ridurre l'esterno all'interno

//#define MAX_DIFF_FOR_STRENGTH 50              //Dimensione minima di un esterno per essere agganciato

//#define MIN_DIM_FOR_STRENGTH(length) 20*length/100
#define MIN_DIM_FOR_STRENGTH(length) 20

//20set07
#define MIN_DIM_FOR_STRENGTH2(length) 20*length/100

//10mag05
#define PRINT_POLYA

//20feb06
#define MIN_POLYA_DIFF 24

//19dic06
//Se questa macro e' definita, i trascritti associati ad un refseq (NM_*) vengono mantenuti
#define DONT_EXTEND_REFSEQ

//07ott05
//#define MERGE_POLYA   //Unisce trascritti uguali con polyA diverse

//09mag08 Settare questa macro a 1 se si vogliono in input trascritti anche di un solo esone
//#define MIN_EXONS_ACCEPTED_INPUT 2
#define MIN_EXONS_ACCEPTED_INPUT 1

//09mag08
//#define MIN_EXONS_ACCEPTED_OUTPUT 2

//09mag08 Settare questa macro a 1 se si vogliono in output trascritti anche di un solo esone
//#define FIRST_MIN_EXONS_ACCEPTED_OUTPUT 2
#define FIRST_MIN_EXONS_ACCEPTED_OUTPUT 1

#define SECOND_MIN_EXONS_ACCEPTED_OUTPUT 4

#define MIN_CONFIRMED_EST_INPUT 1

//27giu05
#define PRUNE_EXON_COMP
//#define MYERS_PRUNING

//18nov05 Se si definisce questa macro fare attenzione a non aggiornare esoni che provengono da refseq!!!
//#define UPDATE_EXON

//21nov05
#define FILTER_BY_INTRONS

//21nov05
#define EXT_NOT_ON_LAST

//21nov05
//#define IGNORE_POLYA

//30gen06
#define MULTI_FASTA_FORMAT

#define READ_ABS_COORD

/*STRUTTURE DATI*/

/*Struttura identificante un trascritto come sequenza degli indici interi dei suoi esoni*/
struct transcript{

  int exons;            /*Numero di esoni che compongono il trascritto*/
  int exon_list[MAX_EXONS];     /*Lista degli indici interi degli esoni interni che compongono il trascritto in
                                                                                  numero pari a (exons-2)*/
  int left_ext;                 /*Esone esterno sinistro*/
  int right_ext;                        /*Esone esterno destro*/

//09mag08
//Se exons==1, allora rappresentare l'unico esone in left_ext e right_ext (che saranno di uguale valore)

//07ott05
#ifdef MERGE_POLYA
  int number_of_polya;
  int polya_end[MAX_EXONS];
#endif

  int ESTs;

//19dic06                                               /*1: refseq, otherwise 0*/
  char type;
  char RefSeq[20];

//03set04
//char sequence[MAX_EXONS][MAX_NLD];    /*sequenze degli esoni che compongono il trascritto*/

};

/*Struttura identificante un percorso come sequenza di indici di trascritti (nodi)*/
struct path{

  struct node *n;               /*Lista dei nodi del percorso (puntatore alla testa della lista)*/
  struct node *tail;    /*Coda della lista (indirizzo finale)*/
  int end;                      /*indice del nodo finale del percorso*/

  struct path *next;

  struct transcript tr; /*Trascritto corrispondente*/
  int L;                        /*Limite per l'estensione*/
  char visit;                   /*TRUE se il path e' da prendere in considerazione, FALSE altrimenti*/
};

/*Struttura identificante un nodo del grafo (trascritto)*/
struct node{

  int index;                    /*indice del nodo*/

  struct node *next;

};

/*Struttura che implementa una coda di percorsi*/
struct queue{
  struct path *head;

  struct path *tail;
};


/*PROTOTIPI*/

/*Costruisce e riempie la transcript_list leggendo un file passato come argomento*/
static void Get_Transcripts_from_File();

/*Costruisce la matrice di estensione (o di adiacenza). Le struttura transcript_list deve essere gia' stata riempita. I
  valori di overlapping sono rispetto al primo esone del trascritto e non rispetto al primo esone interno*/
static void Build_Extension_Matrix();

/*...*/
//08giu05
//char Extends(struct transcript t1, struct transcript t2, int *L, char only_check_name);
//07ott05
//char Extends(struct transcript t1, struct transcript t2, int *L, char only_check_name, char for_ext);
//08gen08
//char Extends(struct transcript t1, struct transcript t2, int *L, char only_check_name, char for_ext, char force_polya);
static char Extends(struct transcript t1, struct transcript t2, int *L, char only_check_name, char for_ext, char force_polya, char filt_phase);

/*Restituisce 1 se t1 estende t2, 2 se t1 e' incluso in t2 e 0 altrimenti. In L viene memorizzato l'indice (dal left ext)
  dell'elemento di t2 che matcha con il primo elemento di t1 nel caso di return not 0*/
//08giu05
//char Overlap(struct transcript t1, struct transcript t2, int *L, char only_check_name);
//07ott05
//char Overlap(struct transcript t1, struct transcript t2, int *L, char only_check_name, char for_ext);
//08gen08
//char Overlap(struct transcript t1, struct transcript t2, int *L, char only_check_name, char for_ext, char force_polya);
static char Overlap(struct transcript t1, struct transcript t2, int *L, char only_check_name, char for_ext, char force_polya, char filt_phase);

static char Check_L_suffix(int exon1, int exon2, char *matching_strength);

static char Check_R_prefix(int exon1, int exon2, char *matching_strength);

static char Check_exons(int exon1, int exon2);

/*Esegue la copia di un percorso*/
static struct path *Copy_of_Path(struct path *path);

static void Add_Path_List(struct path **path_list1, struct path **path_list2);

/*Aggiunge un percorso in testa ad una lista*/
static void Add_Path(struct path **path_list, struct path **path);

static void Empty_Source_List_of_Paths();

static char Equal_Paths(struct path *path1, struct path *path2);

/*Aggiunge un nodo ad un percorso*/
static void Add_Node(struct path **path, int node, char upd_tr);

/*Trova tutti i percorsi da tutte le sorgenti e tutti i trascritti massimali*/
static void Set_Paths();

/*Trova tutti i percorsi da una sorgente (source_index e' legato all'indice del trascritto tramite source_list)*/
static void Set_Paths_for_Source(int source_index);

/*Inizializza una coda di nodi*/
static void InitializeQueue();

/*Vero se la coda e' vuota*/
static char Queue_Is_Empty();

/*Attacca un nodo in coda*/
static void Enqueue(struct path **path);

/*Estrae un elemento dalla coda*/
static struct path *Dequeue();

/*Restituisce l'indirizzo del path in coda che ha lo stesso trascritto di arg_path*/
static struct path *Get_path_with_the_same_exons(struct path *arg_path);

/*Restituisce TRUE se ai due path e' associato lo stesso trascritto*/
static char Equals_transcripts(struct transcript t1, struct transcript t2);

/*Dati due trascritti t1 e t2 ed un limite L di estensione (posizione dell'esone di t1 da cui inizia l'estensione),
  costruisce e restituisce il trascritto complessivo e lo memorizza nella memoria puntata da extension (che deve essere
  gia' stata allocata; la funzione alloca solo exon_list). Non effettua controlli sulla veridicita' dell'argomento L.
  La liberazione della memoria di extension deve essere gestita completamente al di fuori*/
static void Build_extension(struct transcript t1, struct transcript t2, int L, struct transcript *extension);

/*Copia il trascritto t complessivo e lo memorizza nella memoria puntata da copy (che deve essere
  gia' stata allocata; la funzione alloca solo exon_list).
  La liberazione della memoria di copy deve essere gestita completamente al di fuori*/
static void Copy_transcript(struct transcript t, struct transcript *copy);

/*Determina il trascritto massimale a partire da un percorso sul grafo dei trascritti e lo memorizza in return_transcript
  (gia' allocato; la funzione alloca solo exon_list)*/
/* static void Get_Transcript_From_Path(struct path *tr_path, struct transcript *return_transcript); */

static void Set_Path_Transcripts_for_Source(struct path *path);

static void Set_Path_Transcripts();

/*Cancella tutti i trascritti che sono inclusi in altri*/
static void Filter_Path_Transcripts();

//13ott05
static void Filter_Path_Transcripts_by_Introns();

//static void Filter_Path_Transcripts_by_Introns_OLD();

/*Restituisce un percorso con l'unico nodo specificato dal parametro index*/
static struct path *Create_Source_Path(int index);

static char *Substring(char *string, int left, int right);

/* static int Edit_distance(char *seq1, char *seq2); */

//01set04
/* static void Reduce_external_exons(); */

/* static void Update_exon(int exon1, int exon2); */

static void Graph_reduction();

//27giu05
/* static void Graph_reduction_Myers(); */

static int Get_Opposite_Node_Index(int a_index, int b_index, int initial_index);

static void Partial_Graph_reduction_for_arc(int a_index, int b_index);

//Dati a,b,c definiti nella procedura Get_Opposite_Node_Index, effettuare la riduzione del grafo
static void Partial_Graph_reduction_for_node(int a_index, int b_index, int c_index);

static void Add_Node_to_a_node_list(struct node **node_list, int node);

static void Remove_Node_from_a_node_list(struct node **node_list, int node);

//19gen05
/*Elimina i trascritti che praticamente coincidono a meno di una differenza tra gli esterni*/
static void First_Filtering();

#ifdef READ_ABS_COORD
static int GetAbsoluteStart(int left, int right);
static int GetAbsoluteEnd(int left, int right);
#endif

//21feb06
static char *To_lower(char *str);

/*VARIABILI GLOBALI*/

//27giu05
int count_del;

/*Matrice quadrata di dimensioni pari al numero dei trascritti. La posizione (i,j) fornisce l'estensione "e" del trascritto
  j-esimo rispetto al trascritto i-esimo: e=0 se j non estende i o esiste una relazione di inclusione,
  e>0 e' la posizione dell'esone (che non puo' essere il primo) di i da cui incomincia l'estensione di j.
  Tale matrice e' in sostanza la matrice di adiacenza del grafo dei trascritti. Puo' essere inizializzata a -1. Se trovo che j
  estende i allora settero' un valore positivo in (i,j) e uno zero in (j,i): c'e' un arco che va da i a j*/
int **extension_matrix;

//27giu05
char **remove_extension_matrix;

/*Vettore: in posizione i e' TRUE se il trascritto i (nodo) e' una sorgente,
  ovvero un nodo senza archi entranti*/
char *is_source;

/*Vettore: in posizione i e' TRUE se il trascritto i (nodo) ha fornito un path*/
char *give_path;

/*Vettore di dimensioni pari al numero dei trascritti: in posizione i e' memorizzata l'informazione relativa al trascritto i*/
//struct transcript *transcript_list;
struct transcript transcript_list[MAX_TRANSCRIPTS];

/*Numero dei trascritti presenti in transcript_list*/
int number_of_transcripts;

/*Vettore di dimensioni pari al numero delle sorgenti che fornisce la lista degli indici dei trascritti che sono
  sorgenti*/
int *source_list;

int *in_degree;
int *out_degree;

/*Numero dei trascritti delle sorgenti*/
int number_of_sources;

/*Coda di nodi*/
struct queue q;

/*Puntatore per le aggiunte di un nodo ad un percorso*/
struct node *add=NULL;

/*Numero totale dei percorsi trovati*/
int total_paths;

/*Numero totale dei percorsi trovati per una sorgente*/
int source_total_paths;

/*Vettore che memorizza i trascritti determinati sui percorsi trovati*/
//struct transcript *path_transcript_list;
struct transcript path_transcript_list[MAX_TRANSCRIPTS];
struct path *transcript_list_of_paths[MAX_TRANSCRIPTS];

/*Vettore che memorizza i trascritti determinati sui percorsi trovati per una sorgente*/
//struct transcript *source_path_transcript_list;
struct transcript source_path_transcript_list[MAX_TRANSCRIPTS];
struct path *source_list_of_paths[MAX_TRANSCRIPTS];

/*In posizione i e' TRUE se il trascritto i in path_transcript_list e' stato filtrato*/
char *filtered;

/*Liste degli estremi sinistri e destri degli esoni*/
int *list_of_exon_left;
int *list_of_exon_right;
int *list_of_old_exon_left;
int *list_of_old_exon_right;

//10mag05
char *polya;

//03set04
char **sequences;

//01set04
/*0 if the exon is internal in at least one transcript, 1 if sx ext and 2 if dx ext*/
char *is_internal;
//char *is_external;

/*Lista dei flag di comparsa degli esoni (inizialmente settati a 0)*/
char **exists_list;

//22mar07
//char init_reading[100000];
//char init_reading2[100000];
char init_reading[200000];
char init_reading2[200000];

int number_of_exons;    /*Numero degli esoni nella lista totale*/

//22mar07
#ifndef MULTI_FASTA_FORMAT
int alternative_number;
int *vect_alternative_number;
char ***alternative_forms;
#endif

//22mar07
#ifndef MULTI_FASTA_FORMAT
int ***exon_in_transcripts;
int **exon_in_trans_index;
#endif

/*Numero totale dei tascritti massimali dopo filtraggio*/
int actual_path_number=0;

int *vect_actual_path_number;

struct path *path_to_be_added;

int transcript_counter;

#ifdef READ_ABS_COORD
int gen_start, gen_end;
int strand;
int boundary;
#endif

int
main(){
  INFO("BUILD-FULL-LENGTHS");
  PRINT_LICENSE_INFORMATION;
  PRINT_SYSTEM_INFORMATION;
  int i=0, j=0, k=0, p=0;
  char stop=0;
  //char temp_string[MAX_NLD];
  char temp_string2[50];
  //char temp_string3[50];
  //char temp_exon[7];
  //int counter=0;

//22mar07
#ifndef MULTI_FASTA_FORMAT
  int *alt_counter=NULL;
#endif

  char *stampa=NULL;
  int *trans_order=NULL;
  //int first_exon=-1, second_exon=-1;
//FILE *out=NULL;
  FILE **tr_out=NULL;
  FILE **comp_out=NULL;
  struct path *comp_paths=NULL;
  struct node *nodes=NULL;
  char file_name[50];
  char suff_string[5];

  pmytime pt_tot= MYTIME_create_with_name("Total");

  MYTIME_start(pt_tot);

  vect_actual_path_number=(int *)malloc((SECOND_MIN_EXONS_ACCEPTED_OUTPUT-FIRST_MIN_EXONS_ACCEPTED_OUTPUT+1)*sizeof(int));
  if(vect_actual_path_number == NULL){
         fprintf(stderr, "Problem1 in allocating vect_actual_path_number!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  trans_order=(int *)malloc((SECOND_MIN_EXONS_ACCEPTED_OUTPUT-FIRST_MIN_EXONS_ACCEPTED_OUTPUT+1)*sizeof(int));
  if(trans_order == NULL){
         fprintf(stderr, "Problem in allocating trans_order!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  stampa=(char *)malloc((SECOND_MIN_EXONS_ACCEPTED_OUTPUT-FIRST_MIN_EXONS_ACCEPTED_OUTPUT+1)*sizeof(char));
  if(stampa == NULL){
         fprintf(stderr, "Problem in allocating stampa!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  for(i=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; i <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; i++){
         vect_actual_path_number[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=0;
         trans_order[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=0;
         stampa[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=0;
  }

  Get_Transcripts_from_File();

/*for(i=0; i<number_of_transcripts; i++){
  fprintf(stdout, "Transcripts %d (conf %d)*******************\n", i, transcript_list[i].ESTs);
  fprintf(stdout, "TR %d-%d(%d)", list_of_exon_left[transcript_list[i].left_ext], list_of_exon_right[transcript_list[i].left_ext], transcript_list[i].left_ext);
  for(k=0; k<transcript_list[i].exons-2; k++)
  fprintf(stdout, " %d-%d(%d)", list_of_exon_left[transcript_list[i].exon_list[k]], list_of_exon_right[transcript_list[i].exon_list[k]], transcript_list[i].exon_list[k]);
  fprintf(stdout, " %d-%d(%d)", list_of_exon_left[transcript_list[i].right_ext], list_of_exon_right[transcript_list[i].right_ext], transcript_list[i].right_ext);

  fprintf(stdout, "Type %d (polyA %d) *******************\n", transcript_list[i].type, polya[transcript_list[i].right_ext]);
  }*/
//exit(0);

/*for(i=0; i < number_of_exons; i++){
  fprintf(stdout, "%d) %d:%d:%d is int %d\n", i, list_of_exon_left[i], list_of_exon_right[i], polya[i], is_internal[i]);
  }*/
//exit(0);

/*fprintf(stderr, "EXT 33 %d\n", is_internal[33]);
  fprintf(stderr, "EXT 37 %d\n", is_internal[37]);
  fprintf(stderr, "EXT 115 %d\n", is_internal[115]);
  fprintf(stderr, "EXT 137 %d\n", is_internal[137]);
  fprintf(stderr, "HERE 206 %d %d\n", list_of_exon_left[206], list_of_exon_right[206]);*/

//01set04
/*Procedura di aumento (diminuzione) del left (right) end di un esone esterno sinistro (destro) per portarlo
  a coincidere con un taglio confermato se la differenza rientra in certi limiti*/
//19dic06
//Reduce_external_exons();
#ifndef DONT_EXTEND_REFSEQ
  Reduce_external_exons();
#endif

/*for(i=0; i<number_of_transcripts; i++){
  fprintf(stdout, "Transcripts %d (conf %d)*******************\n", i, transcript_list[i].ESTs);
  fprintf(stdout, "TR %d-%d(%d)", list_of_exon_left[transcript_list[i].left_ext], list_of_exon_right[transcript_list[i].left_ext], transcript_list[i].left_ext);
  fprintf(stdout, "Seq %s\n", sequences[transcript_list[i].left_ext]);
  for(k=0; k<transcript_list[i].exons-2; k++){
  fprintf(stdout, " %d-%d(%d)", list_of_exon_left[transcript_list[i].exon_list[k]], list_of_exon_right[transcript_list[i].exon_list[k]], transcript_list[i].exon_list[k]);
  fprintf(stdout, "Seq %s\n", sequences[transcript_list[i].exon_list[k]]);
  }
  fprintf(stdout, " %d-%d(%d)", list_of_exon_left[transcript_list[i].right_ext], list_of_exon_right[transcript_list[i].right_ext], transcript_list[i].right_ext);
  fprintf(stdout, "Seq %s\n", sequences[transcript_list[i].right_ext]);
  fprintf(stdout, "Type %d (polyA %d) *******************\n", transcript_list[i].type, polya[transcript_list[i].right_ext]);
  }*/
//exit(0);

/*fprintf(stderr, "HERE 14 %d %d\n", list_of_exon_left[14], list_of_exon_right[14]);
  fprintf(stderr, "HERE 19 %d %d\n", list_of_exon_left[19], list_of_exon_right[19]);
  fprintf(stderr, "HERE 102 %d %d\n", list_of_exon_left[102], list_of_exon_right[102]);
  fprintf(stderr, "HERE 41 %d %d\n", list_of_exon_left[41], list_of_exon_right[41]);
  exit(0);*/


/*fprintf(stderr, "HERE 922 %d %d\n", list_of_exon_left[922], list_of_exon_right[922]);
  fprintf(stderr, "HERE 925 %d %d\n", list_of_exon_left[926], list_of_exon_right[926]);
  fprintf(stderr, "HERE 38 %d %d\n", list_of_exon_left[38], list_of_exon_right[38]);
  fprintf(stderr, "HERE 65 %d %d\n", list_of_exon_left[65], list_of_exon_right[65]);
  fprintf(stderr, "HERE 120 %d %d\n", list_of_exon_left[120], list_of_exon_right[120]);
  fprintf(stderr, "HERE 132 %d %d\n", list_of_exon_left[132], list_of_exon_right[132]);
  fprintf(stderr, "HERE 145 %d %d\n", list_of_exon_left[145], list_of_exon_right[145]);*/

/*for(i=0; i<number_of_transcripts; i++){
  fprintf(stdout, "Transcripts %d (conf %d)*******************\n", i, transcript_list[i].ESTs);
  fprintf(stdout, "TR %d-%d(%d)", list_of_exon_left[transcript_list[i].left_ext], list_of_exon_right[transcript_list[i].left_ext], transcript_list[i].left_ext);
  fprintf(stdout, "Seq %s\n", sequences[transcript_list[i].left_ext]);
  for(k=0; k<transcript_list[i].exons-2; k++){
  fprintf(stdout, " %d-%d(%d)", list_of_exon_left[transcript_list[i].exon_list[k]], list_of_exon_right[transcript_list[i].exon_list[k]], transcript_list[i].exon_list[k]);
  fprintf(stdout, "Seq %s\n", sequences[transcript_list[i].exon_list[k]]);
  }
  fprintf(stdout, " %d-%d(%d)", list_of_exon_left[transcript_list[i].right_ext], list_of_exon_right[transcript_list[i].right_ext], transcript_list[i].right_ext);
  fprintf(stdout, "Seq %s\n", sequences[transcript_list[i].right_ext]);
  fprintf(stdout, "Type %d (polyA %d) *******************\n", transcript_list[i].type, polya[transcript_list[i].right_ext]);
  }
  exit(0);*/

//19gen05
  First_Filtering();

  Build_Extension_Matrix();

 //27giu05
#ifdef MYERS_PRUNING
  Graph_reduction_Myers();
#else
  Graph_reduction();
#endif

//Conteggio delle sorgenti
  number_of_sources=0;
  for(i=0; i<number_of_transcripts; i++){
         if(in_degree[i] == 0){
                number_of_sources+=1;
                is_source[i]=1;
         }
         else
                is_source[i]=0;
  }

  source_list=(int *)malloc(number_of_sources*sizeof(int));
  if(source_list== NULL){
         fprintf(stderr, "Problem4 of memory allocation in Build_Extension_Matrix!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  j=0;
  for(i=0; i<number_of_transcripts; i++){
         if(is_source[i]){
                source_list[j]=i;
                j++;
         }
  }

/*fprintf(stderr, "HERE 20 %d %d\n", list_of_exon_left[20], list_of_exon_right[20]);
  fprintf(stderr, "HERE 17 %d %d\n", list_of_exon_left[17], list_of_exon_right[17]);
  fprintf(stderr, "HERE 38 %d %d\n", list_of_exon_left[38], list_of_exon_right[38]);
  fprintf(stderr, "HERE 65 %d %d\n", list_of_exon_left[65], list_of_exon_right[65]);
  fprintf(stderr, "HERE 120 %d %d\n", list_of_exon_left[120], list_of_exon_right[120]);
  fprintf(stderr, "HERE 132 %d %d\n", list_of_exon_left[132], list_of_exon_right[132]);
  fprintf(stderr, "HERE 145 %d %d\n", list_of_exon_left[145], list_of_exon_right[145]);

  fprintf(stderr, "HERE 137 %d %d\n", list_of_exon_left[137], list_of_exon_right[137]);*/

/******************************/

/*      for(i=0; i<number_of_transcripts; i++){
        fprintf(stdout, "Transcripts %d (conf %d)*******************\n", i, transcript_list[i].ESTs);
        fprintf(stdout, "TR %d-%d(%d)", list_of_exon_left[transcript_list[i].left_ext], list_of_exon_right[transcript_list[i].left_ext], transcript_list[i].left_ext);
        for(k=0; k<transcript_list[i].exons-2; k++)
        fprintf(stdout, " %d-%d(%d)", list_of_exon_left[transcript_list[i].exon_list[k]], list_of_exon_right[transcript_list[i].exon_list[k]], transcript_list[i].exon_list[k]);
        fprintf(stdout, " %d-%d(%d)", list_of_exon_left[transcript_list[i].right_ext], list_of_exon_right[transcript_list[i].right_ext], transcript_list[i].right_ext);
        fprintf(stdout, "*******************\n");
        }*/
//exit(0);

/*for(i=0; i < number_of_exons; i++){
  fprintf(stdout, "%d) %d:%d:%d is int %d\n", i, list_of_exon_left[i], list_of_exon_right[i], polya[i], is_internal[i]);
  }
  exit(0);*/


/*      for(i=0; i<number_of_transcripts; i++){
        for(j=0; j<number_of_transcripts; j++){
        if(extension_matrix[i][j] != 0){
        fprintf(stdout, "Matrix %d (%d,%d)*******************\n", extension_matrix[i][j], i, j);
        fprintf(stdout, "TR %d", transcript_list[i].left_ext);
        for(k=0; k<transcript_list[i].exons-2; k++)
        fprintf(stdout, " %d", transcript_list[i].exon_list[k]);
        fprintf(stdout, " %d\n", transcript_list[i].right_ext);
        fprintf(stdout, "TR %d", transcript_list[j].left_ext);
        for(k=0; k<transcript_list[j].exons-2; k++)
        fprintf(stdout, " %d", transcript_list[j].exon_list[k]);
        fprintf(stdout, " %d\n", transcript_list[j].right_ext);
        fprintf(stdout, "*******************\n");
        }
        }
        }
        exit(0);*/
/******************************/

  Set_Paths();

//27giu05
#ifndef MYERS_PRUNING
  k=0;
  while(k < number_of_transcripts){
         if(!give_path[k])
                fprintf(stderr, "Transcript %d not in path!\n", k);
         k++;
  }
#endif

/*      i=0;
        while(i<total_paths){
        if(i==41){
        fprintf(stdout, "TR%d type %d\n %d(%d-%d)\n", i, path_transcript_list[i].type, path_transcript_list[i].left_ext, list_of_exon_left[path_transcript_list[i].left_ext], list_of_exon_right[path_transcript_list[i].left_ext]);
        for(k=0; k<path_transcript_list[i].exons-2; k++)
        fprintf(stdout, "       %d(%d-%d)\n", path_transcript_list[i].exon_list[k], list_of_exon_left[path_transcript_list[i].exon_list[k]], list_of_exon_right[path_transcript_list[i].exon_list[k]]);
        fprintf(stdout, "       %d(%d-%d)\n\n", path_transcript_list[i].right_ext, list_of_exon_left[path_transcript_list[i].right_ext], list_of_exon_right[path_transcript_list[i].right_ext]);
//fprintf(stdout, "*******************\n");
}
i++;
}
exit(0);*/

  Filter_Path_Transcripts();

/*      i=0;
        while(i<total_paths){
        fprintf(stdout, "TR(%d filt %d) %d", i, filtered[i], path_transcript_list[i].left_ext);
        for(k=0; k<path_transcript_list[i].exons-2; k++)
        fprintf(stdout, " %d", path_transcript_list[i].exon_list[k]);
        fprintf(stdout, " %d\n", path_transcript_list[i].right_ext);
        fprintf(stdout, "*******************\n");

        i++;
        }*/


//13ott05
//21nov05
#ifdef FILTER_BY_INTRONS
  Filter_Path_Transcripts_by_Introns();
#endif

/*      i=0;
        while(i<total_paths){
        fprintf(stdout, "TR(%d filt %d) %d", i, filtered[i], path_transcript_list[i].left_ext);
        for(k=0; k<path_transcript_list[i].exons-2; k++)
        fprintf(stdout, " %d", path_transcript_list[i].exon_list[k]);
        fprintf(stdout, " %d\n", path_transcript_list[i].right_ext);
        fprintf(stdout, "*******************\n");

        i++;
        }
        exit(0);*/

//fprintf(stderr, "Final paths %d\n", actual_path_number);

/*      i=0;
        while(i<total_paths){
        if(!filtered[i] && i==41){
        fprintf(stdout, "TR%d type %d\n %d(%d-%d)\n", i, path_transcript_list[i].type, path_transcript_list[i].left_ext, list_of_exon_left[path_transcript_list[i].left_ext], list_of_exon_right[path_transcript_list[i].left_ext]);
        for(k=0; k<path_transcript_list[i].exons-2; k++)
        fprintf(stdout, "       %d(%d-%d)\n", path_transcript_list[i].exon_list[k], list_of_exon_left[path_transcript_list[i].exon_list[k]], list_of_exon_right[path_transcript_list[i].exon_list[k]]);
        fprintf(stdout, "       %d(%d-%d)\n\n", path_transcript_list[i].right_ext, list_of_exon_left[path_transcript_list[i].right_ext], list_of_exon_right[path_transcript_list[i].right_ext]);
//fprintf(stdout, "*******************\n");
}
i++;
}
exit(0);*/

/*for(i=0; i < number_of_exons; i++){
  fprintf(stdout, "%d) %d:%d L %d seq length %d\n", i, list_of_exon_left[i], list_of_exon_right[i], (list_of_exon_right[i]-list_of_exon_left[i]+1), strlen(sequences[i]));
  }
  exit(0);*/

//01set04
  for(i=0; i<number_of_exons; i++){
         sprintf(temp_string2, "%d:%d", list_of_exon_left[i], list_of_exon_right[i]);
#ifndef MULTI_FASTA_FORMAT
         strcat(init_reading, temp_string2);
#endif
         strcat(init_reading2, temp_string2);

//10mag05
#ifdef PRINT_POLYA
         sprintf(temp_string2, ":%d", polya[i]);
#ifndef MULTI_FASTA_FORMAT
         strcat(init_reading, temp_string2);
#endif
#endif

         sprintf(temp_string2, ";%d:%d", list_of_old_exon_left[i], list_of_old_exon_right[i]);
         strcat(init_reading2, temp_string2);

//10mag05
#ifdef PRINT_POLYA
         sprintf(temp_string2, ":%d", polya[i]);
         strcat(init_reading2, temp_string2);
#endif

#ifndef MULTI_FASTA_FORMAT
         strcat(init_reading, "\n");
#endif
         strcat(init_reading2, "\n");
  }

/*strcpy(file_name, "F_TEMP_COMPOSITION_TRANS");
  sprintf(suff_string, "%d_%d", MIN_CONFIRMED_EST_INPUT, MIN_EXONS_ACCEPTED_OUTPUT);
  strcat(file_name, suff_string);
  strcat(file_name, ".txt");
//out=fopen("TEMP_COMPOSITION_TRANS.txt", "w");
out=fopen(file_name, "w");
if(out == NULL){
fprintf(stderr, "Problem in opening file!\n");
#ifdef HALT_EXIT_MODE
exit(1);
#else
exit(EXIT_FAILURE);
#endif
}*/

  tr_out=(FILE **)malloc((SECOND_MIN_EXONS_ACCEPTED_OUTPUT-FIRST_MIN_EXONS_ACCEPTED_OUTPUT+1)*sizeof(FILE *));
  comp_out=(FILE **)malloc((SECOND_MIN_EXONS_ACCEPTED_OUTPUT-FIRST_MIN_EXONS_ACCEPTED_OUTPUT+1)*sizeof(FILE *));
  if(tr_out == NULL || comp_out == NULL){
         fprintf(stderr, "Problem1 in opening file!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  for(i=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; i <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; i++){
         tr_out[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=NULL;
         comp_out[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=NULL;
  }

  for(i=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; i <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; i++){
         strcpy(file_name, "TRANSCRIPTS");
         sprintf(suff_string, "%d_%d", MIN_CONFIRMED_EST_INPUT, i);
         strcat(file_name, suff_string);
         strcat(file_name, ".txt");
         tr_out[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=fopen(file_name, "w");
         if(tr_out[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT] == NULL){
                fprintf(stderr, "Problem2 in opening file!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }
  }

  for(i=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; i <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; i++){
         strcpy(file_name, "TEMP_COMPOSITION_TRANS");
         sprintf(suff_string, "%d_%d", MIN_CONFIRMED_EST_INPUT, i);
         strcat(file_name, suff_string);
         strcat(file_name, ".txt");
         comp_out[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=fopen(file_name, "w");
         if(comp_out[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT] == NULL){
                fprintf(stderr, "Problem3 in opening file!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }
  }

//14feb06
  for(i=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; i <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; i++)
         vect_actual_path_number[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=0;
  for(i=0; i<total_paths; i++){
         for(p=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; p <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; p++){
                if(!filtered[i] && path_transcript_list[i].exons >= p){
                  vect_actual_path_number[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]++;
                }
         }
  }

  for(i=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; i <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; i++){
         fprintf(tr_out[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%d\n", vect_actual_path_number[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]);
         fprintf(tr_out[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%s", init_reading);
         fprintf(comp_out[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%d\n", vect_actual_path_number[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]);
         fprintf(comp_out[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%s", init_reading2);
  }

/*printf("%d\n", actual_path_number);
  printf("%s", init_reading);*/
/*fprintf(out, "%d\n", actual_path_number);
  fprintf(out, "%s", init_reading2);*/

//Scrittura dei trascritti massimali e settaggio della struttura exon_in_transcripts, exon_in_trans_index
  for(i=0; i<total_paths; i++){
         for(p=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; p <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; p++){
                if(!filtered[i] && path_transcript_list[i].exons >= p){
                  trans_order[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]++;

#ifdef MULTI_FASTA_FORMAT
//19dic06
//fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], ">%d:%d\n", trans_order[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], path_transcript_list[i].exons);
                  fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], ">%d:%d", trans_order[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], path_transcript_list[i].exons);
#ifdef DONT_EXTEND_REFSEQ
                  if(path_transcript_list[i].type == 1)
                         fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], ":%s", path_transcript_list[i].RefSeq);
#endif
                  fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "\n");
#else
                  fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], ".15\n");
#endif

//if(p == FIRST_MIN_EXONS_ACCEPTED_OUTPUT)
//      fprintf(stdout, ".15(%d)\n", i+1);

                  fprintf(comp_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "..\n");

#ifdef MULTI_FASTA_FORMAT
#ifdef READ_ABS_COORD
                  fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%d:%d:", GetAbsoluteStart(list_of_exon_left[path_transcript_list[i].left_ext], list_of_exon_right[path_transcript_list[i].left_ext]), GetAbsoluteEnd(list_of_exon_left[path_transcript_list[i].left_ext], list_of_exon_right[path_transcript_list[i].left_ext]));
#endif
                  fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%d:%d:%d\n", list_of_exon_left[path_transcript_list[i].left_ext], list_of_exon_right[path_transcript_list[i].left_ext], polya[path_transcript_list[i].left_ext]);
                  fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%s\n", sequences[path_transcript_list[i].left_ext]);
#else
                  fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%d", path_transcript_list[i].left_ext);
#endif

//if(p == FIRST_MIN_EXONS_ACCEPTED_OUTPUT)
//      fprintf(stdout, "%d-%d", list_of_exon_left[path_transcript_list[i].left_ext], list_of_exon_right[path_transcript_list[i].left_ext]);

                  fprintf(comp_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%d", path_transcript_list[i].left_ext);

//22mar07
#ifndef MULTI_FASTA_FORMAT
                  exon_in_transcripts[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][path_transcript_list[i].left_ext][exon_in_trans_index[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][path_transcript_list[i].left_ext]]=trans_order[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT];
                  exon_in_trans_index[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][path_transcript_list[i].left_ext]+=1;
#endif

                  exists_list[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][path_transcript_list[i].left_ext]=1;
                  for(j=0; j<path_transcript_list[i].exons-2; j++){
#ifdef MULTI_FASTA_FORMAT
#ifdef READ_ABS_COORD
                         fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%d:%d:", GetAbsoluteStart(list_of_exon_left[path_transcript_list[i].exon_list[j]], list_of_exon_right[path_transcript_list[i].exon_list[j]]), GetAbsoluteEnd(list_of_exon_left[path_transcript_list[i].exon_list[j]], list_of_exon_right[path_transcript_list[i].exon_list[j]]));
#endif
                         fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%d:%d:%d\n", list_of_exon_left[path_transcript_list[i].exon_list[j]], list_of_exon_right[path_transcript_list[i].exon_list[j]], polya[path_transcript_list[i].exon_list[j]]);
                         fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%s\n", sequences[path_transcript_list[i].exon_list[j]]);
#else
                         fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], ".%d", path_transcript_list[i].exon_list[j]);
#endif

//if(p == FIRST_MIN_EXONS_ACCEPTED_OUTPUT)
//      fprintf(stdout, ";%d-%d", list_of_exon_left[path_transcript_list[i].exon_list[j]], list_of_exon_right[path_transcript_list[i].exon_list[j]]);

                         fprintf(comp_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], ".%d", path_transcript_list[i].exon_list[j]);

//22mar07
#ifndef MULTI_FASTA_FORMAT
                         exon_in_transcripts[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][path_transcript_list[i].exon_list[j]][exon_in_trans_index[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][path_transcript_list[i].exon_list[j]]]=trans_order[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT];
                         exon_in_trans_index[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][path_transcript_list[i].exon_list[j]]+=1;
#endif

                         exists_list[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][path_transcript_list[i].exon_list[j]]=1;
                  }

                  if(path_transcript_list[i].exons >= 2){
#ifdef MULTI_FASTA_FORMAT
#ifdef READ_ABS_COORD
                         fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%d:%d:", GetAbsoluteStart(list_of_exon_left[path_transcript_list[i].right_ext], list_of_exon_right[path_transcript_list[i].right_ext]), GetAbsoluteEnd(list_of_exon_left[path_transcript_list[i].right_ext], list_of_exon_right[path_transcript_list[i].right_ext]));
#endif
                         fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%d:%d:%d\n", list_of_exon_left[path_transcript_list[i].right_ext], list_of_exon_right[path_transcript_list[i].right_ext], polya[path_transcript_list[i].right_ext]);
                         fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%s\n", sequences[path_transcript_list[i].right_ext]);
#else
                         fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], ".%d", path_transcript_list[i].right_ext);
#endif

//if(p == FIRST_MIN_EXONS_ACCEPTED_OUTPUT)
//      fprintf(stdout, ";%d-%d", list_of_exon_left[path_transcript_list[i].right_ext], list_of_exon_right[path_transcript_list[i].right_ext]);

                         fprintf(comp_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], ".%d", path_transcript_list[i].right_ext);

                         exists_list[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][path_transcript_list[i].right_ext]=1;

//22mar07
#ifndef MULTI_FASTA_FORMAT
                         exon_in_transcripts[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][path_transcript_list[i].right_ext][exon_in_trans_index[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][path_transcript_list[i].right_ext]]=trans_order[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT];
                         exon_in_trans_index[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][path_transcript_list[i].right_ext]+=1;
#endif
                  }

//07ott05
#ifdef MERGE_POLYA
#ifndef MULTI_FASTA_FORMAT
                  for(j=0; j<path_transcript_list[i].number_of_polya; j++){
                         fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], ";%d", path_transcript_list[i].polya_end[j]);
                  }
                  if(polya[path_transcript_list[i].right_ext])
                         fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], ";%d", path_transcript_list[i].right_ext);
#endif
#endif
//fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], ";");

#ifndef MULTI_FASTA_FORMAT
                  fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "\n");
#endif
//if(p == FIRST_MIN_EXONS_ACCEPTED_OUTPUT)
//      fprintf(stdout, "\n");

                  fprintf(comp_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "\n");

#ifndef MULTI_FASTA_FORMAT
                  fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%s\n", sequences[path_transcript_list[i].left_ext]);
#endif

                  fprintf(comp_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%s", sequences[path_transcript_list[i].left_ext]);
                  for(j=0; j<path_transcript_list[i].exons-2; j++){
#ifndef MULTI_FASTA_FORMAT
                         fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%s\n", sequences[path_transcript_list[i].exon_list[j]]);
#endif
                         fprintf(comp_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%s", sequences[path_transcript_list[i].exon_list[j]]);
                  }

                  if(path_transcript_list[i].exons >= 2){
#ifndef MULTI_FASTA_FORMAT
                         fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%s\n", sequences[path_transcript_list[i].right_ext]);
#endif
                         fprintf(comp_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%s", sequences[path_transcript_list[i].right_ext]);
                  }

                  fprintf(comp_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "\n");

//Stampare composizione

                  comp_paths=transcript_list_of_paths[i];
                  while(comp_paths != NULL){
                         nodes=comp_paths->n;
                         while(nodes != NULL){
                                fprintf(comp_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], ".%d\n", transcript_list[nodes->index].ESTs);
                                fprintf(comp_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%d", transcript_list[nodes->index].left_ext);
                                for(j=0; j<transcript_list[nodes->index].exons-2; j++){
                                  fprintf(comp_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], ".%d", transcript_list[nodes->index].exon_list[j]);
                                }
                                if(transcript_list[nodes->index].exons >= 2){
                                  fprintf(comp_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], ".%d", transcript_list[nodes->index].right_ext);
                                }
                                fprintf(comp_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "\n");
                                nodes=nodes->next;
                         }
                         fprintf(comp_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "*\n");
                         comp_paths=comp_paths->next;
                  }
                }
         }
  }

/*DA QUI*/
/*fprintf(stderr, "Exon 12 in %d transcripts\n", exon_in_trans_index[12]);
  fprintf(stderr, "Exon 19 in %d transcripts\n", exon_in_trans_index[19]);
  fprintf(stderr, "Exon 52 in %d transcripts\n", exon_in_trans_index[52]);
  fprintf(stderr, "Exon 53 in %d transcripts\n", exon_in_trans_index[53]);
  for(i=0; i<exon_in_trans_index[53]; i++){
  fprintf(stderr, "     in %d\n", exon_in_transcripts[53][i]);
  }*/

  stop=0;
  for(i=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; i <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; i++){
         fprintf(tr_out[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "#\n");
         fprintf(comp_out[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "#\n");
  }

  for(i=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; i <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; i++){
         fclose(comp_out[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]);
  }

//fclose(out);

//09mag08 Da qui in poi non viene piu' utilizzato nel nuovo formato multifasta
//22mar07
#ifndef MULTI_FASTA_FORMAT
//Lettura del numero di forme alternative
  scanf("%s\n", temp_string);
  alternative_number=atoi(temp_string);

  vect_alternative_number=(int *)calloc((SECOND_MIN_EXONS_ACCEPTED_OUTPUT-FIRST_MIN_EXONS_ACCEPTED_OUTPUT+1),sizeof(int));
  if(vect_alternative_number == NULL){
         fprintf(stderr, "Problem in MAIN!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }
  for(p=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; p <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; p++)
         vect_alternative_number[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=alternative_number;

  alt_counter=(int *)calloc((SECOND_MIN_EXONS_ACCEPTED_OUTPUT-FIRST_MIN_EXONS_ACCEPTED_OUTPUT+1),sizeof(int));
  if(alt_counter == NULL){
         fprintf(stderr, "Problem in MAIN!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }
  for(p=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; p <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; p++)
         alt_counter[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=0;

  alternative_forms=(char ***)calloc((SECOND_MIN_EXONS_ACCEPTED_OUTPUT-FIRST_MIN_EXONS_ACCEPTED_OUTPUT+1),sizeof(char **));
  if(alternative_forms == NULL){
         fprintf(stderr, "Problem in MAIN!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  for(p=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; p <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; p++){
         alternative_forms[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=(char **)calloc((alternative_number*3),sizeof(char *));
         for(i=0; i<alternative_number*3; i++){
//05dic05
                alternative_forms[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][i]=(char *)malloc(2000*sizeof(char));
                if(alternative_forms[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][i] == NULL){
                  fprintf(stderr, "Problem in MAIN!\n");
#ifdef HALT_EXIT_MODE
                  exit(1);
#else
                  exit(EXIT_FAILURE);
#endif
                }
         }
  }

  while(!stop){
         scanf("%s\n", temp_string);
         if(!strcmp(temp_string, "*")){
                stop=1;
         }
         else{
                for(p=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; p <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; p++)
                  stampa[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=0;
                i=0;
                counter=0;

                while(i < (int)strlen(temp_string)){
                  j=0;
                  while(temp_string[i] != ',' && i < (int)strlen(temp_string)){
                         temp_exon[j]=temp_string[i];
                         i++;
                         j++;
                  }
                  temp_exon[j]='\0';
                  if(j > 0){
                         if(counter==0){
                                first_exon=atoi(temp_exon);
                                for(p=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; p <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; p++){
                                  if(exists_list[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][first_exon])
                                         stampa[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=1;
                                }
                         }
                         else{
                                if(counter==1){
                                  second_exon=atoi(temp_exon);
                                  for(p=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; p <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; p++){
                                         if(exists_list[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][second_exon])
                                                stampa[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=(stampa[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]==0)?(stampa[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]):(1);
                                         else
                                                stampa[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=0;
                                  }
                                }
                         }
                         counter++;
                  }
                  i++;
                }

                for(p=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; p <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; p++){
                  if(stampa[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]){
                         strcpy(alternative_forms[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][alt_counter[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]], temp_string);
                         alt_counter[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]++;
                  }
                }

                scanf("%s\n", temp_string);

                for(p=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; p <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; p++){
                  if(stampa[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]){
                         strcpy(temp_string, "");
                         for(k=0; k<exon_in_trans_index[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][first_exon]; k++){
                                sprintf(temp_string2, "%d", exon_in_transcripts[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][first_exon][k]);
                                strcat(temp_string, temp_string2);
                                if(k != exon_in_trans_index[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][first_exon]-1)
                                  strcat(temp_string, ",");
                         }

                         strcpy(alternative_forms[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][alt_counter[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]], temp_string);
                         alt_counter[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]++;
                  }
                }

                scanf("%s\n", temp_string);

                for(p=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; p <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; p++){
                  if(stampa[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]){
                         strcpy(temp_string, "");
                         for(k=0; k<exon_in_trans_index[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][second_exon]; k++){
                                sprintf(temp_string2, "%d", exon_in_transcripts[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][second_exon][k]);
                                strcat(temp_string, temp_string2);
                                if(k != exon_in_trans_index[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][second_exon]-1)
                                  strcat(temp_string, ",");
                         }

                         strcpy(alternative_forms[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][alt_counter[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]], temp_string);
                         alt_counter[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]++;
                  }
                }

                for(p=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; p <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; p++){
                  if(!stampa[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]){
                         vect_alternative_number[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]--;
                  }
                }
         }
  }

//22mar07
//#ifndef MULTI_FASTA_FORMAT

  for(p=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; p <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; p++){
         fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%d\n", vect_alternative_number[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]);
         for(i=0; i<vect_alternative_number[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]*3; i++){
                fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "%s\n", alternative_forms[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][i]);
         }

         fprintf(tr_out[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT], "*\n");
  }
#endif

  for(i=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; i <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; i++){
         fclose(tr_out[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]);
  }

  MYTIME_stop(pt_tot);
  MYTIME_LOG(INFO, pt_tot);

  MYTIME_destroy(pt_tot);

  INFO("End");
  resource_usage_log();
}

/*FUNZIONI*/

void Get_Transcripts_from_File(){
  char temp_string[MAX_NLD];
  char temp_exon[7];
  char temp_conf[7];

//19dic06
//  char temp_type[7];
  char temp_refseq[20];

  char temp_coord[20];
  int exon_index=0;
  int i=0, j=0, k=0, p=0;
  int exons1=0, exons2=0;
  int counter=0;
  int temp_exon_list[MAX_EXONS];
  char stop=0, stop2=0;
  int count_exons=0;
  int coord_counter=0;
  int confirming=0;

//19dic06
  int type=0;

//  int max_lgth=0, lgth=0;

#ifdef READ_ABS_COORD
  scanf("%s\n", temp_string);
  gen_start=atoi(temp_string);
  scanf("%s\n", temp_string);
  gen_end=atoi(temp_string);
  scanf("%s\n", temp_string);
  strand=atoi(temp_string);
  scanf("%s\n", temp_string);
  boundary=atoi(temp_string);
#endif

//22mar07
  strcpy(init_reading, "");
  strcpy(init_reading2, "");

//Lettura delle prime tre righe
  for(i=0; i<3; i++){
         scanf("%s\n", temp_string);

         if(i==1){
                number_of_exons=atoi(temp_string);
         }

         if(i > 0){
#ifdef MULTI_FASTA_FORMAT
                if(i == 2){
#endif
                  strcat(init_reading, temp_string);
                  strcat(init_reading, "\n");
#ifdef MULTI_FASTA_FORMAT
                }
#endif
                strcat(init_reading2, temp_string);
                strcat(init_reading2, "\n");
         }
  }

  list_of_exon_left=(int *)malloc(number_of_exons*sizeof(int));
  list_of_exon_right=(int *)malloc(number_of_exons*sizeof(int));
  if(list_of_exon_left == NULL || list_of_exon_right == NULL){
         fprintf(stderr, "Problem11 of memory allocation in Get_Transcripts_from_File!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  list_of_old_exon_left=(int *)malloc(number_of_exons*sizeof(int));
  list_of_old_exon_right=(int *)malloc(number_of_exons*sizeof(int));
  if(list_of_old_exon_left == NULL || list_of_old_exon_right == NULL){
         fprintf(stderr, "Problem12 of memory allocation in Get_Transcripts_from_File!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

//10mag05
  polya=(char *)malloc(number_of_exons*sizeof(char));
  if(polya == NULL){
         fprintf(stderr, "Problem13 of memory allocation in Get_Transcripts_from_File!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  exists_list=(char **)malloc((SECOND_MIN_EXONS_ACCEPTED_OUTPUT-FIRST_MIN_EXONS_ACCEPTED_OUTPUT+1)*sizeof(char *));
  if(exists_list == NULL){
         fprintf(stderr, "Problem3 of memory allocation in Get_Transcripts_from_File!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }
  for(p=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; p <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; p++){
         exists_list[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=(char *)malloc(number_of_exons*sizeof(char));
         if(exists_list[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT] == NULL){
                fprintf(stderr, "Problem2 of memory allocation in Get_Transcripts_from_File!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }
         for(i=0; i<number_of_exons; i++)
                exists_list[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][i]=0;
  }

//03set04
  sequences=(char **)malloc(number_of_exons*sizeof(char *));
  if(sequences == NULL){
         fprintf(stderr, "Problem21 of memory allocation in Get_Transcripts_from_File!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

//22mar07
/*for(i=0; i<number_of_exons; i++){
  sequences[i]=(char *)malloc(MAX_NLD*sizeof(char));
  if(sequences[i] == NULL){
  fprintf(stderr, "Problem22 of memory allocation in Get_Transcripts_from_File!\n");
  #ifdef HALT_EXIT_MODE
  exit(1);
  #else
  exit(EXIT_FAILURE);
  #endif
  }
  }*/

//01set04
/*is_external=(char *)malloc(number_of_exons*sizeof(char));
  if(is_external == NULL){
  fprintf(stderr, "Problem21 of memory allocation in Get_Transcripts_from_File!\n");
  #ifdef HALT_EXIT_MODE
  exit(1);
  #else
  exit(EXIT_FAILURE);
  #endif
  }
  for(i=0; i<number_of_exons; i++)
  is_external[i]=0;*/
  is_internal=(char *)malloc(number_of_exons*sizeof(char));
  if(is_internal == NULL){
         fprintf(stderr, "Problem21 of memory allocation in Get_Transcripts_from_File!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }
  for(i=0; i<number_of_exons; i++)
         is_internal[i]=0;

//22mar07
#ifndef MULTI_FASTA_FORMAT

  exon_in_transcripts=(int ***)malloc((SECOND_MIN_EXONS_ACCEPTED_OUTPUT-FIRST_MIN_EXONS_ACCEPTED_OUTPUT+1)*sizeof(int **));
  if(exon_in_transcripts == NULL){
         fprintf(stderr, "Problem3 of memory allocation in Get_Transcripts_from_File!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }
  for(p=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; p <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; p++){
         exon_in_transcripts[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=(int **)malloc(number_of_exons*sizeof(int *));
         if(exon_in_transcripts[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT] == NULL){
                fprintf(stderr, "Problem3 of memory allocation in Get_Transcripts_from_File!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }
         for(i=0; i<number_of_exons; i++){
                exon_in_transcripts[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][i]=(int *)malloc(MAX_EXON_IN_TRANS*sizeof(int));
                if(exon_in_transcripts[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][i] == NULL){
                  fprintf(stderr, "Problem4 of memory allocation in Get_Transcripts_from_File!\n");
#ifdef HALT_EXIT_MODE
                  exit(1);
#else
                  exit(EXIT_FAILURE);
#endif
                }
                for(k=0; k<MAX_EXON_IN_TRANS; k++)
                  exon_in_transcripts[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][i][k]=0;
         }
  }

  exon_in_trans_index=(int **)malloc((SECOND_MIN_EXONS_ACCEPTED_OUTPUT-FIRST_MIN_EXONS_ACCEPTED_OUTPUT+1)*sizeof(int *));
  if(exon_in_trans_index == NULL){
         fprintf(stderr, "Problem3 of memory allocation in Get_Transcripts_from_File!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }
  for(p=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; p <= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; p++){
         exon_in_trans_index[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=(int *)malloc(number_of_exons*sizeof(int));
         if(exon_in_trans_index[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT] == NULL){
                fprintf(stderr, "Problem5 of memory allocation in Get_Transcripts_from_File!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }
         for(k=0; k<number_of_exons; k++)
                exon_in_trans_index[p-FIRST_MIN_EXONS_ACCEPTED_OUTPUT][k]=0;
  }
//22mar07
#endif

//06mar06
  if(number_of_exons == 0)
         stop=1;

//Lettura e scrittura della lista iniziale degli esoni
  while(!stop){
         scanf("%s\n", temp_string);
         if(temp_string[0] == '.'){
                stop=1;
         }
         else{
                if(count_exons > number_of_exons-1){
                  fprintf(stderr, "Problem3 of memory allocation in Get_Transcripts_from_File!\n");
#ifdef HALT_EXIT_MODE
                  exit(1);
#else
                  exit(EXIT_FAILURE);
#endif
                }

                i=0;
                coord_counter=0;

                while(i < (int)strlen(temp_string)){
                  j=0;
                  while(temp_string[i] != ':' && i < (int)strlen(temp_string)){
                         temp_coord[j]=temp_string[i];
                         i++;
                         j++;
                  }
                  temp_coord[j]='\0';
                  if(j > 0){
                         if(coord_counter==0){
                                list_of_exon_left[count_exons]=atoi(temp_coord);
                                list_of_old_exon_left[count_exons]=list_of_exon_left[count_exons];
                         }
                         else{
                                if(coord_counter==1){
                                  list_of_exon_right[count_exons]=atoi(temp_coord);
                                  list_of_old_exon_right[count_exons]=list_of_exon_right[count_exons];
                                }
                                else{
//10mag05
/*fprintf(stderr, "Problem4 of memory allocation in Get_Transcripts_from_File!\n");
  #ifdef HALT_EXIT_MODE
  exit(1);
  #else
  exit(EXIT_FAILURE);
  #endif*/
                                  if(coord_counter==2){
//21nov05
#ifdef IGNORE_POLYA
                                         polya[count_exons]=0;
#else
                                         polya[count_exons]=(char)atoi(temp_coord);
#endif
                                  }
                                  else{
//10mag05
                                         fprintf(stderr, "Problem4 of memory allocation in Get_Transcripts_from_File!\n");
#ifdef HALT_EXIT_MODE
                                         exit(1);
#else
                                         exit(EXIT_FAILURE);
#endif
                                  }
                                }
                         }
                         coord_counter++;
                  }
                  i++;
                }

                count_exons++;
         }
  }

/*for(p=0; p<number_of_exons; p++){
  if(list_of_exon_right[p]-list_of_exon_left[p]+1 > MAX_NLD){
  fprintf(stderr, "Problem with MAX_NLD!\n");
  exit(0);
  }
  }*/

//fprintf(stderr, "MAX EXON LENGTH: %d\n", max_lgth);


//Lettura dei trascritti
//06mar06
  if(number_of_exons != 0)
         stop=0;

  counter=0;
  while(!stop){
         stop2=0;
         exons1=0;
         exons2=0;

         if(counter > MAX_TRANSCRIPTS-1){
                fprintf(stderr, "Too many transcripts in input!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }

         i=1;
//19dic06
//while(i < strlen(temp_string)){
         while(temp_string[i] != '.' && i < (int)strlen(temp_string)){
                temp_conf[i-1]=temp_string[i];
                i++;
         }
         temp_conf[i-1]='\0';
         confirming=atoi(temp_conf);

//19dic06
         if(i < (int)strlen(temp_string))
                i++;
         while(i < (int)strlen(temp_string)){
                temp_refseq[i-strlen(temp_conf)-2]=temp_string[i];
                i++;
         }

         if(i > (int)strlen(temp_conf)+2){
                temp_refseq[i-strlen(temp_conf)-2]='\0';
                type=1;
         }
         else{
                strcpy(temp_refseq, "");
                type=0;
         }

//19dic06
#ifndef DONT_EXTEND_REFSEQ
         type=0;
         strcpy(temp_refseq, "");
#endif

         scanf("%s\n", temp_string);

         i=0;

         while(i < (int)strlen(temp_string)){
                j=0;
                while(temp_string[i] != '.' && i < (int)strlen(temp_string)){
                  temp_exon[j]=temp_string[i];
                  i++;
                  j++;
                }
                temp_exon[j]='\0';
                if(j > 0){
                  exon_index=atoi(temp_exon);
                  temp_exon_list[exons1]=exon_index;
                  exons1++;
                }
                i++;
         }

         while(!stop2){
                scanf("%s\n", temp_string);
                if(temp_string[0] == '.' || temp_string[0] == '#')
                  stop2=1;
                else{
                  if(exons2 > MAX_EXONS-1){
                         fprintf(stderr, "Too many exons in input transcript number %d!\n", (counter+1));
#ifdef HALT_EXIT_MODE
                         exit(1);
#else
                         exit(EXIT_FAILURE);
#endif
                  }
//03set04
//strcpy(transcript_list[counter].sequence[exons2], temp_string);

//22mar07
                  if(sequences[temp_exon_list[exons2]] == NULL){
//23gen08
//sequences[temp_exon_list[exons2]]=(char *)malloc((strlen(temp_string)+20)*sizeof(char));
//18feb08
//sequences[temp_exon_list[exons2]]=(char *)malloc((strlen(temp_string)+1000)*sizeof(char));
                         sequences[temp_exon_list[exons2]]=(char *)malloc((strlen(temp_string)+5000)*sizeof(char));
                         if(sequences[temp_exon_list[exons2]] == NULL){
                                fprintf(stderr, "Problem221 of memory allocation in Get_Transcripts_from_File!\n");
#ifdef HALT_EXIT_MODE
                                exit(1);
#else
                                exit(EXIT_FAILURE);
#endif
                         }
                         strcpy(sequences[temp_exon_list[exons2]], temp_string);
                  }

//22mar07
//strcpy(sequences[temp_exon_list[exons2]], temp_string);

                  exons2++;
                }
         }

         if(exons1 == 0 || exons2 == 0 || exons1 != exons2){
                fprintf(stderr, "Invalid transcript in input file (%d)!\n", (counter+1));
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }

//09mag08
//if(exons1 >= MIN_EXONS_ACCEPTED_INPUT && confirming >= MIN_CONFIRMED_EST_INPUT){
         if((exons1 >= MIN_EXONS_ACCEPTED_INPUT && confirming >= MIN_CONFIRMED_EST_INPUT) && !(exons1 == 1 && type != 1)){

                transcript_list[counter].exons=exons1;
                transcript_list[counter].ESTs=confirming;

//19dic06
                transcript_list[counter].type=type;
                strcpy(transcript_list[counter].RefSeq, temp_refseq);

                transcript_list[counter].left_ext=temp_exon_list[0];

//16feb06
                if(polya[transcript_list[counter].left_ext] == 1){
                  polya[transcript_list[counter].left_ext]=0;
                }

//01set04
                if(is_internal[temp_exon_list[0]] != 1){
                  if(is_internal[temp_exon_list[0]] == 0){
//09mag08
//is_internal[temp_exon_list[0]]=-1;
                         if(exons1 == 1){
                                is_internal[temp_exon_list[0]]=-3;
                         }
                         else
                                is_internal[temp_exon_list[0]]=-1;
                  }
                  else{
                         if(is_internal[temp_exon_list[0]] == -2){
//09mag08
//is_internal[temp_exon_list[0]]=1;
                                if(exons1 > 1)
                                  is_internal[temp_exon_list[0]]=1;
                         }
//09mag08
                         else{
                                if(is_internal[temp_exon_list[0]] != -1){
                                  if(exons1 > 1)
                                         is_internal[temp_exon_list[0]]=-1;
                                }
                         }
                  }
                }

                for(k=1; k<exons1-1; k++){
//01set04
                  is_internal[temp_exon_list[k]]=1;
                  transcript_list[counter].exon_list[k-1]=temp_exon_list[k];

//16feb06
                  if(polya[transcript_list[counter].exon_list[k-1]] == 1){
                         polya[transcript_list[counter].exon_list[k-1]]=0;
                  }
                }

                transcript_list[counter].right_ext=temp_exon_list[exons1-1];
//01set04
//09mag08
//if(is_internal[temp_exon_list[exons1-1]] != 1){
                if(exons1 > 1 && is_internal[temp_exon_list[exons1-1]] != 1){
                  if(is_internal[temp_exon_list[exons1-1]] == 0){
//10mag05
//is_internal[temp_exon_list[exons1-1]]=-2;
                         if(polya[temp_exon_list[exons1-1]] == 1)
                                is_internal[temp_exon_list[exons1-1]]=1;
                         else
                                is_internal[temp_exon_list[exons1-1]]=-2;
                  }
                  else{
                         if(is_internal[temp_exon_list[exons1-1]] == -1)
                                is_internal[temp_exon_list[exons1-1]]=1;
                  }
                }

                counter++;
         }
         if(temp_string[0] == '#')
                stop=1;
  }

  number_of_transcripts=counter;
}

void Build_Extension_Matrix(){

  int i=0, j=0;
  int extension_limit=0;
  char exists_extension=0;

  extension_matrix=(int **)malloc(number_of_transcripts*sizeof(int *));
  if(extension_matrix == NULL){
         fprintf(stderr, "Problem1 of memory allocation in Build_Extension_Matrix!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  for(i=0; i<number_of_transcripts; i++){
         extension_matrix[i]=(int *)malloc(number_of_transcripts*sizeof(int));
         if(extension_matrix[i] == NULL){
                fprintf(stderr, "Problem2 of memory allocation in Build_Extension_Matrix!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }
  }

//27giu05
  remove_extension_matrix=(char **)malloc(number_of_transcripts*sizeof(char *));
  if(remove_extension_matrix == NULL){
         fprintf(stderr, "Problem11 of memory allocation in Build_Extension_Matrix!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

//27giu05
  for(i=0; i<number_of_transcripts; i++){
         remove_extension_matrix[i]=(char *)malloc(number_of_transcripts*sizeof(char));
         if(remove_extension_matrix[i] == NULL){
                fprintf(stderr, "Problem2 of memory allocation in Build_Extension_Matrix!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }
  }

  is_source=(char *)malloc(number_of_transcripts*sizeof(char));
  if(is_source == NULL){
         fprintf(stderr, "Problem3 of memory allocation in Build_Extension_Matrix!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  for(i=0; i<number_of_transcripts; i++)
         for(j=0; j<number_of_transcripts; j++){
//27giu05
                remove_extension_matrix[i][j]=0;

                if(i == j)
                  extension_matrix[i][j]=0;
                else
                  extension_matrix[i][j]=-1;
         }

  for(i=0; i<number_of_transcripts; i++)
         is_source[i]=1;

  in_degree=(int *)malloc(number_of_transcripts*sizeof(int));
  out_degree=(int *)malloc(number_of_transcripts*sizeof(int));
  if(in_degree == NULL || out_degree == NULL){
         fprintf(stderr, "Problem30 of memory allocation in Build_Extension_Matrix!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  for(i=0; i<number_of_transcripts; i++){
         in_degree[i]=0;
         out_degree[i]=0;
  }

  give_path=(char *)malloc(number_of_transcripts*sizeof(char));
  if(give_path == NULL){
         fprintf(stderr, "Problem31 of memory allocation in Build_Extension_Matrix!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  for(i=0; i<number_of_transcripts; i++)
         give_path[i]=0;

  for(i=0; i<number_of_transcripts; i++){
         for(j=i+1; j<number_of_transcripts; j++)
                if(extension_matrix[i][j] == -1){
/*if(i==19 && j==32){
  fprintf(stdout, "******************************************* %d %d\n", i, j);
  fprintf(stdout, "HERE prima 38 %d %d\n", list_of_exon_left[38], list_of_exon_right[38]);
  fprintf(stdout, "HERE 37 %d %d\n", list_of_exon_left[37], list_of_exon_right[37]);
  }*/

//08giu05
//exists_extension=Extends(transcript_list[i], transcript_list[j], &extension_limit, 0);

//07ott05
//exists_extension=Extends(transcript_list[i], transcript_list[j], &extension_limit, 0, 1);
//19dic06
//exists_extension=Extends(transcript_list[i], transcript_list[j], &extension_limit, 0, 1, 1);
#ifdef DONT_EXTEND_REFSEQ
                  if(transcript_list[i].type == 1){
                         if(transcript_list[j].type == 0){
//19feb07
//exists_extension=Extends(transcript_list[i], transcript_list[j], &extension_limit, 0, 1, 1);
                                exists_extension=0;
                         }
                         else
                                exists_extension=0;
                  }
                  else{
                         if(transcript_list[j].type == 1){
//19feb07
/*exists_extension=Extends(transcript_list[j], transcript_list[i], &extension_limit, 0, 1, 1);
  if(exists_extension == 1)
  exists_extension=-1;*/
                                exists_extension=0;
                         }
                         else
//08gen08
//exists_extension=Extends(transcript_list[i], transcript_list[j], &extension_limit, 0, 1, 1);
                                exists_extension=Extends(transcript_list[i], transcript_list[j], &extension_limit, 0, 1, 1, 0);
                  }
#else
//08gen08
//exists_extension=Extends(transcript_list[i], transcript_list[j], &extension_limit, 0, 1, 1);
                  exists_extension=Extends(transcript_list[i], transcript_list[j], &extension_limit, 0, 1, 1, 0);
#endif

/*if(i==19 && j==32){
  fprintf(stderr, "******************************************* %d %d\n", i, j);
  fprintf(stderr, "EXT %d\n", exists_extension);
  fprintf(stderr, "*******************************************\n");
  }*/
                  extension_matrix[i][j]=0;
                  extension_matrix[j][i]=0;
                  if(exists_extension == 1 || exists_extension == -1){
                         if(exists_extension == 1){
                                extension_matrix[i][j]=extension_limit;
//is_source[j]=0;
                                out_degree[i]=out_degree[i]+1;
                                in_degree[j]=in_degree[j]+1;
                         }
                         else{
                                extension_matrix[j][i]=extension_limit;
                                out_degree[j]=out_degree[j]+1;
                                in_degree[i]=in_degree[i]+1;
//is_source[i]=0;
                         }
                  }
                }
  }

/*      //Conteggio delle sorgenti
        number_of_sources=0;
        for(i=0; i<number_of_transcripts; i++){
        if(is_source[i])
        number_of_sources+=1;
        }

        source_list=(int *)malloc(number_of_sources*sizeof(int));
        if(source_list== NULL){
        fprintf(stderr, "Problem4 of memory allocation in Build_Extension_Matrix!\n");
        #ifdef HALT_EXIT_MODE
        exit(1);
        #else
        exit(EXIT_FAILURE);
        #endif
        }

        j=0;
        for(i=0; i<number_of_transcripts; i++){
        if(is_source[i]){
        source_list[j]=i;
        j++;
        }
        }*/
}

//CONTROLLARE
struct path *Copy_of_Path(struct path *path){
  struct path *copy_path=NULL;
  struct node *head=NULL;

  copy_path=(struct path *)malloc(sizeof(struct path));
  if(copy_path == NULL){
         fprintf(stderr, "Problem1 of memory allocation in Copy_of_Path!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  copy_path->next=NULL;
  copy_path->n=NULL;
  copy_path->tail=NULL;

  head=path->n;

  while(head != NULL){
         Add_Node(&copy_path, head->index, 0);
         head=head->next;
  }

  Copy_transcript(path->tr, &(copy_path->tr));
  copy_path->L=path->L;
  copy_path->visit=path->visit;

  return copy_path;
}

//CONTROLLARE
//upd_tr=TRUE, se il corrispondente trascritto deve essere aggiornato
void Add_Node(struct path **path, int node, char upd_tr){
  struct node *nds=NULL;
  int counter=0;
  struct transcript extension;
  char first_node=0;

  nds=(*path)->n;
  while(nds != NULL){
         if(nds->index == node){
                fprintf(stderr, "Cycle detected!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }
         nds=nds->next;
         counter++;
  }

//07apr06
  if(counter == 40){
         fprintf(stderr, "Too many nodes!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  add=(struct node *)malloc(sizeof(struct node));
  if(add == NULL){
         fprintf(stderr, "Problem1 of memory allocation in Add_Node!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  add->index=node;
  add->next=NULL;

  if((*path)->tail != NULL){
         (*path)->tail->next=add;
  }
  else{
         (*path)->n=add;
         (*path)->L=0;
         first_node=1;
         Copy_transcript(transcript_list[node], &((*path)->tr));
  }

  if(upd_tr && !first_node){
         (*path)->L+=extension_matrix[(*path)->end][node];
         Build_extension((*path)->tr, transcript_list[node], (*path)->L, &extension);
         Copy_transcript(extension, &((*path)->tr));
  }

  (*path)->tail=add;
  (*path)->end=node;
}

void Set_Paths(){
  int i=0;

  total_paths=0;

  for(i=0; i<number_of_sources; i++){
         Set_Paths_for_Source(i);
         total_paths=total_paths+source_total_paths;
         Set_Path_Transcripts();
  }

  filtered=(char *)malloc(total_paths*sizeof(char));
  if(filtered== NULL){
         fprintf(stderr, "Problem of memory allocation in Set_Paths!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }
  for(i=0; i<total_paths; i++)
         filtered[i]=0;
}

void Set_Paths_for_Source(int source_index){
  int i=0;
  struct path *enq_path=NULL;
  struct path *source_p=NULL;
  struct path *copy=NULL;
  char no_edge=0;

  struct node *p_free=NULL;
  struct node *help_p_free=NULL;

//  struct node *head=NULL;
//27giu05
#ifdef PRUNE_EXON_COMP
  struct path *same_exons_path=NULL;
#endif

//Svuota source_list_of_paths
  Empty_Source_List_of_Paths();

  source_total_paths=0;

  InitializeQueue();

//source_list[source_index] e' il vero indice di tale sorgente come trascritto e come nodo
  source_p=Create_Source_Path(source_list[source_index]);
  source_p->visit=1;

  give_path[source_list[source_index]]=1;

  Enqueue(&source_p);

  while(!Queue_Is_Empty()){
         enq_path=Dequeue();

/*head=enq_path->n;
  fprintf(stdout, "Path********************************\n");
  while(head != NULL){
  fprintf(stdout, "%d->", head->index);
  head=head->next;
  }
  fprintf(stdout, "\nPath********************************\n");*/

         no_edge=1;

//Considero enq_path solo se il suo visit e' a 1

//27giu05
#ifdef PRUNE_EXON_COMP
         if(enq_path->visit == 1){
#endif

                for(i=0; i<number_of_transcripts; i++){
//Esiste un arco da u a i
                  if(extension_matrix[enq_path->end][i] != 0){
//27giu05
#ifdef PRUNE_EXON_COMP
                         same_exons_path=NULL;
#endif

                         no_edge=0;

                         copy=Copy_of_Path(enq_path);
                         Add_Node(&copy, i, 1);

                         give_path[i]=1;

//27giu05
#ifdef PRUNE_EXON_COMP
                         same_exons_path=Get_path_with_the_same_exons(copy);

                         if(same_exons_path != NULL){
                                if(out_degree[copy->end] > out_degree[same_exons_path->end]){
                                  same_exons_path->visit=0;
                                  Enqueue(&copy);
                                }
                         }
                         else{
#endif
                                Enqueue(&copy);
//27giu05
#ifdef PRUNE_EXON_COMP
                         }
#endif
                  }
                }

//Se non ci sono archi uscenti dal nodo enq_path->end, enq_path e' terminato
                if(no_edge){
                  Set_Path_Transcripts_for_Source(enq_path);

/*head=enq_path->n;
  fprintf(stdout, "\nPath number %d********************************\n", source_total_paths);
  while(head != NULL){
  fprintf(stdout, "             %d,", head->index);
  head=head->next;
  }
  head=enq_path->n;
  while(head != NULL){
  fprintf(stdout, "\n   Transcript %d\n", head->index);
  fprintf(stdout, "             %d", transcript_list[head->index].left_ext);
  for(k=0; k<transcript_list[head->index].exons-2; k++)
  fprintf(stdout, ".%d", transcript_list[head->index].exon_list[k]);
  fprintf(stdout, ".%d\n", transcript_list[head->index].right_ext);
  fprintf(stdout, "\n");
  head=head->next;
  }
  fprintf(stdout, "\n   TOTAL %d\n", source_total_paths);*/
                }
//27giu05
#ifdef PRUNE_EXON_COMP
         }
#endif

//Libero la memoria di enq_path
         p_free=enq_path->n;
         while(p_free != NULL){
                help_p_free=p_free->next;
                free(p_free);
                p_free=help_p_free;
         }
         free(enq_path);
  }
}

void InitializeQueue(){
  struct path *head=q.head;
  struct path *help=NULL;

  while(head != NULL){
         help=head->next;
         head=help;
  }

  q.head=NULL;
  q.tail=NULL;
}

char Queue_Is_Empty(){
  return (q.head == NULL);
}

void Enqueue(struct path **path){

  if(q.tail == NULL)
         q.head=*path;
  else
         q.tail->next=*path;

  q.tail=*path;
}

struct path *Dequeue(){
  struct path *value=NULL;
  struct path *help=NULL;

  if(q.head != NULL){
         value=q.head;
         help=q.head;
         q.head=q.head->next;

         if(q.head == NULL)
                q.tail=NULL;
  }

  return value;
}

struct path *Get_path_with_the_same_exons(struct path *arg_path){
  struct path *head=q.head;
  char found=0;

  while(head != NULL && !found){
         if(Equals_transcripts(arg_path->tr, head->tr))
                found=1;
         else
                head=head->next;
  }

  return head;
}

char Equals_transcripts(struct transcript t1, struct transcript t2){
  int k=0;
  char stop=0;

  if(t1.exons != t2.exons)
         return 0;

  if(t1.left_ext != t2.left_ext || t1.right_ext != t2.right_ext)
         return 0;

  while(k < t1.exons-2 && !stop){
         if(t1.exon_list[k] == t2.exon_list[k])
                k++;
         else
                stop=1;
  }

  if(!stop)
         return 1;
  else
         return 0;
}

void Build_extension(struct transcript t1, struct transcript t2, int L, struct transcript *extension){
  int i=0;
//  int k=0;

  extension->exons=t2.exons+L;

  if(extension->exons > MAX_EXONS){
         fprintf(stderr, "Too many exons in Extension!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

/*fprintf(stdout, "\n\nBUILD1 L %d:\n", L);
  for(k=0; k<t1.exons; k++){
  if(k == 0)
  fprintf(stdout, "             %d", t1.left_ext);
  else{
  if(k == t1.exons-1)
  fprintf(stdout, ".%d", t1.right_ext);
  else
  fprintf(stdout, ".%d", t1.exon_list[k-1]);
  }
  }
  fprintf(stdout, "\nBUILD2:\n");
  for(k=0; k<t2.exons; k++){
  if(k == 0)
  fprintf(stdout, "             %d", t2.left_ext);
  else{
  if(k == t2.exons-1)
  fprintf(stdout, ".%d", t2.right_ext);
  else
  fprintf(stdout, ".%d", t2.exon_list[k-1]);
  }
  }*/

  extension->left_ext=t1.left_ext;

//Riempio, con gli elementi di t1 i primi elementi di extension
  for(i=0; i<t1.exons-2; i++){
         extension->exon_list[i]=t1.exon_list[i];
  }

//21nov05
  if(i-L < 0){
         if(is_internal[t1.right_ext] == 1 || is_internal[t2.left_ext] != 1){
                extension->exon_list[i]=t1.right_ext;
         }
         else{
                extension->exon_list[i]=t2.left_ext;
         }
         i++;
  }

//Riempio, con gli ultimi elementi di t2 (quelli dopo l'overlap con t1) i rimanenti elementi di extension
  for(; i<extension->exons-2; i++){
         extension->exon_list[i]=t2.exon_list[i-L];
  }

  extension->right_ext=t2.right_ext;

//19dic06. Necessariamente type e' 0 in quanto i refseq non sono estendibili
  extension->type=0;
  strcpy(extension->RefSeq, "");

/*fprintf(stdout, "\nBUILDEXT:\n");
  for(k=0; k<extension->exons; k++){
  if(k == 0)
  fprintf(stdout, "             %d", extension->left_ext);
  else{
  if(k == extension->exons-1)
  fprintf(stdout, ".%d", extension->right_ext);
  else
  fprintf(stdout, ".%d", extension->exon_list[k-1]);
  }
  }
  fprintf(stdout, "\nEND!!!\n");*/

}

void Copy_transcript(struct transcript t, struct transcript *copy){
  int i=0;

  copy->exons=t.exons;

  copy->left_ext=t.left_ext;

//19dic06
  copy->type=t.type;
  strcpy(copy->RefSeq, t.RefSeq);

//Riempio, con gli elementi di t i gli elementi di copy
  for(i=0; i<t.exons-2; i++){
         copy->exon_list[i]=t.exon_list[i];
  }

  copy->right_ext=t.right_ext;
}

/* //MODIFICA: fare in modo da non usarla piu' */
/* void Get_Transcript_From_Path(struct path *tr_path, struct transcript *return_transcript){ */
/*   struct node *head=tr_path->n; */
/*   struct transcript extension; */
/*   struct transcript copy; */
/*   int start=0, end=0, L=0; */
/*   char exists_path=1; */

/*   int k=0; */

/*   if(head != NULL){ */
/*       start=head->index; */
/*       Copy_transcript(transcript_list[start], &copy); */


/* /\*fprintf(stdout, "START:\n"); */
/*   for(k=0; k<copy.exons; k++){ */
/*   if(k == 0) */
/*   fprintf(stdout, "          %d", copy.left_ext); */
/*   else{ */
/*   if(k == copy.exons-1) */
/*   fprintf(stdout, ".%d", copy.right_ext); */
/*   else */
/*   fprintf(stdout, ".%d", copy.exon_list[k-1]); */
/*   } */
/*   } */
/*   fprintf(stdout, "\nSTART:\n");*\/ */

/*   } */
/*   else */
/*       exists_path=0; */

/*   while(head != NULL){ */
/*       start=head->index; */
/*       head=head->next; */
/*       if(head != NULL){ */
/*              end=head->index; */
/*              L+=extension_matrix[start][end]; */
/*              Build_extension(copy, transcript_list[end], L, &extension); */
/*              Copy_transcript(extension, &copy); */

/* /\*fprintf(stdout, "\nWITH:\n"); */
/*   for(k=0; k<transcript_list[end].exons; k++){ */
/*   if(k == 0) */
/*   fprintf(stdout, "          %d", transcript_list[end].left_ext); */
/*   else{ */
/*   if(k == transcript_list[end].exons-1) */
/*   fprintf(stdout, ".%d", transcript_list[end].right_ext); */
/*   else */
/*   fprintf(stdout, ".%d", transcript_list[end].exon_list[k-1]); */
/*   } */
/*   } */
/*   fprintf(stdout, "\nWITH:\n"); */

/*   fprintf(stdout, "\nEXT:\n"); */
/*   for(k=0; k<copy.exons; k++){ */
/*   if(k == 0) */
/*   fprintf(stdout, "          %d", copy.left_ext); */
/*   else{ */
/*   if(k == copy.exons-1) */
/*   fprintf(stdout, ".%d", copy.right_ext); */
/*   else */
/*   fprintf(stdout, ".%d", copy.exon_list[k-1]); */
/*   } */
/*   } */
/*   fprintf(stdout, "\nEXT:\n");*\/ */


/*       } */
/*   } */

/*   if(exists_path) */
/*       Copy_transcript(copy, return_transcript); */
/* } */

void Set_Path_Transcripts_for_Source(struct path *path){
  int i=0;
//struct transcript return_transcript;
  char included=0;
  int type=0;
  char stop=0;
//  struct node *head=NULL;

//DEBUGGING
/*fprintf(stdout, "\n**************************Path:\n");
  head=path->n;
  while(head != NULL){
  fprintf(stdout, "             %d", head->index);
  for(k=0; k < transcript_list[head->index].exons; k++){
  if(k == 0)
  fprintf(stdout, "             %d", transcript_list[head->index].left_ext);
  else{
  if(k == transcript_list[head->index].exons-1)
  fprintf(stdout, ".%d", transcript_list[head->index].right_ext);
  else
  fprintf(stdout, ".%d", transcript_list[head->index].exon_list[k-1]);
  }
  }
  fprintf(stdout, "\n");
  head=head->next;
  }
  fprintf(stdout, "\n");
  fprintf(stdout, "     TR:\n");
  for(k=0; k < path->tr.exons; k++){
  if(k == 0)
  fprintf(stdout, "             %d", path->tr.left_ext);
  else{
  if(k == path->tr.exons-1)
  fprintf(stdout, ".%d", path->tr.right_ext);
  else
  fprintf(stdout, ".%d", path->tr.exon_list[k-1]);
  }
  }
  fprintf(stdout, "\n**************************\n");*/


//MODIFICA
//Get_Transcript_From_Path(path, &return_transcript);

/*fprintf(stdout, "             MAX TR ");
  for(k=0; k < return_transcript.exons; k++){
  if(k == 0)
  fprintf(stdout, "             %d", return_transcript.left_ext);
  else{
  if(k == return_transcript.exons-1)
  fprintf(stdout, ".%d", return_transcript.right_ext);
  else
  fprintf(stdout, ".%d", return_transcript.exon_list[k-1]);
  }
  }
  fprintf(stdout, "\n");*/

  i=0;
  while(i < source_total_paths && !stop){
//MODIFICA
//included=Extends(source_path_transcript_list[i], return_transcript, &type, 0);
//08giu05
//included=Extends(source_path_transcript_list[i], path->tr, &type, 0);

//07ott05
//included=Extends(source_path_transcript_list[i], path->tr, &type, 0, 0);
//08gen08
//included=Extends(source_path_transcript_list[i], path->tr, &type, 0, 0, 1);
         included=Extends(source_path_transcript_list[i], path->tr, &type, 0, 0, 1, 0);

         if(included == 2 || included == -2){

/*fprintf(stdout, "     TR before modified:\n");
  for(k=0; k < path->tr.exons; k++){
  if(k == 0)
  fprintf(stdout, "             %d", path->tr.left_ext);
  else{
  if(k == path->tr.exons-1)
  fprintf(stdout, ".%d", path->tr.right_ext);
  else
  fprintf(stdout, ".%d", path->tr.exon_list[k-1]);
  }
  }
  fprintf(stdout, "\n**************************\n");*/

//25lug07
//Modificare per allungo esoni
                if(included == 2){
//Modificare source_path_transcript_list[i] con esoni piu' lunghi
//23gen08
                  if(list_of_exon_right[source_path_transcript_list[i].left_ext] == list_of_exon_right[path->tr.left_ext]){
                         if(type == 0){
//06nov07
//if(is_internal[source_path_transcript_list[i].left_ext == -1]){
                                if(is_internal[source_path_transcript_list[i].left_ext] == -1){
                                  if(is_internal[path->tr.left_ext] == 1){
                                         source_path_transcript_list[i].left_ext=path->tr.left_ext;
                                  }
                                  else{
                                         if(is_internal[path->tr.left_ext] == -1){
                                                if(list_of_exon_left[path->tr.left_ext] < list_of_exon_left[source_path_transcript_list[i].left_ext]){
                                                  source_path_transcript_list[i].left_ext=path->tr.left_ext;
                                                }
                                         }
                                  }
                                }
                         }
//23gen08
                  }

//23gen08
                  if(list_of_exon_left[source_path_transcript_list[i].right_ext] == list_of_exon_left[path->tr.right_ext]){
                         if(type+path->tr.exons == source_path_transcript_list[i].exons){
                                if(is_internal[source_path_transcript_list[i].right_ext == -2]){
                                  if(is_internal[path->tr.right_ext] == 1){
                                         source_path_transcript_list[i].right_ext=path->tr.right_ext;
                                  }
                                  else{
                                         if(is_internal[path->tr.right_ext] == -2){
                                                if(list_of_exon_right[path->tr.right_ext] > list_of_exon_right[source_path_transcript_list[i].right_ext]){
                                                  source_path_transcript_list[i].right_ext=path->tr.right_ext;
                                                }
                                         }
                                  }
                                }
                         }
//23gen08
                  }
                }
                else{
//Modificare path con esoni piu' lunghi
//23gen08
                  if(list_of_exon_right[source_path_transcript_list[i].left_ext] == list_of_exon_right[path->tr.left_ext]){
                         if(type == 0){
//06nov07
//if(is_internal[path->tr.left_ext == -1]){
                                if(is_internal[path->tr.left_ext] == -1){
                                  if(is_internal[source_path_transcript_list[i].left_ext] == 1){
                                         path->tr.left_ext=source_path_transcript_list[i].left_ext;
                                  }
                                  else{
                                         if(is_internal[source_path_transcript_list[i].left_ext] == -1){
                                                if(list_of_exon_left[source_path_transcript_list[i].left_ext] < list_of_exon_left[path->tr.left_ext]){
                                                  path->tr.left_ext=source_path_transcript_list[i].left_ext;
                                                }
                                         }
                                  }
                                }
                         }
//23gen08
                  }

//23gen08
                  if(list_of_exon_left[source_path_transcript_list[i].right_ext] == list_of_exon_left[path->tr.right_ext]){
                         if(type+source_path_transcript_list[i].exons == path->tr.exons){
                                if(is_internal[path->tr.right_ext == -2]){
                                  if(is_internal[source_path_transcript_list[i].right_ext] == 1){
                                         path->tr.right_ext=source_path_transcript_list[i].right_ext;
                                  }
                                  else{
                                         if(is_internal[source_path_transcript_list[i].right_ext] == -2){
                                                if(list_of_exon_right[source_path_transcript_list[i].right_ext] > list_of_exon_right[path->tr.right_ext]){
                                                  path->tr.right_ext=source_path_transcript_list[i].right_ext;
                                                }
                                         }
                                  }
                                }
                         }
//23gen08
                  }
                }

                stop=1;
                if(included == -2){


/*fprintf(stdout, "                     FILTERING ");
  for(k=0; k < source_path_transcript_list[i].exons; k++){
  if(k == 0)
  fprintf(stdout, "             %d", source_path_transcript_list[i].left_ext);
  else{
  if(k == source_path_transcript_list[i].exons-1)
  fprintf(stdout, ".%d", source_path_transcript_list[i].right_ext);
  else
  fprintf(stdout, ".%d", source_path_transcript_list[i].exon_list[k-1]);
  }
  }
  fprintf(stdout, "\n");
  fprintf(stdout, "     TR modified:\n");
  for(k=0; k < path->tr.exons; k++){
  if(k == 0)
  fprintf(stdout, "             %d", path->tr.left_ext);
  else{
  if(k == path->tr.exons-1)
  fprintf(stdout, ".%d", path->tr.right_ext);
  else
  fprintf(stdout, ".%d", path->tr.exon_list[k-1]);
  }
  }
  fprintf(stdout, "\n**************************\n");*/

//MODIFICA
//Copy_transcript(return_transcript, &source_path_transcript_list[i]);
                  Copy_transcript(path->tr, &source_path_transcript_list[i]);
//AGGIUNGERE IL PERCORSO "path" ALLA LISTA DEI PERCORSI in source_list_of_paths[i]
                  Add_Path(&source_list_of_paths[i], &path);
                }
         }
         else{
                i++;
         }
  }

  if(!stop){
         if(source_total_paths > MAX_TRANSCRIPTS-1){
                fprintf(stderr, "Too many path in source transcript graph!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }
//MODIFICA
//Copy_transcript(return_transcript, &source_path_transcript_list[source_total_paths]);
         Copy_transcript(path->tr, &source_path_transcript_list[source_total_paths]);
//AGGIUNGERE IL PERCORSO "path" ALLA LISTA DEI PERCORSI in source_list_of_paths[source_total_paths]
         Add_Path(&source_list_of_paths[source_total_paths], &path);
         source_total_paths=source_total_paths+1;
  }
}

void Set_Path_Transcripts(){
  int i=0;

  if(total_paths > MAX_TRANSCRIPTS){
         fprintf(stderr, "Too many path in transcript graph!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  for(i=0; i<source_total_paths; i++){
         Copy_transcript(source_path_transcript_list[i], &path_transcript_list[transcript_counter]);
//Metti i paths in source_list_of_paths[i] in transcript_list_of_paths[transcript_counter]
         Add_Path_List(&transcript_list_of_paths[transcript_counter], &source_list_of_paths[i]);
         transcript_counter++;
  }

  if(transcript_counter != total_paths){
         fprintf(stderr, "Problem2 of memory allocation in Set_Path_Transcripts!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }
}

void Filter_Path_Transcripts(){
  int i=0, j=0, k=0;
  char included=0;
  char stop=0;
  int type=0;

  actual_path_number=total_paths;

  while(i<total_paths){

/*fprintf(stdout, "FILTERING %d ", i);
  for(k=0; k < path_transcript_list[i].exons; k++){
  if(k == 0)
  fprintf(stdout, "             %d", path_transcript_list[i].left_ext);
  else{
  if(k == path_transcript_list[i].exons-1)
  fprintf(stdout, ".%d", path_transcript_list[i].right_ext);
  else
  fprintf(stdout, ".%d", path_transcript_list[i].exon_list[k-1]);
  }
  }
  fprintf(stdout, "\n");*/

         if(!filtered[i]){
/*if(path_transcript_list[i].exons == 1)
  actual_path_number--;*/
                j=i+1;
                stop=0;
                while(j<total_paths && !stop){
                  if(!filtered[j]){


/*fprintf(stdout, "             with %d ", j);
  for(k=0; k < path_transcript_list[j].exons; k++){
  if(k == 0)
  fprintf(stdout, "             %d", path_transcript_list[j].left_ext);
  else{
  if(k == path_transcript_list[j].exons-1)
  fprintf(stdout, ".%d", path_transcript_list[j].right_ext);
  else
  fprintf(stdout, ".%d", path_transcript_list[j].exon_list[k-1]);
  }
  }
  fprintf(stdout, "\n");*/

//07ott05
#ifdef MERGE_POLYA
//19dic06
//included=Extends(path_transcript_list[i], path_transcript_list[j], &type, 0, 0, 0);
#ifdef DONT_EXTEND_REFSEQ
                         if(path_transcript_list[i].type == 1){
                                if(path_transcript_list[j].type == 0){
//08gen08
//included=Overlap(path_transcript_list[j], path_transcript_list[i], &type, 0, 0, 0);
                                  included=Overlap(path_transcript_list[j], path_transcript_list[i], &type, 0, 0, 0, 1);
                                }
                                else{
//08gen08
                                  included=Overlap(path_transcript_list[i], path_transcript_list[j], &type, 0, 0, 0, 1);
                                  included=0;
                                }
                         }
                         else{
                                if(path_transcript_list[j].type == 1){
//08gen08
//included=Overlap(path_transcript_list[i], path_transcript_list[j], &type, 0, 0, 0);
                                  included=Overlap(path_transcript_list[i], path_transcript_list[j], &type, 0, 0, 0, 1);
                                  if(included == 2)
                                         included=-2;
                                }
                                else
//08gen08
//included=Extends(path_transcript_list[i], path_transcript_list[j], &type, 0, 0, 0);
                                  included=Extends(path_transcript_list[i], path_transcript_list[j], &type, 0, 0, 0, 1);
                         }
#else
//08gen08
//included=Extends(path_transcript_list[i], path_transcript_list[j], &type, 0, 0, 0);
                         included=Extends(path_transcript_list[i], path_transcript_list[j], &type, 0, 0, 0, 1);
#endif
#else
//19dic06
//included=Extends(path_transcript_list[i], path_transcript_list[j], &type, 0, 0, 1);
#ifdef DONT_EXTEND_REFSEQ
                         if(path_transcript_list[i].type == 1){
                                if(path_transcript_list[j].type == 0){
//08gen08
//included=Overlap(path_transcript_list[j], path_transcript_list[i], &type, 0, 0, 1);
                                  included=Overlap(path_transcript_list[j], path_transcript_list[i], &type, 0, 0, 1, 1);
                                }
                                else{
//08gen08
                                  included=Overlap(path_transcript_list[i], path_transcript_list[j], &type, 0, 0, 1, 1);
                                  included=0;
                                }
                         }
                         else{
                                if(path_transcript_list[j].type == 1){
//08gen08
//included=Overlap(path_transcript_list[i], path_transcript_list[j], &type, 0, 0, 1);
                                  included=Overlap(path_transcript_list[i], path_transcript_list[j], &type, 0, 0, 1, 1);
                                  if(included == 2)
                                         included=-2;
                                }
                                else{
//08gen08
//included=Extends(path_transcript_list[i], path_transcript_list[j], &type, 0, 0, 1);
                                  included=Extends(path_transcript_list[i], path_transcript_list[j], &type, 0, 0, 1, 1);
                                }
                         }
#else
//08gen08
//included=Extends(path_transcript_list[i], path_transcript_list[j], &type, 0, 0, 1);
                         included=Extends(path_transcript_list[i], path_transcript_list[j], &type, 0, 0, 1, 1);
#endif
#endif

//fprintf(stdout, "             Included %d\n", included);

//19dic06
                         if(included == 2 || included == -2){
//if((included == 2 && path_transcript_list[j].type == 0) || (included == -2 && path_transcript_list[i].type == 0)){
                                actual_path_number--;
                                if(included == 2){
//25lug07
//Modificare path_transcript_list[i] con allungo esoni solo se non e' refseq
#ifdef DONT_EXTEND_REFSEQ
                                  if(path_transcript_list[i].type != 1){
#endif
//23gen08
                                         if(list_of_exon_right[path_transcript_list[i].left_ext] == list_of_exon_right[path_transcript_list[j].left_ext]){
                                                if(type == 0){
//06nov07
//if(is_internal[path_transcript_list[i].left_ext == -1]){
                                                  if(is_internal[path_transcript_list[i].left_ext] == -1){
                                                         if(is_internal[path_transcript_list[j].left_ext] == 1){
                                                                path_transcript_list[i].left_ext=path_transcript_list[j].left_ext;
                                                         }
                                                         else{
                                                                if(is_internal[path_transcript_list[j].left_ext] == -1){
                                                                  if(list_of_exon_left[path_transcript_list[j].left_ext] < list_of_exon_left[path_transcript_list[i].left_ext]){
                                                                         path_transcript_list[i].left_ext=path_transcript_list[j].left_ext;
                                                                  }
                                                                }
                                                         }
                                                  }
                                                }
//23gen08
                                         }

//23gen08
                                         if(list_of_exon_left[path_transcript_list[i].right_ext] == list_of_exon_left[path_transcript_list[j].right_ext]){
                                                if(type+path_transcript_list[j].exons == path_transcript_list[i].exons){
//06nov07
//if(is_internal[path_transcript_list[i].right_ext == -2]){
                                                  if(is_internal[path_transcript_list[i].right_ext] == -2){
                                                         if(is_internal[path_transcript_list[j].right_ext] == 1){
                                                                path_transcript_list[i].right_ext=path_transcript_list[j].right_ext;
                                                         }
                                                         else{
                                                                if(is_internal[path_transcript_list[j].right_ext] == -2){
                                                                  if(list_of_exon_right[path_transcript_list[j].right_ext] > list_of_exon_right[path_transcript_list[i].right_ext]){
                                                                         path_transcript_list[i].right_ext=path_transcript_list[j].right_ext;
                                                                  }
                                                                }
                                                         }
                                                  }
                                                }
//23gen08
                                         }
#ifdef DONT_EXTEND_REFSEQ
                                  }
#endif

                                  filtered[j]=1;
//fprintf(stdout, "FILTj %d!\n", j);

//07ott05
#ifdef MERGE_POLYA
                                  if(polya[path_transcript_list[j].right_ext] == 1 && path_transcript_list[j].right_ext != path_transcript_list[i].right_ext){
                                         if(path_transcript_list[i].number_of_polya+1 > MAX_POLYA_END){
                                                fprintf(stderr, "Too many polya\n");
#ifdef HALT_EXIT_MODE
                                                exit(1);
#else
                                                exit(EXIT_FAILURE);
#endif
                                         }

                                         p=0;
                                         stop2=0;
                                         while(p<path_transcript_list[i].number_of_polya && !stop2){
                                                if(path_transcript_list[i].polya_end[p] == path_transcript_list[j].right_ext)
                                                  stop2=1;
                                                else
                                                  p++;
                                         }
                                         if(!stop2){
                                                path_transcript_list[i].polya_end[path_transcript_list[i].number_of_polya]=path_transcript_list[j].right_ext;
                                                path_transcript_list[i].number_of_polya+=1;
                                         }
                                  }
#endif

                                  Add_Path_List(&transcript_list_of_paths[i], &transcript_list_of_paths[j]);
                                }
                                else{
//25lug07
//Modificare path_transcript_list[j] con allungo esoni solo se non e' refseq
#ifdef DONT_EXTEND_REFSEQ
                                  if(path_transcript_list[j].type != 1){
#endif
//23gen08
                                         if(list_of_exon_right[path_transcript_list[i].left_ext] == list_of_exon_right[path_transcript_list[j].left_ext]){
                                                if(type == 0){
//06nov07
//if(is_internal[path_transcript_list[j].left_ext == -1]){
                                                  if(is_internal[path_transcript_list[j].left_ext] == -1){
                                                         if(is_internal[path_transcript_list[i].left_ext] == 1){
                                                                path_transcript_list[j].left_ext=path_transcript_list[i].left_ext;
                                                         }
                                                         else{
                                                                if(is_internal[path_transcript_list[i].left_ext] == -1){
                                                                  if(list_of_exon_left[path_transcript_list[i].left_ext] < list_of_exon_left[path_transcript_list[j].left_ext]){
                                                                         path_transcript_list[j].left_ext=path_transcript_list[i].left_ext;
                                                                  }
                                                                }
                                                         }
                                                  }
                                                }
//23gen08
                                         }

//23gen08
                                         if(list_of_exon_left[path_transcript_list[i].right_ext] == list_of_exon_left[path_transcript_list[j].right_ext]){
                                                if(type+path_transcript_list[i].exons == path_transcript_list[j].exons){
//06nov07
//if(is_internal[path_transcript_list[j].right_ext == -2]){
                                                  if(is_internal[path_transcript_list[j].right_ext] == -2){
                                                         if(is_internal[path_transcript_list[i].right_ext] == 1){
                                                                path_transcript_list[j].right_ext=path_transcript_list[i].right_ext;
                                                         }
                                                         else{
                                                                if(is_internal[path_transcript_list[i].right_ext] == -2){
                                                                  if(list_of_exon_right[path_transcript_list[i].right_ext] > list_of_exon_right[path_transcript_list[j].right_ext]){
                                                                         path_transcript_list[j].right_ext=path_transcript_list[i].right_ext;
                                                                  }
                                                                }
                                                         }
                                                  }
                                                }
//23gen08
                                         }
#ifdef DONT_EXTEND_REFSEQ
                                  }
#endif

//fprintf(stdout, "FILTi %d!\n", i);
                                  filtered[i]=1;

//07ott05
#ifdef MERGE_POLYA
                                  if(polya[path_transcript_list[i].right_ext] == 1 && path_transcript_list[i].right_ext != path_transcript_list[j].right_ext){
                                         if(path_transcript_list[j].number_of_polya+1 > MAX_POLYA_END){
                                                fprintf(stderr, "Too many polya\n");
#ifdef HALT_EXIT_MODE
                                                exit(1);
#else
                                                exit(EXIT_FAILURE);
#endif
                                         }

                                         p=0;
                                         stop2=0;
                                         while(p<path_transcript_list[j].number_of_polya && !stop2){
                                                if(path_transcript_list[j].polya_end[p] == path_transcript_list[i].right_ext)
                                                  stop2=1;
                                                else
                                                  p++;
                                         }

                                         if(!stop2){
                                                path_transcript_list[j].polya_end[path_transcript_list[j].number_of_polya]=path_transcript_list[i].right_ext;
                                                path_transcript_list[j].number_of_polya+=1;
                                         }
                                  }
#endif

                                  Add_Path_List(&transcript_list_of_paths[j], &transcript_list_of_paths[i]);
                                  stop=1;
                                }
                         }
                  }
                  j++;
                }
         }
         i++;
  }

  for(i=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; i<= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; i++){
         vect_actual_path_number[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=actual_path_number;
  }

  i=0;
  while(i<total_paths){
         if(!filtered[i]){
                for(k=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; k<= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; k++){
                  if(path_transcript_list[i].exons < k){
                         vect_actual_path_number[k-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]--;
                  }
                }
         }
         i++;
  }
}

void Filter_Path_Transcripts_by_Introns(){
  int i=0, j=0, k=0, q=0;
  FILE *in=NULL;
  char *tmp_string= NULL;
  size_t tmp_string_l = 0;
  int bytes_read;
//  char tmp_string[10000000];
  int left[500], right[500];
  float donor_err[500], acceptor_err[500];
  int EST_conf[500];
  char **EST_ids;
  char pt_5[500][5], pt_3[500][5];
  char pt[10];
  int donor_exon=0, accept_exon=0;
  int intron_start=0, intron_end=0;
  char found=0;
  int count_introns=0;
  char isRefSeq=0;

  in=fopen("predicted-introns.txt", "r");
  if(in == NULL){
         fprintf(stderr, "predicted-introns.txt not opened!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  EST_ids=(char **)malloc(500*sizeof(char *));
  for(i=0; i<500; i++){
         //EST_ids[i]=(char *)malloc(200000*sizeof(char));
         EST_ids[i]=(char *)malloc(1000000*sizeof(char));
  }

  while(!feof(in)){
	 my_assert(count_introns < 500);
	 bytes_read = my_getline (&tmp_string, &tmp_string_l, in);
	 if (bytes_read == -1) {
		/*fprintf(stderr, "Error reading from predicted-introns.txt!\n");
#ifdef HALT_EXIT_MODE
		exit(1);
#else
		exit(EXIT_FAILURE);
#endif*/
		 DEBUG("no introns from predicted-introns.txt!");

	 } else {
		DEBUG("line >%s<", tmp_string);
		sscanf(tmp_string,
				 "%d %d %*d %*d %*d %d %s %f %f %*f %*f %*f %*d %*d %s %*s %*s %*s %*s %*s",
				 &left[count_introns],
				 &right[count_introns],
				 (EST_conf+count_introns),
				 EST_ids[count_introns],
				 donor_err+count_introns,
				 acceptor_err+count_introns,
				 pt);

		DEBUG("Read intron %d-%d.", left[count_introns], right[count_introns]);
		pt[4]='\0';
		pt_5[count_introns][0]=pt[0];
		pt_5[count_introns][1]=pt[1];
		pt_5[count_introns][2]='\0';

		pt_3[count_introns][0]=pt[2];
		pt_3[count_introns][1]=pt[3];
		pt_3[count_introns][2]='\0';

		EST_ids[count_introns][strlen(EST_ids[count_introns])-1]='\0';
		count_introns++;
	 }
  }
  fclose(in);

  i=0;
  while(i<total_paths){
	 if(!filtered[i]){
		j=0;
		while(j < path_transcript_list[i].exons-1){
		  donor_exon=(j == 0)?(path_transcript_list[i].left_ext):(path_transcript_list[i].exon_list[j-1]);
		  accept_exon=(j == path_transcript_list[i].exons-2)?(path_transcript_list[i].right_ext):(path_transcript_list[i].exon_list[j]);
		  intron_start=list_of_exon_right[donor_exon]+1;
		  intron_end=list_of_exon_left[accept_exon]-1;
		  found=0;
		  k=0;
		  DEBUG("Looking for intron %d-%d.", intron_start, intron_end);
		  while(k < count_introns && !found){
			 DEBUG("  --> encountered (%d) %d-%d.", k, left[k], right[k]);
			 if (intron_start == left[k] && intron_end == right[k]) {
				found=1;
			 } else {
				k++;
			 }
		  }
		  if (found) {
			 if(EST_conf[k] < 2){
				isRefSeq=0;
				q=0;
				while (q < (int)strlen(EST_ids[k])-1 && isRefSeq == 0){
				  if(EST_ids[k][q] == 'N' && EST_ids[k][q+1]){
					 if(q == 0 || EST_ids[k][q-1] == ','){
						if(q < (int)strlen(EST_ids[k])-2 && EST_ids[k][q+2] == '_')
						  isRefSeq=1;
					 }
				  }
				  q++;
				}

				if(isRefSeq == 0){
				  if ((strcmp(To_lower(pt_5[k]), "gt") ||
						 strcmp(To_lower(pt_3[k]), "ag")
						 ) || (donor_err[k]+acceptor_err[k] > 10.00)) {
					 filtered[i]=1;
					 actual_path_number--;
				  }
				}
			 }
		  } else {
			 fprintf(stderr, "Intron not found!\n");
			 filtered[i]=1;
			 actual_path_number--;
/*
#ifdef HALT_EXIT_MODE
			 exit(1);
#else
			 exit(EXIT_FAILURE);
#endif
*/
		  }
		  j++;
		}
	 }
	 i++;
  }

  for(i=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; i<= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; i++){
         vect_actual_path_number[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=actual_path_number;
  }

  i=0;
  while(i<total_paths){
         if(!filtered[i]){
                for(k=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; k<= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; k++){
                  if(path_transcript_list[i].exons < k){
                         vect_actual_path_number[k-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]--;
                  }
                }
         }
         i++;
  }

  for(i=0; i<500; i++)
         free(EST_ids[i]);
  free(EST_ids);
}


/*void Filter_Path_Transcripts_by_Introns_OLD(){
  int i=0, j=0, k=0, q=0;
  //char included=0;
//  char stop=0, stop2=0;
//  int type=0;
  FILE *in=NULL;

//21feb06
  FILE *in2=NULL;

  char tmp_string[400000];

  char ok=0;
  int left[500], right[500];

//21feb06
  float donor_err[500], acceptor_err[500];

  char *pt=NULL;

  int EST_conf[500];

//03mar07
//char EST_ids[500][100000];
  char **EST_ids;

//21feb06
  char pt_5[500][5], pt_3[500][5];

  int donor_exon=0, accept_exon=0;
  int intron_start=0, intron_end=0;
  char found=0;
//  int count_newline=0;
  char *help_string1=NULL, *help_string2=NULL;
  int count_introns=0;

//19dic06
  char isRefSeq=0;

  in=fopen("OUT", "r");
  if(in == NULL){
         fprintf(stderr, "OUT not opened!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

//03mar07
  EST_ids=(char **)malloc(500*sizeof(char *));
  for(i=0; i<500; i++){
//29nov07
//EST_ids[i]=(char *)malloc(100000*sizeof(char));
         EST_ids[i]=(char *)malloc(200000*sizeof(char));
  }

  while(!feof(in)){
         fgets(tmp_string, 400000, in);

         tmp_string[strlen(tmp_string)-1]='\0';

         help_string1=Substring(tmp_string, 0, 7);
         if(!strcmp(help_string1, "Introns:")){
                ok=1;
         }

         help_string2=Substring(tmp_string, 0, 26);
         if(!strcmp(help_string2, "Introns for graphical view:"))
                ok=0;

         if(ok && strcmp(tmp_string, "") && strcmp(help_string1, "Introns:") && tmp_string[0] != '\t'){
                if((int)strlen(tmp_string) > 3){
                  sscanf(tmp_string, "%s %*d %d %d %d", EST_ids[count_introns], (left+count_introns), (right+count_introns), (EST_conf+count_introns));

//sscanf(tmp_string, "%*s %*d %d %d %d", (left+count_introns), (right+count_introns), (EST_conf+count_introns));
                  EST_ids[count_introns][strlen(EST_ids[count_introns])-1]='\0';

                  count_introns++;
                }
         }
         free(help_string1);
         free(help_string2);
  }
//21feb06
  fclose(in);

  if(count_introns > 0){
//21feb06
         in2=fopen("./intron_align.txt", "r");
         if(in2 == NULL){
                fprintf(stderr, "OUT not opened!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }
         count_introns=0;
         while(!feof(in2)){
                fgets(tmp_string, 400000, in2);
                tmp_string[strlen(tmp_string)-1]='\0';

                fscanf(in2, "%f\n", donor_err+count_introns);
                fscanf(in2, "%f\n", acceptor_err+count_introns);

                fgets(tmp_string, 400000, in2);
                fgets(tmp_string, 400000, in2);
                fgets(tmp_string, 400000, in2);

//25mag07
                fgets(tmp_string, 400000, in2);
                fgets(tmp_string, 400000, in2);
                fgets(tmp_string, 400000, in2);
                fgets(tmp_string, 400000, in2);
                fgets(tmp_string, 400000, in2);

//16mag08
                fgets(tmp_string, 400000, in2);

                fgets(tmp_string, 400000, in2);
                tmp_string[strlen(tmp_string)-1]='\0';

                pt=Substring(tmp_string, 0, 1);
                strcpy(pt_5[count_introns], pt);
                free(pt);

                fgets(tmp_string, 400000, in2);
                tmp_string[strlen(tmp_string)-1]='\0';

                pt=Substring(tmp_string, strlen(tmp_string)-2, strlen(tmp_string)-1);
                strcpy(pt_3[count_introns], pt);
                free(pt);

                count_introns++;

                do{
                  fgets(tmp_string, 400000, in2);
                }while(tmp_string[0] != '>');
         }
         fclose(in2);
  }

  i=0;
  while(i<total_paths){
         if(!filtered[i]){
                j=0;
//fprintf(stdout, "%d Transcript %d exons %d\n", total_paths, i, path_transcript_list[i].exons);
                while(j < path_transcript_list[i].exons-1){
                  donor_exon=(j == 0)?(path_transcript_list[i].left_ext):(path_transcript_list[i].exon_list[j-1]);
                  accept_exon=(j == path_transcript_list[i].exons-2)?(path_transcript_list[i].right_ext):(path_transcript_list[i].exon_list[j]);
//fprintf(stdout, "     Ex %d.%d\n", donor_exon, accept_exon);
                  intron_start=list_of_exon_right[donor_exon]+1;
                  intron_end=list_of_exon_left[accept_exon]-1;
//fprintf(stdout, "     Intron from %d to %d\n", intron_start, intron_end);
                  found=0;
                  k=0;
                  while(k < count_introns && !found){
//fprintf(stderr, "             Intron from %d(%d) to %d(%d)\n", left[k], intron_start, right[k], intron_end);
//fprintf(stderr, "             %d == %d TRUE/FALSE -> %d\n", intron_start, left[k], (intron_start == left[k]));
//fprintf(stderr, "             %d == %d TRUE/FALSE -> %d\n", intron_end, right[k], (intron_end == right[k]));
                         if(intron_start == left[k] && intron_end == right[k]){
//fprintf(stderr, "                     FOUND!\n");
                                found=1;
                         }
                         else
                                k++;
                  }
                  if(found){
//fprintf(stderr, "%s\n", EST_ids[k]);
                         if(EST_conf[k] < 2){

//19dic06
                                isRefSeq=0;
                                q=0;
                                while(q < (int)strlen(EST_ids[k])-1 && isRefSeq == 0){
                                  if(EST_ids[k][q] == 'N' && EST_ids[k][q+1]){
                                         if(q == 0 || EST_ids[k][q-1] == ','){
                                                if(q < (int)strlen(EST_ids[k])-2 && EST_ids[k][q+2] == '_')
                                                  isRefSeq=1;
                                         }
                                  }
                                  q++;
                                }

//19dic06
//if(isRefSeq)
//      fprintf(stderr, "Intron RefSeq %d-%d\n", intron_start, intron_end);

//19dic06
//if(EST_ids[k][0] != 'N' ||  EST_ids[k][1] != 'M'){
                                if(isRefSeq == 0){
//21feb06
//18ott06
//if((strcmp(To_lower(pt_5[k]), "gt") || strcmp(To_lower(pt_3[k]), "ag")) || (donor_err[k] != 0.00 || acceptor_err[k] != 0.00)){
                                  if((strcmp(To_lower(pt_5[k]), "gt") || strcmp(To_lower(pt_3[k]), "ag")) || (donor_err[k]+acceptor_err[k] > 10.00)){
                                         filtered[i]=1;
                                         actual_path_number--;
                                  }
                                }
                         }
                  }
                  else{
                         fprintf(stderr, "Intron not found!\n");
#ifdef HALT_EXIT_MODE
                         exit(1);
#else
                         exit(EXIT_FAILURE);
#endif
                  }
                  j++;
                }
         }
         i++;
  }

  for(i=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; i<= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; i++){
         vect_actual_path_number[i-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]=actual_path_number;
  }

  i=0;
  while(i<total_paths){
         if(!filtered[i]){
                for(k=FIRST_MIN_EXONS_ACCEPTED_OUTPUT; k<= SECOND_MIN_EXONS_ACCEPTED_OUTPUT; k++){
                  if(path_transcript_list[i].exons < k){
                         vect_actual_path_number[k-FIRST_MIN_EXONS_ACCEPTED_OUTPUT]--;
                  }
                }
         }
         i++;
  }

//03mar07
  for(i=0; i<500; i++)
         free(EST_ids[i]);
  free(EST_ids);
}*/

struct path *Create_Source_Path(int index){
  struct path *path=NULL;

  path=(struct path *)malloc(sizeof(struct path));
  if(path == NULL){
         fprintf(stderr, "Problem1 of memory allocation in Create_Source_Path!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  path->next=NULL;
  path->n=NULL;
  path->tail=NULL;

  Add_Node(&path, index, 1);

  return path;
}

char *Substring(char *string, int left, int right){
  int i=0;
  char *substring=NULL;

  if(left > right){
         left=1;
         right=0;
  }

  substring=(char *)malloc((right-left+2)*sizeof(char));
  if(substring==NULL){
         fprintf(stderr, "Memory problem in the substring procedure!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  for(i=left; i<=right+1; i++){
         if(i==right+1) {
                substring[i-left]='\0';
         } else {
                substring[i-left]=string[i];
         }
  }

  return substring;
}

/* int Edit_distance(char *seq1, char *seq2){ */
/*   int length1=strlen(seq1); */
/*   int length2=strlen(seq2); */
/*   char *temp_string=NULL; */
/*   int diff=0; */
/*   int k1=0, k2=0; */
/*   int **matrix=NULL; */
/*   int i=0, j=0, k=0, m=0; */
/*   int result=0; */

/* /\*If the two input sequences are the same sequence*\/ */
/*   if(!strcmp(seq1,seq2)){ */
/*       return 0; */
/*   } */

/*   if(length1 < length2){ */
/*       temp_string=seq1; */
/*       seq1=seq2; */
/*       seq2=temp_string; */
/*   } */
/*   length1=strlen(seq1); */
/*   length2=strlen(seq2); */

/*   matrix=(int **)malloc((length2 + 1)*sizeof(int*)); */
/*   if(matrix==NULL){ */
/*       fprintf(stderr, "Edit distance matrix not stored!\n"); */
/* #ifdef HALT_EXIT_MODE */
/*       exit(1); */
/* #else */
/*       exit(EXIT_FAILURE); */
/* #endif */
/*   } */

/*   for(i=0; i<=length2; i++){ */
/*       matrix[i]=(int*)malloc((length1 + 1)*sizeof(int)); */
/*       if(matrix[i]==NULL){ */
/*              fprintf(stderr, "Row %d of the edit distance matrix not stored!\n", (i+1)); */
/* #ifdef HALT_EXIT_MODE */
/*              exit(1); */
/* #else */
/*              exit(EXIT_FAILURE); */
/* #endif */
/*       } */
/*   } */

/*   for(i=0; i<=length2; i++){ */
/*       for(j=0; j<=length1; j++){ */
/*              matrix[i][j]=1; */
/*       } */
/*   } */

/*   for(i=0; i<=length1; i++){ */
/*       matrix[0][i]=i; */
/*   } */

/*   for(i=1; i<=length2; i++){ */
/*       matrix[i][0]=i; */
/*   } */

/*   for(i=1; i<=length2; i++){ */
/*       for(j=1; j<=length1; j++){ */
/*              matrix[i][j]=matrix[i-1][j-1]; */
/*              if(seq1[j-1] != seq2[i-1]) */
/*                matrix[i][j]+=1; */
/*              if(matrix[i][j] > (matrix[i-1][j]+1)) */
/*                matrix[i][j]=matrix[i-1][j]+1; */
/*              if(matrix[i][j] > (matrix[i][j-1]+1)) */
/*                matrix[i][j]=matrix[i][j-1]+1; */
/*       } */
/*   } */

/*   result=matrix[length2][length1]; */

/*   for(i=0; i<=length2; i++){ */
/*       free(matrix[i]); */
/*   } */
/*   free(matrix); */

/*   return result; */
/* } */

/*Return:
  - -1 se t1 estende t2
  - -2 se t1 e' incluso in t2
  - 1 se t2 estende t1
  - 2 se t2 e' incluso in t1
  - 0 altrimenti

  in L e' memorizzato l'indice (da 0- primo esone) dell'esone di t2 (t1 per valori di return > 0) da cui parte l'estensione o l'inlcusione di t1
  (t2 per valori di return > 0)
  only_check_name=1 se si vuole solo il confronto nominale degli esoni*/
//08giu05
//char Extends(struct transcript t1, struct transcript t2, int *L, char only_check_name){
//07ott05
//char Extends(struct transcript t1, struct transcript t2, int *L, char only_check_name, char for_ext){
//08gen08
//char Extends(struct transcript t1, struct transcript t2, int *L, char only_check_name, char for_ext, char force_polya){
char Extends(struct transcript t1, struct transcript t2, int *L, char only_check_name, char for_ext, char force_polya, char filt_phase){
  char overlap1=0;
  char overlap2=0;

//08giu05
//overlap1=Overlap(t1, t2, L, only_check_name);
//07ott05
//overlap1=Overlap(t1, t2, L, only_check_name, for_ext);
//08gen08
//overlap1=Overlap(t1, t2, L, only_check_name, for_ext, force_polya);
  overlap1=Overlap(t1, t2, L, only_check_name, for_ext, force_polya, filt_phase);

//fprintf(stdout, "Overlap1 %d\n", overlap1);

  if(overlap1 == 1){
         return -1;
  }

  if(overlap1 == 2){
         return -2;
  }

//08giu05
//overlap2=Overlap(t2, t1, L, only_check_name);
//07ott05
//overlap2=Overlap(t2, t1, L, only_check_name, for_ext);
//08gen08
//overlap2=Overlap(t2, t1, L, only_check_name, for_ext, force_polya);
  overlap2=Overlap(t2, t1, L, only_check_name, for_ext, force_polya, filt_phase);

//fprintf(stdout, "Overlap2 %d\n", overlap2);

  if(overlap2 == 1){
         return 1;
  }

  if(overlap2 == 2){
         return 2;
  }

  return 0;
}

/*Return:
  - 1 se t1 estende t2
  - 2 se t1 e' incluso in t2
  - 0 altrimenti

  in L e' memorizzato l'indice (da 0- primo esone) dell'esone di t2 da cui parte l'estensione o l'inlcusione di t1
  only_check_name=1 se si vuole solo il confronto nominale degli esoni*/
//08giu05
//char Overlap(struct transcript t1, struct transcript t2, int *L, char only_check_name){
//07ott05
//char Overlap(struct transcript t1, struct transcript t2, int *L, char only_check_name, char for_ext){
//08gen08
//char Overlap(struct transcript t1, struct transcript t2, int *L, char only_check_name, char for_ext, char force_polya){
char Overlap(struct transcript t1, struct transcript t2, int *L, char only_check_name, char for_ext, char force_polya, char filt_phase){
  char found=0, match=0, int_match=0;
  int k=0;
  int first_exon1=t1.left_ext, current_exon2=-1;
  int last_exon1=-1, last_exon2=-1;
  char stop=0;
  int l=0, j=0;
  char strength_l=0, strength_r=0;
//20feb06
  char *temp_string=NULL;

//10mag05
//08giu05
//07ott05
//if(polya[t2.right_ext] == 1 && for_ext)
  if(force_polya && (polya[t2.right_ext] == 1 && for_ext))
         return 0;

/*fprintf(stdout, "\nT1:");
  for(p=0; p < t1.exons; p++){
  if(p == 0)
  fprintf(stdout, "     %d", t1.left_ext);
  else{
  if(p == t1.exons-1)
  fprintf(stdout, ".%d", t1.right_ext);
  else
  fprintf(stdout, ".%d", t1.exon_list[p-1]);
  }
  }
  fprintf(stdout, "     T2:");
  for(p=0; p < t2.exons; p++){
  if(p == 0)
  fprintf(stdout, "     %d", t2.left_ext);
  else{
  if(p == t2.exons-1)
  fprintf(stdout, ".%d", t2.right_ext);
  else
  fprintf(stdout, ".%d", t2.exon_list[p-1]);
  }
  }*/

//Ricerca dell'esone di t2 che matcha con il primo esone in t1
  while(!found && k < t2.exons){
         current_exon2=(k == 0)?(t2.left_ext):((k == t2.exons-1)?(t2.right_ext):(t2.exon_list[k-1]));

         found=(first_exon1 == current_exon2)?(1):(0);

         if(found)
                strength_l=1;

//Controllo del match nominale
         if(!found && !only_check_name){
                found=Check_L_suffix(first_exon1, current_exon2, &strength_l);
         }

         if(!found)
                k++;
  }

/*if(found){
  fprintf(stdout, "     FOUND k %d exons2 %d\n", k, t2.exons);
  }*/

  if(found){
         if(t1.exons == 1){
                *L=k;
                return 2;       //Inclusion of t1 in t2
         }

         if(t2.exons == 1){
                return 0;       //Inclusion of t2 in t1
         }

         l=k+1;
         j=1;

//Aggancio tra primo esone di t1 e ultimo di t2
         if(l == t2.exons){
//21nov05
//31ott07
//#ifdef EXT_NOT_ON_LAST
#ifdef STRONG_FIRST_LAST_MATCH
                return 0;
#else
                if(strength_l == 1){
                  *L=k;
//18nov05
//31ott07
/*#ifdef UPDATE_EXON
  Update_exon(first_exon1, current_exon2);
  #endif*/
                  Update_exon(first_exon1, current_exon2);
                  return 1;     //Extension of t1
                }
                else{
                  return 0;     //Too weak match
                }
#endif
         }

         int_match=0;

         while((l < t2.exons-1 && j < t1.exons-1) && !stop){
                int_match=(t1.exon_list[j-1] == t2.exon_list[l-1])?(1):(0);
                if(!int_match && !only_check_name)
                  int_match=Check_exons(t1.exon_list[j-1], t2.exon_list[l-1]);

                if(int_match){
                  l++;
                  j++;
                }
                else
                  stop=1;
         }

         if(!stop){
//Matching sugli ultimi due esoni di t1 e t2
                if(l == t2.exons-1 && j == t1.exons-1){
                  last_exon1=t1.right_ext;
                  last_exon2=t2.right_ext;

                  match=(last_exon1 == last_exon2)?(1):(0);

                  if(match)
                         strength_r=1;

                  if(!match && !only_check_name){
                         match=Check_R_prefix(last_exon1, last_exon2, &strength_r);
                  }

                  if(match){
                         if(int_match || (strength_l == 1 && strength_r == 1)){
                                *L=k;

//20feb06
//08gen08
/*if(polya[last_exon1] == 1 && polya[last_exon2] == 1){
//13mar07
//if(list_of_exon_right[last_exon1] > list_of_exon_right[last_exon2]){
#ifdef DONT_EXTEND_REFSEQ
if(list_of_exon_right[last_exon1] > list_of_exon_right[last_exon2] && t2.type != 1){
#else
if(list_of_exon_right[last_exon1] > list_of_exon_right[last_exon2]){
#endif

temp_string=Substring(sequences[last_exon1], strlen(sequences[last_exon1])-list_of_exon_right[last_exon1]+list_of_exon_right[last_exon2], strlen(sequences[last_exon1])-1);
list_of_exon_right[last_exon2]=list_of_exon_right[last_exon1];
strcat(sequences[last_exon2], temp_string);
free(temp_string);
}
}*/
                                if(filt_phase){
                                  if(polya[last_exon1] == 1 || polya[last_exon2] == 1){
//15gen08
//if(list_of_exon_right[last_exon1] > list_of_exon_right[last_exon2] && t2.type != 1){
                                         if(!(is_internal[last_exon2] == 1 && polya[last_exon2] == 0) && list_of_exon_right[last_exon1] > list_of_exon_right[last_exon2] && t2.type != 1){
//16gen08
//temp_string=Substring(sequences[last_exon1], strlen(sequences[last_exon1])-list_of_exon_right[last_exon1]+list_of_exon_right[last_exon2], strlen(sequences[last_exon1])-1);
                                                temp_string=Substring(sequences[last_exon1], strlen(sequences[last_exon1])-list_of_exon_right[last_exon1]+list_of_exon_right[last_exon2]+(list_of_exon_left[last_exon1]-list_of_exon_left[last_exon2]), strlen(sequences[last_exon1])-1);
                                                list_of_exon_right[last_exon2]=list_of_exon_right[last_exon1];
                                                strcat(sequences[last_exon2], temp_string);
                                                free(temp_string);
//if(t2.type == 1)
//Refseq extended!
                                         }
//15gen08
//polya[last_exon2]=1;
                                         if(!(is_internal[last_exon2] == 1 && polya[last_exon2] == 0))
                                                polya[last_exon2]=1;
                                  }
                                  else{
                                         if(is_internal[last_exon2] != 1 && t2.type != 1){
                                                if(list_of_exon_right[last_exon1] > list_of_exon_right[last_exon2] && list_of_exon_right[last_exon1]-list_of_exon_right[last_exon2] <= 50){
//16gen08
//temp_string=Substring(sequences[last_exon1], strlen(sequences[last_exon1])-list_of_exon_right[last_exon1]+list_of_exon_right[last_exon2], strlen(sequences[last_exon1])-1);
                                                  temp_string=Substring(sequences[last_exon1], strlen(sequences[last_exon1])-list_of_exon_right[last_exon1]+list_of_exon_right[last_exon2]+(list_of_exon_left[last_exon1]-list_of_exon_left[last_exon2]), strlen(sequences[last_exon1])-1);
                                                  list_of_exon_right[last_exon2]=list_of_exon_right[last_exon1];
                                                  strcat(sequences[last_exon2], temp_string);
                                                  free(temp_string);
//if(t2.type == 1)
//Refseq extended!
                                                }
                                         }
                                  }
//Il primo esone di t1 matcha con il primo di t2
                                  if(k == 0){
                                         if(is_internal[current_exon2] != 1 && t2.type != 1){
                                                if(list_of_exon_left[first_exon1] < list_of_exon_left[current_exon2] && list_of_exon_left[current_exon2]-list_of_exon_left[first_exon1] <= 50){
//15gen08
//temp_string=Substring(sequences[first_exon1], 0, list_of_exon_left[current_exon2]+list_of_exon_left[first_exon1]-1);
                                                  temp_string=Substring(sequences[first_exon1], 0, list_of_exon_left[current_exon2]-list_of_exon_left[first_exon1]-1);
                                                  list_of_exon_left[current_exon2]=list_of_exon_left[first_exon1];
                                                  strcat(sequences[current_exon2], temp_string);
                                                  free(temp_string);
//if(t2.type == 1)
//Refseq extended!
                                                }
                                         }
                                  }
                                }

//18nov05
#ifdef UPDATE_EXON
                                Update_exon(first_exon1, current_exon2);
                                Update_exon(last_exon1, last_exon2);
#endif

                                return 2; //Inclusion of t1 in t2
                         }
                         else{
                                return 0;
                         }
                  }
                  else
                         return 0;
                }
                else{
                  if(l == t2.exons-1){
                         last_exon1=t1.exon_list[j-1];
                         last_exon2=t2.right_ext;

                         match=(last_exon1 == last_exon2)?(1):(0);

                         if(match)
                                strength_r=1;

                         if(!match && !only_check_name){
//02set04
                                match=Check_R_prefix(last_exon1, last_exon2, &strength_r);
                         }

                         if(match){
                                if(k == 0){
                                  return 0; //Inclusion of t2 in t1
                                }
                                else{
                                  if(int_match || (strength_l == 1 && strength_r == 1)){
                                         *L=k;
//18nov05
#ifdef UPDATE_EXON
                                         Update_exon(first_exon1, current_exon2);
                                         Update_exon(last_exon1, last_exon2);
#endif
                                         return 1;      //Extension of t1
                                  }
                                  else
                                         return 0;
                                }
                         }
                         else
                                return 0;
                  }

                  if(j == t1.exons-1){
                         last_exon1=t1.right_ext;
                         last_exon2=t2.exon_list[l-1];

                         match=(last_exon1 == last_exon2)?(1):(0);

                         if(match)
                                strength_r=1;

                         if(!match && !only_check_name){
                                match=Check_R_prefix(last_exon1, last_exon2, &strength_r);
                         }

                         if(match){
//08giu05
//if(int_match || (strength_l == 1 && strength_r == 1)){
//07ott05
//if((polya[last_exon1] == 0) && (int_match || (strength_l == 1 && strength_r == 1))){
                                if((polya[last_exon1] == 0 || !force_polya) && (int_match || (strength_l == 1 && strength_r == 1))){
                                  *L=k;

//08gen08
//Il primo esone di t1 matcha con il primo di t2
                                  if(filt_phase && k == 0){
                                         if(is_internal[current_exon2] != 1 && t2.type != 1){
                                                if(list_of_exon_left[first_exon1] < list_of_exon_left[current_exon2] && list_of_exon_left[current_exon2]-list_of_exon_left[first_exon1] <= 50){
//15gen08
//temp_string=Substring(sequences[first_exon1], 0, list_of_exon_left[current_exon2]+list_of_exon_left[first_exon1]-1);
                                                  temp_string=Substring(sequences[first_exon1], 0, list_of_exon_left[current_exon2]-list_of_exon_left[first_exon1]-1);
                                                  list_of_exon_left[current_exon2]=list_of_exon_left[first_exon1];
                                                  strcat(sequences[current_exon2], temp_string);
                                                  free(temp_string);
//if(t2.type == 1)
//Refseq extended!
                                                }
                                         }
                                  }

//18nov05
#ifdef UPDATE_EXON
                                  Update_exon(first_exon1, current_exon2);
                                  Update_exon(last_exon1, last_exon2);
#endif
                                  return 2; //Inclusion
                                }
                                else
                                  return 0;
                         }
                         else{
                                return 0;
                         }
                  }
                }
         }
         else{
                return 0;
         }
  }
  else
         return 0;
// by Yuri: Aggiunte le prossime istruzioni perche' il programma
// originale non prevedeva questo flusso di esecuzione
  fprintf(stderr, "ERROR - An impossible thing is happened!\n");
  exit(EXIT_FAILURE);
  return 0; // to make the compiler happy
}

//Determina se e' possibile ridurre exon1 a exon2
char Check_L_suffix(int exon1, int exon2, char *matching_strength){
  int right_gap=0;
  int left_gap=0;
//03set04
  //char *temp_string=NULL;
//  char help_string[MAX_NLD];
  int ref_length=0;

//02set04
//Se uno dei due e' esterno dx
  if(is_internal[exon1] == -2){
         fprintf(stderr, "Problem in Check_L_suffix!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  *matching_strength=1;

  right_gap=list_of_exon_right[exon2]-list_of_exon_right[exon1];

//02set04
/*if(right_gap > 2 || right_gap < -2)
  return 0;*/

  left_gap=list_of_exon_left[exon2]-list_of_exon_left[exon1];

//02set04
//Se entrambi sono interni
  if(is_internal[exon1] == 1 && is_internal[exon2] == 1){

         if(right_gap > 2 || right_gap < -2)
                return 0;

         if(left_gap > 2 || left_gap < -2)
                return 0;

         return 1;
  }
//Se exon2 e' interno e exon1 e' ext sx
  if(is_internal[exon2] == 1){
         if(right_gap > 2 || right_gap < -2)
                return 0;

         if(left_gap > MAX_DIFF_FOR_REDUCING)
                return 0;

//MOD12
         ref_length=list_of_exon_right[exon2]-list_of_exon_left[exon2]+1;
//20set07
//if(list_of_exon_right[exon1]-list_of_exon_left[exon1]+1 < MIN_DIM_FOR_STRENGTH(ref_length))
         if(list_of_exon_right[exon1]-list_of_exon_left[exon1]+1 < MIN_DIM_FOR_STRENGTH2(ref_length))
//if(left_gap < -MAX_DIFF_FOR_STRENGTH)
//21set05
//*matching_strength=0;
                return 0;

/*list_of_exon_left[exon1]=list_of_exon_left[exon2];
//UPD SEQ (OK)
if(right_gap >= 0){
temp_string=Substring(sequences[exon2], 0, strlen(sequences[exon2])-right_gap-1);
strcpy(sequences[exon1], temp_string);
}
else{
temp_string=Substring(sequences[exon2], 0, strlen(sequences[exon2])-1);
strcpy(help_string, temp_string);
free(temp_string);
temp_string=Substring(sequences[exon1], strlen(sequences[exon1])+right_gap, strlen(sequences[exon1])-1);
strcat(help_string, temp_string);
strcpy(sequences[exon1], help_string);
}
free(temp_string);*/

         return 1;
  }
//Se solo exon1 e' interno
  if(is_internal[exon1] == 1){
         if(is_internal[exon2] == -1){

                if(right_gap > 2 || right_gap < -2)
                  return 0;

//07lug07
//if(left_gap < -MAX_DIFF_FOR_REDUCING)
                if(left_gap < -MAX_DIFF_FOR_REDUCING || left_gap > MAX_DIFF_FOR_REDUCING)
                  return 0;

//MOD12
                ref_length=list_of_exon_right[exon1]-list_of_exon_left[exon1]+1;
                if(list_of_exon_right[exon2]-list_of_exon_left[exon2]+1 < MIN_DIM_FOR_STRENGTH(ref_length))
//if(left_gap > MAX_DIFF_FOR_STRENGTH)
//21set05
//*matching_strength=0;
                  return 0;

/*list_of_exon_left[exon2]=list_of_exon_left[exon1];
//UPD SEQ (OK)
if(right_gap <= 0){
temp_string=Substring(sequences[exon1], 0, strlen(sequences[exon1])+right_gap-1);
strcpy(sequences[exon2], temp_string);
}
else{
temp_string=Substring(sequences[exon1], 0, strlen(sequences[exon1])-1);
strcpy(help_string, temp_string);
free(temp_string);
temp_string=Substring(sequences[exon2], strlen(sequences[exon2])-right_gap, strlen(sequences[exon2])-1);
strcat(help_string, temp_string);
strcpy(sequences[exon2], help_string);
}
free(temp_string);*/

                return 1;
         }
         else{
                if(left_gap > 2 || left_gap < -2)
                  return 0;

//07lug07
//if(right_gap > MAX_DIFF_FOR_REDUCING)
                if(right_gap > MAX_DIFF_FOR_REDUCING || right_gap < -MAX_DIFF_FOR_REDUCING)
                  return 0;

//MOD12
                ref_length=list_of_exon_right[exon1]-list_of_exon_left[exon1]+1;
                if(list_of_exon_right[exon2]-list_of_exon_left[exon2]+1 < MIN_DIM_FOR_STRENGTH(ref_length))
//if(right_gap < -MAX_DIFF_FOR_STRENGTH)
//21set05
//*matching_strength=0;
                  return 0;

/*list_of_exon_right[exon2]=list_of_exon_right[exon1];
//UPD SEQ (OK)
if(left_gap >= 0){
temp_string=Substring(sequences[exon1], left_gap, strlen(sequences[exon1])-1);
strcpy(sequences[exon2], temp_string);
}
else{
temp_string=Substring(sequences[exon2], 0, -(left_gap+1));
strcpy(help_string, temp_string);
free(temp_string);
temp_string=Substring(sequences[exon1], 0, strlen(sequences[exon1])-1);
strcat(help_string, temp_string);
strcpy(sequences[exon2], help_string);
}
free(temp_string);*/

                return 1;
         }
  }
//Sono entrambi esterni sx
  if(is_internal[exon2] == -1){
         if(right_gap > 2 || right_gap < -2)
                return 0;

         if(list_of_exon_left[exon2] < list_of_exon_left[exon1]){
//MOD12
                ref_length=list_of_exon_right[exon2]-list_of_exon_left[exon2]+1;
                if(list_of_exon_right[exon1]-list_of_exon_left[exon1]+1 < MIN_DIM_FOR_STRENGTH(ref_length))
//if(list_of_exon_left[exon1]-list_of_exon_left[exon2] > MAX_DIFF_FOR_STRENGTH)
                  *matching_strength=0;

//list_of_exon_left[exon1]=list_of_exon_left[exon2];
//UPD SEQ (OK)
//if(right_gap >= 0){
//      temp_string=Substring(sequences[exon2], 0, strlen(sequences[exon2])-right_gap-1);
//      strcpy(sequences[exon1], temp_string);
//}
//else{
//      temp_string=Substring(sequences[exon2], 0, strlen(sequences[exon2])-1);
//      strcpy(help_string, temp_string);
//      free(temp_string);
//      temp_string=Substring(sequences[exon1], strlen(sequences[exon1])+right_gap, strlen(sequences[exon1])-1);
//      strcat(help_string, temp_string);
//      strcpy(sequences[exon1], help_string);
//}
//free(temp_string);
         }
         else{
//MOD12
                ref_length=list_of_exon_right[exon1]-list_of_exon_left[exon1]+1;
                if(list_of_exon_right[exon2]-list_of_exon_left[exon2]+1 < MIN_DIM_FOR_STRENGTH(ref_length))
//if(list_of_exon_left[exon2]-list_of_exon_left[exon1] > MAX_DIFF_FOR_STRENGTH)
                  *matching_strength=0;

//list_of_exon_left[exon2]=list_of_exon_left[exon1];
//UPD SEQ (OK)
//if(right_gap <= 0){
//      temp_string=Substring(sequences[exon1], 0, strlen(sequences[exon1])+right_gap-1);
//      strcpy(sequences[exon2], temp_string);
//}
//else{
//      temp_string=Substring(sequences[exon1], 0, strlen(sequences[exon1])-1);
//      strcpy(help_string, temp_string);
//      free(temp_string);
//      temp_string=Substring(sequences[exon2], strlen(sequences[exon2])-right_gap, strlen(sequences[exon2])-1);
//      strcat(help_string, temp_string);
//      strcpy(sequences[exon2], help_string);
//}
//free(temp_string);
         }
  }
//Exon2 e' ext dx
  else{
//fprintf(stderr, "%d-%d %d-%d\n", list_of_exon_left[exon1], list_of_exon_right[exon1], list_of_exon_left[exon2], list_of_exon_right[exon2]);
//31ott07 Modifica per sperimentazione trascritti
#ifdef STRONG_FIRST_LAST_MATCH
//PROVA
//if(left_gap > MAX_DIFF_FOR_REDUCING || left_gap < -MAX_DIFF_FOR_REDUCING)
         if(left_gap > 2 || left_gap < -2)
                return 0;

//PROVA
//if(right_gap > MAX_DIFF_FOR_REDUCING || right_gap < -MAX_DIFF_FOR_REDUCING)
         if(right_gap > 2 || right_gap < -2)
                return 0;
//31ott07
#else
         if(left_gap > 2)
                return 0;

         if(right_gap > 2)
                return 0;
#endif

//21set01
//14mar07
//if(list_of_exon_right[exon2]-list_of_exon_left[exon1])
//      return 0;
         if(list_of_exon_left[exon2] < list_of_exon_left[exon1]){
                if(list_of_exon_right[exon2] < list_of_exon_right[exon1]){
                  if(list_of_exon_right[exon2]-list_of_exon_left[exon1]+1 < MIN_DIM_FOR_STRENGTH(ref_length))
                         return 0;
                }
                else{
                  if(list_of_exon_right[exon1]-list_of_exon_left[exon1]+1 < MIN_DIM_FOR_STRENGTH(ref_length))
                         return 0;
                }
         }
         else{
                if(list_of_exon_right[exon2] < list_of_exon_right[exon1]){
                  if(list_of_exon_right[exon2]-list_of_exon_left[exon2]+1 < MIN_DIM_FOR_STRENGTH(ref_length))
                         return 0;
                }
                else{
                  if(list_of_exon_right[exon1]-list_of_exon_left[exon2]+1 < MIN_DIM_FOR_STRENGTH(ref_length))
                         return 0;
                }
         }

//fprintf(stderr, "     %d-%d %d-%d\n", list_of_exon_left[exon1], list_of_exon_right[exon1], list_of_exon_left[exon2], list_of_exon_right[exon2]);

//31ott07 Modifica per sperimentazione trascritti
//*matching_strength=0;
#ifdef STRONG_FIRST_LAST_MATCH
         *matching_strength=0;
#else
         *matching_strength=1;
#endif

/*list_of_exon_left[exon1]=list_of_exon_left[exon2];
  list_of_exon_right[exon2]=list_of_exon_right[exon1];
//UPD SEQ (OK)
if(right_gap >= 0){
temp_string=Substring(sequences[exon2], 0, strlen(sequences[exon2])-right_gap-1);
strcpy(sequences[exon1], temp_string);
strcpy(sequences[exon2], temp_string);
}
else{
temp_string=Substring(sequences[exon2], 0, strlen(sequences[exon2])-1);
strcpy(help_string, temp_string);
free(temp_string);
temp_string=Substring(sequences[exon1], strlen(sequences[exon1])+right_gap, strlen(sequences[exon1])-1);
strcat(help_string, temp_string);
strcpy(sequences[exon1], help_string);
strcpy(sequences[exon2], help_string);
}
free(temp_string);*/
  }

  return 1;
}

//02set04
char Check_R_prefix(int exon1, int exon2, char *matching_strength){
  int right_gap=0;
  int left_gap=0;
//20feb06
  int threshold=0;

  //char *temp_string=NULL;
//  char help_string[MAX_NLD];
  int ref_length=0;

/*fprintf(stdout, "Exon1 %d Exon2 %d\n", exon1, exon2);
  fprintf(stdout, "Ex1 %d Ex2 %d\n", is_internal[exon1], is_internal[exon2]);*/

//02set04
//Se uno dei due e' esterno sx
  if(is_internal[exon1] == -1 || is_internal[exon2] == -1){
         fprintf(stderr, "Problem in Check_R_prefix!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  *matching_strength=1;

  left_gap=list_of_exon_left[exon2]-list_of_exon_left[exon1];

  if(left_gap > 2 || left_gap < -2)
         return 0;

  right_gap=list_of_exon_right[exon2]-list_of_exon_right[exon1];

//02set04
//Se entrambi sono interni
  if(is_internal[exon1] == 1 && is_internal[exon2] == 1){

//20feb06
         if(polya[exon1] == 1 && polya[exon2])
                threshold=MIN_POLYA_DIFF;
         else
                threshold=2;
         if(right_gap > threshold || right_gap < -threshold){
                return 0;
         }

         return 1;
  }
//Se solo exon2 e' interno
  if(is_internal[exon2] == 1){
         if(right_gap < -MAX_DIFF_FOR_REDUCING)
                return 0;

//MOD12
         ref_length=list_of_exon_right[exon2]-list_of_exon_left[exon2]+1;
         if(list_of_exon_right[exon1]-list_of_exon_left[exon1]+1 < MIN_DIM_FOR_STRENGTH(ref_length))
//if(right_gap > MAX_DIFF_FOR_STRENGTH)
//21set05
//*matching_strength=0;
                return 0;

/*list_of_exon_right[exon1]=list_of_exon_right[exon2];
//UPD SEQ (OK)
if(left_gap <= 0){
temp_string=Substring(sequences[exon2], -left_gap, strlen(sequences[exon2])-1);
strcpy(sequences[exon1], temp_string);
}
else{
temp_string=Substring(sequences[exon1], 0, left_gap-1);
strcpy(help_string, temp_string);
free(temp_string);
temp_string=Substring(sequences[exon2], 0, strlen(sequences[exon2])-1);
strcat(help_string, temp_string);
strcpy(sequences[exon1], help_string);
}
free(temp_string);*/

         return 1;
  }
//Se solo exon1 e' interno
  if(is_internal[exon1] == 1){
//07lug07
//if(right_gap > MAX_DIFF_FOR_REDUCING)
         if(right_gap > MAX_DIFF_FOR_REDUCING || right_gap < -MAX_DIFF_FOR_REDUCING)
                return 0;

//MOD12
         ref_length=list_of_exon_right[exon1]-list_of_exon_left[exon1]+1;
         if(list_of_exon_right[exon2]-list_of_exon_left[exon2]+1 < MIN_DIM_FOR_STRENGTH(ref_length))
//if(right_gap < -MAX_DIFF_FOR_STRENGTH)
//21set05
//*matching_strength=0;
                return 0;

/*list_of_exon_right[exon2]=list_of_exon_right[exon1];
//UPD SEQ (OK)
if(left_gap >= 0){
temp_string=Substring(sequences[exon1], left_gap, strlen(sequences[exon1])-1);
strcpy(sequences[exon2], temp_string);
}
else{
temp_string=Substring(sequences[exon2], 0, -(left_gap+1));
strcpy(help_string, temp_string);
free(temp_string);
temp_string=Substring(sequences[exon1], 0, strlen(sequences[exon1])-1);
strcat(help_string, temp_string);
strcpy(sequences[exon2], help_string);
}
free(temp_string);*/

         return 1;
  }

//Sono entrambi esterni dx
  if(list_of_exon_right[exon2] > list_of_exon_right[exon1]){
//MOD12
         ref_length=list_of_exon_right[exon2]-list_of_exon_left[exon2]+1;
         if(list_of_exon_right[exon1]-list_of_exon_left[exon1]+1 < MIN_DIM_FOR_STRENGTH(ref_length))
//if(list_of_exon_right[exon2]-list_of_exon_right[exon1] > MAX_DIFF_FOR_STRENGTH)
                *matching_strength=0;

//list_of_exon_right[exon1]=list_of_exon_right[exon2];
//UPD SEQ (OK)
//if(left_gap <= 0){
//      temp_string=Substring(sequences[exon2], -left_gap, strlen(sequences[exon2])-1);
//      strcpy(sequences[exon1], temp_string);
//}
//else{
//      temp_string=Substring(sequences[exon1], 0, left_gap-1);
//      strcpy(help_string, temp_string);
//      free(temp_string);
//      temp_string=Substring(sequences[exon2], 0, strlen(sequences[exon2])-1);
//      strcat(help_string, temp_string);
//      strcpy(sequences[exon1], help_string);
//}
//free(temp_string);
  }
  else{
//MOD12
         ref_length=list_of_exon_right[exon1]-list_of_exon_left[exon1]+1;
         if(list_of_exon_right[exon2]-list_of_exon_left[exon2]+1 < MIN_DIM_FOR_STRENGTH(ref_length))
//if(list_of_exon_right[exon1]-list_of_exon_right[exon2] > MAX_DIFF_FOR_STRENGTH)
                *matching_strength=0;

//list_of_exon_right[exon2]=list_of_exon_right[exon1];
//UPD SEQ (OK)
//if(left_gap >= 0){
//      temp_string=Substring(sequences[exon1], left_gap, strlen(sequences[exon1])-1);
//      strcpy(sequences[exon2], temp_string);
//}
//else{
//      temp_string=Substring(sequences[exon2], 0, -(left_gap+1));
//      strcpy(help_string, temp_string);
//      free(temp_string);
//      temp_string=Substring(sequences[exon1], 0, strlen(sequences[exon1])-1);
//      strcat(help_string, temp_string);
//      strcpy(sequences[exon2], help_string);
//}
//free(temp_string);
  }

  return 1;
}

char Check_exons(int exon1, int exon2){
  int right_gap=0;
  int left_gap=0;

  left_gap=list_of_exon_left[exon2]-list_of_exon_left[exon1];

  if(left_gap > 2 || left_gap < -2)
         return 0;

  right_gap=list_of_exon_right[exon2]-list_of_exon_right[exon1];

  if(right_gap > 2 || right_gap < -2)
         return 0;

  return 1;
}

//01set04
/*Procedura di aumento (diminuzione) del left (right) end di un esone esterno sinistro (destro) per portarlo
  a coincidere con un taglio confermato se la differenza rientra in certi limiti*/
/* void Reduce_external_exons(){ */
/*   int i=0, p=0; */
/*   char stop=0; */

/*   char *temp_string=NULL; */
/*   char help_string[MAX_NLD]; */

/*   for(i=0; i<number_of_exons; i++){ */
/*       stop=0; */
/* //esterno sx */
/*       if(is_internal[i] == -1){ */
/*              p=i-1; */
/*              while(p >= 0 && list_of_exon_right[p] > list_of_exon_left[i] && !stop){ */
/*                if(list_of_exon_left[p] > list_of_exon_left[i] && list_of_exon_left[p]-list_of_exon_left[i] <= 2) */
/*                       if(is_internal[p] == 1) */
/*                              stop=1; */
/*                if(list_of_exon_left[p] < list_of_exon_left[i] && list_of_exon_left[i]-list_of_exon_left[p] <= 2) */
/*                       if(is_internal[p] == 1) */
/*                              stop=1; */
/*                if(stop){ */
/* //UPD SEQ */
/*                       if(list_of_exon_left[p]-list_of_exon_left[i] >= 0){ */
/*                              temp_string=Substring(sequences[i], list_of_exon_left[p]-list_of_exon_left[i], strlen(sequences[i])-1); */
/*                              strcpy(sequences[i], temp_string); */
/*                       } */
/*                       else{ */
/*                              temp_string=Substring(sequences[p], 0, list_of_exon_left[i]-list_of_exon_left[p]-1); */
/*                              strcpy(help_string, temp_string); */
/*                              free(temp_string); */
/*                              temp_string=Substring(sequences[i], 0, strlen(sequences[i])-1); */
/*                              strcat(help_string, temp_string); */
/*                              strcpy(sequences[i], help_string); */
/*                       } */
/*                       free(temp_string); */
/*                       list_of_exon_left[i]=list_of_exon_left[p]; */
/*                } */
/*                else */
/*                       p--; */
/*              } */
/*              p=i+1; */
/*              while(p < number_of_exons && list_of_exon_left[p] < list_of_exon_right[i] && !stop){ */
/*                if(list_of_exon_left[p] > list_of_exon_left[i] && list_of_exon_left[p]-list_of_exon_left[i] <= 2) */
/*                       if(is_internal[p] == 1) */
/*                              stop=1; */
/*                if(list_of_exon_left[p] < list_of_exon_left[i] && list_of_exon_left[i]-list_of_exon_left[p] <= 2) */
/*                       if(is_internal[p] == 1) */
/*                              stop=1; */
/*                if(stop){ */
/* //UPD SEQ */
/*                       if(list_of_exon_left[p]-list_of_exon_left[i] >= 0){ */
/*                              temp_string=Substring(sequences[i], list_of_exon_left[p]-list_of_exon_left[i], strlen(sequences[i])-1); */
/*                              strcpy(sequences[i], temp_string); */
/*                       } */
/*                       else{ */
/*                              temp_string=Substring(sequences[p], 0, list_of_exon_left[i]-list_of_exon_left[p]-1); */
/*                              strcpy(help_string, temp_string); */
/*                              free(temp_string); */
/*                              temp_string=Substring(sequences[i], 0, strlen(sequences[i])-1); */
/*                              strcat(help_string, temp_string); */
/*                              strcpy(sequences[i], help_string); */
/*                       } */
/*                       free(temp_string); */
/*                       list_of_exon_left[i]=list_of_exon_left[p]; */
/*                } */
/*                else */
/*                       p++; */
/*              } */
/*       } */

/* //esterno dx */
/*       if(is_internal[i] == -2){ */
/*              p=i-1; */
/*              while(p >= 0 && list_of_exon_right[p] > list_of_exon_left[i] && !stop){ */
/*                if(list_of_exon_right[p] > list_of_exon_right[i] && list_of_exon_right[p]-list_of_exon_right[i] <= 2) */
/*                       if(is_internal[p] == 1) */
/*                              stop=1; */
/*                if(list_of_exon_right[p] < list_of_exon_right[i] && list_of_exon_right[i]-list_of_exon_right[p] <= 2) */
/*                       if(is_internal[p] == 1) */
/*                              stop=1; */
/*                if(stop){ */
/* //UPD SEQ */
/*                       if(list_of_exon_right[i]-list_of_exon_right[p] >= 0){ */
/*                              temp_string=Substring(sequences[i], 0, strlen(sequences[i])+list_of_exon_right[p]-list_of_exon_right[i]-1); */
/*                              strcpy(sequences[i], temp_string); */
/*                       } */
/*                       else{ */
/*                              temp_string=Substring(sequences[i], 0, strlen(sequences[i])-1); */
/*                              strcpy(help_string, temp_string); */
/*                              free(temp_string); */
/*                              temp_string=Substring(sequences[p], strlen(sequences[p])+list_of_exon_right[i]-list_of_exon_right[p], strlen(sequences[p])-1); */
/*                              strcat(help_string, temp_string); */
/*                              strcpy(sequences[i], help_string); */
/*                       } */
/*                       free(temp_string); */
/*                       list_of_exon_right[i]=list_of_exon_right[p]; */
/*                } */
/*                else */
/*                       p--; */
/*              } */
/*              p=i+1; */
/*              while(p < number_of_exons && list_of_exon_left[p] < list_of_exon_right[i] && !stop){ */
/*                if(list_of_exon_right[p] > list_of_exon_right[i] && list_of_exon_right[p]-list_of_exon_right[i] <= 2) */
/*                       if(is_internal[p] == 1) */
/*                              stop=1; */
/*                if(list_of_exon_right[p] < list_of_exon_right[i] && list_of_exon_right[i]-list_of_exon_right[p] <= 2) */
/*                       if(is_internal[p] == 1) */
/*                              stop=1; */
/*                if(stop){ */
/* //UPD SEQ */
/*                       if(list_of_exon_right[i]-list_of_exon_right[p] >= 0){ */
/*                              temp_string=Substring(sequences[i], 0, strlen(sequences[i])+list_of_exon_right[p]-list_of_exon_right[i]-1); */
/*                              strcpy(sequences[i], temp_string); */
/*                       } */
/*                       else{ */
/*                              temp_string=Substring(sequences[i], 0, strlen(sequences[i])-1); */
/*                              strcpy(help_string, temp_string); */
/*                              free(temp_string); */
/*                              temp_string=Substring(sequences[p], strlen(sequences[p])+list_of_exon_right[i]-list_of_exon_right[p], strlen(sequences[p])-1); */
/*                              strcat(help_string, temp_string); */
/*                              strcpy(sequences[i], help_string); */
/*                       } */
/*                       free(temp_string); */
/*                       list_of_exon_right[i]=list_of_exon_right[p]; */
/*                } */
/*                else */
/*                       p++; */
/*              } */
/*       } */
/*   } */
/* } */

/* void Update_exon(int exon1, int exon2){ */
/*   int right_gap=0, left_gap=0; */
/*   char *temp_string=NULL; */
/*   char help_string[MAX_NLD]; */

/*   if(is_internal[exon1] == 1 && is_internal[exon2] == 1) */
/*       return; */

/*   if(is_internal[exon2] == 1){ */
/*       if(is_internal[exon1] == -1){ */
/* //20dic06 */
/*              is_internal[exon1]=1; */

/*              list_of_exon_left[exon1]=list_of_exon_left[exon2]; */
/*              right_gap=list_of_exon_right[exon2]-list_of_exon_right[exon1]; */
/*              if(right_gap >= 0){ */
/*                temp_string=Substring(sequences[exon2], 0, strlen(sequences[exon2])-right_gap-1); */
/*                strcpy(sequences[exon1], temp_string); */
/*              } */
/*              else{ */
/*                temp_string=Substring(sequences[exon2], 0, strlen(sequences[exon2])-1); */
/*                strcpy(help_string, temp_string); */
/*                free(temp_string); */
/*                temp_string=Substring(sequences[exon1], strlen(sequences[exon1])+right_gap, strlen(sequences[exon1])-1); */
/*                strcat(help_string, temp_string); */
/*                strcpy(sequences[exon1], help_string); */
/*              } */
/*              free(temp_string); */
/*              return; */
/*       } */
/*       else{ */
/* //20dic06 */
/*              is_internal[exon1]=1; */

/*              list_of_exon_right[exon1]=list_of_exon_right[exon2]; */
/*              left_gap=list_of_exon_left[exon2]-list_of_exon_left[exon1]; */
/*              if(left_gap <= 0){ */
/*                temp_string=Substring(sequences[exon2], -left_gap, strlen(sequences[exon2])-1); */
/*                strcpy(sequences[exon1], temp_string); */
/*              } */
/*              else{ */
/*                temp_string=Substring(sequences[exon1], 0, left_gap-1); */
/*                strcpy(help_string, temp_string); */
/*                free(temp_string); */
/*                temp_string=Substring(sequences[exon2], 0, strlen(sequences[exon2])-1); */
/*                strcat(help_string, temp_string); */
/*                strcpy(sequences[exon1], help_string); */
/*              } */
/*              free(temp_string); */
/*              return; */
/*       } */
/*   } */

/*   if(is_internal[exon1] == 1){ */
/*       if(is_internal[exon2] == -1){ */
/* //20dic06 */
/*              is_internal[exon2]=1; */

/*              list_of_exon_left[exon2]=list_of_exon_left[exon1]; */
/*              right_gap=list_of_exon_right[exon2]-list_of_exon_right[exon1]; */
/*              if(right_gap <= 0){ */
/*                temp_string=Substring(sequences[exon1], 0, strlen(sequences[exon1])+right_gap-1); */
/*                strcpy(sequences[exon2], temp_string); */
/*              } */
/*              else{ */
/*                temp_string=Substring(sequences[exon1], 0, strlen(sequences[exon1])-1); */
/*                strcpy(help_string, temp_string); */
/*                free(temp_string); */
/*                temp_string=Substring(sequences[exon2], strlen(sequences[exon2])-right_gap, strlen(sequences[exon2])-1); */
/*                strcat(help_string, temp_string); */
/*                strcpy(sequences[exon2], help_string); */
/*              } */
/*              free(temp_string); */
/*              return; */
/*       } */
/*       else{ */
/* //20dic06 */
/*              is_internal[exon2]=1; */

/*              list_of_exon_right[exon2]=list_of_exon_right[exon1]; */
/*              left_gap=list_of_exon_left[exon2]-list_of_exon_left[exon1]; */
/*              if(left_gap >= 0){ */
/*                temp_string=Substring(sequences[exon1], left_gap, strlen(sequences[exon1])-1); */
/*                strcpy(sequences[exon2], temp_string); */
/*              } */
/*              else{ */
/*                temp_string=Substring(sequences[exon2], 0, -(left_gap+1)); */
/*                strcpy(help_string, temp_string); */
/*                free(temp_string); */
/*                temp_string=Substring(sequences[exon1], 0, strlen(sequences[exon1])-1); */
/*                strcat(help_string, temp_string); */
/*                strcpy(sequences[exon2], help_string); */
/*              } */
/*              free(temp_string); */
/*              return; */
/*       } */
/*   } */

/*   if(is_internal[exon1] == -1){ */
/*       if(is_internal[exon2] == -1){ */
/*              if(list_of_exon_left[exon2] < list_of_exon_left[exon1]){ */
/*                list_of_exon_left[exon1]=list_of_exon_left[exon2]; */
/*                right_gap=list_of_exon_right[exon2]-list_of_exon_right[exon1]; */
/*                if(right_gap >= 0){ */
/*                       temp_string=Substring(sequences[exon2], 0, strlen(sequences[exon2])-right_gap-1); */
/*                       strcpy(sequences[exon1], temp_string); */
/*                } */
/*                else{ */
/*                       temp_string=Substring(sequences[exon2], 0, strlen(sequences[exon2])-1); */
/*                       strcpy(help_string, temp_string); */
/*                       free(temp_string); */
/*                       temp_string=Substring(sequences[exon1], strlen(sequences[exon1])+right_gap, strlen(sequences[exon1])-1); */
/*                       strcat(help_string, temp_string); */
/*                       strcpy(sequences[exon1], help_string); */
/*                } */
/*                free(temp_string); */
/*                return; */
/*              } */
/*              else{ */
/*                list_of_exon_left[exon2]=list_of_exon_left[exon1]; */
/*                right_gap=list_of_exon_right[exon2]-list_of_exon_right[exon1]; */
/*                if(right_gap <= 0){ */
/*                       temp_string=Substring(sequences[exon1], 0, strlen(sequences[exon1])+right_gap-1); */
/*                       strcpy(sequences[exon2], temp_string); */
/*                } */
/*                else{ */
/*                       temp_string=Substring(sequences[exon1], 0, strlen(sequences[exon1])-1); */
/*                       strcpy(help_string, temp_string); */
/*                       free(temp_string); */
/*                       temp_string=Substring(sequences[exon2], strlen(sequences[exon2])-right_gap, strlen(sequences[exon2])-1); */
/*                       strcat(help_string, temp_string); */
/*                       strcpy(sequences[exon2], help_string); */
/*                } */
/*                free(temp_string); */
/*                return; */
/*              } */
/*       } */
/*       else{ */
/*              list_of_exon_left[exon1]=list_of_exon_left[exon2]; */
/*              list_of_exon_right[exon2]=list_of_exon_right[exon1]; */
/*              right_gap=list_of_exon_right[exon2]-list_of_exon_right[exon1]; */
/*              if(right_gap >= 0){ */
/*                temp_string=Substring(sequences[exon2], 0, strlen(sequences[exon2])-right_gap-1); */
/*                strcpy(sequences[exon1], temp_string); */
/*                strcpy(sequences[exon2], temp_string); */
/*              } */
/*              else{ */
/*                temp_string=Substring(sequences[exon2], 0, strlen(sequences[exon2])-1); */
/*                strcpy(help_string, temp_string); */
/*                free(temp_string); */
/*                temp_string=Substring(sequences[exon1], strlen(sequences[exon1])+right_gap, strlen(sequences[exon1])-1); */
/*                strcat(help_string, temp_string); */
/*                strcpy(sequences[exon1], help_string); */
/*                strcpy(sequences[exon2], help_string); */
/*              } */
/*              free(temp_string); */
/*              return; */
/*       } */
/*   } */

/*   if(is_internal[exon2] == -2){ */
/*       if(list_of_exon_right[exon2] > list_of_exon_right[exon1]){ */
/*              list_of_exon_right[exon1]=list_of_exon_right[exon2]; */
/*              left_gap=list_of_exon_left[exon2]-list_of_exon_left[exon1]; */
/*              if(left_gap <= 0){ */
/*                temp_string=Substring(sequences[exon2], -left_gap, strlen(sequences[exon2])-1); */
/*                strcpy(sequences[exon1], temp_string); */
/*              } */
/*              else{ */
/*                temp_string=Substring(sequences[exon1], 0, left_gap-1); */
/*                strcpy(help_string, temp_string); */
/*                free(temp_string); */
/*                temp_string=Substring(sequences[exon2], 0, strlen(sequences[exon2])-1); */
/*                strcat(help_string, temp_string); */
/*                strcpy(sequences[exon1], help_string); */
/*              } */
/*              free(temp_string); */
/*              return; */
/*       } */
/*       else{ */
/*              list_of_exon_right[exon2]=list_of_exon_right[exon1]; */
/*              left_gap=list_of_exon_left[exon2]-list_of_exon_left[exon1]; */
/*              if(left_gap >= 0){ */
/*                temp_string=Substring(sequences[exon1], left_gap, strlen(sequences[exon1])-1); */
/*                strcpy(sequences[exon2], temp_string); */
/*              } */
/*              else{ */
/*                temp_string=Substring(sequences[exon2], 0, -(left_gap+1)); */
/*                strcpy(help_string, temp_string); */
/*                free(temp_string); */
/*                temp_string=Substring(sequences[exon1], 0, strlen(sequences[exon1])-1); */
/*                strcat(help_string, temp_string); */
/*                strcpy(sequences[exon2], help_string); */
/*              } */
/*              free(temp_string); */
/*              return; */
/*       } */
/*   } */
/*   else{ */
/*       list_of_exon_left[exon2]=list_of_exon_left[exon1]; */
/*       list_of_exon_right[exon1]=list_of_exon_right[exon2]; */
/*       right_gap=list_of_exon_right[exon2]-list_of_exon_right[exon1]; */
/*       if(right_gap <= 0){ */
/*              temp_string=Substring(sequences[exon1], 0, strlen(sequences[exon1])+right_gap-1); */
/*              strcpy(sequences[exon1], temp_string); */
/*              strcpy(sequences[exon2], temp_string); */
/*       } */
/*       else{ */
/*              temp_string=Substring(sequences[exon1], 0, strlen(sequences[exon1])-1); */
/*              strcpy(help_string, temp_string); */
/*              free(temp_string); */
/*              temp_string=Substring(sequences[exon2], strlen(sequences[exon2])-right_gap, strlen(sequences[exon2])-1); */
/*              strcat(help_string, temp_string); */
/*              strcpy(sequences[exon1], help_string); */
/*              strcpy(sequences[exon2], help_string); */
/*       } */
/*       free(temp_string); */
/*       return; */
/*   } */
/* } */

//Funzione che aggiunge una lista di percorsi ad una lista di percorsi
void Add_Path_List(struct path **path_list1, struct path **path_list2){
  struct path *head=*path_list2;

  while(head != NULL){
         Add_Path(path_list1, &head);
         head=head->next;
  }
}

//Fare pero' una copia di path prima di aggiungere (variabile globale)
void Add_Path(struct path **path_list, struct path **path){
  struct path *head=*path_list;
  char stop=0;

  //struct node *nd=NULL;
//  struct path *head2=NULL;

  path_to_be_added=Copy_of_Path(*path);


/*fprintf(stdout, "     Path to be added:");
  nd=path_to_be_added->n;
  while(nd != NULL){
  fprintf(stdout, " %d", nd->index);
  nd=nd->next;
  }
  fprintf(stdout, "\n");
  head2=*path_list;
  while(head2 != NULL){
  nd=head2->n;
  while(nd != NULL){
  fprintf(stdout, " %d", nd->index);
  nd=nd->next;
  }
  fprintf(stdout, "----");
  head2=head2->next;
  }
  fprintf(stdout, "\n");*/

  while(head != NULL && !stop){
         if(Equal_Paths(head, path_to_be_added))
                stop=1;
         else
                head=head->next;
  }

  if(head != NULL)
         return;

  path_to_be_added->next=NULL;

  if(*path_list == NULL)
         *path_list=path_to_be_added;
  else{
         path_to_be_added->next=*path_list;
         *path_list=path_to_be_added;
  }
}

//Fare funzione che svuota una lista di percorsi
void Empty_Source_List_of_Paths(){
  struct path *head=NULL;
  struct path *help=NULL;
  int i=0;

  for(i=0; i< source_total_paths; i++){
         head=source_list_of_paths[i];
         while(head != NULL){
                help=head->next;
                free(head);
                head=help;
         }
         source_list_of_paths[i]=NULL;
  }
}

char Equal_Paths(struct path *path1, struct path *path2){
  char stop=0;
  struct node *nodes1=NULL, *nodes2=NULL;

  if(path1->end != path2->end)
         return 0;

  nodes1=path1->n;
  nodes2=path2->n;

  while(nodes1 != NULL && nodes2 != NULL && !stop){
         if(nodes1->index != nodes2->index)
                stop=1;
         else{
                nodes1=nodes1->next;
                nodes2=nodes2->next;
         }
  }

  if(stop)
         return 0;
  else{
         if(nodes1 == NULL && nodes2 == NULL)
                return 1;
         else
                return 0;
  }
}

void Graph_reduction(){
  int i=0, j=0;

  for(i=0; i<number_of_transcripts; i++){
         for(j=0; j<number_of_transcripts; j++){
                if(extension_matrix[i][j] != 0)
                  Partial_Graph_reduction_for_arc(i, j);
         }
  }
}

//27giu05
/* void Graph_reduction_Myers(){ */
/*   int i=0, j=0, k=0; */
/*   char out=0, in=0, stop=0; */

/*   while(k < number_of_transcripts){ */
/*       i=0; */
/*       while(i < number_of_transcripts){ */
/*              if(extension_matrix[k][i] != 0){ */
/*                j=0; */
/*                while(j < number_of_transcripts){ */
/*                       if(extension_matrix[j][k] != 0){ */
/*                              if(extension_matrix[j][i] != 0){ */
/*                                remove_extension_matrix[j][k]=1; */
/*                                remove_extension_matrix[k][i]=1; */
/*                              } */
/*                       } */
/*                       j++; */
/*                } */
/*              } */
/*              i++; */
/*       } */
/*       k++; */
/*   } */

/*   count_del=0; */

/*   k=0; */
/*   while(k < number_of_transcripts){ */
/*       i=0; */
/*       out=0; */
/*       while(i < number_of_transcripts && !stop){ */
/*              if(remove_extension_matrix[k][i] == 0 && extension_matrix[k][i] != 0){ */
/*                out=1; */
/*                stop=1; */
/*              } */
/*              else */
/*                i++; */
/*       } */
/*       i=0; */
/*       in=0; */
/*       stop=0; */
/*       while(i < number_of_transcripts && !stop){ */
/*              if(remove_extension_matrix[i][k] == 0 && extension_matrix[i][k] != 0){ */
/*                in=1; */
/*                stop=1; */
/*              } */
/*              else */
/*                i++; */
/*       } */
/* //30giu05 */
/*       if((in_degree[k] > 0 && out_degree[k] > 0) && (in == 0 && out == 0)){ */
/*              i=0; */
/*              while(i < number_of_transcripts){ */
/*                if(extension_matrix[k][i] != 0){ */
/*                       extension_matrix[k][i]=0; */
/*                       count_del++; */
/*                } */
/*                i++; */
/*              } */
/*              i=0; */
/*              while(i < number_of_transcripts){ */
/*                if(extension_matrix[i][k] != 0){ */
/*                       extension_matrix[i][k]=0; */
/*                       count_del++; */
/*                } */
/*                i++; */
/*              } */
/*       } */
/*       k++; */
/*   } */

/*   k=0; */
/*   while(k < number_of_transcripts){ */
/*       i=0; */
/*       while(i < number_of_transcripts){ */
/*              if(extension_matrix[k][i] != 0){ */
/*                j=0; */
/*                while(j < number_of_transcripts){ */
/*                       if(extension_matrix[j][k] != 0){ */
/*                              if(extension_matrix[j][i] != 0){ */
/*                                remove_extension_matrix[j][i]=2; */
/*                              } */
/*                       } */
/*                       j++; */
/*                } */
/*              } */
/*              i++; */
/*       } */
/*       k++; */
/*   } */

/*   k=0; */
/*   while(k < number_of_transcripts){ */
/*       i=0; */
/*       while(i < number_of_transcripts){ */
/*              if(remove_extension_matrix[k][i] == 2){ */
/*                extension_matrix[k][i]=0; */
/*                count_del++; */
/*              } */
/*              i++; */
/*       } */
/*       k++; */
/*   } */
/* } */

//Dato un arco (a,b) restituisce il nodo c tale che esistano gli archi (a,c) e (c,b). Se restituisce -1, significa
//che o (a,b) non e' un arco, o che non esite un siffatto nodo c, initial_index e' l'indice minimo da cui partire nella
//ricerca
int Get_Opposite_Node_Index(int a_index, int b_index, int initial_index){
  char stop=0;
  int c_index=-1, i=0;

  if(extension_matrix[a_index][b_index] == 0)
         return c_index;

  i=initial_index;
  while(i < number_of_transcripts && !stop){
         if(extension_matrix[a_index][i] != 0 && extension_matrix[i][b_index] != 0){
                stop=1;
                c_index=i;
         }
         else
                i++;
  }

  return c_index;
}

//Dati a,b, definiti nella procedura Get_Opposite_Node_Index, effettuare la riduzione del grafo
void Partial_Graph_reduction_for_arc(int a_index, int b_index){
  int c_index=-1;
  int initial_index=0;

  do{
         c_index=Get_Opposite_Node_Index(a_index, b_index, initial_index);
         if(c_index != -1){
                Partial_Graph_reduction_for_node(a_index, b_index, c_index);
                initial_index=c_index+1;
         }
  }while(c_index != -1);
}

//10gen05
//Dati a,b,c definiti nella procedura Get_Opposite_Node_Index, effettuare la riduzione del grafo
void Partial_Graph_reduction_for_node(int a_index, int b_index, int c_index){
  struct node *node_list=NULL, *out_node_list=NULL;
  struct node *head=NULL, *help_head=NULL, *help=NULL, *out_head=NULL;
  struct node *help_node_list=NULL;
  int i=0;
  char changed=1, stop=0;
  char no_outcoming=1, attached=1;;

//c non ha archi uscenti?
  i=0;
  while(i < number_of_transcripts){
         if(extension_matrix[c_index][i] != 0 && i != b_index){
                no_outcoming=0;
                Add_Node_to_a_node_list(&out_node_list, i);
         }
         i++;
  }

//Considero tutti i nodi che entrano in quello relativo a c_index (eccetto il nodo a_index)
  i=0;
  while(i < number_of_transcripts){
         if(extension_matrix[i][c_index] != 0 && i != a_index){
                Add_Node_to_a_node_list(&node_list, i);
         }
         i++;
  }

//Per tutti gli n in node list tale che esiste l'arco (n,a) o (n,b), cancello l'arco (n,c)
  head=node_list;
  while(head != NULL){
//09mar07
         help=head->next;

         if(extension_matrix[head->index][a_index] != 0){
                extension_matrix[head->index][c_index]=0;
                out_degree[head->index]=out_degree[head->index]-1;
                in_degree[c_index]=in_degree[c_index]-1;
                Add_Node_to_a_node_list(&help_node_list, head->index);
                Remove_Node_from_a_node_list(&node_list, head->index);
         }
         else{
                if(extension_matrix[head->index][b_index] != 0){
                  if(no_outcoming){
                         extension_matrix[head->index][c_index]=0;
                         out_degree[head->index]=out_degree[head->index]-1;
                         in_degree[c_index]=in_degree[c_index]-1;
                         Add_Node_to_a_node_list(&help_node_list, head->index);
                         Remove_Node_from_a_node_list(&node_list, head->index);
                  }
                  else{
                         attached=1;
                         out_head=out_node_list;
                         while(out_head != NULL && attached){
                                if(extension_matrix[head->index][out_head->index] == 0)
                                  attached=0;

                                out_head=out_head->next;
                         }
                         if(attached){
                                extension_matrix[head->index][c_index]=0;
                                out_degree[head->index]=out_degree[head->index]-1;
                                in_degree[c_index]=in_degree[c_index]-1;
                                Add_Node_to_a_node_list(&help_node_list, head->index);
                                Remove_Node_from_a_node_list(&node_list, head->index);
                         }
                  }
                }
         }
//09mar07
//head=head->next;
         head=help;
  }

  while(changed){
         changed=0;
         head=node_list;
         while(head != NULL){
                stop=0;
                help_head=help_node_list;
                while(help_head != NULL && !stop){
                  if(extension_matrix[head->index][help_head->index] != 0){
                         stop=1;
                         changed=1;
                         extension_matrix[head->index][c_index]=0;
                         out_degree[head->index]=out_degree[head->index]-1;
                         in_degree[c_index]=in_degree[c_index]-1;
                         Add_Node_to_a_node_list(&help_node_list, head->index);
                         help=head->next;
                         Remove_Node_from_a_node_list(&node_list, head->index);
                         head=help;
                  }
                  else
                         help_head=help_head->next;
                }
                if(!stop)
                  head=head->next;
         }
  }

//Rimuovo anche l'arco (c,b)
  if(node_list == NULL){
         extension_matrix[c_index][b_index]=0;
         out_degree[c_index]=out_degree[c_index]-1;
         in_degree[b_index]=in_degree[b_index]-1;
  }
}

//Dati a,b,c definiti nella procedura Get_Opposite_Node_Index, effettuare la riduzione del grafo
/*void Partial_Graph_reduction_for_node(int a_index, int b_index, int c_index){
  struct node *node_list=NULL;
  struct node *head=NULL, *help_head=NULL, *help=NULL;
  struct node *help_node_list=NULL;
  int i=0;
  char changed=1, stop=0;

//Considero tutti i nodi che entrano in quello relativo a c_index (eccetto il nodo a_index)
i=0;
while(i < number_of_transcripts){
if(extension_matrix[i][c_index] != 0 && i != a_index){
Add_Node_to_a_node_list(&node_list, i);
}
i++;
}

//Per tutti gli n in node list tale che esiste l'arco (n,a) o (n,b), cancello l'arco (n,c)
head=node_list;
while(head != NULL){
if(extension_matrix[head->index][a_index] != 0 || extension_matrix[head->index][b_index] != 0){
extension_matrix[head->index][c_index]=0;
out_degree[head->index]=out_degree[head->index]-1;
in_degree[c_index]=in_degree[c_index]-1;
Add_Node_to_a_node_list(&help_node_list, head->index);
Remove_Node_from_a_node_list(&node_list, head->index);
}
head=head->next;
}

while(changed){
changed=0;
head=node_list;
while(head != NULL){
stop=0;
help_head=help_node_list;
while(help_head != NULL && !stop){
if(extension_matrix[head->index][help_head->index] != 0){
stop=1;
changed=1;
extension_matrix[head->index][c_index]=0;
out_degree[head->index]=out_degree[head->index]-1;
in_degree[c_index]=in_degree[c_index]-1;
Add_Node_to_a_node_list(&help_node_list, head->index);
help=head->next;
Remove_Node_from_a_node_list(&node_list, head->index);
head=help;
}
else
help_head=help_head->next;
}
if(!stop)
head=head->next;
}
}

//Rimuovo anche l'arco (c,b)
if(node_list == NULL){
extension_matrix[c_index][b_index]=0;
out_degree[c_index]=out_degree[c_index]-1;
in_degree[b_index]=in_degree[b_index]-1;
}
}*/

void Add_Node_to_a_node_list(struct node **node_list, int node){
  struct node *nds=NULL;
  int counter=0;

  nds=*node_list;
  while(nds != NULL){
         if(nds->index == node){
                fprintf(stdout, "Cycle detected 2!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }
         nds=nds->next;
         counter++;
  }

  add=(struct node *)malloc(sizeof(struct node));
  if(add == NULL){
         fprintf(stderr, "Problem1 of memory allocation in Add_Node_to_a_node_list!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  add->index=node;
  add->next=*node_list;

  *node_list=add;
}

void Remove_Node_from_a_node_list(struct node **node_list, int node){
  struct node *nds=NULL;
  struct node *prev=NULL, *next_one=NULL;
  char stop=0;

  nds=*node_list;
  while(nds != NULL && !stop){
         if(nds->index == node){
                stop=1;
         }
         else{
                prev=nds;
                nds=nds->next;
                next_one=nds->next;
         }
  }

  if(prev == NULL)
         *node_list=next_one;
  else
         prev->next=next_one;

  free(nds);
}

//19gen05
/*Elimina i trascritti che praticamente coincidono a meno di una differenza tra gli esterni*/
void First_Filtering(){
  int i=0, j=0, pos=0;
  char *contained=NULL;
  char stop=0;
  int limit=0, inclusion=0;

  contained=(char *)malloc(number_of_transcripts*sizeof(char));
  if(contained == NULL){
         fprintf(stdout, "Contained not stored!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  for(i=0; i<number_of_transcripts; i++){
         contained[i]=0;
  }

  i=0;
  while(i<number_of_transcripts){
         if(!contained[i]){

                j=i+1;
                stop=0;

                while(j<number_of_transcripts && !stop){

//08giu05
//inclusion=Extends(transcript_list[i], transcript_list[j], &limit, 0);

//07ott05
//inclusion=Extends(transcript_list[i], transcript_list[j], &limit, 0, 0);
//19dic06
//inclusion=Extends(transcript_list[i], transcript_list[j], &limit, 0, 0, 1);
#ifdef DONT_EXTEND_REFSEQ
                  if(transcript_list[i].type == 1){
                         if(transcript_list[j].type == 0){
//08gen08
//inclusion=Overlap(transcript_list[j], transcript_list[i], &limit, 0, 0, 1);
                                inclusion=Overlap(transcript_list[j], transcript_list[i], &limit, 0, 0, 1, 1);
                         }
                         else
                                inclusion=0;
                  }
                  else{
                         if(transcript_list[j].type == 1){
//08gen08
//inclusion=Overlap(transcript_list[i], transcript_list[j], &limit, 0, 0, 1);
                                inclusion=Overlap(transcript_list[i], transcript_list[j], &limit, 0, 0, 1, 1);
                                if(inclusion == 2)
                                  inclusion=-2;
                         }
                         else{
//08gen08
//inclusion=Extends(transcript_list[i], transcript_list[j], &limit, 0, 0, 1);
                                inclusion=Extends(transcript_list[i], transcript_list[j], &limit, 0, 0, 1, 1);
                         }
                  }
#else
//08gen08
//inclusion=Extends(transcript_list[i], transcript_list[j], &limit, 0, 0, 1);
                  inclusion=Extends(transcript_list[i], transcript_list[j], &limit, 0, 0, 1, 1);
#endif

/*if(i == 0){
  if(inclusion == 2){
  fprintf(stderr, "%d included in 0 (%d)\n", j, transcript_list[j].type);
  }
  if(inclusion == -2)
  fprintf(stderr, "\n\nNot possible %d!!\n", j);
  }*/

                  if(inclusion == -2 || inclusion == 2){
//if((inclusion == -2 && transcript_list[i].type == 0) || (inclusion == 2 && transcript_list[j].type == 0)){

                         if(limit == 0){
                                if(transcript_list[i].exons == transcript_list[j].exons){
                                  if(inclusion == -2){
//25lug07
//Modificare transcript_list[j] con allungo esoni solo se non e' refseq
#ifdef DONT_EXTEND_REFSEQ
                                         if(transcript_list[j].type != 1){
#endif

//23gen08
                                                if(list_of_exon_right[transcript_list[j].left_ext] == list_of_exon_right[transcript_list[i].left_ext]){
                                                  if(limit == 0){
//06nov07
//if(is_internal[transcript_list[j].left_ext == -1]){
                                                         if(is_internal[transcript_list[j].left_ext] == -1){
                                                                if(is_internal[transcript_list[i].left_ext] == 1){
                                                                  transcript_list[j].left_ext=transcript_list[i].left_ext;
                                                                }
                                                                else{
                                                                  if(is_internal[transcript_list[i].left_ext] == -1){
                                                                         if(list_of_exon_left[transcript_list[i].left_ext] < list_of_exon_left[transcript_list[j].left_ext]){
                                                                                transcript_list[j].left_ext=transcript_list[i].left_ext;
                                                                         }
                                                                  }
                                                                }
                                                         }
                                                  }
//23gen08
                                                }

//23gen08
                                                if(list_of_exon_left[transcript_list[j].right_ext] == list_of_exon_left[transcript_list[i].right_ext]){
                                                  if(limit+transcript_list[i].exons == transcript_list[j].exons){

//06nov07
//if(is_internal[transcript_list[j].right_ext == -2]){
                                                         if(is_internal[transcript_list[j].right_ext] == -2){

                                                                if(is_internal[transcript_list[i].right_ext] == 1){
                                                                  transcript_list[j].right_ext=transcript_list[i].right_ext;
                                                                }
                                                                else{
                                                                  if(is_internal[transcript_list[i].right_ext] == -2){
                                                                         if(list_of_exon_right[transcript_list[i].right_ext] > list_of_exon_right[transcript_list[j].right_ext]){
                                                                                transcript_list[j].right_ext=transcript_list[i].right_ext;
                                                                         }
                                                                  }
                                                                }
                                                         }
                                                  }
//23gen08
                                                }
#ifdef DONT_EXTEND_REFSEQ
                                         }
#endif

                                         contained[i]=1;

                                         transcript_list[j].ESTs+=transcript_list[i].ESTs;
                                         stop=1;
                                  }
                                  else{
//25lug07
//Modificare transcript_list[i] con allungo esoni solo se non e' refseq
#ifdef DONT_EXTEND_REFSEQ
                                         if(transcript_list[i].type != 1){
#endif
//23gen08
                                                if(list_of_exon_right[transcript_list[j].left_ext] == list_of_exon_right[transcript_list[i].left_ext]){
                                                  if(limit == 0){
//06nov07
//if(is_internal[transcript_list[i].left_ext == -1]){
                                                         if(is_internal[transcript_list[i].left_ext] == -1){
                                                                if(is_internal[transcript_list[j].left_ext] == 1){
                                                                  transcript_list[i].left_ext=transcript_list[j].left_ext;
                                                                }
                                                                else{
                                                                  if(is_internal[transcript_list[j].left_ext] == -1){
                                                                         if(list_of_exon_left[transcript_list[j].left_ext] < list_of_exon_left[transcript_list[i].left_ext]){
                                                                                transcript_list[i].left_ext=transcript_list[j].left_ext;
                                                                         }
                                                                  }
                                                                }
                                                         }
                                                  }
//23gen08
                                                }

//23gen08
                                                if(list_of_exon_left[transcript_list[j].right_ext] == list_of_exon_left[transcript_list[i].right_ext]){
                                                  if(limit+transcript_list[j].exons == transcript_list[i].exons){
//06nov07
//if(is_internal[transcript_list[i].right_ext == -2]){
                                                         if(is_internal[transcript_list[i].right_ext] == -2){
                                                                if(is_internal[transcript_list[j].right_ext] == 1){
                                                                  transcript_list[i].right_ext=transcript_list[j].right_ext;
                                                                }
                                                                else{
                                                                  if(is_internal[transcript_list[j].right_ext] == -2){
                                                                         if(list_of_exon_right[transcript_list[j].right_ext] > list_of_exon_right[transcript_list[i].right_ext]){
                                                                                transcript_list[i].right_ext=transcript_list[j].right_ext;
                                                                         }
                                                                  }
                                                                }
                                                         }
                                                  }
//23gen08
                                                }
#ifdef DONT_EXTEND_REFSEQ
                                         }
#endif

                                         contained[j]=1;
                                         transcript_list[i].ESTs+=transcript_list[j].ESTs;
                                  }
/*if(i == 4){
  fprintf(stdout, "Transcripts %d and %d (%d) *******************\n", i, j, inclusion);
  fprintf(stdout, "TR1 %d", transcript_list[i].left_ext);
  for(k=0; k<transcript_list[i].exons-2; k++)
  fprintf(stdout, " %d", transcript_list[i].exon_list[k]);
  fprintf(stdout, " %d\n", transcript_list[i].right_ext);

  fprintf(stdout, "TR2 %d", transcript_list[j].left_ext);
  for(k=0; k<transcript_list[j].exons-2; k++)
  fprintf(stdout, " %d", transcript_list[j].exon_list[k]);
  fprintf(stdout, " %d\n", transcript_list[j].right_ext);
  fprintf(stdout, "*******************\n");
  fprintf(stdout, "i cont %d j cont %d\n", contained[i], contained[j]);
  }*/
                                }
                         }
                  }
                  j++;
                }
         }
         i++;
  }

  i=0;
  pos=0;
  while(i<number_of_transcripts){
         if(contained[i]){
                pos--;
         }
         else{
                transcript_list[i+pos].left_ext=transcript_list[i].left_ext;
                for(j=0; j<transcript_list[i].exons-2; j++)
                  transcript_list[i+pos].exon_list[j]=transcript_list[i].exon_list[j];
                transcript_list[i+pos].right_ext=transcript_list[i].right_ext;
                transcript_list[i+pos].exons=transcript_list[i].exons;
                transcript_list[i+pos].ESTs=transcript_list[i].ESTs;

//19dic06(aggiunto il 10gen07)
                transcript_list[i+pos].type=transcript_list[i].type;
                if(transcript_list[i+pos].type == 1)
                  strcpy(transcript_list[i+pos].RefSeq, transcript_list[i].RefSeq);
         }
         i++;
  }

  number_of_transcripts+=pos;

/*for(q=0; q<number_of_transcripts; q++){
  fprintf(stdout, "Transcripts %d (conf %d)*******************\n", q, transcript_list[q].ESTs);
  fprintf(stdout, "TR %d-%d(%d)", list_of_exon_left[transcript_list[q].left_ext], list_of_exon_right[transcript_list[q].left_ext], transcript_list[q].left_ext);
  fprintf(stdout, "Seq %s\n", sequences[transcript_list[q].left_ext]);
  for(k=0; k<transcript_list[q].exons-2; k++){
  fprintf(stdout, " %d-%d(%d)", list_of_exon_left[transcript_list[q].exon_list[k]], list_of_exon_right[transcript_list[q].exon_list[k]], transcript_list[q].exon_list[k]);
  fprintf(stdout, "Seq %s\n", sequences[transcript_list[q].exon_list[k]]);
  }
  fprintf(stdout, " %d-%d(%d)", list_of_exon_left[transcript_list[q].right_ext], list_of_exon_right[transcript_list[q].right_ext], transcript_list[q].right_ext);
  fprintf(stdout, "Seq %s\n", sequences[transcript_list[q].right_ext]);
  fprintf(stdout, "Type %d (polyA %d) *******************\n", transcript_list[q].type, polya[transcript_list[q].right_ext]);
  }*/

  free(contained);
}

#ifdef READ_ABS_COORD
int GetAbsoluteStart(int left, int right){
  if(strand == 1)
         return gen_start+left-(boundary+1);
  else
         return gen_end-right+(boundary+1);
}

int GetAbsoluteEnd(int left, int right){
  if(strand == 1)
         return gen_start+right-(boundary+1);
  else
         return gen_end-left+(boundary+1);
}
#endif

//21feb06
char *To_lower(char *str){
  size_t length=strlen(str);
  int i=0;

  if(isupper(str[0])){
         for(i=0; i < (int)length; i++){
                str[i]=tolower(str[i]);
         }
  }

  return str;
}
