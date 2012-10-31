/**
 *
 *
 *                              PIntron
 *
 * A novel pipeline for computational gene-structure prediction based on
 * spliced alignment of expressed sequences (ESTs and mRNAs).
 *
 * Copyright (C) 2010,2011  Yuri Pirola, Raffaella Rizzi
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "my_time.h"
#include "log.h"
#include "util.h"
#include "log-build-info.h"

#define MAX_NLD 50000           //Massimo numero di nucleotidi gestibili in una sequenza
#define MAX_EXONS 80            //Massimo numero di esoni gestibili per un trascritto
#define FROM_ONE                        //Se definita, allora le coord abs di ASPIC sono da 1 altrimenti da 0
#define TR_IND 33
#define PRINT_GTF_VARIANT
#define PRINT_FRAME

struct cds{

//Coordinate assolute
  int *cds_from;                /*left end degli esoni*/
  int *cds_to;                  /*right end degli esoni*/

  int exons;                            /*Numero di esoni nel cds*/
};

struct annotated_cds{

  char RefSeq[20];                      //ID del trascritto
//Coordinate assolute
  int rel_start;                //start relativo della CDS
  int rel_end;                  //end relativo della CDS
  int exons;
  char *RefSeq_sequence;                //Sequenza nt
};

struct intron{

  int left;             /*left end degli introni*/
  int right;            /*right end degli introni*/

  char **IDs;                   /*Lista degli est id*/

  int ESTs;

  char type;                    /*0: mRNA; 1: confermato da almeno 2 EST; 2: altrimenti*/

  char **RefSeq;
  int RefSeqNum;
};

struct region{

  struct region *prev;                  /*pointer to the previous region*/

  struct region *next;                  /*pointer to the next region*/

  int left;                             /*position of the first character on the genomic sequence*/

  int right;                            /*position of the last character on the genomic sequence*/

  char phase1;
  char phase2;

  char type;
};

struct genomic_exon{
  int rel_left;
  int rel_right;

  char *sequence;

  struct genomic_exon *prev;
  struct genomic_exon *next;
};

struct exon{

//Coordinate assolute (su 5'3' per strand=1 e su 3'5' per strand=-1)

  int left;                     /*left end degli esoni*/
  int right;                    /*right end degli esoni*/

  int rel_left;                 /*relative left end degli esoni*/
  int rel_right;                        /*relative right end degli esoni*/

//first UTR e' 5'UTR per strand=1 e 3'UTR per strand=-1; second UTR e' 3'UTR per strand=1 e 5'UTR per strand=-1
  char pos_flag_from;           /*0: il left e' nel CDS; 1: e' in first UTR; 2: e' in second UTR*/
  char pos_flag_to;                     /*0: il right e' nel CDS; 1: e' in first UTR; 2: e' in second UTR*/

  char polyA;
  char is_int;                  /*1: entrambi gli estremi confermati; 2: solo estr dx; 3: solo estr sx*/

  int matrix_index;

  char covered_exon;
  int cover_index;

  int variant_label_first;
  int variant_label_second;

  struct exon *prev;
  struct exon *next;

  char *sequence;
};

struct transcript{

  char is_annotated;            /*1 se e' stata recuperata l'annotazione del CDS, 0 altrimenti (vale solo per i RefSeq*/

  int *tr_from;         /*left end degli esoni (temporaneo solo come supporto iniziale)*/
  int *tr_to;                   /*right end degli esoni (temporaneo solo come supporto iniziale)*/

  char **temp_sequences;                /*sequenza degli esoni (temporaneo solo come supporto iniziale)*/

  int *exon_index;                      /*indici degli esoni*/

  int exons;                            /*Numero di esoni nel trascritto*/

  int length;                   /*Lunghezza in basi*/

  int first_ORF_index;          /*Indice (in exon_index) dell'esone in cui inizia l'ORF*/
  int second_ORF_index;         /*Indice (in exon_index) dell'esone in cui finisce l'ORF*/

  int abs_ORF_start;
  int abs_ORF_end;

  int ORF_start;
  int ORF_end;
  char start_cons;
  char end_cons;
  char start_c[4];
  char stop_c[4];

  char type;                            /*0: mRNA; 1: confermato da almeno 2 EST; 2: altrimenti*/
  char *RefSeq;

  char has_stop;
  char no_ATG;

//              Allineamenti EST-geomica per ogni esoni (sse il trascritto e' refseq type=0)
  char **EST_exon_alignments;
  char **GEN_exon_alignments;
};

struct new_label_for_exon{
  struct new_label_for_exon *next;
  struct new_label_for_exon *prev;

  int left, right;

  char *representation;
  //char character;
};

int align_dim;  //Dimensione della matrice di allineamento

char AlignEST[40000];
char AlignGenomic[40000];

struct new_label_for_exon **list_of_new_labels;
struct new_label_for_exon *computed_newlabel_copy;

struct transcript *trs;                 /*Struttura che contiene tutti i trascritti*/
int number_of_transcripts;

int *order_index;

struct intron *introns;
int number_of_introns;

int number_of_cds;
struct annotated_cds *a_cds;

char strand;

struct exon *help_exons=NULL;
struct exon *exons=NULL;

struct genomic_exon *gen_exons=NULL;
struct genomic_exon *gen_exon_copy=NULL;

int number_of_cover_exons;
int number_of_exons;
struct exon *exon_copy=NULL;

struct region *computed_region_copy=NULL;
struct exon *computed_exon_copy=NULL;

char gen_length_str[20];
char *out_path;
char *current_path;

char *in_path;

char *gene;

static void exit_with_problem(const char* problem) {
  fprintf(stderr, "%s\n", problem);
#ifdef HALT_EXIT_MODE
  exit(1);
#else
  exit(EXIT_FAILURE);
#endif
}

static void exit_with_problem_if(const bool condition,
											const char* problem) {
  if (condition) {
	 exit_with_problem(problem);
  }
}

static int GetCDSStart(struct cds cds_inst);
static int GetCDSEnd(struct cds cds_inst);
static void MarkExonEndpoints(struct cds cds_to_be_mapped);
static void GetIntronList(char *fileName);
static void MarkIntronType();
static void MarkTranscriptType(struct transcript *tr);

static void PrintTABOutput(int ref, struct cds cds_for_gene);


static void Set_cover_exon();
static void SetAltSplMatrixIndexes();


static void getCompetingLabels(int index, int ref, char *label);

static void getEXInitTermSkipNewLabels(int index, int ref, char *label);
static void getnewIRLabels(int index, int ref, char *label);

static void GetLocalization(char *local, int index, int exon);

static void GetLongestORF(int ref, int i, int min_length);
static char getContext(int ATG_start, char *tr_seq);

static void GetCDSAnnotations(char *fileName);

static char GetCDSAnnotationForRefSeq_2(int i);
static void CheckStartEndWRTref(int ref, int i);

static bool Check_start_codon(int pos, char *tr_seq);
static char Check_stop_codon(int pos, char *tr_seq);
static void Get_Transcripts_from_File_FASTA_format(char *fileName);
static struct exon *Insert_exon_into_a_exon_list(struct exon *arg_exon_list, int left, int right, int rel_left, int rel_right, char polyA, char *sequence, int *incr);

static struct genomic_exon *Insert_genexon_into_a_genexon_list(struct genomic_exon *arg_genexon_list, int rel_left, int rel_right, char *sequence);

static char isInFrame(int index, int ref);
static void PrintOutputFile(int ref);
static void SetPrintOrder(int ref);

static int SetREFToLongestTranscript();
static char GetLongestORFforCCDS(struct cds *cds_for_gene, int i, pmytime p);

static struct new_label_for_exon *Insert_newlabel_into_a_newlabel_list(struct new_label_for_exon *arg_newlabel_list, int left, int right, char **representation);

static void ComputeAlignment(char *EST_exon, char *genomic_exon);
static void ComputeAlignMatrix(char *EST_exon, char *genomic_exon, int n, int m, char **Mdir);
static void TracebackAlignment(char *EST_exon, char *genomic_exon, char **Mdir, int i, int j);
static void GetExonAlignments();
static void GetGenomicExons(char *fileName);
static char *GetGENexonSequence(int rel_left, int rel_right);

char* int2alpha(unsigned int num);

int Tcds;       //Lunghezza minima per le ORF

char* pathToDirectory(char* path){

  int i = strlen(path) - 1;
  while( path[i] != '/' ){
         i--;
         if(i < 0) break;
  }

  if(i < 0) return "";
  char* newDir = (char*)malloc((i+1) * sizeof(char));

  return (char*)memcpy(newDir, path, (i+1));
}

int main(int argc, char *argv[]){
  INFO("CDS-ANNOTATION");
  PRINT_LICENSE_INFORMATION;
  PRINT_SYSTEM_INFORMATION;
  int trs_length=0;

  if(argc != 5)
         {
                printf("Error!!\nUsage: ./CCDS [in_files_path] [out_files_path] [gene_name] [organism]\n\n"); return -3;
         }
  out_path = argv[2];

  in_path = argv[1];

  if( (current_path = pathToDirectory(argv[0])) == NULL ){
         fprintf(stderr, "\nCannot get current_path\n");
         return -1;
  }

  gene = argv[3];

  char *organism;
  organism=argv[4];

  struct cds cds_for_gene;
  int i=0, j=0;
  int ref=-1;

  char *temp = (char *) malloc(255*sizeof(char)); temp[0] = '\0';

  pmytime pt_tot= MYTIME_create_with_name("Total");

  MYTIME_start(pt_tot);

  sprintf(temp,"%scds",in_path);
  GetCDSAnnotations(temp);

  sprintf(temp,"%sisoforms.txt",out_path);
  Get_Transcripts_from_File_FASTA_format(temp);

  sprintf(temp,"%spredicted-introns.txt",out_path);

  GetIntronList(temp);

  sprintf(temp,"%sgenomic-exonforCCDS.txt",out_path);

  GetGenomicExons(temp);
  free(temp);

  GetExonAlignments();

  MarkIntronType();

  for(i=0; i<number_of_transcripts; i++){
         trs_length=0;
         for(j=0; j<trs[i].exons; j++){
                trs_length+=strlen(exons[trs[i].exon_index[j]].sequence);
         }
         trs[i].length=trs_length;
  }

  for(i=0; i<number_of_transcripts; i++){
         MarkTranscriptType(&trs[i]);
  }

  //ref=SetREFToLongestTranscript();

  i=0;
  Tcds=100;     //Lunghezza minima delle ORF
  while(i < number_of_transcripts){
         if(trs[i].type == 0){
                if(GetCDSAnnotationForRefSeq_2(i)){
                  trs[i].is_annotated=1;
                }
                else{
                  trs[i].is_annotated=0;
                }
         }
         i++;
  }

  i=0;
  while(i < number_of_transcripts){
         if(trs[i].type != 0 || trs[i].is_annotated == 0){
                GetLongestORF(ref, i, Tcds);
         }
         i++;
  }

  ref=SetREFToLongestTranscript();

  i=0;
  while(i < number_of_transcripts){
         CheckStartEndWRTref(ref, i);
         i++;
  }

  if(number_of_transcripts > 0){
         GetLongestORFforCCDS(&cds_for_gene, ref, pt_tot);
  }
  else{
          cds_for_gene.exons=0;
          cds_for_gene.cds_from=NULL;
          cds_for_gene.cds_to=NULL;
  }

  if(number_of_transcripts > 0)
         MarkExonEndpoints(cds_for_gene);

  Set_cover_exon();

  SetAltSplMatrixIndexes();

  if(number_of_transcripts > 0)
         SetPrintOrder(ref);

  PrintTABOutput(ref, cds_for_gene);

  PrintOutputFile(ref);

  free(order_index);

  MYTIME_stop(pt_tot);
  MYTIME_LOG(INFO, pt_tot);

  MYTIME_destroy(pt_tot);

  INFO("End");
  resource_usage_log();

}

int GetCDSStart(struct cds cds_inst){
  if(cds_inst.cds_from == NULL)
         return -1;

  return cds_inst.cds_from[0];
}

int GetCDSEnd(struct cds cds_inst){
  if(cds_inst.cds_from == NULL)
         return -1;

  return cds_inst.cds_to[cds_inst.exons-1];
}

void MarkExonEndpoints(struct cds cds_to_be_mapped){
  int i=0;
  int cds_start=0;
  int cds_end=0;

  if(cds_to_be_mapped.exons > 0){
         cds_start=GetCDSStart(cds_to_be_mapped);
         cds_end=GetCDSEnd(cds_to_be_mapped);
  }
  else{
         fprintf(stderr, "ERROR: CCDS not set 1!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  for(i=0; i<number_of_exons; i++){
         if(exons[i].left >= cds_start && exons[i].left <= cds_end)
                exons[i].pos_flag_from=0;
         else{
                if(exons[i].left < cds_start){
                  exons[i].pos_flag_from=1;
                }
                else{
                  exons[i].pos_flag_from=2;
                }
         }

         if(exons[i].right >= cds_start && exons[i].right <= cds_end)
                exons[i].pos_flag_to=0;
         else{
                if(exons[i].right < cds_start){
                  exons[i].pos_flag_to=1;
                }
                else{
                  exons[i].pos_flag_to=2;
                }
         }
  }
}

void Get_Transcripts_from_File_FASTA_format(char *fileName){
  FILE *in=NULL;
  char temp_string[MAX_NLD];
  char is_int=0;
  int i=0, j=0, k=0, p=0, z=0;
  char stop=0;
  int coord_counter=0;
  char temp_coord[20];
  int counter=0;
  int exons1=0;

  char temp_exons1[10];

  struct exon *head=NULL, *help=NULL;
  int tr_ID=0;

  char temp_tr_ID[20];
  char temp_refseq[20];

  struct exon *temp_exons=NULL;
  int left=0, right=0, rel_left=0, rel_right=0;
  char tmp_polyA=0;
  int incr=0;

  number_of_exons=0;

  in=fopen(fileName, "r");
  if(in == NULL){
         fprintf(stderr, "Error3!\n");
         exit(0);
  }

//Lettura delle prime due righe
  for(i=0; i<2; i++){
         fscanf(in, "%s\n", temp_string);

         if(i==0)
                number_of_transcripts=atoi(temp_string);
         else{
                strcpy(gen_length_str, temp_string);
         }
  }

  trs=(struct transcript *)malloc(number_of_transcripts*sizeof(struct transcript));
  exit_with_problem_if(trs == NULL, "Error4!");

  for(i=0; i<number_of_transcripts; i++){
	 fscanf(in, "%s\n", temp_string);

	 exit_with_problem_if(temp_string[0] != '>', "Invalid format!");

	 j=1;
	 while(j < (int)strlen(temp_string) && temp_string[j] != ':'){
		temp_tr_ID[j-1]=temp_string[j];
		j++;
	 }
	 temp_tr_ID[j-1]='\0';
	 tr_ID=atoi(temp_tr_ID);

	 if(j < (int)strlen(temp_string))
		j++;
	 while(j < (int)strlen(temp_string) && temp_string[j] != ':'){
		temp_exons1[j-strlen(temp_tr_ID)-2]=temp_string[j];
		j++;
	 }
	 temp_exons1[j-strlen(temp_tr_ID)-2]='\0';
	 exons1=atoi(temp_exons1);
	 if(j < (int)strlen(temp_string))
		j++;
	 while(j < (int)strlen(temp_string)){
		temp_refseq[j-strlen(temp_tr_ID)-strlen(temp_exons1)-3]=temp_string[j];
		j++;
	 }
	 if(j > (int)strlen(temp_tr_ID)+(int)strlen(temp_exons1)+3){
		temp_refseq[j-strlen(temp_tr_ID)-strlen(temp_exons1)-3]='\0';
	 }
	 else{
		strcpy(temp_refseq, "");
	 }

	 exons1=atoi(temp_exons1);

	 trs[i].exons=exons1;

	 if(strcmp(temp_refseq, "")){
		trs[i].type=0;
		trs[i].RefSeq=(char *)malloc((strlen(temp_refseq)+1)*sizeof(char));
		exit_with_problem_if(trs[i].RefSeq == NULL,
									"Problem10 of memory allocation in Get_Transcripts_from_File_FASTA_format!");
		strcpy(trs[i].RefSeq, temp_refseq);
	 } else {
		trs[i].type= -1;
		trs[i].RefSeq= NULL;
	 }

         trs[i].exon_index=(int *)malloc(exons1*sizeof(int));
         if(trs[i].exon_index == NULL){
                fprintf(stderr, "Problem11 of memory allocation in Get_Transcripts_from_File_FASTA_format!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }
         trs[i].tr_from=(int *)malloc(trs[i].exons*sizeof(int));
         trs[i].tr_to=(int *)malloc(trs[i].exons*sizeof(int));
         if(trs[i].tr_from == NULL || trs[i].tr_to == NULL){
                fprintf(stderr, "Error5!\n");
                exit(0);
         }

         trs[i].temp_sequences=(char **)malloc(trs[i].exons*sizeof(char *));
         if(trs[i].temp_sequences == NULL){
                fprintf(stderr, "Error51!\n");
                exit(0);
         }
         for(j=0; j<trs[i].exons; j++){
                trs[i].temp_sequences[j]=(char *)malloc(MAX_NLD*sizeof(char));
                if(trs[i].temp_sequences[j] == NULL){
                  fprintf(stderr, "Error511!\n");
                  exit(0);
                }
         }

         for(j=0; j<exons1; j++){
                fscanf(in, "%s\n", temp_string);

                k=0;
                coord_counter=0;
                while(k < (int)strlen(temp_string)){
                  p=0;
                  while(temp_string[k] != ':' && k < (int)strlen(temp_string)){
                         temp_coord[p]=temp_string[k];
                         k++;
                         p++;
                  }
                  temp_coord[p]='\0';
                  if(p > 0){
                         if(coord_counter==0){
                                left=atoi(temp_coord);
                         }
                         else{
                                if(coord_counter==1){
                                  right=atoi(temp_coord);
                                }
                                else{
                                  if(coord_counter==2){
                                         rel_left=atoi(temp_coord);
                                  }
                                  else{
                                         if(coord_counter==3)
                                                rel_right=atoi(temp_coord);
                                         else
                                                tmp_polyA=atoi(temp_coord);
                                  }
                                }
                         }
                         coord_counter++;
                  }
                  k++;
                }

                fscanf(in, "%s\n", temp_string);

                temp_exons=Insert_exon_into_a_exon_list(temp_exons, left, right, rel_left, rel_right, tmp_polyA, temp_string, &incr);
                number_of_exons+=incr;

                trs[i].tr_from[j]=left;
                trs[i].tr_to[j]=right;

                strcpy(trs[i].temp_sequences[j], temp_string);
         }
  }

  //CORREZIONE PER JOB 288 (gene TBCC) - se ogni trascritto ha un solo esone, non va bene... :(((
  /*if(number_of_transcripts != 0){
	  	  int p=0;
	  	  while(p < number_of_transcripts && trs[p].exons == 1){
	  		  p++;
	  	  }
	  	  if(p == number_of_transcripts){
	  		fprintf(stderr, "The strand is not retrieved!\n");
	  		#ifdef HALT_EXIT_MODE
	  		         exit(1);
	  		#else
	  		         exit(EXIT_FAILURE);
	  		#endif
	  	  }else{
	  		  if(trs[p].tr_from[0] > trs[p].tr_from[trs[p].exons-1])
	  			  strand=-1;
	  		  else
	  			  strand=1;
	  	  }
  }*/
  char *temp = (char *) malloc(255*sizeof(char)); temp[0] = '\0';
  sprintf(temp,"%sgenomic.txt",out_path);
  FILE *in_gen=fopen(temp, "r");
  if(in_gen == NULL){
         fprintf(stderr, "Error genomic file!\n");
         exit(0);
  }
  
  char *tmp_line= NULL;
  size_t tmp_string_l = 0;
  int bytes_read;
  bytes_read=my_getline(&tmp_line, &tmp_string_l, in_gen);
  if (bytes_read == -1) {
		 DEBUG("Empty genomic file!");
		 strand=1;
  }
  else{
	 char *occurrence=strrchr(tmp_line, ':');
	 if(occurrence == NULL){
	 	DEBUG("Strand not found in genomic file!");
	 	strand=1;
	 }else{
	 	strand=atoi(occurrence+1);
		DEBUG("Genomic strand: %d", strand);
	 }
  }
  if (tmp_line!=NULL) {
	 pfree(tmp_line);
  }

  //fscanf(in_gen, ">chr%*d:%*d:%*d:%d\n", &strand);
  fclose(in_gen);
  free(temp);

  exons=(struct exon *)malloc(number_of_exons*sizeof(struct exon));
  if(exons == NULL){
         fprintf(stderr, "Problem12 of memory allocation in Get_Transcripts_from_File_FASTA_format!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  counter=0;
  head=temp_exons;
  while(head != NULL){
         exons[counter].left=head->left;
         exons[counter].right=head->right;
         exons[counter].rel_left=head->rel_left;
         exons[counter].rel_right=head->rel_right;
         exons[counter].polyA=head->polyA;
         exons[counter].sequence=(char *)malloc((strlen(head->sequence)+1)*sizeof(char));
         if(exons[counter].sequence == NULL){
                fprintf(stderr, "Memory problem2 in the Get_Transcripts... procedure!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }
         strcpy(exons[counter].sequence, head->sequence);
         counter++;
         head=head->next;
  }

  head=temp_exons;
  while(head != NULL){
         help=head->next;
         free(head->sequence);
         free(head);
         head=help;
  }

  if(strand == -1){
         for(i=0; i<number_of_transcripts; i++){
                k=0;
                p=trs[i].exons-1;
                while(k < trs[i].exons){
                  left=trs[i].tr_from[p];
                  right=trs[i].tr_to[p];

                  z=0;
                  stop=0;
                  while(z < number_of_exons && !stop){
                         if(exons[z].left == left && exons[z].right == right && !strcmp(exons[z].sequence, trs[i].temp_sequences[p]))
                                stop=1;
                         else
                                z++;
                  }
                  if(stop)
                         trs[i].exon_index[k]=z;
                  else{
                         fprintf(stderr, "Problem in exon in Get_Transcripts_from_File_FASTA_format!\n");
#ifdef HALT_EXIT_MODE
                         exit(1);
#else
                         exit(EXIT_FAILURE);
#endif
                  }

                  if(exons[z].polyA == 1)
                         is_int=1;
                  else
                         is_int=(p != 0 && p != trs[i].exons-1)?(1):((p == 0)?(3):(2));

                  exons[z].is_int=is_int;

                  k++;
                  p--;
                }
         }
  }
  else{
         for(i=0; i<number_of_transcripts; i++){
                k=0;
                while(k < trs[i].exons){
                  left=trs[i].tr_from[k];
                  right=trs[i].tr_to[k];

                  z=0;
                  stop=0;
                  while(z < number_of_exons && !stop){
                         if(exons[z].left == left && exons[z].right == right && !strcmp(exons[z].sequence, trs[i].temp_sequences[k]))
                                stop=1;
                         else
                                z++;
                  }
                  if(stop)
                         trs[i].exon_index[k]=z;
                  else{
                         fprintf(stderr, "Problem in exon in Get_Transcripts_from_File_FASTA_format!\n");
#ifdef HALT_EXIT_MODE
                         exit(1);
#else
                         exit(EXIT_FAILURE);
#endif
                  }

                  if(exons[z].polyA == 1)
                         is_int=1;
                  else
                         is_int=(k != 0 && k != trs[i].exons-1)?(1):((k == 0)?(3):(2));

                  exons[z].is_int=is_int;

                  k++;
                }
         }
  }

  for(i=0; i<number_of_transcripts; i++){
         free(trs[i].tr_from);
         free(trs[i].tr_to);

         for(j=0; j<trs[i].exons; j++)
                free(trs[i].temp_sequences[j]);
         free(trs[i].temp_sequences);
  }

  fclose(in);
}

void GetIntronList(char *fileName){
  FILE *in=NULL;
  //char temp_string[400000];
  const int temp_string_l= 1000000;
  char temp_string[temp_string_l];
  //char temp_string2[10000000];
  char temp_ID[20];
  int conf_EST=0;
  int left=0, right=0;
  int counter=0;
  int i=0, j=0;
  int ID_counter=0;

  in=fopen(fileName, "r");
  if(in == NULL){
         fprintf(stderr, "Error6!\n");
         exit(0);
  }

  number_of_introns=0;
  while(!feof(in)){
          //fgets(temp_string, 400000, in);
          fgets(temp_string, temp_string_l, in);
          if(strcmp(temp_string,  "")){
                  number_of_introns++;
          }
  }

  fclose(in);

  if(number_of_introns == 0)
          return;

  in=fopen(fileName, "r");
  if(in == NULL){
         fprintf(stderr, "Error6!\n");
         exit(0);
  }

  introns=(struct intron *)calloc(number_of_introns,sizeof(struct intron));
  if(introns == NULL){
         fprintf(stderr, "Error7!\n");
         exit(0);
  }

  while(!feof(in)){
         //fgets(temp_string, 400000, in);
         fgets(temp_string, temp_string_l, in);
         temp_string[strlen(temp_string)-1]='\0';

         sscanf(temp_string, "%*d %*d %d %d %*d %d %s %*f %*f %*f %*f %*f %*d %*d %*s %*s %*s %*s %*s %*s\n", &left, &right, &conf_EST, temp_string);
         introns[counter].left=left;
         introns[counter].right=right;
         introns[counter].ESTs=conf_EST;

         introns[counter].IDs=(char **)malloc(conf_EST*sizeof(char *));
         if(introns[counter].IDs == NULL){
                fprintf(stderr, "Error8!\n");
                exit(0);
         }

         i=0;
         ID_counter=0;
         while(i < (int)strlen(temp_string)){
                j=0;
                while(temp_string[i] != ',' && i < (int)strlen(temp_string)){
                  temp_ID[j]=temp_string[i];
                  i++;
                  j++;
                }
                temp_ID[j]='\0';
                if(j > 0){
                  introns[counter].IDs[ID_counter]=(char *)malloc(strlen(temp_ID)+1*sizeof(char));
                  if(introns[counter].IDs[ID_counter] == NULL){
                         fprintf(stderr, "Error9!\n");
                         exit(0);
                  }
                  strcpy(introns[counter].IDs[ID_counter], temp_ID);
                  ID_counter++;
                }
                i++;
         }

         counter++;
  }

  fclose(in);
}

void GetCDSAnnotations(char *fileName){
  FILE *in=NULL;
  char temp_ID[20];
  int rel_start=0, rel_end=0;
  int exons=0;
  int counter=0;
  int length=0;

  in=fopen(fileName, "r");
  if(in == NULL)
         fprintf(stderr, "WARNING: CDS annotation %s file does not exist!\n", fileName);
  else
         fprintf(stderr, "CDS for RefSeq obtained from %s!\n", fileName);

  if(in == NULL)
         number_of_cds=0;
  else
         fscanf(in, "%d\n", &number_of_cds);

  if(in == NULL)
         return;

   a_cds=(struct annotated_cds *)malloc(number_of_cds*sizeof(struct annotated_cds));
  if(a_cds == NULL){
         fprintf(stderr, "Error1 in memory allocation GetCDSAnnotations!\n");
         exit(0);
  }

  while(!feof(in)){
	 fscanf(in, "%d\n", &length);

	 if(length > 0){
		a_cds[counter].RefSeq_sequence=(char *)malloc((length+1)*sizeof(char));
		if (a_cds == NULL) {
		  fprintf(stderr, "Error2 in memory allocation GetCDSAnnotations!\n");
		  exit(0);
		}

		DEBUG("Reading RefSeq no. %d...", counter+1);
		int fscanf_res= fscanf(in, "%s\t%d\t%d\t%d\t%s\n",
									  temp_ID, &rel_start, &rel_end, &exons,
									  a_cds[counter].RefSeq_sequence);
		my_assert(fscanf_res == 5);
		strcpy(a_cds[counter].RefSeq, temp_ID);
		a_cds[counter].rel_start=rel_start;
		a_cds[counter].rel_end=rel_end;
		a_cds[counter].exons=exons;

		counter++;
	 } else {
		fprintf(stderr, "WARNING: CDS annotation %s file not correct!\n", fileName);
		fscanf(in, "%s\t%d\t%d\t%d\n", temp_ID, &rel_start, &rel_end, &exons);
		fprintf(stderr, "\tRefSeq %s has null length!\n", temp_ID);
	 }
  }

  fclose(in);
}

void CheckStartEndWRTref(int ref, int i){
  if(ref != -1){
         trs[i].start_cons=0;
         trs[i].end_cons=0;

         if(trs[ref].abs_ORF_start != -1 && trs[ref].abs_ORF_end != -1){
                if(i == ref){
                  trs[i].start_cons=1;
                  trs[i].end_cons=1;
                }
                else{
                  if(trs[i].abs_ORF_start == trs[ref].abs_ORF_start){
                         if(strand == 1)
                                trs[i].start_cons=1;
                         else
                                trs[i].end_cons=1;
                  }
                  if(trs[i].abs_ORF_end == trs[ref].abs_ORF_end){
                         if(strand == 1)
                                trs[i].end_cons=1;
                         else
                                trs[i].start_cons=1;
                  }
                }
         }
  }
}

char GetCDSAnnotationForRefSeq_2(int i){
  char stop=0, found=0;
  char *tr_seq=NULL;
  int j=0, z=0, p=0, k=0, r_index;
  int tmp_ORF_start=0, tmp_ORF_end=0;
  int length=0;
  int start_align_index=0, end_align_index=0;
  char ORF_found=0;
  char *EST_temp, *GEN_temp;
  int cfr_length=0;

  if(trs[i].type != 0)
         return 0;

  while(j < number_of_cds && !stop){
         if(!strcmp(a_cds[j].RefSeq, trs[i].RefSeq))
                stop=1;
         else
                j++;
  }

  if(!stop)
         return 0;

  trs[i].ORF_start=-1;
  trs[i].ORF_end=-1;

  tr_seq=(char *)malloc((trs[i].length+1)*sizeof(char));

  if(tr_seq == NULL){
         fprintf(stderr, "Memory problem in GetCDSAnnotationForRefSeq\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  strcpy(tr_seq, "");

  if(strand == 1){
         for(j=0; j<trs[i].exons; j++){
                strcat(tr_seq, exons[trs[i].exon_index[j]].sequence);
         }
  }
  else{
         for(j=trs[i].exons-1; j>=0; j--){
                strcat(tr_seq, exons[trs[i].exon_index[j]].sequence);
         }
  }

  trs[i].no_ATG=0;
  ORF_found=1;

//Confrontare la sequenza trs[i] con quella annotata per ricavare trs[i].ORF_start e trs[i].ORF_end
//ed effettuare i controlli di congruenza, altrimenti return 0
//Effettuare il controllo sul multiplo di 3, altrimenti return 0

//Trovo il refseq nel file delle annotazioni
  r_index=0;
  stop=0;
  while(r_index < number_of_cds && !stop){
         if(!strcmp(a_cds[r_index].RefSeq, trs[i].RefSeq))
                stop=1;
         else
                r_index++;
  }
  if(!stop)
         return 0;

 //Trovo l'ORF sul trascritto in questione che sia uguale a quello annotato
  z=0;
  while(z < (int)strlen(tr_seq) && !found){
         found=0;
	//Non necessariamente il primo codine dell'annotazione deve essere ATG (issue #31)
         //if(Check_start_codon(z, tr_seq)){
                k=a_cds[r_index].rel_start-1;
                p=z;
                stop=0;
                while(k < a_cds[r_index].rel_end && p < (int)strlen(tr_seq) && !stop){
                  if(tolower(a_cds[r_index].RefSeq_sequence[k]) != tolower(tr_seq[p]))
                         stop=1;
                  else{
                         k++;
                         p++;
                  }
                }
                if(!stop && k == a_cds[r_index].rel_end)
                  found=1;
         //}
         if(!found)
                z++;
  }

  if(!found)
         return 0;

  trs[i].ORF_start=z+1;
  trs[i].ORF_end=p;

  if((trs[i].ORF_end-trs[i].ORF_start+1) %3 != 0)
         return 0;

  if(Tcds > trs[i].ORF_end-trs[i].ORF_start+1)
         Tcds=trs[i].ORF_end-trs[i].ORF_start+1;

  for(z=0; z<3; z++){
         trs[i].start_c[z]=tr_seq[trs[i].ORF_start+z-1];
  }
  trs[i].start_c[z]='\0';

  //Controllo presenza ATG (issue #31)
   if(!(!strcmp(trs[i].start_c, "atg") || !strcmp(trs[i].start_c, "ATG")))
         trs[i].no_ATG=1;

  for(z=0; z<3; z++){
         trs[i].stop_c[z]=tr_seq[trs[i].ORF_end+z-3];
  }
  trs[i].stop_c[z]='\0';
  if((!strcmp(trs[i].stop_c, "tga") || !strcmp(trs[i].stop_c, "TGA")) || (!strcmp(trs[i].stop_c, "tag") || !strcmp(trs[i].stop_c, "TAG")) || (!strcmp(trs[i].stop_c, "taa") || !strcmp(trs[i].stop_c, "TAA")))
         trs[i].has_stop=1;

  if(strand == -1){
         tmp_ORF_start=trs[i].length-trs[i].ORF_end+1;
         tmp_ORF_end=trs[i].length-trs[i].ORF_start+1;
  }
  else{
         tmp_ORF_start=trs[i].ORF_start;
         tmp_ORF_end=trs[i].ORF_end;
  }

  p=0;
  length=0;
  stop=0;
  while(p < trs[i].exons && !stop){
         cfr_length=strlen(exons[trs[i].exon_index[p]].sequence);

         if(tmp_ORF_start <= length+cfr_length){
                trs[i].first_ORF_index=p;
                stop=1;
         }
         else{
                length+=strlen(exons[trs[i].exon_index[p]].sequence);

                p++;
         }
  }

  EST_temp=trs[i].EST_exon_alignments[trs[i].first_ORF_index];
  GEN_temp=trs[i].GEN_exon_alignments[trs[i].first_ORF_index];

  if(strand == 1){
         k=0;
         start_align_index=0;
         while(k < tmp_ORF_start-length){
                if(EST_temp[start_align_index] != '-')
                  k++;
                start_align_index++;
         }
         start_align_index--;

         k=0;
         while(start_align_index >= 0){
                if(GEN_temp[start_align_index] != '-')
                  k++;
                start_align_index--;
         }
         trs[i].abs_ORF_start=exons[trs[i].exon_index[trs[i].first_ORF_index]].left+k-1;
  }
  else{
         k=0;
         start_align_index=strlen(trs[i].EST_exon_alignments[trs[i].first_ORF_index])-1;
         while(k < tmp_ORF_start-length){
                if(EST_temp[start_align_index] != '-')
                  k++;
                start_align_index--;
         }
         start_align_index++;
         k=0;
         while(start_align_index < (int)strlen(trs[i].GEN_exon_alignments[trs[i].first_ORF_index])){
                if(GEN_temp[start_align_index] != '-')
                  k++;
                start_align_index++;
         }
         trs[i].abs_ORF_start=exons[trs[i].exon_index[trs[i].first_ORF_index]].left+k-1;
  }

  length=0;
  p=0;
  stop=0;
  while(p < trs[i].exons && !stop){
         cfr_length=strlen(exons[trs[i].exon_index[p]].sequence);

         if(tmp_ORF_end <= length+cfr_length){
                trs[i].second_ORF_index=p;
                stop=1;
         }
         else{
                length+=strlen(exons[trs[i].exon_index[p]].sequence);

                p++;
         }
  }

  EST_temp=trs[i].EST_exon_alignments[trs[i].second_ORF_index];
  GEN_temp=trs[i].GEN_exon_alignments[trs[i].second_ORF_index];
  if(strand == 1){
         k=0;
         end_align_index=0;
         while(k < tmp_ORF_end-length){
                if(EST_temp[end_align_index] != '-')
                  k++;
                end_align_index++;
         }
         end_align_index--;

         k=0;
         while(end_align_index >= 0){
                if(GEN_temp[end_align_index] != '-')
                  k++;
                end_align_index--;
         }
         trs[i].abs_ORF_end=exons[trs[i].exon_index[trs[i].second_ORF_index]].left+k-1;
  }
  else{
         k=0;
         end_align_index=strlen(trs[i].EST_exon_alignments[trs[i].second_ORF_index])-1;
         while(k < tmp_ORF_end-length){
                if(EST_temp[end_align_index] != '-')
                  k++;
                end_align_index--;
         }
         end_align_index++;
         k=0;
         while(end_align_index < (int)strlen(trs[i].GEN_exon_alignments[trs[i].second_ORF_index])){
                if(GEN_temp[end_align_index] != '-')
                  k++;
                end_align_index++;
         }
         trs[i].abs_ORF_end=exons[trs[i].exon_index[trs[i].second_ORF_index]].left+k-1;
  }

  free(tr_seq);

  return 1;
}

void MarkIntronType(){
  int i=0, j=0;
  int refseqnum=0;

  for(i=0; i<number_of_introns; i++){
         j=0;
         refseqnum=0;
         while(j<introns[i].ESTs){
                if(introns[i].IDs[j][0] == 'N' && introns[i].IDs[j][1] == 'M'){
                  introns[i].type=0;
                  refseqnum++;
                }
                j++;
         }

         introns[i].RefSeqNum=refseqnum;

         if(refseqnum == 0){
                if(introns[i].ESTs > 1)
                  introns[i].type=1;
                else
                  introns[i].type=2;
         }
         else{
                introns[i].RefSeq=(char **)malloc(refseqnum*sizeof(char *));
                if(introns[i].RefSeq == NULL){
                  fprintf(stderr, "Error in MarkIntronType!\n");
                  exit(0);
                }
                j=0;
                refseqnum=0;
                while(j<introns[i].ESTs){
                  if(introns[i].IDs[j][0] == 'N' && introns[i].IDs[j][1] == 'M'){
                         introns[i].RefSeq[refseqnum]=(char *)malloc((strlen(introns[i].IDs[j])+1)*sizeof(char));
                         if(introns[i].RefSeq[refseqnum] == NULL){
                                fprintf(stderr, "Error in MarkIntronType!\n");
                                exit(0);
                         }
                         strcpy(introns[i].RefSeq[refseqnum], introns[i].IDs[j]);
                         refseqnum++;
                  }
                  j++;
                }
         }
  }
}

void MarkTranscriptType(struct transcript *tr){
  int k=0;
  int intron_left=0, intron_right=0;
  char stop=0;
  char conf_at_least2=1;

  if(tr->type == -1){
         intron_left=exons[tr->exon_index[0]].right+1;
         intron_right=exons[tr->exon_index[1]].left-1;

         k=0;
         stop=0;
         while(k<number_of_introns && !stop){
                if(intron_left == introns[k].left && intron_right == introns[k].right){
                  stop=1;
                  if(introns[k].type != 1)
                         conf_at_least2=0;
                }
                else
                  k++;
         }

         if(conf_at_least2)
                tr->type=1;
         else
                tr->type=2;
  }
}

void PrintTABOutput(int ref, struct cds cds_for_gene){
  int i=0, order=0;
  int prev_cursor=0, cursor=0;

  char comp_label[10000], new_label[10000];
  char ir_label[10000];
  int ORF_start=0, ORF_end=0;
  char isINF=0;
  char start_cons=0, end_cons=0;
  int print_counter=0;

#ifdef PRINT_GTF_VARIANT
  FILE *gtf=NULL;
  char *temp = (char *) malloc(255*sizeof(char)); temp[0] = '\0';
  sprintf(temp,"%sVariantGTF.txt",out_path);
  gtf=fopen(temp, "w");
  free(temp);
#endif

  fprintf(stderr, "Transcript");
  fprintf(stderr, "\tExons");
  fprintf(stderr, "\tL (nt)");
  fprintf(stderr, "\tCDS");
  fprintf(stderr, "\tReference CDS start/stop");
  fprintf(stderr, "\tProt L (aa)");
#ifdef PRINT_FRAME
  fprintf(stderr, "\tReference Frame");
#endif
  fprintf(stderr, "\tVariant Type\n");

  if(ref != -1){
         list_of_new_labels=(struct new_label_for_exon **)malloc((trs[ref].exons+2)*sizeof(struct new_label_for_exon *));
         if(list_of_new_labels == NULL){
                fprintf(stderr, "Memory problem in the PrintHTMLOutput!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }
         for(i=0; i<=trs[ref].exons+1; i++){
                list_of_new_labels[i]=NULL;
         }
  }

  order=0;

  while(order < number_of_transcripts){
         prev_cursor=0;
         cursor=0;

         i=order_index[order];

         print_counter++;

         if(i == ref)
                fprintf(stderr, "1\t%s.Ref", gene);
         else
                fprintf(stderr, "0\t%s.tr%d", gene, print_counter);

#ifdef PRINT_GTF_VARIANT
         fprintf(gtf, "variant_isoform#%d", print_counter);
#endif

         fprintf(stderr, "\t%d", trs[i].exons);

#ifdef PRINT_GTF_VARIANT
         fprintf(gtf, " /nex=%d", trs[i].exons);
#endif

         ORF_start=trs[i].ORF_start;
         ORF_end=trs[i].ORF_end;
         start_cons=trs[i].start_cons;
         end_cons=trs[i].end_cons;

         fprintf(stderr, "\t%d", trs[i].length);

#ifdef PRINT_GTF_VARIANT
         fprintf(gtf, " /L=%d", trs[i].length);
#endif

         if(ORF_start != -1 && ORF_end != -1){
                fprintf(stderr, "\t%s%d..%d%s", (trs[i].no_ATG)?("<"):(""), ORF_start, ORF_end, (trs[i].has_stop)?(""):(">"));

#ifdef PRINT_GTF_VARIANT
                fprintf(gtf, " /CDS=%s%d..%d%s", (trs[i].no_ATG)?("<"):(""), ORF_start, ORF_end, (trs[i].has_stop)?(""):(">"));
#endif
         }
         else{
                fprintf(stderr, "\t..");
#ifdef PRINT_GTF_VARIANT
                fprintf(gtf, " /CDS=..");
#endif
         }

         if(i == ref){
                fprintf(stderr, "\tyes/yes");

#ifdef PRINT_GTF_VARIANT
                fprintf(gtf, " /RefSeq=%s", trs[i].RefSeq);
#endif
         }
         else{
                if(ORF_start != -1 && ORF_end != -1){
                  if(trs[i].no_ATG && !trs[i].has_stop)
                         fprintf(stderr, "\t..");
                  else
                         fprintf(stderr, "\t%s/%s", (trs[i].no_ATG == 0)?((start_cons == 1)?("yes"):("no")):(".."), (trs[i].has_stop)?((end_cons == 1)?("yes"):("no")):(".."));
                }
                else{
                  fprintf(stderr, "\t..");
                }

#ifdef PRINT_GTF_VARIANT
                if(!trs[i].has_stop)
                  fprintf(gtf, " /RefSeq=%s", ((trs[i].RefSeq == NULL)?(""):(trs[i].RefSeq)));
                else
                  fprintf(gtf, " /RefSeq=%s(%s%s)", ((trs[i].RefSeq == NULL)?(""):(trs[i].RefSeq)), (start_cons == 1)?("Y"):("N"), (end_cons == 1)?("Y"):("N"));
#endif
         }

         if(ORF_start != -1 && ORF_end != -1){
                fprintf(stderr, "\t%s%d", (trs[i].no_ATG == 1 || !trs[i].has_stop)?(">"):(""), (ORF_end-ORF_start+1)/3-1);
#ifdef PRINT_GTF_VARIANT
                fprintf(gtf, " /ProtL=%s%d", (trs[i].no_ATG == 1 || !trs[i].has_stop)?(">"):(""), (ORF_end-ORF_start+1)/3-1);
#endif
         }
         else{
                fprintf(stderr, "\t..");
#ifdef PRINT_GTF_VARIANT
                fprintf(gtf, " /ProtL=..");
#endif
         }

         if(i != ref){
#ifdef PRINT_FRAME
                if(!trs[i].has_stop){
                  fprintf(stderr, "\t..");
#ifdef PRINT_GTF_VARIANT
                  fprintf(gtf, " /Frame=..");
#endif
                }
                else{
                  isINF=isInFrame(i, ref);
                  if(isINF == 0){
                         fprintf(stderr, "\tno");
#ifdef PRINT_GTF_VARIANT
                         fprintf(gtf, " /Frame=no");
#endif
                  }
                  else{
                         fprintf(stderr, "\tyes");
#ifdef PRINT_GTF_VARIANT
                         fprintf(gtf, " /Frame=yes");
#endif
                  }
                }
#endif
         }

         if(i == ref){
#ifdef PRINT_FRAME
                fprintf(stderr, "\tyes");
#endif
                if(trs[i].RefSeq == NULL)
                  fprintf(stderr, "\tReference TR");
                else
                  fprintf(stderr, "\t%s (Reference TR)", trs[i].RefSeq);
#ifdef PRINT_GTF_VARIANT
                fprintf(gtf, " /Type=Ref");
                if(print_counter < number_of_transcripts)
                  fprintf(gtf, "\n");
#endif
         } else {
			  getCompetingLabels(i, ref, comp_label);
			  getnewIRLabels(i, ref, ir_label);

			  getEXInitTermSkipNewLabels(i, ref, new_label);

			  if (trs[i].RefSeq == NULL)
				 fprintf(stderr, "\t%s%s%s", comp_label, ir_label, new_label);
			  else
				 fprintf(stderr, "\t%s (%s%s%s)", trs[i].RefSeq, comp_label, ir_label, new_label);

#ifdef PRINT_GTF_VARIANT
                fprintf(gtf, " /Type=%s%s%s", comp_label, ir_label, new_label);
                if(print_counter < number_of_transcripts)
                  fprintf(gtf, "\n");
#endif
         }

         fprintf(stderr, "\n");

         order++;
  }

#ifdef PRINT_GTF_VARIANT
  fclose(gtf);
#endif
}

void Set_cover_exon(){
  char stop=0;
  int i=0, j=0;
  int k=0;

  while(i<number_of_exons){
         exons[i].covered_exon=0;
         exons[i].cover_index=-1;
         exons[i].variant_label_first=0;
         exons[i].variant_label_second=0;
         i++;
  }

  i=0;

  while(i < number_of_exons){
         if(exons[i].covered_exon == 0){
                j=i+1;
                stop=0;

                while(j < number_of_exons && !stop){
                  if(exons[i].left >= exons[j].left && exons[i].right <= exons[j].right && !stop){
                         exons[i].covered_exon=1;
                         exons[i].cover_index=j;
                         stop=1;
                  }
                  else{
                         if(exons[j].left >= exons[i].left && exons[j].right <= exons[i].right){
                                exons[j].covered_exon=1;
                                exons[j].cover_index=i;
                         }
                  }
                  j++;
                }
         }
         i++;
  }

  i=0;
  while(i < number_of_exons){
         if(exons[i].covered_exon){
                k=i;

                do{
                  j=exons[k].cover_index;
                  k=j;

                }while(exons[j].covered_exon != 0);

                exons[i].cover_index=j;
         }
         i++;
  }

  i=0;
  while(i<number_of_exons){
         if(exons[i].covered_exon == 0){
                number_of_cover_exons++;
         }
         i++;
  }
}

void SetAltSplMatrixIndexes(){
  int i=0;
  int index=0;

  while(i<number_of_exons){
         if(exons[i].covered_exon == 0){
                exons[i].matrix_index=index;
                index++;
         }
         i++;
  }

  i=0;
  while(i<number_of_exons){
         if(exons[i].covered_exon == 1){
                if(exons[exons[i].cover_index].covered_exon != 0){
                  fprintf(stderr, "Problem2!\n");
#ifdef HALT_EXIT_MODE
                  exit(1);
#else
                  exit(EXIT_FAILURE);
#endif
                }
                exons[i].matrix_index=exons[exons[i].cover_index].matrix_index;
         }
         i++;
  }

  if(index != number_of_cover_exons){
         fprintf(stderr, "Problem3!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }
}

void getnewIRLabels(int index, int ref, char *label){
  int i=0, j=0;
  char ir_label[20];
  i=0;
  char localize_str[20];

  strcpy(label, "");

  if(ref == -1)
         return;

//IR+ introne in RefSeq
  i=0;
  j=0;
  while(i < trs[index].exons){
         while(j < trs[ref].exons-1 && exons[trs[ref].exon_index[j]].right < exons[trs[index].exon_index[i]].left){
                j++;
         }
         while(j < trs[ref].exons-1 && exons[trs[ref].exon_index[j]].right <= exons[trs[index].exon_index[i]].right){
                if(exons[trs[ref].exon_index[j+1]].left <= exons[trs[index].exon_index[i]].right){

                  if(strand == 1)
                         sprintf(ir_label, "IR+(I%d),", j+1);
                  else
                         sprintf(ir_label, "IR+(I%d),", trs[ref].exons-j-1);

                  GetLocalization(localize_str, ref, j);
                  strcat(label, ir_label);
                  strcat(label, localize_str);

                  strcat(label, "; ");

                }
                j++;
         }
         i++;
  }

//IR- introne in trascritto
  i=0;
  j=0;
  while(i < trs[ref].exons){
         while(j < trs[index].exons-1 && exons[trs[index].exon_index[j]].right < exons[trs[ref].exon_index[i]].left){
                j++;
         }
         while(j < trs[index].exons-1 && exons[trs[index].exon_index[j]].right <= exons[trs[ref].exon_index[i]].right){
                if(exons[trs[index].exon_index[j+1]].left <= exons[trs[ref].exon_index[i]].right){

                  if(strand == 1)
                         sprintf(ir_label, "IR-(E%d),", i+1);
                  else
                         sprintf(ir_label, "IR-(E%d),", trs[ref].exons-i);

                  GetLocalization(localize_str, ref, i);
                  strcat(label, ir_label);
                  strcat(label, localize_str);

                  strcat(label, "; ");
                }
                j++;
         }
         i++;
  }
}

void getCompetingLabels(int index, int ref, char *label){
  int i=0, j=0;
  char add[20];
  char comp_label[20];
  char overlap=1;

  i=0;
  strcpy(label, "");

  if(ref == -1)
         return;

  i=0;
  while(i < trs[index].exons-1){
         j=0;
         while(j < trs[ref].exons && exons[trs[index].exon_index[i]].matrix_index != exons[trs[ref].exon_index[j]].matrix_index){
                j++;
         }

         do{
                if(j+1 < trs[ref].exons && exons[trs[index].exon_index[i+1]].matrix_index == exons[trs[ref].exon_index[j+1]].matrix_index){
                  overlap=1;
                  if(exons[trs[index].exon_index[i]].left > exons[trs[ref].exon_index[j]].right || exons[trs[index].exon_index[i]].right < exons[trs[ref].exon_index[j]].left)
                         overlap=0;
                  if(exons[trs[index].exon_index[i+1]].left > exons[trs[ref].exon_index[j+1]].right || exons[trs[index].exon_index[i+1]].right < exons[trs[ref].exon_index[j+1]].left)
                         overlap=0;

                  if(exons[trs[index].exon_index[i]].right != exons[trs[ref].exon_index[j]].right && overlap){
                         if(strand == 1)
                                sprintf(comp_label, "A5E (I%d, ", j+1);
                         else
                                sprintf(comp_label, "A3E (I%d, ", trs[ref].exons-j-1);

                         strcat(label, comp_label);
                         sprintf(add, "%s%d nt), ", (exons[trs[ref].exon_index[j]].right-exons[trs[index].exon_index[i]].right < 0)?(""):("+"), exons[trs[ref].exon_index[j]].right-exons[trs[index].exon_index[i]].right);
                         strcat(label, add);

                         if(exons[trs[ref].exon_index[j]].pos_flag_to == 0){
                                if(exons[trs[index].exon_index[i]].pos_flag_to == 0)
                                  strcat(label, "CDS");
                                else{
                                  if(exons[trs[index].exon_index[i]].pos_flag_to == 1){
                                         if(strand == 1)
                                                strcat(label, "5UTR_CDS");
                                         else
                                                strcat(label, "CDS_3UTR");
                                  }
                                  else{
                                         if(strand == 1)
                                                strcat(label, "CDS_3UTR");
                                         else
                                                strcat(label, "5UTR_CDS");
                                  }
                                }
                         }
                         else{
                                if(exons[trs[ref].exon_index[j]].pos_flag_to == 1){
                                  if(exons[trs[index].exon_index[i]].pos_flag_to == 1){
                                         if(strand == 1)
                                                strcat(label, "5UTR");
                                         else
                                                strcat(label, "3UTR");
                                  }
                                  else{
                                         if(exons[trs[index].exon_index[i]].pos_flag_to == 0){
                                                if(strand == 1)
                                                  strcat(label, "5UTR_CDS");
                                                else
                                                  strcat(label, "CDS_3UTR");
                                         }
                                         else{
                                                strcat(label, "5UTR_3UTR");
                                         }
                                  }
                                }
                                else{
                                  if(exons[trs[index].exon_index[i]].pos_flag_to == 2){
                                         if(strand == 1)
                                                strcat(label, "3UTR");
                                         else
                                                strcat(label, "5UTR");
                                  }
                                  else{
                                         if(exons[trs[index].exon_index[i]].pos_flag_to == 0){
                                                if(strand == 1)
                                                  strcat(label, "CDS_3UTR");
                                                else
                                                  strcat(label, "5UTR_CDS");
                                         }
                                         else{
                                                strcat(label, "5UTR_3UTR");
                                         }
                                  }
                                }
                         }

                         strcat(label, "; ");
                  }

                  if(exons[trs[index].exon_index[i+1]].left != exons[trs[ref].exon_index[j+1]].left && overlap){
                         if(strand == 1)
                                sprintf(comp_label, "A3E (I%d, ", j+1);
                         else
                                sprintf(comp_label, "A5E (I%d, ", trs[ref].exons-j-1);

                         strcat(label, comp_label);

                         sprintf(add, "%s%d nt), ", (exons[trs[index].exon_index[i+1]].left-exons[trs[ref].exon_index[j+1]].left < 0)?(""):("+"), exons[trs[index].exon_index[i+1]].left-exons[trs[ref].exon_index[j+1]].left);
                         strcat(label, add);

                         if(exons[trs[ref].exon_index[j+1]].pos_flag_from == 0){
                                if(exons[trs[index].exon_index[i+1]].pos_flag_from == 0)
                                  strcat(label, "CDS");
                                else{
                                  if(exons[trs[index].exon_index[i+1]].pos_flag_from == 1){
                                         if(strand == 1)
                                                strcat(label, "5UTR_CDS");
                                         else
                                                strcat(label, "CDS_3UTR");
                                  }
                                  else{
                                         strcat(label, "CDS");
                                  }
                                }
                         }
                         else{
                                if(exons[trs[ref].exon_index[j+1]].pos_flag_from == 1){
                                  if(exons[trs[index].exon_index[i+1]].pos_flag_from == 1){
                                         if(strand == 1)
                                                strcat(label, "5UTR");
                                         else
                                                strcat(label, "3UTR");
                                  }
                                  else{
                                         if(exons[trs[index].exon_index[i+1]].pos_flag_from == 0){
                                                if(strand == 1)
                                                  strcat(label, "5UTR_CDS");
                                                else
                                                  strcat(label, "CDS_3UTR");
                                         }
                                         else{
                                                strcat(label, "5UTR_3UTR");
                                         }
                                  }
                                }
                                else{
                                  if(exons[trs[index].exon_index[i+1]].pos_flag_from == 2){
                                         if(strand == 1)
                                                strcat(label, "3UTR");
                                         else
                                                strcat(label, "5UTR");
                                  }
                                  else{
                                         if(exons[trs[index].exon_index[i+1]].pos_flag_from == 0){
                                                if(strand == 1)
                                                  strcat(label, "CDS_3UTR");
                                                else
                                                  strcat(label, "5UTR_CDS");
                                         }
                                         else{
                                                strcat(label, "5UTR_3UTR");
                                         }
                                  }
                                }
                         }

                         strcat(label, "; ");
                  }
                }
                j++;
         }while(j < trs[ref].exons && exons[trs[index].exon_index[i]].matrix_index == exons[trs[ref].exon_index[j]].matrix_index);

         i++;
  }
}

void getEXInitTermSkipNewLabels(int index, int ref, char *label){
  int i=0, j=0, k=0, q=0;

  int p=0;

  char add[20];
  char extr_variant=1;

  //char label_char=0;
  char *representation;

  char localize_str[20];
  char first_time=1;

  char stop=0;

  strcpy(label, "");

  if(ref == -1)
         return;

//INIT per strand 1, altrimenti TERM
  if(exons[trs[ref].exon_index[0]].right == exons[trs[index].exon_index[0]].right){
         if(exons[trs[ref].exon_index[0]].left == exons[trs[index].exon_index[0]].left){
                extr_variant=0;
         }
         else{
                if(exons[trs[ref].exon_index[0]].left > exons[trs[index].exon_index[0]].left){
                  if(exons[trs[ref].exon_index[0]].polyA != 1 || exons[trs[ref].exon_index[0]].left-exons[trs[index].exon_index[0]].left <= 20)
                         extr_variant=0;
                }
                else{
                  if(exons[trs[index].exon_index[0]].polyA != 1 || exons[trs[index].exon_index[0]].left-exons[trs[ref].exon_index[0]].left <= 20)
                         extr_variant=0;
                }
         }
  }

//Solo per variante INIT
  if(extr_variant == 1 && exons[trs[index].exon_index[0]].polyA != 1){
         stop=0;
         p=1;
         while(p < trs[ref].exons && !stop){
                if(exons[trs[ref].exon_index[p]].left == exons[trs[index].exon_index[0]].left && exons[trs[ref].exon_index[p]].right == exons[trs[index].exon_index[0]].right)
                  stop=1;
                else
                  p++;
         }
         if(stop)
                extr_variant=0;
  }

//TERM e' sempre riferito all'ultimo esone
  i=1;
  if(extr_variant == 1){
         GetLocalization(localize_str, ref, 0);

         int r_index=1;
         if(exons[trs[index].exon_index[0]].left < exons[trs[ref].exon_index[0]].left){
                r_index=0;
         }
         list_of_new_labels[r_index]=Insert_newlabel_into_a_newlabel_list(list_of_new_labels[r_index], exons[trs[index].exon_index[0]].left, exons[trs[index].exon_index[0]].right, &representation);

         if(strand == 1){
                sprintf(label, "init(E%d%s),", r_index, representation);
         }
         else{
                if(r_index == 1)
                  sprintf(label, "term(E%d%s),", trs[ref].exons, representation);
                else
                  sprintf(label, "term(%da%s),", trs[ref].exons, representation);
         }

         strcat(label, localize_str);

         strcat(label, "; ");

         first_time=0;

         while(i < trs[index].exons && exons[trs[index].exon_index[i]].right < exons[trs[ref].exon_index[0]].left){
                list_of_new_labels[0]=Insert_newlabel_into_a_newlabel_list(list_of_new_labels[0], exons[trs[index].exon_index[i]].left, exons[trs[index].exon_index[i]].right, &representation);

                if(strand == 1){
                  sprintf(add, "init(E0%s),", representation);
                }
                else{
                  sprintf(add, "term(%da%s),", trs[ref].exons, representation);
                }

                strcat(label, add);
                strcat(label, localize_str);
                strcat(label, "; ");

                i++;
         }
//Da i si parte per la variante NEW
  }

  extr_variant=1;
//TERM per strand 1, altrimenti INIT
  if(exons[trs[ref].exon_index[trs[ref].exons-1]].left == exons[trs[index].exon_index[trs[index].exons-1]].left){
         if(exons[trs[ref].exon_index[trs[ref].exons-1]].right == exons[trs[index].exon_index[trs[index].exons-1]].right){
                extr_variant=0;
         }
         else{
                if(exons[trs[ref].exon_index[trs[ref].exons-1]].right < exons[trs[index].exon_index[trs[index].exons-1]].right){
                  if(exons[trs[ref].exon_index[trs[ref].exons-1]].polyA != 1 || exons[trs[index].exon_index[trs[index].exons-1]].right-exons[trs[ref].exon_index[trs[ref].exons-1]].right <= 20)
                         extr_variant=0;
                }
                else{
                  if(exons[trs[index].exon_index[trs[index].exons-1]].polyA != 1 || exons[trs[ref].exon_index[trs[ref].exons-1]].right-exons[trs[index].exon_index[trs[index].exons-1]].right <= 20)
                         extr_variant=0;
                }
         }
  }


//sostituzione di i con p
  if(extr_variant == 1 && exons[trs[index].exon_index[trs[index].exons-1]].polyA != 1){
         stop=0;
         p=trs[ref].exons-2;
         while(p >= 0 && !stop){
                if(exons[trs[ref].exon_index[p]].left == exons[trs[index].exon_index[trs[index].exons-1]].left && exons[trs[ref].exon_index[p]].right == exons[trs[index].exon_index[trs[index].exons-1]].right)
                  stop=1;
                else
                  p--;
         }
         if(stop)
                extr_variant=0;
  }

  j=trs[index].exons-2;
  if(extr_variant == 1){
         GetLocalization(localize_str, ref, trs[ref].exons-1);

         int r_index=trs[ref].exons;
         if(exons[trs[index].exon_index[trs[index].exons-1]].right > exons[trs[ref].exon_index[trs[ref].exons-1]].right){
                r_index=trs[ref].exons+1;
         }
         list_of_new_labels[r_index]=Insert_newlabel_into_a_newlabel_list(list_of_new_labels[r_index], exons[trs[index].exon_index[0]].left, exons[trs[index].exon_index[0]].right, &representation);

         if(strand == 1){
                if(r_index == trs[ref].exons)
                  sprintf(add, "term(E%d%s),", trs[ref].exons, representation);
                else
                  sprintf(add, "term(%da%s),", trs[ref].exons, representation);
         }
         else{
                sprintf(add, "init(E%d%s),", (trs[ref].exons-r_index+1), representation);
         }

         strcat(label, add);

         strcat(label, localize_str);
         strcat(label, "; ");

         while(j >= 0 && exons[trs[index].exon_index[j]].left > exons[trs[ref].exon_index[trs[ref].exons-1]].right){

                list_of_new_labels[trs[ref].exons+1]=Insert_newlabel_into_a_newlabel_list(list_of_new_labels[trs[ref].exons+1], exons[trs[index].exon_index[j]].left, exons[trs[index].exon_index[j]].right, &representation);

                if(strand == 1){
                  sprintf(add, "term(%da%s),", trs[ref].exons, representation);
                }
                else{
                  sprintf(add, "init(E0%s),", representation);
                }

                strcat(label, add);
                strcat(label, localize_str);

                strcat(label, "; ");

                j--;
         }
//In j si finisce per la variante NEW
  }

  q=0;
  k=i;
  while(k <= j){
         while(q < trs[ref].exons && exons[trs[ref].exon_index[q]].right < exons[trs[index].exon_index[k]].left)
                q++;

         if(q < trs[ref].exons && exons[trs[ref].exon_index[q]].left > exons[trs[index].exon_index[k]].right){
                GetLocalization(localize_str, ref, q-1);
                list_of_new_labels[q-1]=Insert_newlabel_into_a_newlabel_list(list_of_new_labels[q-1], exons[trs[index].exon_index[k]].left, exons[trs[index].exon_index[k]].right, &representation);

                sprintf(add, "new(E%d%s),", (strand == 1)?(q):(trs[ref].exons-q), representation);

                strcat(label, add);
                strcat(label, localize_str);

                strcat(label, "; ");
         }
         k++;
  }

  i=1;
  while(i < trs[ref].exons-1 && exons[trs[ref].exon_index[i]].left <= exons[trs[index].exon_index[0]].right){
         i++;
  }

  q=0;
  while(i < trs[ref].exons-1){
         while(q < trs[index].exons && exons[trs[index].exon_index[q]].right < exons[trs[ref].exon_index[i]].left)
                q++;

         if(q < trs[index].exons && exons[trs[index].exon_index[q]].left > exons[trs[ref].exon_index[i]].right){
                GetLocalization(localize_str, ref, i);

                sprintf(add, "skip(E%d),", (strand == 1)?(i+1):(trs[ref].exons-i));

                strcat(label, add);
                strcat(label, localize_str);

                strcat(label, "; ");
         }
         i++;
  }
}

void GetLocalization(char *local, int index, int exon){

  strcpy(local, "");

  if(exons[trs[index].exon_index[exon]].pos_flag_from == 1){
         if(exons[trs[index].exon_index[exon]].pos_flag_to == 1){
                if(strand == 1)
                  strcat(local, "5UTR");
                else
                  strcat(local, "3UTR");
         }
         else{
                if(exons[trs[index].exon_index[exon]].pos_flag_to == 0){
                  if(strand == 1)
                         strcat(local, "5UTR_CDS");
                  else
                         strcat(local, "CDS_3UTR");
                }
                else{
                  strcat(local, "5UTR_3UTR");
                }
         }
  }
  else{
         if(exons[trs[index].exon_index[exon]].pos_flag_from == 2){
                if(strand == 1)
                  strcat(local, "3UTR");
                else
                  strcat(local, "5UTR");
         }
         else{
                if(exons[trs[index].exon_index[exon]].pos_flag_to == 0){
                  strcat(local, "CDS");
                }
                else{
                  if(strand == 1)
                         strcat(local, "CDS_3UTR");
                  else
                         strcat(local, "5UTR_CDS");
                }
         }
  }
}

void GetLongestORF(int ref, int i, int min_length){
  char *tr_seq= NULL;
  int j=0, p=0, k=0;
  int tmp_ORF_start=0, tmp_ORF_end=0;
  int length=0;
  char stop=0;
  int start_align_index=0, end_align_index=0;

  char *EST_temp, *GEN_temp;
  int cfr_length=0;

  DEBUG("Looking for an ORF for transcript %d (%dbp long)", i, trs[i].length);

  tr_seq=(char *)malloc((trs[i].length+1)*sizeof(char));
  exit_with_problem_if(tr_seq == NULL, "Memory problem");
  tr_seq[0]= '\0';
  if (strand == 1) {
	 for (j=0; j<trs[i].exons; j++) {
		strcat(tr_seq, exons[trs[i].exon_index[j]].sequence);
	 }
  } else {
	 for (j=trs[i].exons-1; j>=0; j--) {
		strcat(tr_seq, exons[trs[i].exon_index[j]].sequence);
	 }
  }

  DEBUG("Transcript sequence: %s", tr_seq);

  trs[i].has_stop= 0;
  trs[i].no_ATG= 0;
  trs[i].ORF_start= -1;
  trs[i].ORF_end= -1;

  int z= 0;
  const int ccds_end= strlen(tr_seq)-3;
  bool ORF_found= false;
  int ORF_length= 0;

  for (int frame= 0; frame<3; ++frame) {
	 DEBUG("Considering frame %d...", frame);
	 z= frame;

	 while (z <= ccds_end) {

		if (Check_start_codon(z, tr_seq)) {
		  DEBUG("  Found a start codon at position %d.", z);
		  j= z+3;
		  while (j <= ccds_end && !Check_stop_codon(j, tr_seq)) {
			 j+=3;
		  }

		  if (j <= ccds_end) {
			 const int this_ORF_length= j-z+3;
			 DEBUG("  Found a in-frame stop codon at position %d (ORF length=%d)",
					 j, this_ORF_length);
			 if (this_ORF_length >= min_length) {
// An ORF is annotated if (see issue #20):
// - it has a context and the previous one has not
// OR
// - ( * it is longer than the previous AND
//     * it has a context if the previous one has it )
				const bool has_context= getContext(z, tr_seq)>0;
				DEBUG("    New ORF found (start=%d, end=%d, length=%d, ORF_found=%s).",
						z+1, j+3, this_ORF_length, has_context?"YES":"NO");
				if ( (!ORF_found && has_context) ||
					  ( (this_ORF_length > ORF_length) &&
						 (!ORF_found || has_context) ) ) {
				  DEBUG("    The new ORF is better than the previous one. Saved.");
				  ORF_length= this_ORF_length;
				  trs[i].ORF_start= z+1;
				  trs[i].ORF_end= j+3;
				  ORF_found= has_context;
				} else {
				  DEBUG("    The new ORF is NOT better than the previous one. Discarded.");
				}
			 } else {
				DEBUG("    The ORF is shorter than the minimum allowed length. Discarded.");
			 }
		  } else {
			 DEBUG("  No in-frame stop codon found.");
		  }
		  z= j+3;
		} else {
		  z+= 3;
		}
	 }
  }

  if(trs[i].ORF_start != -1 && trs[i].ORF_end != -1){
         for(z=0; z<3; z++){
                trs[i].start_c[z]=tr_seq[trs[i].ORF_start+z-1];
         }
         trs[i].start_c[z]='\0';

         for(z=0; z<3; z++){
                trs[i].stop_c[z]=tr_seq[trs[i].ORF_end+z-3];
         }
         trs[i].stop_c[z]='\0';

         if((!strcmp(trs[i].stop_c, "tga") || !strcmp(trs[i].stop_c, "TGA")) || (!strcmp(trs[i].stop_c, "tag") || !strcmp(trs[i].stop_c, "TAG")) || (!strcmp(trs[i].stop_c, "taa") || !strcmp(trs[i].stop_c, "TAA")))
                trs[i].has_stop=1;
         else{
                fprintf(stderr, "Stop problem\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }

         if(!trs[i].has_stop){
                trs[i].ORF_end=trs[i].length;
         }

         if(strand == -1){
                tmp_ORF_start=trs[i].length-trs[i].ORF_end+1;
                tmp_ORF_end=trs[i].length-trs[i].ORF_start+1;
         }
         else{
                tmp_ORF_start=trs[i].ORF_start;
                tmp_ORF_end=trs[i].ORF_end;
         }

         p=0;
         length=0;
         stop=0;
         while(p < trs[i].exons && !stop){
                if(trs[i].type == 0)
                  cfr_length=strlen(exons[trs[i].exon_index[p]].sequence);
                else
                  cfr_length=exons[trs[i].exon_index[p]].right-exons[trs[i].exon_index[p]].left+1;

                if(tmp_ORF_start <= length+cfr_length){
                  trs[i].first_ORF_index=p;
                  stop=1;
                }
                else{
                  if(trs[i].type == 0)
                         length+=strlen(exons[trs[i].exon_index[p]].sequence);
                  else
                         length+=exons[trs[i].exon_index[p]].right-exons[trs[i].exon_index[p]].left+1;
                  p++;
                }
         }

         if(trs[i].type == 0){
                EST_temp=trs[i].EST_exon_alignments[trs[i].first_ORF_index];
                GEN_temp=trs[i].GEN_exon_alignments[trs[i].first_ORF_index];

                if(strand == 1){
                  k=0;
                  start_align_index=0;
                  while(k < tmp_ORF_start-length){
                         if(EST_temp[start_align_index] != '-')
                                k++;
                         start_align_index++;
                  }
                  start_align_index--;

                  k=0;
                  while(start_align_index >= 0){
                         if(GEN_temp[start_align_index] != '-')
                                k++;
                         start_align_index--;
                  }
                  trs[i].abs_ORF_start=exons[trs[i].exon_index[trs[i].first_ORF_index]].left+k-1;
                }
                else{
                  k=0;
                  start_align_index=strlen(trs[i].EST_exon_alignments[trs[i].first_ORF_index])-1;
                  while(k < tmp_ORF_start-length){
                         if(EST_temp[start_align_index] != '-')
                                k++;
                         start_align_index--;
                  }
                  start_align_index++;
                  k=0;
                  while(start_align_index < (int)strlen(trs[i].GEN_exon_alignments[trs[i].first_ORF_index])){
                         if(GEN_temp[start_align_index] != '-')
                                k++;
                         start_align_index++;
                  }
                  trs[i].abs_ORF_start=exons[trs[i].exon_index[trs[i].first_ORF_index]].left+k-1;
                }
         }
         else{
                trs[i].abs_ORF_start=tmp_ORF_start-length+exons[trs[i].exon_index[trs[i].first_ORF_index]].left-1;
         }

         length=0;
         p=0;
         stop=0;
         while(p < trs[i].exons && !stop){
                if(trs[i].type == 0)
                  cfr_length=strlen(exons[trs[i].exon_index[p]].sequence);
                else
                  cfr_length=exons[trs[i].exon_index[p]].right-exons[trs[i].exon_index[p]].left+1;

                if(tmp_ORF_end <= length+cfr_length){
                  trs[i].second_ORF_index=p;
                  stop=1;
                }
                else{
                  if(trs[i].type == 0)
                         length+=strlen(exons[trs[i].exon_index[p]].sequence);
                  else
                         length+=exons[trs[i].exon_index[p]].right-exons[trs[i].exon_index[p]].left+1;

                  p++;
                }
         }

         if(trs[i].type == 0){
                EST_temp=trs[i].EST_exon_alignments[trs[i].second_ORF_index];
                GEN_temp=trs[i].GEN_exon_alignments[trs[i].second_ORF_index];
                if(strand == 1){
                  k=0;
                  end_align_index=0;
                  while(k < tmp_ORF_end-length){
                         if(EST_temp[end_align_index] != '-')
                                k++;
                         end_align_index++;
                  }
                  end_align_index--;

                  k=0;
                  while(end_align_index >= 0){
                         if(GEN_temp[end_align_index] != '-')
                                k++;
                         end_align_index--;
                  }
                  trs[i].abs_ORF_end=exons[trs[i].exon_index[trs[i].second_ORF_index]].left+k-1;
                }
                else{
                  k=0;
                  end_align_index=strlen(trs[i].EST_exon_alignments[trs[i].second_ORF_index])-1;
                  while(k < tmp_ORF_end-length){
                         if(EST_temp[end_align_index] != '-')
                                k++;
                         end_align_index--;
                  }
                  end_align_index++;
                  k=0;
                  while(end_align_index < (int)strlen(trs[i].GEN_exon_alignments[trs[i].second_ORF_index])){
                         if(GEN_temp[end_align_index] != '-')
                                k++;
                         end_align_index++;
                  }
                  trs[i].abs_ORF_end=exons[trs[i].exon_index[trs[i].second_ORF_index]].left+k-1;
                }
         }
         else{
                trs[i].abs_ORF_end=tmp_ORF_end-length+exons[trs[i].exon_index[trs[i].second_ORF_index]].left-1;
         }
  }
  else{
         trs[i].abs_ORF_start=-1;
         trs[i].first_ORF_index=-1;
         trs[i].abs_ORF_end=-1;
         trs[i].second_ORF_index=-1;
  }

  free(tr_seq);
}

char getContext(int ATG_start, char *tr_seq){

  if (!Check_start_codon(ATG_start, tr_seq)) {
	 exit_with_problem("ATG problem");
  }

  char context=2; //0=weak, 1=medium, 2=strong

  if ( (ATG_start-3 < 0) ||
		 ( (tr_seq[ATG_start-3] != 'a' && tr_seq[ATG_start-3] != 'A') &&
			(tr_seq[ATG_start-3] != 'g' && tr_seq[ATG_start-3] != 'G') ) ) {
	 context--;
  }

  if ( (ATG_start+3 >= (int)strlen(tr_seq)) ||
		 ( (tr_seq[ATG_start+3] != 'a' && tr_seq[ATG_start+3] != 'A') &&
			(tr_seq[ATG_start+3] != 'g' && tr_seq[ATG_start+3] != 'G') ) ) {
	 context--;
  }

  return context;
}

bool Check_start_codon(int pos, char *tr_seq){
  return ((tr_seq[pos] == 'a' && tr_seq[pos+1] == 't' && tr_seq[pos+2] == 'g') ||
			 (tr_seq[pos] == 'A' && tr_seq[pos+1] == 'T' && tr_seq[pos+2] == 'G'));
}

char Check_stop_codon(int pos, char *tr_seq){
  return ( ( (tr_seq[pos] == 't' && tr_seq[pos+1] == 'a' && tr_seq[pos+2] == 'a') ||
				 (tr_seq[pos] == 'T' && tr_seq[pos+1] == 'A' && tr_seq[pos+2] == 'A') ) ||
			  ( (tr_seq[pos] == 't' && tr_seq[pos+1] == 'a' && tr_seq[pos+2] == 'g') ||
				 (tr_seq[pos] == 'T' && tr_seq[pos+1] == 'A' && tr_seq[pos+2] == 'G') ) ||
			  ( (tr_seq[pos] == 't' && tr_seq[pos+1] == 'g' && tr_seq[pos+2] == 'a') ||
				 (tr_seq[pos] == 'T' && tr_seq[pos+1] == 'G' && tr_seq[pos+2] == 'A') ) );
}

struct exon *Insert_exon_into_a_exon_list(struct exon *arg_exon_list, int left, int right, int rel_left, int rel_right, char polyA, char *sequence, int *incr){
  struct exon *head=arg_exon_list;
  struct exon *y=NULL;
  char add=0;

  char found=0;

  *incr=0;

  head=arg_exon_list;
  while(head!=NULL && !(left <= head->left)){
         y=head;
         head=head->next;
  }

  if(head != NULL && left == head->left){

  while(head != NULL && (left == head->left && right < head->right)){
                y=head;
                head=head->next;
         }
         if(head != NULL && (left == head->left && right == head->right)){
                found=0;
                while(head != NULL && (left == head->left && right == head->right) && !found){
                  if(!strcmp(sequence, head->sequence))
                         found=1;
                  else{
                         y=head;
                         head=head->next;
                  }
                }

                if(!found){
                  add=1;
                }
         }
         else{
                add=1;
         }
  }
  else{
         add=1;
  }

  if(add){
         *incr=1;
         computed_exon_copy=(struct exon *)malloc(sizeof(struct exon));
         if(computed_exon_copy == NULL){
                fprintf(stderr, "Memory problem in the Insert_exon_into_a_exon_list procedure!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }

         computed_exon_copy->left=left;
         computed_exon_copy->right=right;
         computed_exon_copy->rel_left=rel_left;
         computed_exon_copy->rel_right=rel_right;
         computed_exon_copy->polyA=polyA;
         computed_exon_copy->sequence=(char *)malloc((strlen(sequence)+1)*sizeof(char));
         if(computed_exon_copy->sequence == NULL){
                fprintf(stderr, "Memory problem2 in the Insert_exon_into_a_exon_list procedure!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }
         strcpy(computed_exon_copy->sequence, sequence);

         computed_exon_copy->next=head;
         computed_exon_copy->prev=y;

         if(head!=NULL){
                head->prev=computed_exon_copy;
         }

         if(y != NULL){
                y->next=computed_exon_copy;
         }
         else{
                arg_exon_list=computed_exon_copy;
         }
  }

  return arg_exon_list;
}

struct genomic_exon *Insert_genexon_into_a_genexon_list(struct genomic_exon *arg_genexon_list, int rel_left, int rel_right, char *sequence){
  struct genomic_exon *head=arg_genexon_list;
  struct genomic_exon *y=NULL;
  char add=0;

  head=arg_genexon_list;
  while(head!=NULL && !(rel_left <= head->rel_left)){
         y=head;
         head=head->next;
  }

  if(head != NULL && rel_left == head->rel_left){

         while(head != NULL && (rel_left == head->rel_left && rel_right < head->rel_right)){
                y=head;
                head=head->next;
         }
         if(head != NULL && (rel_left == head->rel_left && rel_right == head->rel_right)){
                add=0;
         }
         else{
                add=1;
         }
  }
  else{
         add=1;
  }

  if(add){
         gen_exon_copy=(struct genomic_exon *)malloc(sizeof(struct genomic_exon));
         if(gen_exon_copy == NULL){
                fprintf(stderr, "Memory problem in the Insert_genexon_into_a_genexon_list procedure!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }

         gen_exon_copy->rel_left=rel_left;
         gen_exon_copy->rel_right=rel_right;
         gen_exon_copy->sequence=(char *)malloc((strlen(sequence)+1)*sizeof(char));
         if(gen_exon_copy->sequence == NULL){
                fprintf(stderr, "Memory problem2 in the Insert_genexon_into_a_genexon_list procedure!\n");
#ifdef HALT_EXIT_MODE
                exit(1);
#else
                exit(EXIT_FAILURE);
#endif
         }
         strcpy(gen_exon_copy->sequence, sequence);

         gen_exon_copy->next=head;
         gen_exon_copy->prev=y;

         if(head!=NULL){
                head->prev=gen_exon_copy;
         }

         if(y != NULL){
                y->next=gen_exon_copy;
         }
         else{
                arg_genexon_list=gen_exon_copy;
         }
  }

  return arg_genexon_list;
}

//Inizio e fine degli ORF sono rispetto alla sequenza del tr e non in senso assoluto
char isInFrame(int index, int ref){
  char stop=0;
  int i=0, j=0;
  int ref_left=0, ref_right=0, left=0, right=0;
  int region_left=0, region_right=0, region_length=0, tr_length=0;
  int ref_partial_length=0, partial_length=0;
  char phase1=0, phase2=0;
  int f_cds_i=-1, s_cds_i=-1;

  if(ref == -1)
         return 2;

  if(trs[index].abs_ORF_start == -1){
         return 0;
  }

  if(trs[index].no_ATG || !trs[index].has_stop){
         return 0;
  }

  if(trs[ref].abs_ORF_end < trs[index].abs_ORF_start || trs[index].abs_ORF_end < trs[ref].abs_ORF_start){
         return 0;
  }

  f_cds_i=trs[ref].first_ORF_index;
  s_cds_i=trs[ref].second_ORF_index;

  region_length=0;

  if(strand == -1){
         i=s_cds_i;
         stop=0;
         while(i >= f_cds_i && !stop){
                ref_left=(i == s_cds_i)?(exons[trs[ref].exon_index[i]].rel_left+(exons[trs[ref].exon_index[i]].right-trs[ref].abs_ORF_end)):(exons[trs[ref].exon_index[i]].rel_left);
                ref_right=(i == f_cds_i)?(exons[trs[ref].exon_index[i]].rel_right-(trs[ref].abs_ORF_start-exons[trs[ref].exon_index[i]].left)):(exons[trs[ref].exon_index[i]].rel_right);

                j=trs[index].second_ORF_index;
                left=(j == trs[index].second_ORF_index)?(exons[trs[index].exon_index[j]].rel_left+(exons[trs[index].exon_index[j]].right-trs[index].abs_ORF_end)):(exons[trs[index].exon_index[j]].rel_left);
                right=(j == trs[index].first_ORF_index)?(exons[trs[index].exon_index[j]].rel_right-(trs[index].abs_ORF_start-exons[trs[index].exon_index[j]].left)):(exons[trs[index].exon_index[j]].rel_right);

                partial_length=0;

                while(j >= trs[index].first_ORF_index && left <= ref_right && !stop){

                  if(right >= ref_left){
                         region_left=(left >= ref_left)?(left):(ref_left);
                         region_right=(right <= ref_right)?(right):(ref_right);
                         region_length+=region_right-region_left+1;

                         phase1=(region_left-ref_left+ref_partial_length)%3;
                         phase2=(region_left-left+partial_length)%3;

                         if(phase1 != phase2){
                                stop=1;
                         }
                  }

                  if(!stop){
                         partial_length+=right-left+1;
                         j--;
                         if(j >= 0){
                                left=(j == trs[index].second_ORF_index)?(exons[trs[index].exon_index[j]].rel_left+(exons[trs[index].exon_index[j]].right-trs[index].abs_ORF_end)):(exons[trs[index].exon_index[j]].rel_left);
                                right=(j == trs[index].first_ORF_index)?(exons[trs[index].exon_index[j]].rel_right-(trs[index].abs_ORF_start-exons[trs[index].exon_index[j]].left)):(exons[trs[index].exon_index[j]].rel_right);
                         }
                  }
                }

                if(!stop){
                  ref_partial_length+=ref_right-ref_left+1;
                  i--;
                }
         }
  }
  else{
         i=f_cds_i;
         stop=0;
         while(i <= s_cds_i && !stop){
                ref_left=(i == f_cds_i)?(trs[ref].abs_ORF_start):(exons[trs[ref].exon_index[i]].left);
                ref_right=(i == s_cds_i)?(trs[ref].abs_ORF_end):(exons[trs[ref].exon_index[i]].right);

                j=trs[index].first_ORF_index;
                left=(j == trs[index].first_ORF_index)?(trs[index].abs_ORF_start):(exons[trs[index].exon_index[j]].left);
                right=(j == trs[index].second_ORF_index)?(trs[index].abs_ORF_end):(exons[trs[index].exon_index[j]].right);
                partial_length=0;

                while(j <= trs[index].second_ORF_index && left <= ref_right && !stop){
                  if(right >= ref_left){
                         region_left=(left >= ref_left)?(left):(ref_left);
                         region_right=(right <= ref_right)?(right):(ref_right);
                         region_length+=region_right-region_left+1;

                         phase1=(region_left-ref_left+ref_partial_length)%3;
                         phase2=(region_left-left+partial_length)%3;

                         if(phase1 != phase2){
                                stop=1;
                         }
                  }

                  if(!stop){
                         partial_length+=right-left+1;
                         j++;
                         if(j < trs[index].exons){
                                left=(j == trs[index].first_ORF_index)?(trs[index].abs_ORF_start):(exons[trs[index].exon_index[j]].left);
                                right=(j == trs[index].second_ORF_index)?(trs[index].abs_ORF_end):(exons[trs[index].exon_index[j]].right);
                         }
                  }
                }

                if(!stop){
                  ref_partial_length+=ref_right-ref_left+1;
                  i++;
                }
         }
  }

  i=trs[index].first_ORF_index;
  tr_length=0;

  while(i <= trs[index].second_ORF_index){
         left=(i == trs[index].first_ORF_index)?(trs[index].abs_ORF_start):(exons[trs[index].exon_index[i]].left);
         right=(i == trs[index].second_ORF_index)?(trs[index].abs_ORF_end):(exons[trs[index].exon_index[i]].right);
         tr_length+=right-left+1;
         i++;
  }

  if((double)(region_length*100/tr_length) < 50.0)
         return 0;

  if(!stop)
         return 1;
  else
         return 0;
}

void PrintOutputFile(int ref){
  FILE *out=NULL;
  int i=0, j=0, order=0;
  char one_color=1;
  int left=0, right=0;
  int first_UTR_length=0, second_UTR_length=0;
  int print_counter=0;

  char *temp = (char *) malloc(255*sizeof(char)); temp[0] = '\0';
  sprintf(temp,"%sCCDS_transcripts.txt",out_path);
  out=fopen(temp, "w");
  free(temp);

  if(out == NULL){
         fprintf(stderr, "File of transcripts not created!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  fprintf(out, "%d\n", number_of_transcripts);
  fprintf(out, "%s\n", gen_length_str);

  order=0;

  while(order < number_of_transcripts){
         print_counter++;

         i=order_index[order];

         fprintf(out, ">%d:%d:%d:%d:", print_counter, trs[i].exons, (i == ref), (trs[i].type == 0)?(1):(0));

         if(!trs[i].has_stop || (trs[i].abs_ORF_start == -1 && trs[i].abs_ORF_end == -1)){
                fprintf(out, "-1\n");
         }
         else{
                if(strand == 1){
                  if(trs[i].second_ORF_index == trs[i].exons-1){
                         fprintf(out, "0\n");
                  }
                  else{
                         if(exons[trs[i].exon_index[trs[i].second_ORF_index]].right-trs[i].abs_ORF_end > 50){
                                fprintf(out, "1\n");
                         }
                         else{
                                fprintf(out, "0\n");
                         }
                  }
                }
                else{
                  if(trs[i].first_ORF_index == 0){
                         fprintf(out, "0\n");
                  }
                  else{
                         if(trs[i].abs_ORF_start-exons[trs[i].exon_index[trs[i].first_ORF_index]].left > 50){
                                fprintf(out, "1\n");
                         }
                         else{
                                fprintf(out, "0\n");
                         }
                  }
                }
         }

         j=0;

         while(j < trs[i].exons){
                left=exons[trs[i].exon_index[j]].left;
                right=exons[trs[i].exon_index[j]].right;

                fprintf(out, "%d:%d:", left, right);
                fprintf(out, "%d:%d:", exons[trs[i].exon_index[j]].rel_left, exons[trs[i].exon_index[j]].rel_right);
                fprintf(out, "%d:", exons[trs[i].exon_index[j]].polyA);

                first_UTR_length=0;
                second_UTR_length=0;
                if(trs[i].abs_ORF_start != -1 && trs[i].abs_ORF_end != -1){
                  one_color=1;
                  if(trs[i].first_ORF_index == j){
                         one_color=0;
                         first_UTR_length=trs[i].abs_ORF_start-left;
                  }
                  if(trs[i].second_ORF_index == j){
                         one_color=0;
                         second_UTR_length=right-trs[i].abs_ORF_end;
                  }

                  if(one_color){
                         if(left > trs[i].abs_ORF_end){
                                second_UTR_length=right-left+1;
                         }
                         else{
                                if(right < trs[i].abs_ORF_start)
                                  first_UTR_length=right-left+1;
                         }
                  }

                  if(strand == 1){
                         fprintf(out, "%d:%d\n", first_UTR_length, second_UTR_length);
                  }
                  else{
                         fprintf(out, "%d:%d\n", second_UTR_length, first_UTR_length);
                  }
                }
                else{
                  fprintf(out, "-1:-1\n");
                }

                fprintf(out, "%s\n", exons[trs[i].exon_index[j]].sequence);

                j++;
         }

         order++;
  }

  fprintf(out, "#\n");

  fclose(out);
}

void SetPrintOrder(int ref){
  int start=0;
  int i=0, pos=0, j=0;
  int help=0;

  order_index=(int *)malloc(number_of_transcripts*sizeof(int));
  if(order_index == NULL){
         fprintf(stderr, "Problem in SetPrintOrder!\n");
#ifdef HALT_EXIT_MODE
         exit(1);
#else
         exit(EXIT_FAILURE);
#endif
  }

  if(ref != -1){
         order_index[0]=ref;
         start=1;
  }
  else
         start=0;

  i=0;
  pos=start;
  while(i<number_of_transcripts){
         if(i != ref){
                order_index[pos]=i;
                pos++;
         }
         i++;
  }

  for(i=start+1; i<number_of_transcripts; i++){
         help=order_index[i];
         j=i-1;
         while(j>=start && trs[help].exons > trs[order_index[j]].exons){
                order_index[j+1]=order_index[j];
                j--;
         }
         order_index[j+1]=help;
  }
}

int SetREFToLongestTranscript(){

  int j=0, k=0;
  char stop=0, first=1;
  int intron_left=0, intron_right=0;
  int product=0;

  int i=0;
  int index=-1;
  int trs_length=0;
  int trs_exons=0;

#ifndef EXON_LONGEST_REF
  int *min_E=NULL;
  min_E=(int *)calloc(number_of_transcripts,sizeof(int));
  if(min_E == NULL){
	 fprintf(stderr, "Memory problem in the SetREFToLongestTranscript procedure!\n");
#ifdef HALT_EXIT_MODE
	 exit(1);
#else
	 exit(EXIT_FAILURE);
#endif
  }
  i=0;

  while(i<number_of_transcripts){

	//Si considerano i soli full-lengths annotati con CDS
	if(trs[i].abs_ORF_start != -1 && trs[i].abs_ORF_end != -1){
	  j=0;
	 first=1;

	 while(j < trs[i].exons-1){
		 //CORREZIONE PER JOB 288 (gene TBCC)
		 if(strand == 1){
			 intron_left=exons[trs[i].exon_index[j]].right+1;
			 intron_right=exons[trs[i].exon_index[j+1]].left-1;
		 }else{
			 intron_right=exons[trs[i].exon_index[j+1]].left-1;
			 intron_left=exons[trs[i].exon_index[j]].right+1;
		 }
		k=0;
		stop=0;
		DEBUG("Looking for intron %d-%d.", intron_left, intron_right);
		while(k < number_of_introns && !stop){
		  DEBUG("  --> considering intron %d-%d.", introns[k].left, introns[k].right);
		  if(intron_left == introns[k].left && intron_right == introns[k].right){
			 stop=1;
		  } else {
			 k++;
		  }
		}

		my_assert(k<number_of_introns);
		if(first){
		  first=0;
		  min_E[i]=introns[k].ESTs;
		}
		else{
		  if(introns[k].ESTs < min_E[i])
			 min_E[i]=introns[k].ESTs;
		}
		j++;
	 }
	}
	 i++;
  }
#endif

//PRIMA SI RICERCA TRA I REFSEQ CHE SONO ANNOTATI nel file cds (vedi issue #32)
i=0;

  trs_length=0;
  trs_exons=0;

#ifndef EXON_LONGEST_REF
  product=0;
#endif

  while(i<number_of_transcripts){
	//Si considerano i soli full-lengths annotati con CDS
	if(trs[i].abs_ORF_start != -1 && trs[i].abs_ORF_end != -1){
#ifdef EXON_LONGEST_REF
	 if((trs[i].type == 0 && trs[i].is_annotated == 1) && (trs[i].exons >= trs_exons && trs[i].length >= trs_length)){
		trs_exons=trs[i].exons;
		trs_length=trs[i].length;
		index=i;
	 }
#else
	 if(trs[i].type == 0 && trs[i].exons*min_E[i] > product) {
		product=trs[i].exons*min_E[i];
		index=i;
	 }
#endif
	}
	 i++;
  }

  if(index != -1){
	 return index;
  }
//FINE DELLA RICERCA TRA I REFSEQ CHE SONO ANNOTATI

  i=0;

  trs_length=0;
  trs_exons=0;

#ifndef EXON_LONGEST_REF
  product=0;
#endif

  while(i<number_of_transcripts){
	//Si considerano i soli full-lengths annotati con CDS
	if(trs[i].abs_ORF_start != -1 && trs[i].abs_ORF_end != -1){
#ifdef EXON_LONGEST_REF
	 if((trs[i].type == 0 && trs[i].is_annotated == 0)  && (trs[i].exons >= trs_exons && trs[i].length >= trs_length)){
		trs_exons=trs[i].exons;
		trs_length=trs[i].length;
		index=i;
	 }
#else
	 if(trs[i].type == 0 && trs[i].exons*min_E[i] > product) {
		product=trs[i].exons*min_E[i];
		index=i;
	 }
#endif
	}
	 i++;
  }

  if(index != -1){
	 return index;
  }

  i=0;

  trs_length=0;
  trs_exons=0;

#ifndef EXON_LONGEST_REF
  product=0;
#endif

  while(i<number_of_transcripts){
	//Si considerano i soli full-lengths annotati con CDS
	if(trs[i].abs_ORF_start != -1 && trs[i].abs_ORF_end != -1){
#ifdef EXON_LONGEST_REF
	 if(trs[i].type == 1 && (trs[i].exons >= trs_exons && trs[i].length >= trs_length)){
		trs_exons=trs[i].exons;
		trs_length=trs[i].length;
		index=i;
	 }
#else
	 if(trs[i].type == 1 && trs[i].exons*min_E[i] > product){
		product=trs[i].exons*min_E[i];
		index=i;
	 }
#endif
	}
	 i++;
  }

  if(index != -1){
	 return index;
  }

  i=0;

  trs_length=0;
  trs_exons=0;

#ifndef EXON_LONGEST_REF
  product=0;
#endif

  while(i<number_of_transcripts){
	//Si considerano i soli full-lengths annotati con CDS
	if(trs[i].abs_ORF_start != -1 && trs[i].abs_ORF_end != -1){
#ifdef EXON_LONGEST_REF
	 if(trs[i].exons >= trs_exons && trs[i].length >= trs_length){
		trs_exons=trs[i].exons;
		trs_length=trs[i].length;
		index=i;
	 }
#else
	 if(trs[i].exons*min_E[i] > product){
		product=trs[i].exons*min_E[i];
		index=i;
	 }
#endif
	}
	 i++;
  }

  //INIZIO 30nov10
  if(index != -1){
  	 return index;
  }

  i=0;
  trs_length=0;
  trs_exons=0;
  int current_type=-1;
  while(i<number_of_transcripts){
	//Si considerano i soli full-lengths annotati con CDS
	if(trs[i].abs_ORF_start != -1 && trs[i].abs_ORF_end != -1){
	 if(current_type != 0){
		 if(trs[i].exons >= trs_exons && trs[i].length >= trs_length){
		   	trs_exons=trs[i].exons;
		   	trs_length=trs[i].length;
		   	current_type=trs[i].type;
		   	index=i;
		 }
	 }
	 else{
		 if(trs[i].type == 0 && (trs[i].exons >= trs_exons && trs[i].length >= trs_length)){
		 	trs_exons=trs[i].exons;
		 	trs_length=trs[i].length;
		 	current_type=trs[i].type;
		 	index=i;
		 }
	 }
	}
  	 i++;
    }
  //FINE 30nov10

  DEBUG("Index %d", index);
  if(index == -1 && number_of_transcripts != 0){
	 fprintf(stderr, "Error!\n");
	 exit(0);
  }

#ifndef EXON_LONGEST_REF
  free(min_E);
#endif

  return index;
}

char GetLongestORFforCCDS(struct cds *cds_for_gene, int i, pmytime pt_tot){
                  int counter=0, j=0;

                  if(trs[i].abs_ORF_start == -1 || trs[i].abs_ORF_end == -1){
                         fprintf(stderr, "ERROR: CCDS not set 2!\n");
                         MYTIME_stop(pt_tot);
                           MYTIME_LOG(INFO, pt_tot);

                           MYTIME_destroy(pt_tot);

                           INFO("End");
                           resource_usage_log();

#ifdef HALT_EXIT_MODE
                         exit(1);
#else
                         exit(EXIT_FAILURE);
#endif
                  }

                  cds_for_gene->exons=trs[i].second_ORF_index-trs[i].first_ORF_index+1;
                  cds_for_gene->cds_from=(int *)malloc(cds_for_gene->exons*sizeof(int));
                  cds_for_gene->cds_to=(int *)malloc(cds_for_gene->exons*sizeof(int));

                  if(cds_for_gene->cds_from == NULL || cds_for_gene->cds_to == NULL){
                         fprintf(stderr, "Error2 here1!\n");
                         exit(0);
                  }

                  counter=0;
                  for(j=trs[i].first_ORF_index; j<=trs[i].second_ORF_index; j++){
                         if(j == trs[i].first_ORF_index)
                                cds_for_gene->cds_from[counter]=trs[i].abs_ORF_start;
                         else
                                cds_for_gene->cds_from[counter]=exons[trs[i].exon_index[j]].left;

                         if(j == trs[i].second_ORF_index)
                                cds_for_gene->cds_to[counter]=trs[i].abs_ORF_end;
                         else
                                cds_for_gene->cds_to[counter]=exons[trs[i].exon_index[j]].right;

                         counter++;
                  }

                  return 1;
                }

              struct new_label_for_exon *Insert_newlabel_into_a_newlabel_list(struct new_label_for_exon *arg_newlabel_list, int left, int right, char **representation){
                  struct new_label_for_exon *head=arg_newlabel_list;
                  struct new_label_for_exon *tail=NULL;
                  unsigned int counter=0;
                  char stop=0;

                  stop=0;
                  while(head != NULL && !stop){
                         if(head->left == left && head->right == right){
                                stop=1;
                         }
                         else{
                                counter++;
                                tail=head;
                                head=head->next;
                         }
                  }

                  if(stop){
	                         *representation=head->representation;
                  }
                  else{
                         computed_newlabel_copy=(struct new_label_for_exon *)malloc(sizeof(struct new_label_for_exon));
                         if(computed_newlabel_copy == NULL){
                                fprintf(stderr, "Memory problem in the Insert_newlabel_into_a_newlabel_list!\n");
#ifdef HALT_EXIT_MODE
                                exit(1);
#else
                                exit(EXIT_FAILURE);
#endif
                         }

                         computed_newlabel_copy->left=left;
                         computed_newlabel_copy->right=right;
                         computed_newlabel_copy->next=NULL;
                         computed_newlabel_copy->prev=tail;

                         //computed_newlabel_copy->character=97+counter;
                         computed_newlabel_copy->representation=int2alpha(counter);

                         //*character=computed_newlabel_copy->character;
                         *representation=computed_newlabel_copy->representation;

                         if(tail == NULL)
                                arg_newlabel_list=computed_newlabel_copy;
                         else
                                tail->next=computed_newlabel_copy;
                  }


                  return arg_newlabel_list;
                }


void ComputeAlignment(char *EST_exon, char *genomic_exon){
  int n=strlen(EST_exon);
  int m=strlen(genomic_exon);

  char **Mdir=NULL;

  int i, j;

  align_dim=0;

  Mdir=(char **)malloc((n+1)*sizeof(char *));
  for(i=0; i<n+1; i++){
         Mdir[i]=(char *)malloc((m+1)*sizeof(char));
         if(Mdir[i] == NULL){
                fprintf(stderr, "Error in ComputeAlignment\n");
                exit(0);
         }
  }

  for(i=0; i<n+1; i++){
         for(j=0; j<m+1; j++){
                Mdir[i][j]=0;
         }
  }

  ComputeAlignMatrix(EST_exon, genomic_exon, n, m, Mdir);
  TracebackAlignment(EST_exon, genomic_exon, Mdir, n, m);

  AlignEST[align_dim]='\0';
  AlignGenomic[align_dim]='\0';

  for(i=0; i<n+1; i++){
         free(Mdir[i]);
  }

  free(Mdir);
}

void ComputeAlignMatrix(char *EST_exon, char *genomic_exon, int n, int m, char **Mdir){
  int **M=NULL;
  int i, j;
  //int i_del;

  M=(int **)malloc((n+1)*sizeof(int *));
  if(M == NULL){
         fprintf(stderr, "Error in ComputeAlignMatrix\n");
         exit(0);
  }

  for(i=0; i<n+1; i++){
         M[i]=(int *)malloc((m+1)*sizeof(int));
         if(M[i] == NULL){
                fprintf(stderr, "Error2 in ComputeAlignMatrix\n");
                exit(0);
         }
  }

  for(i=0; i<n+1; i++){
         for(j=0; j<m+1; j++){
                M[i][j]=0;
         }
  }

//Casi base
  for(i=0; i<n+1; i++){
         M[i][0]=i;
  }
  for(i=0; i<m+1; i++){
         M[0][i]=i;
  }

//Costruzione della matrice M
  for(i=1; i<n+1; i++){
         for(j=1; j<m+1; j++){

                M[i][j]=M[i-1][j-1];

                if(EST_exon[i-1] == genomic_exon[j-1] || EST_exon[i-1] == 'n' || genomic_exon[j-1] == 'n' || EST_exon[i-1] == 'N' || genomic_exon[j-1] == 'N')
                  M[i][j]=M[i][j];//Costo match a 0
                else
                  M[i][j]=M[i][j]+1;    //Costo mismatch a +1

                Mdir[i][j]=0;           //0 per allineamento caratteri in i-1 e j-1


                if(M[i][j] > M[i-1][j]+1){
                  M[i][j]= M[i-1][j]+1; //Costo spazio a +1
                  Mdir[i][j]=1; //1 per cancellazione in genomic_seq; il carattere in i-1 matcha con -
                }

                if(M[i][j] > M[i][j-1]+1){
                  M[i][j]= M[i][j-1]+1; //Costo spazio a +1
                  Mdir[i][j]=2; //2 per cancellazione in EST_seq; il carattere in j-1 matcha con -
                }
         }
  }

  for(i=0; i<n+1; i++){
         free(M[i]);
  }

  free(M);
}

void TracebackAlignment(char *EST_exon, char *genomic_exon, char **Mdir, int i, int j){
  char direction=0;

  if(i > 0 && j > 0){
         direction=Mdir[i][j];

         if (direction == 0){
                TracebackAlignment(EST_exon, genomic_exon, Mdir, i-1, j-1);
                AlignEST[align_dim]=EST_exon[i-1];
                AlignGenomic[align_dim]=genomic_exon[j-1];
                align_dim++;
         }
         else{
                if (direction == 1){
                  TracebackAlignment(EST_exon, genomic_exon, Mdir, i-1, j);
                  AlignEST[align_dim]=EST_exon[i-1];
                  AlignGenomic[align_dim]='-';
                  align_dim++;
                }
                else{
                  TracebackAlignment(EST_exon, genomic_exon, Mdir, i, j-1);
                  AlignEST[align_dim]='-';
                  AlignGenomic[align_dim]=genomic_exon[j-1];
                  align_dim++;
                }
         }
  }
  else{
         if(i > 0){
                TracebackAlignment(EST_exon, genomic_exon, Mdir, i-1, j);
                AlignEST[align_dim]=EST_exon[i-1];
                AlignGenomic[align_dim]='-';
                align_dim++;
         }
         else{
                if(j > 0){
                  TracebackAlignment(EST_exon, genomic_exon, Mdir, i, j-1);
                  AlignEST[align_dim]='-';
                  AlignGenomic[align_dim]=genomic_exon[j-1];
                  align_dim++;
                }
         }
  }
}

void GetExonAlignments(){
  int i=0, j=0;
  int rel_left=0, rel_right=0;

  for(i=0; i<number_of_transcripts; i++){
         if(trs[i].type == 0){
                trs[i].EST_exon_alignments=(char **)malloc(trs[i].exons*sizeof(char *));
                trs[i].GEN_exon_alignments=(char **)malloc(trs[i].exons*sizeof(char *));

                if(trs[i].EST_exon_alignments == NULL || trs[i].GEN_exon_alignments == NULL){
                  fprintf(stderr, "Problem1 in GetExonAlignments!\n");
#ifdef HALT_EXIT_MODE
                  exit(1);
#else
                  exit(EXIT_FAILURE);
#endif
                }

                for(j=0; j<trs[i].exons; j++){
                  rel_left=exons[trs[i].exon_index[j]].rel_left;
                  rel_right=exons[trs[i].exon_index[j]].rel_right;

                  if(strcmp(exons[trs[i].exon_index[j]].sequence, GetGENexonSequence(rel_left, rel_right))){
                         ComputeAlignment(exons[trs[i].exon_index[j]].sequence, GetGENexonSequence(rel_left, rel_right));
                  }
                  else{
                         align_dim=strlen(exons[trs[i].exon_index[j]].sequence);
                         strcpy(AlignEST, exons[trs[i].exon_index[j]].sequence);
                         strcpy(AlignGenomic, exons[trs[i].exon_index[j]].sequence);
                  }
                  trs[i].EST_exon_alignments[j]=(char *)malloc((align_dim+1)*sizeof(char));
                  trs[i].GEN_exon_alignments[j]=(char *)malloc((align_dim+1)*sizeof(char));
                  if(trs[i].EST_exon_alignments[j] == NULL || trs[i].GEN_exon_alignments[j] == NULL){
                         fprintf(stderr, "Problem2 in GetExonAlignments!\n");
#ifdef HALT_EXIT_MODE
                         exit(1);
#else
                         exit(EXIT_FAILURE);
#endif
                  }

                  strcpy(trs[i].EST_exon_alignments[j], AlignEST);
                  strcpy(trs[i].GEN_exon_alignments[j], AlignGenomic);
                }
         }
         else{
                trs[i].EST_exon_alignments=NULL;
                trs[i].GEN_exon_alignments=NULL;
         }
  }
}

void GetGenomicExons(char *fileName){
  FILE *in=NULL;
  char tmp_str[100000];

  int rel_left=0, rel_right=0;

  in=fopen(fileName, "r");
  if(in == NULL){
         fprintf(stderr, "Error11!\n");
         exit(0);
  }

  while(!feof(in)){
         fscanf(in, "%d %d %s\n", &rel_left, &rel_right, tmp_str);
         gen_exons=Insert_genexon_into_a_genexon_list(gen_exons, rel_left, rel_right, tmp_str);
  }

  fclose(in);
}

char *GetGENexonSequence(int rel_left, int rel_right){
  struct genomic_exon *head=gen_exons;
  char found=0;

  while(head!=NULL && !(rel_left <= head->rel_left)){
         head=head->next;
  }

  if(head != NULL && rel_left == head->rel_left){

         while(head != NULL && (rel_left == head->rel_left && rel_right < head->rel_right)){
                head=head->next;
         }
         if(head != NULL && (rel_left == head->rel_left && rel_right == head->rel_right)){
                found=1;
         }
         else{
                found=0;
         }
  }
  else{
         found=0;
  }

  if(found)
         return head->sequence;
  else
         return NULL;
}

char* int2alpha(unsigned int num) {
  unsigned int n_digits= 0;
  unsigned int drift= 0;
  while ((drift+1)*26<=num) {
	 drift= (drift+1)*26;
	 ++n_digits;
  };
  ++n_digits;
  //printf("%6d ", drift);
  char* repr= (char*)calloc(sizeof(char), n_digits+1);
  unsigned int quotient= num-drift;
  unsigned int i= n_digits;
  do {
	 unsigned int remainder= quotient % 26;
	 quotient= quotient / 26;
	 --i;
	 repr[i]= 'a'+remainder;
  } while (i>0);
  return repr;
}
