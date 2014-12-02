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
/**
 *
 * @file types.h
 *
 * Tipi fondamentali condivisi e funzioni di creazione/distruzione.
 *
 **/

#ifndef _TYPES_H_
#define _TYPES_H_

#include "list.h"
#include "util.h"
#include "ext_array.h"
#include "bool_list.h"
#include <limits.h>

typedef struct _pointer* ppointer;

typedef struct _factor* pfactor;

typedef struct _genomic_intron* pgenomic_intron;
typedef struct _genomic_intron_info* pgenomic_intron_info;
typedef struct _burset_frequency* pburset_frequency;
typedef struct _intron* pintron;

typedef struct _subtree_factorizations* psubtree_factorizations;
typedef struct _subtree_embeddings* psubtree_embeddings;

typedef struct _alignment* palignment;
typedef struct _gap_alignment* pgap_alignment;

typedef struct _EST_info* pEST_info;
typedef struct _EST* pEST;
typedef plist pfactorization;
typedef struct _pairing* ppairing;
typedef struct _GEN_ESTS* pGEN_ESTS;
typedef struct _EST_MEG* pEST_MEG;

typedef long int lint;

struct _pointer{
	void *pointer;
};

struct _subtree_factorizations
{
	ppairing root;
	plist factorizations;
};

struct _subtree_embeddings
{
	ppairing root;
	plist embeddings;
};

struct _genomic_intron
{
  int end;
  int start;
  char type;	//0=U12, 1=U2, 2=unclassified
  char *donor_pt;	//pattern in 5'
  char *acceptor_pt;	//pattern in 3'
  int burset_frequency;
  double score5;
  double score3;
  int BPS_position;
  double BPS_score;
  plist info; //List of pgenomic_intro_info
  int supportingESTs;
  bool classified;
  char agree_type;//0: if supported from at least one refseq, 1: from at least one canonical or at-ac U2/U12, 2: other
};

struct _genomic_intron_info
{
	pEST_info info;
	int EST_cut;
};

struct _burset_frequency
{
	int start;	//genomic start
	int end;	//genomic end
	int frequency; //frequency of the burset pattern of the start-end intron
};

struct _intron
{
  pgenomic_intron gen_intron;
  pfactor donor;
  pfactor acceptor;
  bool agreed;
  bool try_agree; //false if the transcript is a RefSeq or is (gt-ag, gc-ag) or is a U12/U2 at-ac
  pEST_info est_info;
  bool isReal;	//false if the intron is not real but is simply a genomic suffix/prefix
  //and donor_EST_start and donor_GEN_start are set to -1 (genomic prefix)
  //or donor_EST_end and donor_GEN_end are set to the gen length (genomic suffix)
  char agree_type; //0: from refseq, 1: canonical or at-ac U2/U12, 2: other
};

struct _factor
{
  int EST_start;
  int EST_end;
  int GEN_start;
  int GEN_end;
};

struct _EST_info
{
  char* original_EST_seq;		//Sequenza della EST 5'3'
  char* EST_seq;				//Sequenza della EST 5'3'
  char* EST_id;					//FASTA header
  char* EST_gb;					//GB id della EST estratto dall'header
										//Unigene (/gb=AA11111)
  char* EST_chr;				//Numero del cromosoma (solo per la genomica)
  char* EST_strand_as_read;
             // Strand estratto dall'header.
             // Per genomica:  "+1/1" per 5'3' e "-1" per 3'5'
             // Per ESTs:      "3" -> concorde con strand genomica,
             //                "5" -> discorde con strand genomica
  int EST_strand;
             // Strand interpretato (o settato di default se mancante dall'header).
             // Per genomica:  1 per 5'3' e -1 per 3'5'
             // Per ESTs:      +1 -> concorde con strand genomica,
             //                -1 -> discorde con strand genomica
  bool fixed_strand;
             // Se true, non cercare di allinearne il reverse & complement

  int abs_start;				//Start rispetto al cromosoma (per genomica)
  int abs_end;					//End rispetto al cromosoma (per genomica)

  int pref_polyA_length;		//Lunghezza del prefisso di polyA tagliato
										//(se -1, non esiste)
  int suff_polyA_length;		//Lunghezza del suffisso di polyA tagliato
										//(se -1, non esiste)
  int pref_polyT_length;		//Lunghezza del prefisso di polyT tagliato
										//(se -1, non esiste)
  int suff_polyT_length;		//Lunghezza del suffisso di polyT tagliato
										//(se -1, non esiste)
  int pref_N_length; // Lunghezza del prefisso di N tagliato
  int suff_N_length; // Lunghezza del suffisso di N tagliato
};

struct _EST
{
  pEST_info info;
  plist factorizations;
  plist bin_factorizations;
  pboollist polyA_signals;
  pboollist polyadenil_signals;
};

struct _pairing {

  unsigned int id;

  int p;
  int t;
  int l;
  plist adjs; //Lista degli adiacenti
  plist incs; //Lista degli incidenti
  bool visited;
  int number_of_visits;
};

/**
 * Definitions of "false" pairings which denote the first
 * and the last pairing (source and sink) of a MEG
 **/

#define SOURCE_PAIRING_LEN (200)
#define SINK_PAIRING_LEN (200)
#define SOURCE_PAIRING_START (INT_MIN)
#define SINK_PAIRING_START ((INT_MAX)-SINK_PAIRING_LEN)


struct _GEN_ESTS
{
  pEST_info gen;
  plist ests;
};

struct _EST_MEG
{
 pEST_info est;
 pext_array meg;
};

struct _alignment
{
  char* EST_alignment;//riga della matrice di allineamento relativa alla EST
  char* GEN_alignment;//riga della matrice di allineamento relativa alla genomica
  int alignment_dim;//Colonne della matrice di allineamento
  int score;
};

struct _gap_alignment
{
  char* EST_gap_alignment;//riga della matrice di allineamento relativa alla EST
  char* GEN_gap_alignment;//riga della matrice di allineamento relativa alla genomica

  int gap_alignment_dim;//Colonne della matrice di allineamento

  int factor_cut;//Left end del fattore di EST in acceptor dopo
				  //allineamento (rispetto alla regione di EST allineata)

  int intron_start;//Left end dell'introne (rispetto alla regione di
					  //genomica allineata)
  int intron_end;//Right end dell'introne (rispetto alla regione di
				  //genomica allineata)

  int intron_start_on_align;//Left end dell'introne nell'allineamento
  int intron_end_on_align;//Right end dell'introne nell'allineamento

  int new_acceptor_factor_left;//Left end del fattore di EST in acceptor dopo
									  //allineamento (rispetto alla EST)
  int new_donor_right_on_gen;//Right end dell'esone in donor dopo
									  //allineamento (rispetto alla genomica)
  int new_acceptor_left_on_gen;//Left end dell'esone in acceptor dopo
									  //allineamento (rispetto alla genomica)
};

ppointer  pointer_create(void);
void pointer_destroy(ppointer);

pEST  EST_create(void);
void EST_destroy(pEST);

palignment alignment_create(size_t);
void alignment_destroy(palignment);

void alignments_destroy(plist);

pgap_alignment gap_alignment_create(size_t);
void gap_alignment_destroy(pgap_alignment);

void gap_alignments_destroy(plist);

/*
 * Non libera la memoria di info
 */
void EST_destroy_just_factorizations(pEST);

pEST_info  EST_info_create(void);
void EST_info_destroy(pEST_info);

pfactor factor_create(void);
void factor_destroy(pfactor);

pgenomic_intron genomic_intron_create(int, int);
void genomic_intron_destroy(pgenomic_intron);

pburset_frequency burset_frequency_create(int, int);
void burset_frequency_destroy(pburset_frequency);

pgenomic_intron_info genomic_intron_info_create(pEST_info, int);
void genomic_intron_info_destroy(pgenomic_intron_info);

pintron intron_create(void);
void intron_destroy(pintron);
//For destroying also the genomic intron
void intron_destroy_free(pintron);
void intron_destroy_free_if_false(pintron);

psubtree_factorizations subtree_factorizations_create(void);
void subtree_factorizations_destroy(psubtree_factorizations);

psubtree_embeddings subtree_embeddings_create(void);
void subtree_embeddings_destroy(psubtree_embeddings);

pfactorization factorization_create(void);
void factorization_destroy(pfactorization);

void factorization_list_destroy(plist);

void intron_factorization_destroy(plist);
void intron_factorization_destroy_free(plist);
void intron_factorization_destroy_free_if_false(plist);

void intron_factorization_list_destroy(plist);
void intron_factorization_list_destroy_free(plist);
void intron_factorization_list_destroy_free_if_false(plist);

plist embedding_create(void);
void embedding_destroy(plist embedding);

ppairing pairing_create(void);
void pairing_destroy(ppairing p);
void pairing_destroy_2(ppairing p);

int pairing_compare(const ppairing* p1, const ppairing* p2);
int genomic_intron_compare(const pgenomic_intron* pgi1, const pgenomic_intron* pgi2);
int burset_frequency_compare(const pburset_frequency* pbf1, const pburset_frequency* pbf2);

ppairing pairing_simple_copy(const ppairing p);

pGEN_ESTS GEN_ESTS_create(void);
void GEN_ESTS_destroy(pGEN_ESTS);

pEST_MEG EST_MEG_create(void);
void EST_MEG_destroy(pEST_MEG);

void vi_destroy(plist vi);

#endif

