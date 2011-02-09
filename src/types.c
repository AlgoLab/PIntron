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
#include "types.h"
#include "list.h"
#include "bit_vector.h"


ppointer  pointer_create(void)
{
  ppointer pp=(ppointer)palloc(sizeof(struct _pointer));
  pp->pointer=NULL;

  return pp;
}

void pointer_destroy(ppointer pp)
{
  pfree(pp);
}

pEST  EST_create(void)
{
  pEST pest=(pEST)palloc(sizeof(struct _EST));
  pest->info=NULL;
  pest->factorizations=NULL;
  pest->bin_factorizations=NULL;
  pest->polyA_signals=NULL;	//List of booleans
  pest->polyadenil_signals=NULL;	//List of booleans

  return pest;
}

void EST_destroy(pEST pest)
{
  if(pest->info!=NULL)EST_info_destroy(pest->info);
  if(pest->factorizations!=NULL) list_destroy(pest->factorizations,(delete_function)factorization_destroy);
  if(pest->bin_factorizations!=NULL) list_destroy(pest->bin_factorizations,(delete_function)BV_destroy);
  if(pest->polyA_signals!=NULL) boollist_destroy(pest->polyA_signals);
  if(pest->polyadenil_signals!=NULL) boollist_destroy(pest->polyadenil_signals);
  pfree(pest);
}

//length: number of columns in the alignment matrix
palignment alignment_create(size_t length)
{
  palignment pal=PALLOC(struct _alignment);
  pal->EST_alignment=c_palloc(length+1);
  pal->GEN_alignment=c_palloc(length+1);
  pal->alignment_dim=0;
  pal->score=-1;

  return pal;
}

//length: number of columns in the alignment matrix
pgap_alignment gap_alignment_create(size_t length)
{
  pgap_alignment pgap=PALLOC(struct _gap_alignment);
  pgap->EST_gap_alignment=c_palloc(length+1);
  pgap->GEN_gap_alignment=c_palloc(length+1);
  pgap->gap_alignment_dim=0;
  pgap->factor_cut=0;
  pgap->intron_start=0;
  pgap->intron_end=0;
  pgap->intron_start_on_align=0;
  pgap->intron_end_on_align=0;
  pgap->new_acceptor_factor_left=0;
  pgap->new_donor_right_on_gen=0;
  pgap->new_donor_right_on_gen=0;

  return pgap;
}

void alignment_destroy(palignment palign)
{
 if (palign!=NULL) {
	 if (palign->EST_alignment!=NULL)
		pfree(palign->EST_alignment);
	 if (palign->GEN_alignment!=NULL)
		pfree(palign->GEN_alignment);
	 pfree(palign);
 }
}

void gap_alignment_destroy(pgap_alignment pg_align)
{
 if (pg_align!=NULL) {
	 if (pg_align->EST_gap_alignment!=NULL)
		pfree(pg_align->EST_gap_alignment);
	 if (pg_align->GEN_gap_alignment!=NULL)
		pfree(pg_align->GEN_gap_alignment);
	 pfree(pg_align);
 }
}

void EST_destroy_just_factorizations(pEST pest)
{
  if(pest->factorizations!=NULL) list_destroy(pest->factorizations,(delete_function)factorization_destroy);
  if(pest->bin_factorizations!=NULL) list_destroy(pest->bin_factorizations,(delete_function)BV_destroy);
  pfree(pest);
}


pEST_info  EST_info_create(void)
{
  pEST_info pesti= PALLOC(struct _EST_info);
  pesti->original_EST_seq= NULL;
  pesti->EST_seq= NULL;
  pesti->EST_id= NULL;
  pesti->EST_gb= NULL;
  pesti->EST_chr= NULL;
  pesti->EST_strand_as_read= NULL;
  pesti->EST_strand= 1;

  return pesti;
}

void EST_info_destroy(pEST_info pest_info)
{
 if (pest_info!=NULL) {
	 if (pest_info->original_EST_seq!=NULL)
		pfree(pest_info->original_EST_seq);
	 if (pest_info->EST_seq!=NULL)
		pfree(pest_info->EST_seq);
	 if (pest_info->EST_id!=NULL)
		pfree(pest_info->EST_id);
	 if (pest_info->EST_gb!=NULL)
		pfree(pest_info->EST_gb);
	 if (pest_info->EST_chr!=NULL)
		pfree(pest_info->EST_chr);
	 if (pest_info->EST_strand_as_read!=NULL)
		pfree(pest_info->EST_strand_as_read);
	 pfree(pest_info);
 }
}

pfactorization factorization_create(void) {
  return (pfactorization)list_create();
}

void factorization_list_destroy(plist factorization_list){
  if(factorization_list != NULL)
	 list_destroy(factorization_list,(delete_function)factorization_destroy);
}

void intron_factorization_list_destroy(plist intron_factorization_list){
  if(intron_factorization_list != NULL)
	 list_destroy(intron_factorization_list,(delete_function)intron_factorization_destroy);
}

void intron_factorization_list_destroy_free(plist intron_factorization_list){
  if(intron_factorization_list != NULL)
	 list_destroy(intron_factorization_list,(delete_function)intron_factorization_destroy_free);
}

void intron_factorization_list_destroy_free_if_false(plist intron_factorization_list){
  if(intron_factorization_list != NULL)
	 list_destroy(intron_factorization_list,(delete_function)intron_factorization_destroy_free_if_false);
}

void factorization_destroy(pfactorization pfact)
{
  if(pfact!=NULL)
	 list_destroy(pfact,(delete_function)factor_destroy);
}

void intron_factorization_destroy(plist pintronfact)
{
  if(pintronfact!=NULL)
	 list_destroy(pintronfact,(delete_function)intron_destroy);
}

void intron_factorization_destroy_free(plist pintronfact)
{
  if(pintronfact!=NULL)
	 list_destroy(pintronfact,(delete_function)intron_destroy_free);
}

void intron_factorization_destroy_free_if_false(plist pintronfact)
{
  if(pintronfact!=NULL)
	 list_destroy(pintronfact,(delete_function)intron_destroy_free_if_false);
}

plist embedding_create(void) {
  return (plist)list_create();
}

void embedding_destroy(plist embedding)
{
  if(embedding!=NULL)
	 list_destroy(embedding,(delete_function)pairing_destroy);
}

void alignments_destroy(plist alignments)
{
  if(alignments!=NULL)
	 list_destroy(alignments,(delete_function)alignment_destroy);
}

void gap_alignments_destroy(plist alignments)
{
  if(alignments!=NULL)
	 list_destroy(alignments,(delete_function)gap_alignment_destroy);
}

pfactor factor_create(void)
{
  pfactor pf=(pfactor)palloc(sizeof(struct _factor));
  return pf;
}

void factor_destroy(pfactor pf)
{
  if(pf!=NULL)pfree(pf);
}

pgenomic_intron genomic_intron_create(int start, int end)
{
  pgenomic_intron pi=(pgenomic_intron)palloc(sizeof(struct _genomic_intron));
  pi->start=start;
  pi->end=end;
  pi->donor_pt=NULL;
  pi->acceptor_pt=NULL;
  pi->burset_frequency=-1;
  pi->info=list_create();
  pi->supportingESTs=0;
  pi->classified=false;
  pi->agree_type=2;
  return pi;
}

pburset_frequency burset_frequency_create(int start, int end)
{
  pburset_frequency pb=(pburset_frequency)palloc(sizeof(struct _burset_frequency));
  pb->start=start;
  pb->end=end;
  pb->frequency=-1;
  return pb;
}

/*
 * The parameter gb must be allocated and set before calling this procedure
 */
pgenomic_intron_info genomic_intron_info_create(pEST_info est, int EST_cut)
{
  pgenomic_intron_info pi=(pgenomic_intron_info)palloc(sizeof(struct _genomic_intron_info));
  pi->info=est;
  pi->EST_cut=EST_cut;
  return pi;
}

void genomic_intron_destroy(pgenomic_intron pi)
{
  if(pi!=NULL){
	if(pi->donor_pt != NULL)pfree(pi->donor_pt);
	if(pi->acceptor_pt != NULL)pfree(pi->acceptor_pt);
	if(pi->info != NULL)
		list_destroy(pi->info, (delete_function) genomic_intron_info_destroy);
	pfree(pi);
  }
}

void burset_frequency_destroy(pburset_frequency pb)
{
  if(pb!=NULL){
	pfree(pb);
  }
}

void genomic_intron_info_destroy(pgenomic_intron_info pi)
{
 if(pi != NULL){
	 pfree(pi);
 }
}
pintron intron_create(void)
{
  pintron pi=(pintron)palloc(sizeof(struct _intron));
  pi->gen_intron=NULL;
  pi->est_info=NULL;
  pi->donor=NULL;
  pi->acceptor=NULL;
  return pi;
}

void intron_destroy(pintron pi)
{
  if(pi!=NULL)pfree(pi);
}

void intron_destroy_free(pintron pi)
{
  if(pi!=NULL){
	  if(pi->gen_intron != NULL)
		  pfree(pi->gen_intron);
	  pfree(pi);
  }
}

void intron_destroy_free_if_false(pintron pi)
{
  if(pi!=NULL){
	  if(pi->isReal == false && pi->gen_intron != NULL)
		  pfree(pi->gen_intron);
	  pfree(pi);
  }
}

psubtree_factorizations subtree_factorizations_create(void)
{
  psubtree_factorizations pf=(psubtree_factorizations)palloc(sizeof(struct _subtree_factorizations));
  return pf;
}

void subtree_factorizations_destroy(psubtree_factorizations pf)
{
  if(pf!=NULL)pfree(pf);
}

psubtree_embeddings subtree_embeddings_create(void)
{
  psubtree_embeddings pe=(psubtree_embeddings)palloc(sizeof(struct _subtree_embeddings));
  return pe;
}

void subtree_embeddings_destroy(psubtree_embeddings pe)
{
  if(pe!=NULL)pfree(pe);
}

ppairing pairing_create(void) {

  static int next_id= 0;

  ppairing p= PALLOC(struct _pairing);

  p->id= next_id;
  ++next_id;

  p->p= 0;
  p->t= 0;
  p->l= 0;
  p->adjs= list_create();
  p->incs= list_create();
  p->visited=false;
  return p;
}

void pairing_destroy(ppairing p) {
  my_assert(p!=NULL);
  list_destroy(p->incs, noop_free);
  list_destroy(p->adjs, noop_free);
  pfree(p);
}

//Senza distruzione delle liste di adiacenti e incidenti
void pairing_destroy_2(ppairing p) {
  my_assert(p!=NULL);
  pfree(p);
}

int pairing_compare(const ppairing* pp1, const ppairing* pp2) {
  my_assert(pp1!=NULL);
  my_assert(pp2!=NULL);
  ppairing p1= *pp1;
  ppairing p2= *pp2;
  my_assert(p1!=NULL);
  my_assert(p2!=NULL);
  if (p1->p < p2->p) {
	 return -1;
  } else if (p1->p > p2->p) {
	 return 1;
  } else {
	 if (p1->t < p2->t) {
		return -1;
	 } else if (p1->t > p2->t) {
		return 1;
	 } else {
		return (p1->l)-(p2->l);
	 }
  }
}

int genomic_intron_compare(const pgenomic_intron* pgi1, const pgenomic_intron* pgi2) {
  my_assert(pgi1!=NULL);
  my_assert(pgi2!=NULL);
  pgenomic_intron p1= *pgi1;
  pgenomic_intron p2= *pgi2;
  my_assert(p1!=NULL);
  my_assert(p2!=NULL);
  if (p1->start < p2->start) {
	 return -1;
  } else if (p1->start > p2->start) {
	 return 1;
  } else {
	 if (p1->end < p2->end) {
		return -1;
	 } else {
		return 1;
	 }
  }
}

int burset_frequency_compare(const pburset_frequency* pbf1, const pburset_frequency* pbf2) {
  my_assert(pbf1!=NULL);
  my_assert(pbf2!=NULL);
  pburset_frequency p1= *pbf1;
  pburset_frequency p2= *pbf2;
  my_assert(p1!=NULL);
  my_assert(p2!=NULL);
  if (p1->frequency > p2->frequency) {
	 return -1;
  } else {
	 return 1;
  }
}

ppairing pairing_simple_copy(const ppairing p) {
  my_assert(p!=NULL);
  ppairing copy= pairing_create();
  copy->p= p->p;
  copy->t= p->t;
  copy->l= p->l;
  return copy;
}

pGEN_ESTS GEN_ESTS_create(void){
  pGEN_ESTS GEN_ESTS = PALLOC(struct _GEN_ESTS);
  GEN_ESTS->gen= NULL;
  GEN_ESTS->ests= list_create();
  return GEN_ESTS;
}

void GEN_ESTS_destroy(pGEN_ESTS GEN_ESTS){
  EST_info_destroy(GEN_ESTS->gen);
  list_destroy(GEN_ESTS->ests, (delete_function) EST_MEG_destroy);
  pfree(GEN_ESTS);
}

pEST_MEG EST_MEG_create(void){
  pEST_MEG EST_MEG = PALLOC(struct _EST_MEG);
  EST_MEG->est = NULL;
  EST_MEG->meg = NULL;
  return EST_MEG;
}

void EST_MEG_destroy(pEST_MEG EST_MEG){
  EST_info_destroy(EST_MEG->est);
  EA_destroy(EST_MEG->meg, (delete_function) vi_destroy);
  pfree(EST_MEG);
}

void vi_destroy(plist vi) {
  list_destroy(vi, (delete_function)pairing_destroy);
}

