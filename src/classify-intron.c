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
#include <ctype.h>
#include <math.h>

#include "est-factorizations.h"
#include "classify-intron.h"
#include "types.h"

#include "log.h"


//#define LOG_THRESHOLD LOG_LEVEL_TRACE

plist classify_genomic_intron_list(char *genomic_sequence, plist gen_intron_list){
	my_assert(genomic_sequence != NULL);
	my_assert(gen_intron_list != NULL);

	double **pwm_matx=LoadPWMMatrices();
	double **CVector=LoadCVPWMMatrices(pwm_matx);
	double **MAXVector=LoadMAXPWMMatrices(pwm_matx);

	plistit list_it_for_gen_intron=list_first(gen_intron_list);
	while(listit_has_next(list_it_for_gen_intron)){
		 pgenomic_intron current_gen_intron=(pgenomic_intron)listit_next(list_it_for_gen_intron);
		 double score5, score3, BPS_score;
		 int BPS_position;
		 current_gen_intron->type=classify_genomic_intron_start_end(genomic_sequence, current_gen_intron->start, current_gen_intron->end, &score5, &score3, &BPS_position, &BPS_score, pwm_matx, CVector, MAXVector);
		 current_gen_intron->score5=score5;
		 current_gen_intron->score3=score3;
		 current_gen_intron->BPS_position=BPS_position;
		 current_gen_intron->BPS_score=BPS_score;
		 current_gen_intron->classified=true;
	}
	listit_destroy(list_it_for_gen_intron);

	FreePWMMatrices(pwm_matx);
	FreeCVPWMMatrices(CVector);
	FreeMAXPWMMatrices(MAXVector);

	return gen_intron_list;
}

pgenomic_intron classify_genomic_intron(char *genomic_sequence, pgenomic_intron gen_intron){
	my_assert(genomic_sequence != NULL);
	my_assert(gen_intron != NULL);

	double **pwm_matx=LoadPWMMatrices();
	double **CVector=LoadCVPWMMatrices(pwm_matx);
	double **MAXVector=LoadMAXPWMMatrices(pwm_matx);

	 double score5, score3, BPS_score;
	 int BPS_position;
	 gen_intron->type=classify_genomic_intron_start_end(genomic_sequence, gen_intron->start, gen_intron->end, &score5, &score3, &BPS_position, &BPS_score, pwm_matx, CVector, MAXVector);
	 gen_intron->score5=score5;
	 gen_intron->score3=score3;
	 gen_intron->BPS_position=BPS_position;
	 gen_intron->BPS_score=BPS_score;
	 gen_intron->classified=true;

	FreePWMMatrices(pwm_matx);
	FreeCVPWMMatrices(CVector);
	FreeMAXPWMMatrices(MAXVector);

	return gen_intron;
}

char classify_genomic_intron_start_end(char *genomic_sequence, int start, int end, double *score5, double *score3, int *BPS_position, double *BPS_score, double **pwm_matx, double **CVector, double **MAXVector){
#ifndef NDEBUG
	my_assert(genomic_sequence != NULL);
	size_t gen_length=strlen(genomic_sequence);

	my_assert(start >= 0 && start < (int)gen_length);
	my_assert(end >= 0 && end < (int)gen_length);

	my_assert(score5 != NULL);
	my_assert(score3 != NULL);
	my_assert(BPS_position != NULL);
	my_assert(BPS_score != NULL);
	my_assert(pwm_matx != NULL);
	my_assert(CVector != NULL);
	my_assert(MAXVector != NULL);
#endif

	//BPS is searched inside a window from the 14th to 30th nucleotide before the intron 3' site
	char *intronSequence=real_substring(start, end-start+1, genomic_sequence);
	*BPS_position=ExistsGoodBPSinIntronSequenceWithMathInspector(intronSequence, pwm_matx[0], pwm_matx[1], CVector[0], CVector[1], MAXVector[0], MAXVector[1], BPS_score, 14, 30);

	char *pt_5=real_substring(0, 2, intronSequence);
	char *pt_3=real_substring(strlen(intronSequence)-2, 2, intronSequence);

	double scoreU12_3=0.0f, scoreU2_3=0.0f;
	double scoreU12_3_2=0.0f, scoreU2_3_2=0.0f;
	double scoreU12_5=0.0f, scoreU2_5=0.0f;
	double scoreU12_5_2=0.0f, scoreU2_5_2=0.0f;

	char pt_type=1;	/*0=gt-ag or gc-ag; 1=altro*/

	if((!strcmp(pt_5, "gt") || !strcmp(pt_5, "GT")) && (!strcmp(pt_3, "ag") || !strcmp(pt_3, "AG"))){
		pt_type=0;
		scoreU12_5=GetScoreOf5PrimeGTAGU12BySS(genomic_sequence, start, pwm_matx[2], CVector[2], MAXVector[2]);
		scoreU2_5=GetScoreOf5PrimeGTAGU2BySS(genomic_sequence, start, pwm_matx[4], CVector[4], MAXVector[4]);

		scoreU12_3=GetScoreOf3PrimeGTAGU12BySS(genomic_sequence, end, pwm_matx[6], CVector[6], MAXVector[6]);
		scoreU2_3=GetScoreOf3PrimeGTAGU2BySS(genomic_sequence, end, pwm_matx[8], CVector[8], MAXVector[8]);
	}
	else{
		if((!strcmp(pt_5, "gc") || !strcmp(pt_5, "GC")) && (!strcmp(pt_3, "ag") || !strcmp(pt_3, "AG"))){
			pt_type=0;

			scoreU2_5=GetScoreOf5PrimeGCAGU2BySS(genomic_sequence, start, pwm_matx[5], CVector[5], MAXVector[5]);
			scoreU2_3=GetScoreOf3PrimeGCAGU2BySS(genomic_sequence, end, pwm_matx[9], CVector[9], MAXVector[9]);
			scoreU12_5=GetScoreOf5PrimeGTAGU12BySS(genomic_sequence, start, pwm_matx[2], CVector[2], MAXVector[2]);

			scoreU12_5_2=GetScoreOf5PrimeATACU12BySS(genomic_sequence, start, pwm_matx[3], CVector[3], MAXVector[3]);
			if(scoreU12_5_2 > scoreU12_5)
				scoreU12_5=scoreU12_5_2;

			scoreU12_3=GetScoreOf3PrimeGTAGU12BySS(genomic_sequence, end, pwm_matx[6], CVector[6], MAXVector[6]);

			scoreU12_3_2=GetScoreOf3PrimeATACU12BySS(genomic_sequence, end, pwm_matx[7], CVector[7], MAXVector[7]);
			if(scoreU12_3_2 > scoreU12_3)
				scoreU12_3=scoreU12_3_2;
		}
		else{
			if((!strcmp(pt_5, "at") || !strcmp(pt_5, "AT")) && (!strcmp(pt_3, "ac") || !strcmp(pt_3, "AC"))){
				scoreU12_5=GetScoreOf5PrimeATACU12BySS(genomic_sequence, start, pwm_matx[3], CVector[3], MAXVector[3]);
				scoreU12_3=GetScoreOf3PrimeATACU12BySS(genomic_sequence, end, pwm_matx[7], CVector[7], MAXVector[7]);

				scoreU2_5=GetScoreOf5PrimeGTAGU2BySS(genomic_sequence, start, pwm_matx[4], CVector[4], MAXVector[4]);
				scoreU2_5_2=GetScoreOf5PrimeGCAGU2BySS(genomic_sequence, start, pwm_matx[5], CVector[5], MAXVector[5]);
				if(scoreU2_5_2 > scoreU2_5)
					scoreU2_5=scoreU2_5_2;

				scoreU2_3=GetScoreOf3PrimeGTAGU2BySS(genomic_sequence, end, pwm_matx[8], CVector[8], MAXVector[8]);
				scoreU2_3_2=GetScoreOf3PrimeGCAGU2BySS(genomic_sequence, end, pwm_matx[9], CVector[9], MAXVector[9]);
				if(scoreU2_3_2 > scoreU2_3)
					scoreU2_3=scoreU2_3_2;
			}
			else{
				scoreU12_5=GetScoreOf5PrimeGTAGU12BySS(genomic_sequence, start, pwm_matx[2], CVector[2], MAXVector[2]);
				scoreU12_5_2=GetScoreOf5PrimeATACU12BySS(genomic_sequence, start, pwm_matx[3], CVector[3], MAXVector[3]);
				if(scoreU12_5_2 > scoreU12_5)
					scoreU12_5=scoreU12_5_2;

				scoreU2_5=GetScoreOf5PrimeGTAGU2BySS(genomic_sequence, start, pwm_matx[4], CVector[4], MAXVector[4]);
				scoreU2_5_2=GetScoreOf5PrimeGCAGU2BySS(genomic_sequence, start, pwm_matx[5], CVector[5], MAXVector[5]);
				if(scoreU2_5_2 > scoreU2_5)
					scoreU2_5=scoreU2_5_2;

				scoreU12_3=GetScoreOf3PrimeGTAGU12BySS(genomic_sequence, end, pwm_matx[6], CVector[6], MAXVector[6]);
				scoreU12_3_2=GetScoreOf3PrimeATACU12BySS(genomic_sequence, end, pwm_matx[7], CVector[7], MAXVector[7]);
				if(scoreU12_3_2 > scoreU12_3)
					scoreU12_3=scoreU12_3_2;

				scoreU2_3=GetScoreOf3PrimeGTAGU2BySS(genomic_sequence, end, pwm_matx[8], CVector[8], MAXVector[8]);
				scoreU2_3_2=GetScoreOf3PrimeGCAGU2BySS(genomic_sequence, end, pwm_matx[9], CVector[9], MAXVector[9]);
				if(scoreU2_3_2 > scoreU2_3)
					scoreU2_3=scoreU2_3_2;
			}
		}
	}

	char type=2;
	if(*BPS_position != -1){
		if(scoreU12_5 > scoreU2_5)
			type=0;
		else
			type=1;
	}
	else{
		//Se e' gt-ag o gc-ag
		if(pt_type == 0){
			type=1;
			*BPS_position=ExistsGoodBPSinIntronSequenceWithMathInspector(intronSequence, pwm_matx[0], pwm_matx[1], CVector[0], CVector[1], MAXVector[0], MAXVector[1], BPS_score, 30, 200);
		}
		else{
			if(scoreU12_5-scoreU2_5 > 0.25 && scoreU12_5 >= 0.75){
				type=0;
				*BPS_position=ExistsGoodBPSinIntronSequenceWithMathInspector(intronSequence, pwm_matx[0], pwm_matx[1], CVector[0], CVector[1], MAXVector[0], MAXVector[1], BPS_score, 30, 200);
			}
		}
	}

	*score5=0.0f;
	*score3=0.0f;

	if(type == 0){
		*score5=scoreU12_5;
		*score3=scoreU12_3;
	}
	else{
		*score5=scoreU2_5;
		*score3=scoreU2_3;
	}

	pfree(pt_5);
	pfree(pt_3);
	pfree(intronSequence);

	return type;
}

double GetScoreOf5PrimeGTAGU12BySS(char *genomic_sequence, int splice5, double *pwm_5PrimeGTAGU12, double *CV_pwm_5PrimeGTAGU12, double *MAXV_pwm_5PrimeGTAGU12){

	my_assert(splice5 >= 0);
	my_assert(genomic_sequence != NULL);
	my_assert(pwm_5PrimeGTAGU12 != NULL);
	my_assert(CV_pwm_5PrimeGTAGU12 != NULL);
	my_assert(MAXV_pwm_5PrimeGTAGU12 != NULL);

	double score=0.0f;

	char *spliceSequence=real_substring(splice5-3, 14, genomic_sequence);

	score=GetScoreOf5PrimeGTAGU12(spliceSequence, pwm_5PrimeGTAGU12, CV_pwm_5PrimeGTAGU12, MAXV_pwm_5PrimeGTAGU12);

	pfree(spliceSequence);

	return score;
}

double GetScoreOf5PrimeGTAGU12(char *sequence, double *pwm_5PrimeGTAGU12, double *CV_pwm_5PrimeGTAGU12, double *MAXV_pwm_5PrimeGTAGU12){

	my_assert(sequence != NULL);
	my_assert(pwm_5PrimeGTAGU12 != NULL);
	my_assert(CV_pwm_5PrimeGTAGU12 != NULL);
	my_assert(MAXV_pwm_5PrimeGTAGU12 != NULL);

	double score=0.0f;

#ifndef NDEBUG
	size_t length=strlen(sequence);
	my_assert(length == 14);
#endif

	score=GetMatInspectorScoreOfaMotif(sequence, pwm_5PrimeGTAGU12, CV_pwm_5PrimeGTAGU12, MAXV_pwm_5PrimeGTAGU12, 14);

	return score;
}

double GetScoreOf5PrimeATACU12BySS(char *genomic_sequence, int splice5, double *pwm_5PrimeATACU12, double *CV_pwm_5PrimeATACU12, double *MAXV_pwm_5PrimeATACU12){

	my_assert(splice5 >= 0);
	my_assert(genomic_sequence != NULL);
	my_assert(pwm_5PrimeATACU12 != NULL);
	my_assert(CV_pwm_5PrimeATACU12 != NULL);
	my_assert(MAXV_pwm_5PrimeATACU12 != NULL);

	double score=0.0f;

	char *spliceSequence=real_substring(splice5-3, 14, genomic_sequence);

	score=GetScoreOf5PrimeATACU12(spliceSequence, pwm_5PrimeATACU12, CV_pwm_5PrimeATACU12, MAXV_pwm_5PrimeATACU12);

	pfree(spliceSequence);

	return score;
}

double GetScoreOf5PrimeATACU12(char *sequence, double *pwm_5PrimeATACU12, double *CV_pwm_5PrimeATACU12, double *MAXV_pwm_5PrimeATACU12){

	my_assert(sequence != NULL);
	my_assert(pwm_5PrimeATACU12 != NULL);
	my_assert(CV_pwm_5PrimeATACU12 != NULL);
	my_assert(MAXV_pwm_5PrimeATACU12 != NULL);

	double score=0.0f;

#ifndef NDEBUG
	size_t length=strlen(sequence);
	my_assert(length == 14);
#endif

	score=GetMatInspectorScoreOfaMotif(sequence, pwm_5PrimeATACU12, CV_pwm_5PrimeATACU12, MAXV_pwm_5PrimeATACU12, 14);

	return score;
}

double GetScoreOf5PrimeGTAGU2BySS(char *genomic_sequence, int splice5, double *pwm_5PrimeGTAGU2, double *CV_pwm_5PrimeGTAGU2, double *MAXV_pwm_5PrimeGTAGU2){

	my_assert(splice5 >= 0);
	my_assert(genomic_sequence != NULL);
	my_assert(pwm_5PrimeGTAGU2 != NULL);
	my_assert(CV_pwm_5PrimeGTAGU2 != NULL);
	my_assert(MAXV_pwm_5PrimeGTAGU2 != NULL);

	double score=0.0f;

	char *spliceSequence=real_substring(splice5-3, 13, genomic_sequence);

	score=GetScoreOf5PrimeGTAGU2(spliceSequence, pwm_5PrimeGTAGU2, CV_pwm_5PrimeGTAGU2, MAXV_pwm_5PrimeGTAGU2);

	pfree(spliceSequence);

	return score;
}

double GetScoreOf5PrimeGTAGU2(char *sequence, double *pwm_5PrimeGTAGU2, double *CV_pwm_5PrimeGTAGU2, double *MAXV_pwm_5PrimeGTAGU2){

	my_assert(sequence != NULL);
	my_assert(pwm_5PrimeGTAGU2 != NULL);
	my_assert(CV_pwm_5PrimeGTAGU2 != NULL);
	my_assert(MAXV_pwm_5PrimeGTAGU2 != NULL);

	double score=0.0f;

#ifndef NDEBUG
	int length=strlen(sequence);
	my_assert(length == 13);
#endif

	score=GetMatInspectorScoreOfaMotif(sequence, pwm_5PrimeGTAGU2, CV_pwm_5PrimeGTAGU2, MAXV_pwm_5PrimeGTAGU2, 13);

	return score;
}

double GetScoreOf5PrimeGCAGU2BySS(char *genomic_sequence, int splice5, double *pwm_5PrimeGCAGU2, double *CV_pwm_5PrimeGCAGU2, double *MAXV_pwm_5PrimeGCAGU2){

	my_assert(splice5 >= 0);
	my_assert(genomic_sequence != NULL);
	my_assert(pwm_5PrimeGCAGU2 != NULL);
	my_assert(CV_pwm_5PrimeGCAGU2 != NULL);
	my_assert(MAXV_pwm_5PrimeGCAGU2 != NULL);

	double score=0.0f;

	char *spliceSequence=real_substring(splice5-3, 14, genomic_sequence);

	score=GetScoreOf5PrimeGCAGU2(spliceSequence, pwm_5PrimeGCAGU2, CV_pwm_5PrimeGCAGU2, MAXV_pwm_5PrimeGCAGU2);

	free(spliceSequence);

	return score;
}

double GetScoreOf5PrimeGCAGU2(char *sequence, double *pwm_5PrimeGCAGU2, double *CV_pwm_5PrimeGCAGU2, double *MAXV_pwm_5PrimeGCAGU2){

	my_assert(sequence != NULL);
	my_assert(pwm_5PrimeGCAGU2 != NULL);
	my_assert(CV_pwm_5PrimeGCAGU2 != NULL);
	my_assert(MAXV_pwm_5PrimeGCAGU2 != NULL);

	double score=0.0f;

#ifndef NDEBUG
	size_t length=strlen(sequence);
	my_assert(length == 14);
#endif

	score=GetMatInspectorScoreOfaMotif(sequence, pwm_5PrimeGCAGU2, CV_pwm_5PrimeGCAGU2, MAXV_pwm_5PrimeGCAGU2, 14);

	return score;
}

double GetScoreOf3PrimeGTAGU12BySS(char *genomic_sequence, int splice3, double *pwm_3PrimeGTAGU12, double *CV_pwm_3PrimeGTAGU12, double *MAXV_pwm_3PrimeGTAGU12){

	my_assert(splice3 >= 0);
	my_assert(genomic_sequence != NULL);
	my_assert(pwm_3PrimeGTAGU12 != NULL);
	my_assert(CV_pwm_3PrimeGTAGU12 != NULL);
	my_assert(MAXV_pwm_3PrimeGTAGU12 != NULL);

	double score=0.0f;

	char *spliceSequence=real_substring(splice3-14+1, 18, genomic_sequence);

	score=GetScoreOf3PrimeGTAGU12(spliceSequence, pwm_3PrimeGTAGU12, CV_pwm_3PrimeGTAGU12, MAXV_pwm_3PrimeGTAGU12);

	pfree(spliceSequence);

	return score;
}

double GetScoreOf3PrimeGTAGU12(char *sequence, double *pwm_3PrimeGTAGU12, double *CV_pwm_3PrimeGTAGU12, double *MAXV_pwm_3PrimeGTAGU12){

	my_assert(sequence != NULL);
	my_assert(pwm_3PrimeGTAGU12 != NULL);
	my_assert(CV_pwm_3PrimeGTAGU12 != NULL);
	my_assert(MAXV_pwm_3PrimeGTAGU12 != NULL);

	double score=0.0f;

#ifndef NDEBUG
	size_t length=strlen(sequence);
	my_assert(length == 18);
#endif

	score=GetMatInspectorScoreOfaMotif(sequence, pwm_3PrimeGTAGU12, CV_pwm_3PrimeGTAGU12, MAXV_pwm_3PrimeGTAGU12, 18);

	return score;
}

double GetScoreOf3PrimeATACU12BySS(char *genomic_sequence, int splice3, double *pwm_3PrimeATACU12, double *CV_pwm_3PrimeATACU12, double *MAXV_pwm_3PrimeATACU12){

	my_assert(splice3 >= 0);
	my_assert(genomic_sequence != NULL);
	my_assert(pwm_3PrimeATACU12 != NULL);
	my_assert(CV_pwm_3PrimeATACU12 != NULL);
	my_assert(MAXV_pwm_3PrimeATACU12 != NULL);

	double score=0.0f;

	char *spliceSequence=real_substring(splice3-14+1, 17, genomic_sequence);

	score=GetScoreOf3PrimeATACU12(spliceSequence, pwm_3PrimeATACU12, CV_pwm_3PrimeATACU12, MAXV_pwm_3PrimeATACU12);

	pfree(spliceSequence);

	return score;
}

double GetScoreOf3PrimeATACU12(char *sequence, double *pwm_3PrimeATACU12, double *CV_pwm_3PrimeATACU12, double *MAXV_pwm_3PrimeATACU12){

	my_assert(sequence != NULL);
	my_assert(pwm_3PrimeATACU12 != NULL);
	my_assert(CV_pwm_3PrimeATACU12 != NULL);
	my_assert(MAXV_pwm_3PrimeATACU12 != NULL);

	double score=0.0f;

#ifndef NDEBUG
	size_t length=strlen(sequence);
	my_assert(length == 17);
#endif

	score=GetMatInspectorScoreOfaMotif(sequence, pwm_3PrimeATACU12, CV_pwm_3PrimeATACU12, MAXV_pwm_3PrimeATACU12, 17);

	return score;
}

double GetScoreOf3PrimeGTAGU2BySS(char *genomic_sequence, int splice3, double *pwm_3PrimeGTAGU2, double *CV_pwm_3PrimeGTAGU2, double *MAXV_pwm_3PrimeGTAGU2){

	my_assert(splice3 >= 0);
	my_assert(genomic_sequence != NULL);
	my_assert(pwm_3PrimeGTAGU2 != NULL);
	my_assert(CV_pwm_3PrimeGTAGU2 != NULL);
	my_assert(MAXV_pwm_3PrimeGTAGU2 != NULL);

	double score=0.0f;

	char *spliceSequence=real_substring(splice3-14+1, 17, genomic_sequence);

	score=GetScoreOf3PrimeGTAGU2(spliceSequence, pwm_3PrimeGTAGU2, CV_pwm_3PrimeGTAGU2, MAXV_pwm_3PrimeGTAGU2);

	pfree(spliceSequence);

	return score;
}

double GetScoreOf3PrimeGTAGU2(char *sequence, double *pwm_3PrimeGTAGU2, double *CV_pwm_3PrimeGTAGU2, double *MAXV_pwm_3PrimeGTAGU2){

	my_assert(sequence != NULL);
	my_assert(pwm_3PrimeGTAGU2 != NULL);
	my_assert(CV_pwm_3PrimeGTAGU2 != NULL);
	my_assert(MAXV_pwm_3PrimeGTAGU2 != NULL);

	double score=0.0f;

#ifndef NDEBUG
	size_t length=strlen(sequence);
	my_assert(length == 17);
#endif

	score=GetMatInspectorScoreOfaMotif(sequence, pwm_3PrimeGTAGU2, CV_pwm_3PrimeGTAGU2, MAXV_pwm_3PrimeGTAGU2, 17);

	return score;
}

double GetScoreOf3PrimeGCAGU2BySS(char *genomic_sequence, int splice3, double *pwm_3PrimeGCAGU2, double *CV_pwm_3PrimeGCAGU2, double *MAXV_pwm_3PrimeGCAGU2){

	my_assert(splice3 >= 0);
	my_assert(genomic_sequence != NULL);
	my_assert(pwm_3PrimeGCAGU2 != NULL);
	my_assert(CV_pwm_3PrimeGCAGU2 != NULL);
	my_assert(MAXV_pwm_3PrimeGCAGU2 != NULL);

	double score=0.0f;

	char *spliceSequence=real_substring(splice3-14+1, 18, genomic_sequence);

	score=GetScoreOf3PrimeGCAGU2(spliceSequence, pwm_3PrimeGCAGU2, CV_pwm_3PrimeGCAGU2, MAXV_pwm_3PrimeGCAGU2);

	pfree(spliceSequence);

	return score;
}

double GetScoreOf3PrimeGCAGU2(char *sequence, double *pwm_3PrimeGCAGU2, double *CV_pwm_3PrimeGCAGU2, double *MAXV_pwm_3PrimeGCAGU2){

	my_assert(sequence != NULL);
	my_assert(pwm_3PrimeGCAGU2 != NULL);
	my_assert(CV_pwm_3PrimeGCAGU2 != NULL);
	my_assert(MAXV_pwm_3PrimeGCAGU2 != NULL);

	double score=0.0f;

#ifndef NDEBUG
	size_t length=strlen(sequence);
	my_assert(length == 18);
#endif

	score=GetMatInspectorScoreOfaMotif(sequence, pwm_3PrimeGCAGU2, CV_pwm_3PrimeGCAGU2, MAXV_pwm_3PrimeGCAGU2, 18);

	return score;
}

int ExistsGoodBPSinIntronSequenceWithMathInspector(char *intronSequence, double *pwm_9, double *pwm_10, double *CV_pwm_9, double *CV_pwm_10, double *MAXV_pwm_9, double *MAXV_pwm_10, double *score, int range_start, int range_end){

	my_assert(intronSequence != NULL);
	my_assert(pwm_9 != NULL);
	my_assert(pwm_10 != NULL);
	my_assert(CV_pwm_9 != NULL);
	my_assert(CV_pwm_10 != NULL);
	my_assert(MAXV_pwm_9 != NULL);
	my_assert(MAXV_pwm_10 != NULL);
	my_assert(score != NULL);

	*score=0.0f;

	if(range_end > (int)strlen(intronSequence)){
		return -1;
	}

	double score_9=0.0f;
	double score_10=0.0f;
	int bps_9=SearchBPSinIntronSequenceWithMathInspector(intronSequence, pwm_9, CV_pwm_9, MAXV_pwm_9, &score_9, range_start, range_end);
	int bps_10=SearchBPSinIntronSequenceWithMathInspector(intronSequence, pwm_10, CV_pwm_10, MAXV_pwm_10, &score_10, range_start, range_end);

	int bps_pos=-1;

	if(score_9 > score_10){
		if(score_9 > 0.75f){
			*score=score_9;
			bps_pos=bps_9;
		}
	}
	else{
		if(score_10 > 0.75f){
			*score=score_10;
			bps_pos=bps_10;
		}
	}

	return bps_pos;
}

int SearchBPSinIntronSequenceWithMathInspector(char *intronSequence, double *score_matx, double *CVector, double *MAXVector, double *score, int range_start, int range_end){

	my_assert(intronSequence != NULL);
	my_assert(score_matx != NULL);
	my_assert(CVector != NULL);
	my_assert(MAXVector != NULL);
	my_assert(score != NULL);

	*score=0.0f;

	size_t length=strlen(intronSequence);
	if(length < (unsigned int)range_start){
		return -1;
	}

	int start_w=(int)length-range_end;
	int end_w=(int)length-range_start;

	if(start_w < 0)
		start_w=0;

	int start_bps=-1;
	char *bps=NULL;
	int i=start_w;
	bool firsttime=true;

	while(i <= end_w){
		bps=real_substring(i, 12, intronSequence);

		double score_bps=GetMatInspectorScoreOfaMotif(bps, score_matx, CVector, MAXVector, 12);

		if(firsttime == true || score_bps >= *score){
			*score=score_bps;
			start_bps=i;
			if(firsttime == true)
				firsttime=false;
		}

		pfree(bps);
		i++;
	}

	return start_bps;
}

double GetMatInspectorScoreOfaMotif(char *sequence, double *pwm, double *CVector, double *MAXVector, int motif_length){

#ifndef NDEBUG
	my_assert(sequence != NULL);
	size_t length=strlen(sequence);
	my_assert(length == (unsigned int)motif_length);

	my_assert(pwm != NULL);
	my_assert(CVector != NULL);
	my_assert(MAXVector != NULL);
#endif

	int i;
	double den=0.0f;
	double num=0.0f;
	for(i=0; i<motif_length; i++){
		int index=-1;

		if(sequence[i] == 'N' || sequence[i] == 'n'){
			index=0;
		}

		if(sequence[i] == 'A' || sequence[i] == 'a'){
			index=0;
		}
		if(sequence[i] == 'C' || sequence[i] == 'c'){
			index=1;
		}
		if(sequence[i] == 'G' || sequence[i] == 'g'){
			index=2;
		}
		if(sequence[i] == 'T' || sequence[i] == 't'){
			index=3;
		}
		my_assert(index != -1);

		num+=CVector[i]*pwm[index*motif_length+i];
		den+=CVector[i]*MAXVector[i];
	}

	double score=num/den;

	return score;
}

double *GetPWMforBPS_9(){
	double *pwm=NULL;

	pwm=NPALLOC(double, 4*12);
	my_assert(pwm != NULL);

	pwm[0*12+0]=0.16;
	pwm[0*12+1]=0.19;
	pwm[0*12+2]=0.08;
	pwm[0*12+3]=0.09;
	pwm[0*12+4]=0.01;
	pwm[0*12+5]=0.01;
	pwm[0*12+6]=0.01;
	pwm[0*12+7]=0.01;
	pwm[0*12+8]=1.00;
	pwm[0*12+9]=0.94;
	pwm[0*12+10]=0.01;
	pwm[0*12+11]=0.28;

	pwm[1*12+0]=0.15;
	pwm[1*12+1]=0.18;
	pwm[1*12+2]=0.12;
	pwm[1*12+3]=0.09;
	pwm[1*12+4]=0.90;
	pwm[1*12+5]=0.90;
	pwm[1*12+6]=0.01;
	pwm[1*12+7]=0.01;
	pwm[1*12+8]=0.00;
	pwm[1*12+9]=0.01;
	pwm[1*12+10]=0.84;
	pwm[1*12+11]=0.16;

	pwm[2*12+0]=0.18;
	pwm[2*12+1]=0.14;
	pwm[2*12+2]=0.13;
	pwm[2*12+3]=0.07;
	pwm[2*12+4]=0.01;
	pwm[2*12+5]=0.01;
	pwm[2*12+6]=0.01;
	pwm[2*12+7]=0.01;
	pwm[2*12+8]=0.00;
	pwm[2*12+9]=0.04;
	pwm[2*12+10]=0.01;
	pwm[2*12+11]=0.04;

	pwm[3*12+0]=0.51;
	pwm[3*12+1]=0.49;
	pwm[3*12+2]=0.67;
	pwm[3*12+3]=0.75;
	pwm[3*12+4]=0.08;
	pwm[3*12+5]=0.08;
	pwm[3*12+6]=0.97;
	pwm[3*12+7]=0.97;
	pwm[3*12+8]=0.00;
	pwm[3*12+9]=0.01;
	pwm[3*12+10]=0.14;
	pwm[3*12+11]=0.52;

	int i;
	for(i=0; i<4; i++){
		int j;
		for(j=0; j<12; j++){
			pwm[i*12+j]=pwm[i*12+j]+0.00001f;
		}
	}

	return pwm;
}
double *GetPWMforBPS_10(){
	double *pwm=NULL;

	pwm=NPALLOC(double, 4*12);
	my_assert(pwm != NULL);

	pwm[0*12+0]=0.12;
	pwm[0*12+1]=0.12;
	pwm[0*12+2]=0.13;
	pwm[0*12+3]=0.09;
	pwm[0*12+4]=0.01;
	pwm[0*12+5]=0.01;
	pwm[0*12+6]=0.01;
	pwm[0*12+7]=0.01;
	pwm[0*12+8]=0.70;
	pwm[0*12+9]=1.00;
	pwm[0*12+10]=0.02;
	pwm[0*12+11]=0.23;

	pwm[1*12+0]=0.15;
	pwm[1*12+1]=0.17;
	pwm[1*12+2]=0.20;
	pwm[1*12+3]=0.10;
	pwm[1*12+4]=0.86;
	pwm[1*12+5]=0.92;
	pwm[1*12+6]=0.03;
	pwm[1*12+7]=0.01;
	pwm[1*12+8]=0.02;
	pwm[1*12+9]=0.00;
	pwm[1*12+10]=0.82;
	pwm[1*12+11]=0.25;

	pwm[2*12+0]=0.21;
	pwm[2*12+1]=0.18;
	pwm[2*12+2]=0.12;
	pwm[2*12+3]=0.05;
	pwm[2*12+4]=0.01;
	pwm[2*12+5]=0.01;
	pwm[2*12+6]=0.01;
	pwm[2*12+7]=0.03;
	pwm[2*12+8]=0.24;
	pwm[2*12+9]=0.00;
	pwm[2*12+10]=0.02;
	pwm[2*12+11]=0.05;

	pwm[3*12+0]=0.52;
	pwm[3*12+1]=0.53;
	pwm[3*12+2]=0.55;
	pwm[3*12+3]=0.76;
	pwm[3*12+4]=0.12;
	pwm[3*12+5]=0.06;
	pwm[3*12+6]=0.95;
	pwm[3*12+7]=0.95;
	pwm[3*12+8]=0.04;
	pwm[3*12+9]=0.00;
	pwm[3*12+10]=0.14;
	pwm[3*12+11]=0.47;

	int i;
	for(i=0; i<4; i++){
		int j;
		for(j=0; j<12; j++){
			pwm[i*12+j]=pwm[i*12+j]+0.00001f;
		}
	}

	return pwm;
}
double *GetPWMfor5PrimeGTAGU12(){
	double *pwm=NULL;

	pwm=NPALLOC(double, 4*14);
	my_assert(pwm != NULL);

	pwm[0*14+0]=0.293478260869565;
	pwm[0*14+1]=0.271739130434783;
	pwm[0*14+2]=0.217391304347826;
	pwm[0*14+3]=0.00;
	pwm[0*14+4]=0.00;
	pwm[0*14+5]=0.983695652173913;
	pwm[0*14+6]=0.00543478260869565;
	pwm[0*14+7]=0.00543478260869565;
	pwm[0*14+8]=0.0217391304347826;
	pwm[0*14+9]=0.0108695652173913;
	pwm[0*14+10]=0.0380434782608696;
	pwm[0*14+11]=0.103260869565217;
	pwm[0*14+12]=0.184782608695652;
	pwm[0*14+13]=0.40;

	pwm[1*14+0]=0.239130434782609;
	pwm[1*14+1]=0.326086956521739;
	pwm[1*14+2]=0.184782608695652;
	pwm[1*14+3]=0.00;
	pwm[1*14+4]=0.00;
	pwm[1*14+5]=0.00543478260869565;
	pwm[1*14+6]=0.00543478260869565;
	pwm[1*14+7]=0.983695652173913;
	pwm[1*14+8]=0.907608695652174;
	pwm[1*14+9]=0.016304347826087;
	pwm[1*14+10]=0.0489130434782609;
	pwm[1*14+11]=0.195652173913043;
	pwm[1*14+12]=0.266304347826087;
	pwm[1*14+13]=0.20;

	pwm[2*14+0]=0.239130434782609;
	pwm[2*14+1]=0.152173913043478;
	pwm[2*14+2]=0.0543478260869565;
	pwm[2*14+3]=1.00;
	pwm[2*14+4]=0.00;
	pwm[2*14+5]=0.00543478260869565;
	pwm[2*14+6]=0.00543478260869565;
	pwm[2*14+7]=0.00543478260869565;
	pwm[2*14+8]=0.00543478260869565;
	pwm[2*14+9]=0.0108695652173913;
	pwm[2*14+10]=0.0489130434782609;
	pwm[2*14+11]=0.0869565217391304;
	pwm[2*14+12]=0.141304347826087;
	pwm[2*14+13]=0.20;

	pwm[3*14+0]=0.228260869565217;
	pwm[3*14+1]=0.25;
	pwm[3*14+2]=0.543478260869565;
	pwm[3*14+3]=0.00;
	pwm[3*14+4]=1.00;
	pwm[3*14+5]=0.00543478260869565;
	pwm[3*14+6]=0.983695652173913;
	pwm[3*14+7]=0.00543478260869565;
	pwm[3*14+8]=0.0652173913043478;
	pwm[3*14+9]=0.96195652173913;
	pwm[3*14+10]=0.864130434782609;
	pwm[3*14+11]=0.614130434782609;
	pwm[3*14+12]=0.407608695652174;
	pwm[3*14+13]=0.20;

	int i;
	for(i=0; i<4; i++){
		int j;
		for(j=0; j<14; j++){
			pwm[i*14+j]=pwm[i*14+j]+0.00001f;
		}
	}

	return pwm;
}
double *GetPWMfor5PrimeATACU12(){
	double *pwm=NULL;

	pwm=NPALLOC(double, 4*14);
	my_assert(pwm != NULL);

	pwm[0*14+0]=0.271028037383178;
	pwm[0*14+1]=0.280373831775701;
	pwm[0*14+2]=0.299065420560748;
	pwm[0*14+3]=1.00;
	pwm[0*14+4]=0.00;
	pwm[0*14+5]=0.97196261682243;
	pwm[0*14+6]=0.00934579439252336;
	pwm[0*14+7]=0.00934579439252336;
	pwm[0*14+8]=0.00934579439252336;
	pwm[0*14+9]=0.0186915887850467;
	pwm[0*14+10]=0.0186915887850467;
	pwm[0*14+11]=0.0467289719626168;
	pwm[0*14+12]=0.177570093457944;
	pwm[0*14+13]=0.40;

	pwm[1*14+0]=0.280373831775701;
	pwm[1*14+1]=0.271028037383178;
	pwm[1*14+2]=0.289719626168224;
	pwm[1*14+3]=0.00;
	pwm[1*14+4]=0.00;
	pwm[1*14+5]=0.00934579439252336;
	pwm[1*14+6]=0.00934579439252336;
	pwm[1*14+7]=0.97196261682243;
	pwm[1*14+8]=0.962616822429907;
	pwm[1*14+9]=0.0186915887850467;
	pwm[1*14+10]=0.0747663551401869;
	pwm[1*14+11]=0.205607476635514;
	pwm[1*14+12]=0.224299065420561;
	pwm[1*14+13]=0.20;

	pwm[2*14+0]=0.224299065420561;
	pwm[2*14+1]=0.149532710280374;
	pwm[2*14+2]=0.299065420560748;
	pwm[2*14+3]=0.00;
	pwm[2*14+4]=0.00;
	pwm[2*14+5]=0.00934579439252336;
	pwm[2*14+6]=0.00934579439252336;
	pwm[2*14+7]=0.00934579439252336;
	pwm[2*14+8]=0.00934579439252336;
	pwm[2*14+9]=0.00934579439252336;
	pwm[2*14+10]=0.0186915887850467;
	pwm[2*14+11]=0.0373831775700935;
	pwm[2*14+12]=0.233644859813084;
	pwm[2*14+13]=0.20;

	pwm[3*14+0]=0.224299065420561;
	pwm[3*14+1]=0.299065420560748;
	pwm[3*14+2]=0.11214953271028;
	pwm[3*14+3]=0.00;
	pwm[3*14+4]=1.00;
	pwm[3*14+5]=0.00934579439252336;
	pwm[3*14+6]=0.97196261682243;
	pwm[3*14+7]=0.00934579439252336;
	pwm[3*14+8]=0.0186915887850467;
	pwm[3*14+9]=0.953271028037383;
	pwm[3*14+10]=0.88785046728972;
	pwm[3*14+11]=0.710280373831776;
	pwm[3*14+12]=0.364485981308411;
	pwm[3*14+13]=0.20;

	int i;
	for(i=0; i<4; i++){
		int j;
		for(j=0; j<14; j++){
			pwm[i*14+j]=pwm[i*14+j]+0.00001f;
		}
	}

	return pwm;
}
double *GetPWMfor5PrimeGTAGU2(){
	double *pwm=NULL;

	pwm=NPALLOC(double, 4*13);
	my_assert(pwm != NULL);

	pwm[0*13+0]=0.341467901547831;
	pwm[0*13+1]=0.660619132199949;
	pwm[0*13+2]=0.0973420451662015;
	pwm[0*13+3]=0.00;
	pwm[0*13+4]=0.00;
	pwm[0*13+5]=0.596517381375285;
	pwm[0*13+6]=0.763251712763258;
	pwm[0*13+7]=0.0781781273788379;
	pwm[0*13+8]=0.183893681806648;
	pwm[0*13+9]=0.297487947221517;
	pwm[0*13+10]=0.224714539456991;
	pwm[0*13+11]=0.222602131438721;
	pwm[0*13+12]=0.225862725196651;

	pwm[1*13+0]=0.375069779243847;
	pwm[1*13+1]=0.0834496320730779;
	pwm[1*13+2]=0.00990230905861456;
	pwm[1*13+3]=0.00;
	pwm[1*13+4]=0.00;
	pwm[1*13+5]=0.0225640700329866;
	pwm[1*13+6]=0.0241055569652372;
	pwm[1*13+7]=0.0350291804110632;
	pwm[1*13+8]=0.144626998223801;
	pwm[1*13+9]=0.189672671910683;
	pwm[1*13+10]=0.247881248414108;
	pwm[1*13+11]=0.259737376300431;
	pwm[1*13+12]=0.23369068764273;

	pwm[2*13+0]=0.183646282669373;
	pwm[2*13+1]=0.10822126363867;
	pwm[2*13+2]=0.840998477543771;
	pwm[2*13+3]=1.00;
	pwm[2*13+4]=0.00;
	pwm[2*13+5]=0.369373255518904;
	pwm[2*13+6]=0.11949378330373;
	pwm[2*13+7]=0.830106571936057;
	pwm[2*13+8]=0.188219994925146;
	pwm[2*13+9]=0.304884547069272;
	pwm[2*13+10]=0.23943161634103;
	pwm[2*13+11]=0.245045673686882;
	pwm[2*13+12]=0.258709718345598;

	pwm[3*13+0]=0.0998160365389495;
	pwm[3*13+1]=0.147709972088302;
	pwm[3*13+2]=0.0517571682314133;
	pwm[3*13+3]=0.00;
	pwm[3*13+4]=1.00;
	pwm[3*13+5]=0.0115452930728242;
	pwm[3*13+6]=0.0931489469677747;
	pwm[3*13+7]=0.0566861202740421;
	pwm[3*13+8]=0.483259325044405;
	pwm[3*13+9]=0.207954833798528;
	pwm[3*13+10]=0.287972595787871;
	pwm[3*13+11]=0.272614818573966;
	pwm[3*13+12]=0.281736868815022;

	int i;
	for(i=0; i<4; i++){
		int j;
		for(j=0; j<13; j++){
			pwm[i*13+j]=pwm[i*13+j]+0.00001f;
		}
	}

	return pwm;
}
double *GetPWMfor5PrimeGCAGU2(){
	double *pwm=NULL;

	pwm=NPALLOC(double, 4*14);
	my_assert(pwm != NULL);

	pwm[0*14+0]=0.402203856749311;
	pwm[0*14+1]=0.873278236914601;
	pwm[0*14+2]=0.0199724517906336;
	pwm[0*14+3]=0.00;
	pwm[0*14+4]=0.00;
	pwm[0*14+5]=0.924931129476584;
	pwm[0*14+6]=0.831955922865014;
	pwm[0*14+7]=0.00619834710743802;
	pwm[0*14+8]=0.075068870523416;
	pwm[0*14+9]=0.330578512396694;
	pwm[0*14+10]=0.192148760330579;
	pwm[0*14+11]=0.194214876033058;
	pwm[0*14+12]=0.245179063360882;
	pwm[0*14+13]=0.40;

	pwm[1*14+0]=0.368457300275482;
	pwm[1*14+1]=0.0130853994490358;
	pwm[1*14+2]=0.00206611570247934;
	pwm[1*14+3]=0.00;
	pwm[1*14+4]=1.00;
	pwm[1*14+5]=0.0137741046831956;
	pwm[1*14+6]=0.0254820936639118;
	pwm[1*14+7]=0.00206611570247934;
	pwm[1*14+8]=0.0867768595041322;
	pwm[1*14+9]=0.158402203856749;
	pwm[1*14+10]=0.272038567493113;
	pwm[1*14+11]=0.305785123966942;
	pwm[1*14+12]=0.214187327823691;
	pwm[1*14+13]=0.20;

	pwm[2*14+0]=0.172176308539945;
	pwm[2*14+1]=0.0433884297520661;
	pwm[2*14+2]=0.974517906336088;
	pwm[2*14+3]=1.00;
	pwm[2*14+4]=0.00;
	pwm[2*14+5]=0.0564738292011019;
	pwm[2*14+6]=0.087465564738292;
	pwm[2*14+7]=0.988292011019284;
	pwm[2*14+8]=0.0929752066115702;
	pwm[2*14+9]=0.34228650137741;
	pwm[2*14+10]=0.213498622589532;
	pwm[2*14+11]=0.221763085399449;
	pwm[2*14+12]=0.260330578512397;
	pwm[2*14+13]=0.20;

	pwm[3*14+0]=0.0571625344352617;
	pwm[3*14+1]=0.0702479338842975;
	pwm[3*14+2]=0.0034435261707989;
	pwm[3*14+3]=0.00;
	pwm[3*14+4]=0.00;
	pwm[3*14+5]=0.00482093663911846;
	pwm[3*14+6]=0.0550964187327824;
	pwm[3*14+7]=0.0034435261707989;
	pwm[3*14+8]=0.745179063360882;
	pwm[3*14+9]=0.168732782369146;
	pwm[3*14+10]=0.322314049586777;
	pwm[3*14+11]=0.278236914600551;
	pwm[3*14+12]=0.28030303030303;
	pwm[3*14+13]=0.20;

	int i;
	for(i=0; i<4; i++){
		int j;
		for(j=0; j<14; j++){
			pwm[i*14+j]=pwm[i*14+j]+0.00001f;
		}
	}

	return pwm;
}
double *GetPWMfor3PrimeGTAGU12(){
	double *pwm=NULL;

	pwm=NPALLOC(double, 4*18);
	my_assert(pwm != NULL);

	pwm[0*18+0]=0.347826086956522;
	pwm[0*18+1]=0.331521739130435;
	pwm[0*18+2]=0.239130434782609;
	pwm[0*18+3]=0.298913043478261;
	pwm[0*18+4]=0.293478260869565;
	pwm[0*18+5]=0.228260869565217;
	pwm[0*18+6]=0.168478260869565;
	pwm[0*18+7]=0.260869565217391;
	pwm[0*18+8]=0.16304347826087;
	pwm[0*18+9]=0.195652173913043;
	pwm[0*18+10]=0.16304347826087;
	pwm[0*18+11]=0.119565217391304;
	pwm[0*18+12]=1.00;
	pwm[0*18+13]=0.00;
	pwm[0*18+14]=0.494565217391304;
	pwm[0*18+15]=0.130434782608696;
	pwm[0*18+16]=0.41304347826087;
	pwm[0*18+17]=0.40;

	pwm[1*18+0]=0.271739130434783;
	pwm[1*18+1]=0.282608695652174;
	pwm[1*18+2]=0.326086956521739;
	pwm[1*18+3]=0.266304347826087;
	pwm[1*18+4]=0.271739130434783;
	pwm[1*18+5]=0.217391304347826;
	pwm[1*18+6]=0.255434782608696;
	pwm[1*18+7]=0.25;
	pwm[1*18+8]=0.358695652173913;
	pwm[1*18+9]=0.309782608695652;
	pwm[1*18+10]=0.28804347826087;
	pwm[1*18+11]=0.652173913043478;
	pwm[1*18+12]=0.00;
	pwm[1*18+13]=0.00;
	pwm[1*18+14]=0.130434782608696;
	pwm[1*18+15]=0.217391304347826;
	pwm[1*18+16]=0.135869565217391;
	pwm[1*18+17]=0.20;

	pwm[2*18+0]=0.114130434782609;
	pwm[2*18+1]=0.0923913043478261;
	pwm[2*18+2]=0.0978260869565217;
	pwm[2*18+3]=0.146739130434783;
	pwm[2*18+4]=0.125;
	pwm[2*18+5]=0.146739130434783;
	pwm[2*18+6]=0.135869565217391;
	pwm[2*18+7]=0.141304347826087;
	pwm[2*18+8]=0.141304347826087;
	pwm[2*18+9]=0.0760869565217391;
	pwm[2*18+10]=0.141304347826087;
	pwm[2*18+11]=0.016304347826087;
	pwm[2*18+12]=0.00;
	pwm[2*18+13]=1.00;
	pwm[2*18+14]=0.168478260869565;
	pwm[2*18+15]=0.146739130434783;
	pwm[2*18+16]=0.201086956521739;
	pwm[2*18+17]=0.20;

	pwm[3*18+0]=0.266304347826087;
	pwm[3*18+1]=0.293478260869565;
	pwm[3*18+2]=0.33695652173913;
	pwm[3*18+3]=0.28804347826087;
	pwm[3*18+4]=0.309782608695652;
	pwm[3*18+5]=0.407608695652174;
	pwm[3*18+6]=0.440217391304348;
	pwm[3*18+7]=0.347826086956522;
	pwm[3*18+8]=0.33695652173913;
	pwm[3*18+9]=0.418478260869565;
	pwm[3*18+10]=0.407608695652174;
	pwm[3*18+11]=0.21195652173913;
	pwm[3*18+12]=0.00;
	pwm[3*18+13]=0.00;
	pwm[3*18+14]=0.206521739130435;
	pwm[3*18+15]=0.505434782608696;
	pwm[3*18+16]=0.25;
	pwm[3*18+17]=0.20;

	int i;
	for(i=0; i<4; i++){
		int j;
		for(j=0; j<18; j++){
			pwm[i*18+j]=pwm[i*18+j]+0.00001f;
		}
	}

	return pwm;
}
double *GetPWMfor3PrimeATACU12(){
	double *pwm=NULL;

	pwm=NPALLOC(double, 4*17);
	my_assert(pwm != NULL);

	pwm[0*17+0]=0.214953271028037;
	pwm[0*17+1]=0.355140186915888;
	pwm[0*17+2]=0.495327102803738;
	pwm[0*17+3]=0.327102803738318;
	pwm[0*17+4]=0.233644859813084;
	pwm[0*17+5]=0.149532710280374;
	pwm[0*17+6]=0.168224299065421;
	pwm[0*17+7]=0.224299065420561;
	pwm[0*17+8]=0.233644859813084;
	pwm[0*17+9]=0.0934579439252336;
	pwm[0*17+10]=0.14018691588785;
	pwm[0*17+11]=0.102803738317757;
	pwm[0*17+12]=1.00;
	pwm[0*17+13]=0.00;
	pwm[0*17+14]=0.392523364485981;
	pwm[0*17+15]=0.0654205607476635;
	pwm[0*17+16]=0.299065420560748;

	pwm[1*17+0]=0.196261682242991;
	pwm[1*17+1]=0.0934579439252336;
	pwm[1*17+2]=0.196261682242991;
	pwm[1*17+3]=0.373831775700935;
	pwm[1*17+4]=0.345794392523364;
	pwm[1*17+5]=0.317757009345794;
	pwm[1*17+6]=0.252336448598131;
	pwm[1*17+7]=0.214953271028037;
	pwm[1*17+8]=0.233644859813084;
	pwm[1*17+9]=0.289719626168224;
	pwm[1*17+10]=0.289719626168224;
	pwm[1*17+11]=0.476635514018692;
	pwm[1*17+12]=0.00;
	pwm[1*17+13]=1.00;
	pwm[1*17+14]=0.252336448598131;
	pwm[1*17+15]=0.149532710280374;
	pwm[1*17+16]=0.168224299065421;

	pwm[2*17+0]=0.0467289719626168;
	pwm[2*17+1]=0.14018691588785;
	pwm[2*17+2]=0.0934579439252336;
	pwm[2*17+3]=0.0654205607476635;
	pwm[2*17+4]=0.0747663551401869;
	pwm[2*17+5]=0.214953271028037;
	pwm[2*17+6]=0.186915887850467;
	pwm[2*17+7]=0.168224299065421;
	pwm[2*17+8]=0.224299065420561;
	pwm[2*17+9]=0.177570093457944;
	pwm[2*17+10]=0.196261682242991;
	pwm[2*17+11]=0.0186915887850467;
	pwm[2*17+12]=0.00;
	pwm[2*17+13]=0.00;
	pwm[2*17+14]=0.149532710280374;
	pwm[2*17+15]=0.196261682242991;
	pwm[2*17+16]=0.224299065420561;

	pwm[3*17+0]=0.542056074766355;
	pwm[3*17+1]=0.411214953271028;
	pwm[3*17+2]=0.214953271028037;
	pwm[3*17+3]=0.233644859813084;
	pwm[3*17+4]=0.345794392523364;
	pwm[3*17+5]=0.317757009345794;
	pwm[3*17+6]=0.392523364485981;
	pwm[3*17+7]=0.392523364485981;
	pwm[3*17+8]=0.308411214953271;
	pwm[3*17+9]=0.439252336448598;
	pwm[3*17+10]=0.373831775700935;
	pwm[3*17+11]=0.401869158878505;
	pwm[3*17+12]=0.00;
	pwm[3*17+13]=0.00;
	pwm[3*17+14]=0.205607476635514;
	pwm[3*17+15]=0.588785046728972;
	pwm[3*17+16]=0.308411214953271;

	int i;
	for(i=0; i<4; i++){
		int j;
		for(j=0; j<17; j++){
			pwm[i*17+j]=pwm[i*17+j]+0.00001f;
		}
	}

	return pwm;
}
double *GetPWMfor3PrimeGTAGU2(){
	double *pwm=NULL;

	pwm=NPALLOC(double, 4*17);
	my_assert(pwm != NULL);

	pwm[0*17+0]=0.113239025628013;
	pwm[0*17+1]=0.101370210606445;
	pwm[0*17+2]=0.0929903577772139;
	pwm[0*17+3]=0.084775437706166;
	pwm[0*17+4]=0.0874714539456991;
	pwm[0*17+5]=0.0980842425780259;
	pwm[0*17+6]=0.108665313372241;
	pwm[0*17+7]=0.113410301953819;
	pwm[0*17+8]=0.0862598325298148;
	pwm[0*17+9]=0.0897931996955088;
	pwm[0*17+10]=0.239241309312357;
	pwm[0*17+11]=0.0591791423496574;
	pwm[0*17+12]=1.00;
	pwm[0*17+13]=0.00;
	pwm[0*17+14]=0.257047703628521;
	pwm[0*17+15]=0.246948743973611;
	pwm[0*17+16]=0.258189545800558;

	pwm[1*17+0]=0.280690180157321;
	pwm[1*17+1]=0.276865008880995;
	pwm[1*17+2]=0.279611773661507;
	pwm[1*17+3]=0.259673940624207;
	pwm[1*17+4]=0.282161887845724;
	pwm[1*17+5]=0.295629281908145;
	pwm[1*17+6]=0.323883532098452;
	pwm[1*17+7]=0.333557472722659;
	pwm[1*17+8]=0.341004821111393;
	pwm[1*17+9]=0.296415884293327;
	pwm[1*17+10]=0.272532352194874;
	pwm[1*17+11]=0.646929713270743;
	pwm[1*17+12]=0.00;
	pwm[1*17+13]=0.00;
	pwm[1*17+14]=0.141861202740421;
	pwm[1*17+15]=0.189552144125856;
	pwm[1*17+16]=0.234914996193859;

	pwm[2*17+0]=0.123471200202994;
	pwm[2*17+1]=0.116512306521188;
	pwm[2*17+2]=0.107288759198173;
	pwm[2*17+3]=0.102309058614565;
	pwm[2*17+4]=0.108855620400913;
	pwm[2*17+5]=0.11354351687389;
	pwm[2*17+6]=0.103615833544786;
	pwm[2*17+7]=0.0917660492260848;
	pwm[2*17+8]=0.0640256280131946;
	pwm[2*17+9]=0.0628457244354225;
	pwm[2*17+10]=0.204516620147171;
	pwm[2*17+11]=0.00315909667597057;
	pwm[2*17+12]=0.00;
	pwm[2*17+13]=1.00;
	pwm[2*17+14]=0.489228622177112;
	pwm[2*17+15]=0.193795990865263;
	pwm[2*17+16]=0.236094899771632;

	pwm[3*17+0]=0.482599594011672;
	pwm[3*17+1]=0.505252473991373;
	pwm[3*17+2]=0.520109109363106;
	pwm[3*17+3]=0.553241563055062;
	pwm[3*17+4]=0.521511037807663;
	pwm[3*17+5]=0.492742958639939;
	pwm[3*17+6]=0.463835320984522;
	pwm[3*17+7]=0.461266176097437;
	pwm[3*17+8]=0.508709718345598;
	pwm[3*17+9]=0.550945191575742;
	pwm[3*17+10]=0.283709718345598;
	pwm[3*17+11]=0.290732047703629;
	pwm[3*17+12]=0.00;
	pwm[3*17+13]=0.00;
	pwm[3*17+14]=0.111862471453946;
	pwm[3*17+15]=0.36970312103527;
	pwm[3*17+16]=0.270800558233951;

	int i;
	for(i=0; i<4; i++){
		int j;
		for(j=0; j<17; j++){
			pwm[i*17+j]=pwm[i*17+j]+0.00001f;
		}
	}

	return pwm;
}

double *GetPWMfor3PrimeGCAGU2(){
	double *pwm=NULL;

	pwm=NPALLOC(double, 4*18);
	my_assert(pwm != NULL);

	pwm[0*18+0]=0.128787878787879;
	pwm[0*18+1]=0.115702479338843;
	pwm[0*18+2]=0.110881542699725;
	pwm[0*18+3]=0.0819559228650138;
	pwm[0*18+4]=0.0943526170798898;
	pwm[0*18+5]=0.0984848484848485;
	pwm[0*18+6]=0.122589531680441;
	pwm[0*18+7]=0.127410468319559;
	pwm[0*18+8]=0.0881542699724518;
	pwm[0*18+9]=0.0902203856749311;
	pwm[0*18+10]=0.25137741046832;
	pwm[0*18+11]=0.0440771349862259;
	pwm[0*18+12]=1.00;
	pwm[0*18+13]=0.00;
	pwm[0*18+14]=0.269972451790634;
	pwm[0*18+15]=0.261019283746556;
	pwm[0*18+16]=0.238292011019284;
	pwm[0*18+17]=0.40;

	pwm[1*18+0]=0.288567493112948;
	pwm[1*18+1]=0.254820936639118;
	pwm[1*18+2]=0.28236914600551;
	pwm[1*18+3]=0.232782369146006;
	pwm[1*18+4]=0.272727272727273;
	pwm[1*18+5]=0.287190082644628;
	pwm[1*18+6]=0.329201101928375;
	pwm[1*18+7]=0.339531680440771;
	pwm[1*18+8]=0.351239669421488;
	pwm[1*18+9]=0.28099173553719;
	pwm[1*18+10]=0.252066115702479;
	pwm[1*18+11]=0.644628099173554;
	pwm[1*18+12]=0.00;
	pwm[1*18+13]=0.00;
	pwm[1*18+14]=0.123278236914601;
	pwm[1*18+15]=0.176308539944904;
	pwm[1*18+16]=0.207988980716253;
	pwm[1*18+17]=0.20;

	pwm[2*18+0]=0.121212121212121;
	pwm[2*18+1]=0.115013774104683;
	pwm[2*18+2]=0.108815426997245;
	pwm[2*18+3]=0.101928374655647;
	pwm[2*18+4]=0.12396694214876;
	pwm[2*18+5]=0.119834710743802;
	pwm[2*18+6]=0.114325068870523;
	pwm[2*18+7]=0.0847107438016529;
	pwm[2*18+8]=0.0626721763085399;
	pwm[2*18+9]=0.0743801652892562;
	pwm[2*18+10]=0.196280991735537;
	pwm[2*18+11]=0.00413223140495868;
	pwm[2*18+12]=0.00;
	pwm[2*18+13]=1.00;
	pwm[2*18+14]=0.520661157024793;
	pwm[2*18+15]=0.176308539944904;
	pwm[2*18+16]=0.25;
	pwm[2*18+17]=0.20;

	pwm[3*18+0]=0.461432506887052;
	pwm[3*18+1]=0.514462809917355;
	pwm[3*18+2]=0.497933884297521;
	pwm[3*18+3]=0.583333333333333;
	pwm[3*18+4]=0.508953168044077;
	pwm[3*18+5]=0.494490358126722;
	pwm[3*18+6]=0.433884297520661;
	pwm[3*18+7]=0.448347107438017;
	pwm[3*18+8]=0.497933884297521;
	pwm[3*18+9]=0.554407713498623;
	pwm[3*18+10]=0.300275482093664;
	pwm[3*18+11]=0.307162534435262;
	pwm[3*18+12]=0.00;
	pwm[3*18+13]=0.00;
	pwm[3*18+14]=0.0860881542699724;
	pwm[3*18+15]=0.386363636363636;
	pwm[3*18+16]=0.303719008264463;
	pwm[3*18+17]=0.20;

	int i;
	for(i=0; i<4; i++){
		int j;
		for(j=0; j<18; j++){
			pwm[i*18+j]=pwm[i*18+j]+0.00001f;
		}
	}

	return pwm;
}
double *GetLODforBPSPWM(double *pwm){
	double *lod=NULL;

	lod=NPALLOC(double, 4*12);
	my_assert(lod != NULL);

	int i;
	for(i=0; i<4; i++){
		int j;
		for(j=0; j<12; j++){
			lod[i*12+j]=(double)log((double)(pwm[i*12+j]/0.25f));
			lod[i*12+j]=(double)(lod[i*12+j]/(double)log(2));
		}
	}

	return lod;
}

double *GetLODforPWM(double *pwm, int motif_length){
	double *lod=NULL;

	lod=NPALLOC(double, 4*motif_length);
	my_assert(lod != NULL);

	int i;
	for(i=0; i<4; i++){
		int j;
		for(j=0; j<motif_length; j++){
			lod[i*motif_length+j]=(double)log((double)(pwm[i*motif_length+j]/0.25f));
			lod[i*motif_length+j]=(double)(lod[i*motif_length+j]/(double)log(2));
		}
	}

	return lod;
}

double *GetCVectorForPWM(double *pwm, int motif_length){
	double *Cvector=NULL;

	Cvector=NPALLOC(double, motif_length);
	my_assert(Cvector != NULL);

	int i;
	for(i=0; i<motif_length; i++){
		Cvector[i]=0;
		int j=0;
		while(j<4){
			Cvector[i]+=pwm[j*motif_length+i]*log(pwm[j*motif_length+i]);
			j++;
		}

		Cvector[i]+=log(5.0f);
		Cvector[i]*=(100.0f/log(5.0f));
	}

	return Cvector;
}

double *GetMAXVectorForPWM(double *pwm, int motif_length){
	double *MAXvector=NULL;

	MAXvector=NPALLOC(double, motif_length);
	my_assert(MAXvector != NULL);

	int i;
	for(i=0; i<motif_length; i++){
		MAXvector[i]=0.0f;
		int j=0;
		while(j<4){
			if(pwm[j*motif_length+i] > MAXvector[i])
				MAXvector[i]=pwm[j*motif_length+i];
			j++;
		}
	}

	return MAXvector;
}

/*
 * Loading of the 10 PWM matrices
 *
 * Index 0 ==> pwm_9
 * Index 1 ==> pwm_10
 * Index 2 ==> pwm_5PrimeGTAGU12
 * Index 3 ==> pwm_5PrimeATACU12
 * Index 4 ==> pwm_5PrimeGTAGU2
 * Index 5 ==> pwm_5PrimeGCAGU2
 * Index 6 ==> pwm_3PrimeGTAGU12
 * Index 7 ==> pwm_3PrimeATACU12
 * Index 8 ==> pwm_3PrimeGTAGU2
 * Index 9 ==> pwm_3PrimeGCAGU2
 */

double **LoadPWMMatrices(){

	double **pwm_matx=NPALLOC(double *, 10);
	my_assert(pwm_matx != NULL);

	pwm_matx[0]=GetPWMforBPS_9();
	pwm_matx[1]=GetPWMforBPS_10();
	pwm_matx[2]=GetPWMfor5PrimeGTAGU12();
	pwm_matx[3]=GetPWMfor5PrimeATACU12();
	pwm_matx[4]=GetPWMfor5PrimeGTAGU2();
	pwm_matx[5]=GetPWMfor5PrimeGCAGU2();
	pwm_matx[6]=GetPWMfor3PrimeGTAGU12();
	pwm_matx[7]=GetPWMfor3PrimeATACU12();
	pwm_matx[8]=GetPWMfor3PrimeGTAGU2();
	pwm_matx[9]=GetPWMfor3PrimeGCAGU2();

	return pwm_matx;
}

void FreePWMMatrices(double **pwm_matx){

	int i;
	for(i=0; i < 10; i++)
		pfree(pwm_matx[i]);

	pfree(pwm_matx);
}

/*
 * Loading of the 10 CV vectors
 *
 * Index 0 ==> CV vector for pwm_9
 * Index 1 ==> CV vector for pwm_10
 * Index 2 ==> CV vector for pwm_5PrimeGTAGU12
 * Index 3 ==> CV vector for pwm_5PrimeATACU12
 * Index 4 ==> CV vector for pwm_5PrimeGTAGU2
 * Index 5 ==> CV vector for pwm_5PrimeGCAGU2
 * Index 6 ==> CV vector for pwm_3PrimeGTAGU12
 * Index 7 ==> CV vector for pwm_3PrimeATACU12
 * Index 8 ==> CV vector for pwm_3PrimeGTAGU2
 * Index 9 ==> CV vector for pwm_3PrimeGCAGU2
 */

double **LoadCVPWMMatrices(double **pwm_matx){

	my_assert(pwm_matx != NULL);

	double **CV_vector=NPALLOC(double *, 10);
	my_assert(CV_vector != NULL);

	CV_vector[0]=GetCVectorForPWM(pwm_matx[0], 12);
	CV_vector[1]=GetCVectorForPWM(pwm_matx[1], 12);
	CV_vector[2]=GetCVectorForPWM(pwm_matx[2], 14);
	CV_vector[3]=GetCVectorForPWM(pwm_matx[3], 14);
	CV_vector[4]=GetCVectorForPWM(pwm_matx[4], 13);
	CV_vector[5]=GetCVectorForPWM(pwm_matx[5], 14);
	CV_vector[6]=GetCVectorForPWM(pwm_matx[6], 18);
	CV_vector[7]=GetCVectorForPWM(pwm_matx[7], 17);
	CV_vector[8]=GetCVectorForPWM(pwm_matx[8], 17);
	CV_vector[9]=GetCVectorForPWM(pwm_matx[9], 18);

	return CV_vector;
}

void FreeCVPWMMatrices(double **CV_vector){

	int i;
	for(i=0; i < 10; i++)
		pfree(CV_vector[i]);

	pfree(CV_vector);
}

/*
 * Loading of the 10 MAX vectors
 *
 * Index 0 ==> MAX vector for pwm_9
 * Index 1 ==> MAX vector for pwm_10
 * Index 2 ==> MAX vector for pwm_5PrimeGTAGU12
 * Index 3 ==> MAX vector for pwm_5PrimeATACU12
 * Index 4 ==> MAX vector for pwm_5PrimeGTAGU2
 * Index 5 ==> MAX vector for pwm_5PrimeGCAGU2
 * Index 6 ==> MAX vector for pwm_3PrimeGTAGU12
 * Index 7 ==> MAX vector for pwm_3PrimeATACU12
 * Index 8 ==> MAX vector for pwm_3PrimeGTAGU2
 * Index 9 ==> MAX vector for pwm_3PrimeGCAGU2
 */

double **LoadMAXPWMMatrices(double **pwm_matx){

	my_assert(pwm_matx != NULL);

	double **MAX_vector=NPALLOC(double *, 10);
	my_assert(MAX_vector != NULL);

	MAX_vector[0]=GetMAXVectorForPWM(pwm_matx[0], 12);
	MAX_vector[1]=GetMAXVectorForPWM(pwm_matx[1], 12);
	MAX_vector[2]=GetMAXVectorForPWM(pwm_matx[2], 14);
	MAX_vector[3]=GetMAXVectorForPWM(pwm_matx[3], 14);
	MAX_vector[4]=GetMAXVectorForPWM(pwm_matx[4], 13);
	MAX_vector[5]=GetMAXVectorForPWM(pwm_matx[5], 14);
	MAX_vector[6]=GetMAXVectorForPWM(pwm_matx[6], 18);
	MAX_vector[7]=GetMAXVectorForPWM(pwm_matx[7], 17);
	MAX_vector[8]=GetMAXVectorForPWM(pwm_matx[8], 17);
	MAX_vector[9]=GetMAXVectorForPWM(pwm_matx[9], 18);

	return MAX_vector;
}

void FreeMAXPWMMatrices(double **MAX_vector){

	int i;
	for(i=0; i < 10; i++)
		pfree(MAX_vector[i]);

	pfree(MAX_vector);
}

pgenomic_intron set_pattern(char *genomic_sequence, pgenomic_intron gen_intron){

#ifndef NDEBUG
	my_assert(gen_intron != NULL);
	my_assert(genomic_sequence != NULL);

	size_t gen_length=strlen(genomic_sequence);

	my_assert(gen_intron->start >= 0 && gen_intron->start < (int)gen_length);
	my_assert(gen_intron->end >= 0 && gen_intron->end < (int)gen_length);
#endif

	gen_intron->donor_pt=real_substring(gen_intron->start, 2, genomic_sequence);
	gen_intron->acceptor_pt=real_substring(gen_intron->end-1, 2, genomic_sequence);

	return gen_intron;
}

char *get_donor_EST_suffix(char *est_sequence, int EST_cut, int suffix_dim){

	char *donor_suffix=real_substring(EST_cut-suffix_dim, suffix_dim, est_sequence);

	return donor_suffix;
}

char *get_acceptor_EST_prefix(char *est_sequence, int EST_cut, int prefix_dim){

	char *acceptor_prefix=real_substring(EST_cut, prefix_dim, est_sequence);

	return acceptor_prefix;
}

char *get_donor_suffix(char *genomic_sequence, pgenomic_intron gen_intron, int suffix_dim){

	char *donor_suffix=real_substring(gen_intron->start-suffix_dim, suffix_dim, genomic_sequence);

	return donor_suffix;
}

char *get_acceptor_prefix(char *genomic_sequence, pgenomic_intron gen_intron, int prefix_dim){

	char *acceptor_prefix=real_substring(gen_intron->end+1, prefix_dim, genomic_sequence);

	return acceptor_prefix;
}

char *get_intron_suffix(char *genomic_sequence, pgenomic_intron gen_intron, int suffix_dim){

	char *intron_suffix=real_substring(gen_intron->end-suffix_dim+1, suffix_dim, genomic_sequence);

	return intron_suffix;
}

char *get_intron_prefix(char *genomic_sequence, pgenomic_intron gen_intron, int prefix_dim){

	char *intron_prefix=real_substring(gen_intron->start, prefix_dim, genomic_sequence);

	return intron_prefix;
}

char *GetRepeatSequence(char *genomic_sequence, int intron_left, int intron_right){
	my_assert(genomic_sequence != NULL);
	my_assert(intron_left < intron_right);

	int i=intron_left-1;
	while(genomic_sequence[i] == genomic_sequence[intron_right-intron_left+i+1]){
		i--;
	}

	//If a repeat sequence does exist in 5'
	char *five_repeat=NULL;
	int f_r_l=0;
	if(intron_left-i-1 > 0){
		five_repeat=real_substring(i+1, intron_left-i-1, genomic_sequence);
		f_r_l=intron_left-i-1;
	}

	i=intron_right+1;
	while(genomic_sequence[i] == genomic_sequence[-intron_right+intron_left+i-1]){
		i++;
	}

	//If a repeat sequence does exist in 3'
	char *three_repeat=NULL;
	int t_r_l=0;
	if(i-intron_right-1 > 0){
		three_repeat=real_substring(intron_right+1, i-intron_right-1, genomic_sequence);
		t_r_l=i-intron_right-1;
	}

	int repeat_length=f_r_l+t_r_l;

	char *repeat=NULL;
	if(repeat_length > 0){
		repeat=c_palloc(repeat_length+1);
		strcpy(repeat, "");
		if(five_repeat != NULL){
			strcat(repeat, five_repeat);
			pfree(five_repeat);
		}
		if(three_repeat != NULL){
			strcat(repeat, three_repeat);
			pfree(three_repeat);
		}
	}

	return repeat;
}
