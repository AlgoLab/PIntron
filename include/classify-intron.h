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
/**
 *
 * @file classifying-intron.h
 *
 * Procedures for classifying an intron.
 *
 **/

#ifndef _CLASSIFY_INTRON_H_
#define _CLASSIFY_INTRON_H_

#include "types.h"
#include "configuration.h"

//Include

/*Classifies an intron from left to right on the genomic sequence.
 * Returns an integer:
 * 0 for U12 introns
 * 1 for U2 introns
 * 2 for unclassified introns.
 *  score5 is the 5' scores, score3 is the 3' score, BPS_position is the
 *  position of the BPS sequence (if it exists) wrt to the intron left (0 index) and
 *  BPS_score is the computed BPS score.
 */
plist classify_genomic_intron_list(char *, plist);
pgenomic_intron classify_genomic_intron(char *, pgenomic_intron);
char classify_genomic_intron_start_end(char *, int, int, double *, double *, int *, double *, double **, double **, double **);

/*Compute the score of the 14-long input sequence wrt to 5' ss of type GT-AG (U12)
==> 3nt inside exon and 11nt inside intron*/

double GetScoreOf5PrimeGTAGU12BySS(char *, int, double *, double *, double *);
double GetScoreOf5PrimeGTAGU12(char *, double *, double *, double *);

/*Compute the score of the 14-long input sequence wrt to 5' ss of type AT-AC (U12)
==> 3nt inside exon and 11nt inside intron*/

double GetScoreOf5PrimeATACU12BySS(char *, int, double *, double *, double *);
double GetScoreOf5PrimeATACU12(char *, double *, double *, double *);

/*Compute the score of the 13-long input sequence wrt to 5' ss of type GT-AG (U2)
==> 3nt inside exon and 10nt inside intron*/

double GetScoreOf5PrimeGTAGU2BySS(char *, int, double *, double *, double *);
double GetScoreOf5PrimeGTAGU2(char *, double *, double *, double *);

/*Compute the score of the 14-long input sequence wrt to 5' ss of type GC-AG (U2)
==> 3nt inside exon and 11nt inside intron*/

double GetScoreOf5PrimeGCAGU2BySS(char *, int, double *, double *, double *);
double GetScoreOf5PrimeGCAGU2(char *, double *, double *, double *);

/*Compute the score of the 18-long input sequence wrt to 3' ss of type GT-AG (U12)
==> 14nt inside intron and 4nt inside exon*/

double GetScoreOf3PrimeGTAGU12BySS(char *, int, double *, double *, double *);
double GetScoreOf3PrimeGTAGU12(char *, double *, double *, double *);

/*Compute the score of the 17-long input sequence wrt to 3' ss of type AT-AC (U12)
==> 14nt inside intron and 3nt inside exon*/

double GetScoreOf3PrimeATACU12BySS(char *, int, double *, double *, double *);
double GetScoreOf3PrimeATACU12(char *, double *, double *, double *);

/*Compute the score of the 17-long input sequence wrt to 3' ss of type GT-AG (U2)
==> 14nt inside intron and 3nt inside exon*/

double GetScoreOf3PrimeGTAGU2BySS(char *, int, double *, double *, double *);
double GetScoreOf3PrimeGTAGU2(char *, double *, double *, double *);

/*Compute the score of the 18-long input sequence wrt to 3' ss of type GC-AG (U2)
==> 14nt inside intron and 4nt inside exon*/

double GetScoreOf3PrimeGCAGU2BySS(char *, int, double *, double *, double *);
double GetScoreOf3PrimeGCAGU2(char *, double *, double *, double *);

/*Search in the input intron sequence the 12-long substring with the best BPS score.
This searching is performed inside a window from the nucleotide in positions (L-range_end)
to the nucleotide in position (L-range_start), where L is the intron length (usually
range_end is 30 bp and range_start is 14 bp).
The return value is the position of the first BPS character wrt the first intron
character (index 0), if this value is -1, then no BPS was found.
The BPS score is returned in the score parameter and is the best (> 0.75) between
BP1 and BP2.
*/

int ExistsGoodBPSinIntronSequenceWithMathInspector(char *, double *, double *, double *, double *, double *, double *, double *, int, int);
int SearchBPSinIntronSequenceWithMathInspector(char *, double *, double *, double *, double *, int, int);
double GetMatInspectorScoreOfaMotif(char *, double *, double *, double *, int);

//Loading of PWM matrices

double *GetPWMforBPS_9();
double *GetPWMforBPS_10();
double *GetPWMfor5PrimeGTAGU12();
double *GetPWMfor5PrimeATACU12();
double *GetPWMfor5PrimeGTAGU2();
double *GetPWMfor5PrimeGCAGU2();
double *GetPWMfor3PrimeGTAGU12();
double *GetPWMfor3PrimeATACU12();
double *GetPWMfor3PrimeGTAGU2();
double *GetPWMfor3PrimeGCAGU2();

//Loading of LOD matrices

double *GetLODforBPSPWM(double *);
double *GetLODforPWM(double *, int);

//Loading of CV vector

double *GetCVectorForPWM(double *, int);

//Loading of MAX vector

double *GetMAXVectorForPWM(double *, int);

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

double **LoadPWMMatrices();
void FreePWMMatrices(double **);

/*
 * Loading of the 10 C vectors
 *
 * Index 0 ==> C vector for pwm_9
 * Index 1 ==> C vector for pwm_10
 * Index 2 ==> C vector for pwm_5PrimeGTAGU12
 * Index 3 ==> C vector for pwm_5PrimeATACU12
 * Index 4 ==> C vector for pwm_5PrimeGTAGU2
 * Index 5 ==> C vector for pwm_5PrimeGCAGU2
 * Index 6 ==> C vector for pwm_3PrimeGTAGU12
 * Index 7 ==> C vector for pwm_3PrimeATACU12
 * Index 8 ==> C vector for pwm_3PrimeGTAGU2
 * Index 9 ==> C vector for pwm_3PrimeGCAGU2
 */

double **LoadCVPWMMatrices(double **);
void FreeCVPWMMatrices(double **);

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

double **LoadMAXPWMMatrices(double **);
void FreeMAXPWMMatrices(double **);

pgenomic_intron set_pattern(char *, pgenomic_intron);

char *get_donor_EST_suffix(char *, int, int);
char *get_acceptor_EST_prefix(char *, int, int);

char *get_donor_suffix(char *, pgenomic_intron, int);
char *get_acceptor_prefix(char *, pgenomic_intron, int);
char *get_intron_suffix(char *, pgenomic_intron, int);
char *get_intron_prefix(char *, pgenomic_intron, int);

char *GetRepeatSequence(char *, int, int);

#endif
