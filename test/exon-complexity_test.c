//gcc exon-complexity_test.c -o exon-complexity_test -l criterion -I '/home/lorenzo/PIntron/include'

#include "list.h"
#include "util.h"
#include <ctype.h>
#include <math.h>
#include "exon-complexity.h"
#include "est-factorizations.h"
#include "log.h"

#include <string.h>
#include "../src/exon-complexity.c"
#include "../src/util.c"
#include <criterion/criterion.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]){
	//for test execution
	struct criterion_test_set *tests = criterion_initialize();
	int result = 0;
	if (criterion_handle_args(argc, argv, true))
	result = !criterion_run_all_tests(tests);
	criterion_finalize(tests);
	return result;
}

/*
	create variables needed for the function dustScoreByLeftAndRight,
	verify that the function return the expected value,
	it was used >= instead of == because the result can be a periodic number
*/
Test(exonComplexityTest,dustScoreByLeftAndRightTest) {
	char *s="aCgTccAAtGTaC";
	int n1=1;
	int n2=11;
	cr_expect(dustScoreByLeftAndRight(s,n1,n2)>=0.101010);
}

/*
	create variable needed for the function dustScore,
	verify that the function return the expected value,
	it was used >= instead of == because the result can be a irrational number
*/
Test(exonComplexityTest,dustScoreTest) {
	char *s="aCgTccAAtGTaC";
	cr_expect(dustScore(s)>=0.139860);
}


/*
	create variables needed for the function getDinucleotideIndex,
	verify that the function return the expected value
*/
Test(exonComplexityTest,getDinucleotideIndexTest) {
	char c1='a';
	char c2='A';
	cr_expect(getDinucleotideIndex(c1,c2)==0);
}

/*
	create variables needed for the function getDinucleotideIndex,
	verify that the function return the expected value
*/
Test(exonComplexityTest,getDinucleotideIndexTest2) {
	char c1='A';
	char c2='c';
	cr_expect(getDinucleotideIndex(c1,c2)==1);
}

/*
	create variables needed for the function getDinucleotideIndex,
	verify that the function return the expected value
*/
Test(exonComplexityTest,getDinucleotideIndexTest3) {
	char c1='a';
	char c2='G';
	cr_expect(getDinucleotideIndex(c1,c2)==2);
}

/*
	create variables needed for the function getDinucleotideIndex,
	verify that the function return the expected value
*/
Test(exonComplexityTest,getDinucleotideIndexTest4) {
	char c1='a';
	char c2='t';
	cr_expect(getDinucleotideIndex(c1,c2)==3);
}

/*
	create variables needed for the function getDinucleotideIndex,
	verify that the function return the expected value
*/
Test(exonComplexityTest,getDinucleotideIndexTest5) {
	char c1='C';
	char c2='A';
	cr_expect(getDinucleotideIndex(c1,c2)==4);
}

/*
	create variables needed for the function getDinucleotideIndex,
	verify that the function return the expected value
*/
Test(exonComplexityTest,getDinucleotideIndexTest6) {
	char c1='c';
	char c2='c';
	cr_expect(getDinucleotideIndex(c1,c2)==5);
}

/*
	create variables needed for the function getDinucleotideIndex,
	verify that the function return the expected value
*/
Test(exonComplexityTest,getDinucleotideIndexTest7) {
	char c1='C';
	char c2='g';
	cr_expect(getDinucleotideIndex(c1,c2)==6);
}

/*
	create variables needed for the function getDinucleotideIndex,
	verify that the function return the expected value
*/
Test(exonComplexityTest,getDinucleotideIndexTest8) {
	char c1='c';
	char c2='T';
	cr_expect(getDinucleotideIndex(c1,c2)==7);
}

/*
	create variables needed for the function getDinucleotideIndex,
	verify that the function return the expected value
*/
Test(exonComplexityTest,getDinucleotideIndexTest9) {
	char c1='G';
	char c2='A';
	cr_expect(getDinucleotideIndex(c1,c2)==8);
}

/*
	create variables needed for the function getDinucleotideIndex,
	verify that the function return the expected value
*/
Test(exonComplexityTest,getDinucleotideIndexTest10) {
	char c1='G';
	char c2='c';
	cr_expect(getDinucleotideIndex(c1,c2)==9);
}

/*
	create variables needed for the function getDinucleotideIndex,
	verify that the function return the expected value
*/
Test(exonComplexityTest,getDinucleotideIndexTest11) {
	char c1='G';
	char c2='G';
	cr_expect(getDinucleotideIndex(c1,c2)==10);
}

/*
	create variables needed for the function getDinucleotideIndex,
	verify that the function return the expected value
*/
Test(exonComplexityTest,getDinucleotideIndexTest12) {
	char c1='G';
	char c2='t';
	cr_expect(getDinucleotideIndex(c1,c2)==11);
}

/*
	create variables needed for the function getDinucleotideIndex,
	verify that the function return the expected value
*/
Test(exonComplexityTest,getDinucleotideIndexTest13) {
	char c1='t';
	char c2='A';
	cr_expect(getDinucleotideIndex(c1,c2)==12);
}

/*
	create variables needed for the function getDinucleotideIndex,
	verify that the function return the expected value
*/
Test(exonComplexityTest,getDinucleotideIndexTest14) {
	char c1='T';
	char c2='c';
	cr_expect(getDinucleotideIndex(c1,c2)==13);
}

/*
	create variables needed for the function getDinucleotideIndex,
	verify that the function return the expected value
*/
Test(exonComplexityTest,getDinucleotideIndexTest15) {
	char c1='T';
	char c2='G';
	cr_expect(getDinucleotideIndex(c1,c2)==14);
}

/*
	create variables needed for the function getDinucleotideIndex,
	verify that the function return the expected value
*/
Test(exonComplexityTest,getDinucleotideIndexTest16) {
	char c1='t';
	char c2='t';
	cr_expect(getDinucleotideIndex(c1,c2)==15);
}

/*
	create variables needed for the function getDinucleotideIndex,
	verify that the function return the expected value
*/
Test(exonComplexityTest,getDinucleotideIndexTest17) {
	char c1='z';
	char c2='k';
	cr_expect(getDinucleotideIndex(c1,c2)==16);
}
