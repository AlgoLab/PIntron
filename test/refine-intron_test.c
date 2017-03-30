//gcc refine-intron_test.c -o refine-intron_test -l criterion -I '/home/lorenzo/PIntron/include'

#define _GNU_SOURCE
#include "../src/refine.c"
#include "../src/refine-intron.c"
#include "../src/util.c"
#include "../src/list.c"
#include "../src/types.c"
#include "../src/bit_vector.c"
#include "../src/bool_list.c"
#include "../src/ext_array.c"
#include <criterion/criterion.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>

#include <ctype.h>
#include "refine-intron.h"
#include "refine.h"
#include "est-factorizations.h"
#include "list.h"
#include "types.h"
#include "log.h"

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
	create variables needed for the function Check_Burset_patterns,
	the Check_Burset_patterns function returns the frequency of the pairs,
	verify that the function return the expected value,
	the expected value is 0 because the pairs in this test case does not exist
*/
Test(refineIntronTest,checkBursetPatternsTest) {
	char *gs="acctgtactac";
	int nl=2;
	int nr=5;
	cr_expect(Check_Burset_patterns(gs,nl,nr)==0);
}

/*
	create variables needed for the function Check_Burset_patterns,
	the Check_Burset_patterns function returns the frequency of the pairs,
	verify that the function return the expected value,
	the expected value is 0 because the pairs in this test case does not exist
*/
Test(refineIntronTest,checkBursetPatternsTest2) {
	char *gs="acctgtaaaaaaaaaactac";
	int nl=9;
	int nr=19;
	cr_expect(Check_Burset_patterns(gs,nl,nr)==0);
}

/*
	create variables needed for the function Check_Burset_patterns,
	the Check_Burset_patterns function returns the frequency of the pairs,
	verify that the function return the expected value,
	the expected value is 0 because the pairs in this test case does not exist
*/
Test(refineIntronTest,checkBursetPatternsTest3) {
	char *gs="acctgtactac";
	int nl=9;
	int nr=2;
	cr_expect(Check_Burset_patterns(gs,nl,nr)==0);
}

/*
	create variables needed for the function Check_Burset_patterns,
	the Check_Burset_patterns function returns the frequency of the pairs,
	verify that the function return the expected value,
	the expected value is 0 because the pairs in this test case does not exist
*/
Test(refineIntronTest,checkBursetPatternsTest4) {
	char *gs="acctgaaagtacagtaaac";
	int nl=1;
	int nr=30;
	cr_expect(Check_Burset_patterns(gs,nl,nr)==0);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with wrong values,
	verify that this will cause a failure
*/
Test(refineIntronTest,getBursetFrequencyErrorTest,.signal = SIGSEGV) {
	char *gsd="acctgaaagtacagtaaac";
	char *gsa="acctgaaagtacagtaaac";
	cr_expect(getBursetFrequency(gsd,gsa));
}

/*
	create variables needed for the function getBursetFrequency,
	leave these variables empty,
	verify that the function with these variables returns 0
*/
Test(refineIntronTest,getBursetFrequencyEmptyTest) {
	char *gsd="";
	char *gsa="";
	cr_expect(getBursetFrequency(gsd,gsa)==0);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with wrong values,
	verify that the function with these variables returns 0
*/
Test(refineIntronTest,getBursetFrequencyWrongTest) {
	char *gsd;
	asprintf(&gsd, "Az");
	char *gsa;
	asprintf(&gsa, "yK");
	cr_expect(getBursetFrequency(gsd,gsa)==0);
}

/*
	create an example string written in uppercase,
	use To_lower function on it,
	compare it with an equivalent string written in lowercase
*/
Test(refineIntronTest,toLowerTest) {
	char s=*"ACGTACTGGA";
	cr_expect(strcmp(To_lower(&s),"acgtactgga"));
}

/*
	create an example string written in lowercase,
	use To_upper function on it,
	compare it with an equivalent string written in uppercase
*/
Test(refineIntronTest,toUpperTest) {
	char s=*"acgtactgga";
	cr_expect(strcmp(To_upper(&s),"ACGTACTGGA"));
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest) {
	char *gsd;
	asprintf(&gsd, "AA");
	char *gsa;
	asprintf(&gsa, "AG");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest2) {
	char *gsd;
	asprintf(&gsd, "AA");
	char *gsa;
	asprintf(&gsa, "AT");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest3) {
	char *gsd;
	asprintf(&gsd, "aa");
	char *gsa;
	asprintf(&gsa, "gt");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest4) {
	char *gsd;
	asprintf(&gsd, "AC");
	char *gsa;
	asprintf(&gsa, "CC");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest5) {
	char *gsd;
	asprintf(&gsd, "AG");
	char *gsa;
	asprintf(&gsa, "Ac");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 5
*/
Test(refineIntronTest,getBursetFrequencyTest6) {
	char *gsd;
	asprintf(&gsd, "aG");
	char *gsa;
	asprintf(&gsa, "AG");
	cr_expect(getBursetFrequency(gsd,gsa)==5);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 2
*/
Test(refineIntronTest,getBursetFrequencyTest7) {
	char *gsd;
	asprintf(&gsd, "Ag");
	char *gsa;
	asprintf(&gsa, "cT");
	cr_expect(getBursetFrequency(gsd,gsa)==2);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest8) {
	char *gsd;
	asprintf(&gsd, "AG");
	char *gsa;
	asprintf(&gsa, "gc");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 2
*/
Test(refineIntronTest,getBursetFrequencyTest9) {
	char *gsd;
	asprintf(&gsd, "aG");
	char *gsa;
	asprintf(&gsa, "Tg");
	cr_expect(getBursetFrequency(gsd,gsa)==2);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest10) {
	char *gsd;
	asprintf(&gsd, "at");
	char *gsa;
	asprintf(&gsa, "AA");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}


/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest11) {
	char *gsd;
	asprintf(&gsd, "At");
	char *gsa;
	asprintf(&gsa, "AA");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 8
*/
Test(refineIntronTest,getBursetFrequencyTest12) {
	char *gsd;
	asprintf(&gsd, "AT");
	char *gsa;
	asprintf(&gsa, "aC");
	cr_expect(getBursetFrequency(gsd,gsa)==8);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 7
*/
Test(refineIntronTest,getBursetFrequencyTest13) {
	char *gsd;
	asprintf(&gsd, "AT");
	char *gsa;
	asprintf(&gsa, "Ag");
	cr_expect(getBursetFrequency(gsd,gsa)==7);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 2
*/
Test(refineIntronTest,getBursetFrequencyTest14) {
	char *gsd;
	asprintf(&gsd, "At");
	char *gsa;
	asprintf(&gsa, "At");
	cr_expect(getBursetFrequency(gsd,gsa)==2);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest15) {
	char *gsd;
	asprintf(&gsd, "At");
	char *gsa;
	asprintf(&gsa, "gc");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest16) {
	char *gsd;
	asprintf(&gsd, "At");
	char *gsa;
	asprintf(&gsa, "Gt");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest17) {
	char *gsd;
	asprintf(&gsd, "CA");
	char *gsa;
	asprintf(&gsa, "Ag");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest18) {
	char *gsd;
	asprintf(&gsd, "CA");
	char *gsa;
	asprintf(&gsa, "TT");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 2
*/
Test(refineIntronTest,getBursetFrequencyTest19) {
	char *gsd;
	asprintf(&gsd, "CC");
	char *gsa;
	asprintf(&gsa, "Ag");
	cr_expect(getBursetFrequency(gsd,gsa)==2);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest20) {
	char *gsd;
	asprintf(&gsd, "CG");
	char *gsa;
	asprintf(&gsa, "Ag");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest21) {
	char *gsd;
	asprintf(&gsd, "CG");
	char *gsa;
	asprintf(&gsa, "CA");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 2
*/
Test(refineIntronTest,getBursetFrequencyTest22) {
	char *gsd;
	asprintf(&gsd, "CT");
	char *gsa;
	asprintf(&gsa, "AC");
	cr_expect(getBursetFrequency(gsd,gsa)==2);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest23) {
	char *gsd;
	asprintf(&gsd, "CT");
	char *gsa;
	asprintf(&gsa, "ca");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 8
*/
Test(refineIntronTest,getBursetFrequencyTest24) {
	char *gsd;
	asprintf(&gsd, "ga");
	char *gsa;
	asprintf(&gsa, "Ag");
	cr_expect(getBursetFrequency(gsd,gsa)==8);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest25) {
	char *gsd;
	asprintf(&gsd, "ga");
	char *gsa;
	asprintf(&gsa, "gt");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest26) {
	char *gsd;
	asprintf(&gsd, "ga");
	char *gsa;
	asprintf(&gsa, "Tc");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest27) {
	char *gsd;
	asprintf(&gsd, "ga");
	char *gsa;
	asprintf(&gsa, "tG");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 126
*/
Test(refineIntronTest,getBursetFrequencyTest28) {
	char *gsd;
	asprintf(&gsd, "gc");
	char *gsa;
	asprintf(&gsa, "aG");
	cr_expect(getBursetFrequency(gsd,gsa)==126);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest29) {
	char *gsd;
	asprintf(&gsd, "ga");
	char *gsa;
	asprintf(&gsa, "gt");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest30) {
	char *gsd;
	asprintf(&gsd, "gC");
	char *gsa;
	asprintf(&gsa, "gG");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest31) {
	char *gsd;
	asprintf(&gsd, "gC");
	char *gsa;
	asprintf(&gsa, "TA");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest32) {
	char *gsd;
	asprintf(&gsd, "gG");
	char *gsa;
	asprintf(&gsa, "AC");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 11
*/
Test(refineIntronTest,getBursetFrequencyTest33) {
	char *gsd;
	asprintf(&gsd, "gG");
	char *gsa;
	asprintf(&gsa, "AG");
	cr_expect(getBursetFrequency(gsd,gsa)==11);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest34) {
	char *gsd;
	asprintf(&gsd, "gG");
	char *gsa;
	asprintf(&gsa, "Ca");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 2
*/
Test(refineIntronTest,getBursetFrequencyTest35) {
	char *gsd;
	asprintf(&gsd, "gG");
	char *gsa;
	asprintf(&gsa, "ga");
	cr_expect(getBursetFrequency(gsd,gsa)==2);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 2
*/
Test(refineIntronTest,getBursetFrequencyTest36) {
	char *gsd;
	asprintf(&gsd, "gG");
	char *gsa;
	asprintf(&gsa, "tc");
	cr_expect(getBursetFrequency(gsd,gsa)==2);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 200
*/
Test(refineIntronTest,getBursetFrequencyTest37) {
	char *gsd;
	asprintf(&gsd, "gt");
	char *gsa;
	asprintf(&gsa, "AG");
	cr_expect(getBursetFrequency(gsd,gsa)==200);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 4
*/
Test(refineIntronTest,getBursetFrequencyTest38) {
	char *gsd;
	asprintf(&gsd, "gT");
	char *gsa;
	asprintf(&gsa, "Ac");
	cr_expect(getBursetFrequency(gsd,gsa)==4);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 2
*/
Test(refineIntronTest,getBursetFrequencyTest39) {
	char *gsd;
	asprintf(&gsd, "gt");
	char *gsa;
	asprintf(&gsa, "At");
	cr_expect(getBursetFrequency(gsd,gsa)==2);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 9
*/
Test(refineIntronTest,getBursetFrequencyTest40) {
	char *gsd;
	asprintf(&gsd, "Gt");
	char *gsa;
	asprintf(&gsa, "ca");
	cr_expect(getBursetFrequency(gsd,gsa)==9);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 4
*/
Test(refineIntronTest,getBursetFrequencyTest41) {
	char *gsd;
	asprintf(&gsd, "Gt");
	char *gsa;
	asprintf(&gsa, "cg");
	cr_expect(getBursetFrequency(gsd,gsa)==4);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 3
*/
Test(refineIntronTest,getBursetFrequencyTest42) {
	char *gsd;
	asprintf(&gsd, "Gt");
	char *gsa;
	asprintf(&gsa, "cT");
	cr_expect(getBursetFrequency(gsd,gsa)==3);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest43) {
	char *gsd;
	asprintf(&gsd, "Gt");
	char *gsa;
	asprintf(&gsa, "GC");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 10
*/
Test(refineIntronTest,getBursetFrequencyTest44) {
	char *gsd;
	asprintf(&gsd, "Gt");
	char *gsa;
	asprintf(&gsa, "GG");
	cr_expect(getBursetFrequency(gsd,gsa)==10);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest45) {
	char *gsd;
	asprintf(&gsd, "Gt");
	char *gsa;
	asprintf(&gsa, "Gt");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 7
*/
Test(refineIntronTest,getBursetFrequencyTest46) {
	char *gsd;
	asprintf(&gsd, "Gt");
	char *gsa;
	asprintf(&gsa, "TA");
	cr_expect(getBursetFrequency(gsd,gsa)==7);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 2
*/
Test(refineIntronTest,getBursetFrequencyTest47) {
	char *gsd;
	asprintf(&gsd, "GT");
	char *gsa;
	asprintf(&gsa, "tC");
	cr_expect(getBursetFrequency(gsd,gsa)==2);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 8
*/
Test(refineIntronTest,getBursetFrequencyTest48) {
	char *gsd;
	asprintf(&gsd, "GT");
	char *gsa;
	asprintf(&gsa, "tG");
	cr_expect(getBursetFrequency(gsd,gsa)==8);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 2
*/
Test(refineIntronTest,getBursetFrequencyTest49) {
	char *gsd;
	asprintf(&gsd, "GT");
	char *gsa;
	asprintf(&gsa, "tt");
	cr_expect(getBursetFrequency(gsd,gsa)==2);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 6
*/
Test(refineIntronTest,getBursetFrequencyTest50) {
	char *gsd;
	asprintf(&gsd, "ta");
	char *gsa;
	asprintf(&gsa, "ag");
	cr_expect(getBursetFrequency(gsd,gsa)==6);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest51) {
	char *gsd;
	asprintf(&gsd, "TA");
	char *gsa;
	asprintf(&gsa, "CG");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest52) {
	char *gsd;
	asprintf(&gsd, "TA");
	char *gsa;
	asprintf(&gsa, "TC");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest53) {
	char *gsd;
	asprintf(&gsd, "TC");
	char *gsa;
	asprintf(&gsa, "AG");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest54) {
	char *gsd;
	asprintf(&gsd, "tc");
	char *gsa;
	asprintf(&gsa, "GG");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest55) {
	char *gsd;
	asprintf(&gsd, "TG");
	char *gsa;
	asprintf(&gsa, "AC");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 7
*/
Test(refineIntronTest,getBursetFrequencyTest56) {
	char *gsd;
	asprintf(&gsd, "TG");
	char *gsa;
	asprintf(&gsa, "AG");
	cr_expect(getBursetFrequency(gsd,gsa)==7);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 2
*/
Test(refineIntronTest,getBursetFrequencyTest57) {
	char *gsd;
	asprintf(&gsd, "TG");
	char *gsa;
	asprintf(&gsa, "GG");
	cr_expect(getBursetFrequency(gsd,gsa)==2);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 5
*/
Test(refineIntronTest,getBursetFrequencyTest58) {
	char *gsd;
	asprintf(&gsd, "TT");
	char *gsa;
	asprintf(&gsa, "ag");
	cr_expect(getBursetFrequency(gsd,gsa)==5);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest59) {
	char *gsd;
	asprintf(&gsd, "TT");
	char *gsa;
	asprintf(&gsa, "aT");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}

/*
	create variables needed for the function getBursetFrequency,
	set these variables with right values,
	verify that the function with these variables returns 1
*/
Test(refineIntronTest,getBursetFrequencyTest60) {
	char *gsd;
	asprintf(&gsd, "TT");
	char *gsa;
	asprintf(&gsa, "GG");
	cr_expect(getBursetFrequency(gsd,gsa)==1);
}
