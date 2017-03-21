//gcc refine-intron_test.c -o refine-intron_test -l criterion -I '/home/lorenzo/PIntron/include'

#include <ctype.h>

#include "refine-intron.h"
#include "refine.h"
#include "est-factorizations.h"
#include "list.h"
#include "types.h"

#include "log.h"

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
