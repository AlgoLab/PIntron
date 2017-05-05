//gcc util_test.c -o util_test -l criterion -I '/home/lorenzo/PIntron/include'

#include "util.h"
#include "log.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>

#include "../src/util.c"
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
	create a string containing testTESTtest,
	extract from s1 TESTtest with the substring function and save it in s2,
	verify that s2 contains TESTtest
*/
Test(utilTest,substringTest) {
	char *s1="testTESTtest";
	char *s2=substring(4,s1);
	cr_expect(!strcmp(s2,"TESTtest"));
}

/*
	create an empty string,
	try to extract something with the substring function from s1,
	verify that this will cause a failure
*/
Test(utilTest,substringEmptyTest,.signal = SIGSEGV) {
	char *s1="";
	char *s2=substring(7,s1);
	cr_expect(!strcmp(s2,""));
}

/*
	create a string with length of 4,
	try to extract something from index 7 with the substring function from s1,
	verify that this will cause a failure
*/
Test(utilTest,substringWrongIndexTest,.signal = SIGSEGV) {
	char *s1="test";
	char *s2=substring(7,s1);
	cr_expect(!strcmp(s2,""));
}

/*
	create a string with length of 4,
	try to extract something from index -7 with the substring function from s1,
	verify that this will cause a failure
*/
Test(utilTest,substringNegativeIndexTest,.signal = SIGSEGV) {
	char *s1="test";
	char *s2=substring(-7,s1);
	cr_expect(!strcmp(s2,""));
}

/*
	create a string containing testTESTtest,
	extract from s1 TESTtest with the real_substring function and save it in s2,
	verify that s2 contains TESTtest
*/
Test(utilTest,realSubstringTest) {
	char *s1="testTESTtest";
	char *s2=real_substring(4,12,s1);
	cr_expect(!strcmp(s2,"TESTtest"));
}
