//gcc ext_array_test.c -o ext_array_test -l criterion -I '/home/lorenzo/PIntron/include'

#include "ext_array.h"
#include "util.h"
#include <stdlib.h>
#include <string.h>

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
	create an empty ext_array,
	verify that the size of the ext_array is 0,
	verify that the ext_array is empty
*/
Test(extArrayTest,isEmptyTest) {
	pext_array v1=EA_create();
	cr_expect(EA_size(v1)==0);
	cr_expect(EA_is_empty(v1)==1);
}

/*
	create an empty ext_array,
	try to use the get function on it,
	verify that this will cause a failure
*/
Test(extArrayTest,getEmptyTest,.signal = SIGSEGV) {
	pext_array v1=EA_create();
	cr_expect(EA_get(v1,0));
}
