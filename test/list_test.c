//gcc list_test.c -o list_test -l criterion -I '/home/lorenzo/PIntron/include'

#include "list.h"
#include "log.h"
#include "util.h"
#include <stdlib.h>

#include "../src/list.c"
#include <criterion/criterion.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]){
	//for test execution
	struct criterion_test_set *tests = criterion_initialize();
	int result = 0;
	if (criterion_handle_args(argc, argv, true))
	result = !criterion_run_all_tests(tests);
	criterion_finalize(tests);
	return result;
}

Test(listTest,isEmptyTest) {
	plist l1=list_create();
	cr_expect(list_is_empty(l1)==1);
}

Test(listTest,listSizeTest) {
	plist l1=list_create();
	cr_expect(list_size(l1)==0);
}
