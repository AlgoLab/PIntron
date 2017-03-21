//gcc BuildTranscripts_test.c -o BuildTranscripts_test -l criterion -I '/home/lorenzo/PIntron/include' -e mainT

#include <string.h>
#include "my_time.h"
#include "log.h"
#include "util.h"
#include "log-build-info.h"
#include <stdlib.h>

#include "../src/BuildTranscripts.c"
#include "../src/my_time.c"
#include "../src/util.c"
#include <criterion/criterion.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int mainT(int argc, char *argv[]){
	//for test execution
	struct criterion_test_set *tests = criterion_initialize();
	int result = 0;
	if (criterion_handle_args(argc, argv, true))
	result = !criterion_run_all_tests(tests);
	criterion_finalize(tests);
	return result;
}

/*
	create a simplelist,
	verify that the size of the list is 0,
	verify that the list is empty
*/
Test(buildTranscriptsTest,isEmptyTest) {
	structSimpleList *sl1=initializeSimpleList();
	cr_expect(lengthSimpleList(sl1)==0);
	cr_expect(isEmptySimpleList(sl1));
}
