//gcc double_list_test.c -o double_list_test -l criterion -I '/home/lorenzo/PIntron/include'


#include "double_list.h"
#include "util.h"
#include <stdlib.h>
#include "log.h"

#include "double_list.c"
#include <criterion/criterion.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char *argv[]){
	//eseguo i test
	struct criterion_test_set *tests = criterion_initialize();
	int result = 0;
	if (criterion_handle_args(argc, argv, true))
	result = !criterion_run_all_tests(tests);
	criterion_finalize(tests);
	return result;
}

Test(doubleListTest,headTailTest) {
	pdoublelist d1=doublelist_create();
	pdoublelist d2=doublelist_create();
	doublelist_add_to_tail(d1,79.82);
	doublelist_add_to_head(d2,79.82);
	double v1=doublelist_head(d1);
	double v2=doublelist_tail(d2);
	cr_expect(v1==v2);
}

Test(doubleListTest,emptyTest) {
	pdoublelist d1=doublelist_create();
	double v1=doublelist_head(d1);
	double v2=doublelist_tail(d1);
	cr_expect(v1==0);
	cr_expect(v2==0);
}

Test(doubleListTest,headSingleValueTest) {
	pdoublelist d1=doublelist_create();
	doublelist_add_to_head(d1,83.87);
	double v1=doublelist_head(d1);
	cr_expect(v1==83.87);
}

Test(doubleListTest,tailSingleValueTest) {
	pdoublelist d1=doublelist_create();
	doublelist_add_to_tail(d1,196.124);
	double v1=doublelist_tail(d1);
	cr_expect(v1==196.124);
}

Test(doubleListTest,multipleValueTest) {
	pdoublelist d1=doublelist_create();
	doublelist_add_to_head(d1,85.21);
	doublelist_add_to_head(d1,72.32);
	doublelist_add_to_tail(d1,83.64);
	doublelist_add_to_tail(d1,79.89);
	double v1=doublelist_head(d1);
	double v2=doublelist_tail(d1);
	cr_expect(v1==72.32);
	cr_expect(v2==79.89);
}
