//gcc int_list_test.c -o int_list_test -l criterion -I '/home/lorenzo/PIntron/include'


#include "int_list.h"
#include "util.h"
#include <stdlib.h>
#include "log.h"

#include "../src/int_list.c"
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

Test(intListTest,headTailTest) {
	pintlist i1=intlist_create();
	pintlist i2=intlist_create();
	intlist_add_to_tail(i1,79);
	intlist_add_to_head(i2,79);
	int v1=intlist_head(i1);
	int v2=intlist_tail(i2);
	cr_expect(v1==v2);
}

Test(intListTest,emptyTest) {
	pintlist i1=intlist_create();
	int v1=intlist_head(i1);
	int v2=intlist_tail(i1);
	cr_expect(v1==0);
	cr_expect(v2==0);
}

Test(intListTest,headSingleValueTest) {
	pintlist i1=intlist_create();
	intlist_add_to_head(i1,83);
	int v1=intlist_head(i1);
	cr_expect(v1==83);
}

Test(intListTest,tailSingleValueTest) {
	pintlist i1=intlist_create();
	intlist_add_to_tail(i1,196);
	int v1=intlist_tail(i1);
	cr_expect(v1==196);
}

Test(intListTest,multipleValueTest) {
	pintlist i1=intlist_create();
	intlist_add_to_head(i1,85);
	intlist_add_to_head(i1,72);
	intlist_add_to_tail(i1,83);
	intlist_add_to_tail(i1,79);
	int v1=intlist_head(i1);
	int v2=intlist_tail(i1);
	cr_expect(v1==72);
	cr_expect(v2==79);
}
