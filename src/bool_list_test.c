//gcc bool_list_test.c -o bool_list_test -l criterion -I '/home/lorenzo/PIntron/include'


#include "bool_list.h"
#include "util.h"
#include <stdlib.h>
#include "log.h"

#include "bool_list.c"
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

Test(boolListTest,headTailTest) {
	pboollist b1=boollist_create();
	pboollist b2=boollist_create();
	boollist_add_to_tail(b1,1);
	boollist_add_to_head(b2,1);
	//boollist_add_to_head(b2,0);
	bool v1=boollist_head(b1);
	bool v2=boollist_tail(b2);
	cr_expect(v1==v2);
}

Test(boolListTest,emptyTest) {
	pboollist b1=boollist_create();
	bool v1=boollist_head(b1);
	bool v2=boollist_tail(b1);
	cr_expect(v1==0);
	cr_expect(v2==0);
}

Test(boolListTest,headSingleValueTest) {
	pboollist b1=boollist_create();
	boollist_add_to_head(b1,1);
	bool v1=boollist_head(b1);
	cr_expect(v1==1);
}

Test(boolListTest,tailSingleValueTest) {
	pboollist b1=boollist_create();
	boollist_add_to_tail(b1,1);
	bool v1=boollist_tail(b1);
	cr_expect(v1==1);
}

Test(boolListTest,multipleValueTest) {
	pboollist b1=boollist_create();
	boollist_add_to_head(b1,0);
	boollist_add_to_head(b1,1);
	boollist_add_to_tail(b1,0);
	boollist_add_to_tail(b1,1);
	bool v1=boollist_head(b1);
	bool v2=boollist_tail(b1);
	cr_expect(v1==1);
	cr_expect(v2==1);
}
