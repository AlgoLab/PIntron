//gcc bool_list_test.c -o bool_list_test -l criterion -I '/home/lorenzo/PIntron/include'


#include "bool_list.h"
#include "util.h"
#include <stdlib.h>
#include "log.h"

#include "../src/bool_list.c"
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

/*
	create two boollists,
	add 1 to the tail of the first list,
	add 1 to the head of the second list,
	extract the head of the first list,
	extract the tail of the second list,
	test that the two extracted elements are equal
*/
Test(boolListTest,headTailTest) {
	pboollist b1=boollist_create();
	pboollist b2=boollist_create();
	boollist_add_to_tail(b1,1);
	boollist_add_to_head(b2,1);
	bool v1=boollist_head(b1);
	bool v2=boollist_tail(b2);
	cr_expect(v1==v2);
}

/*
	create two boollists,
	add 1 to the tail of the first list,
	add 0 to the tail of the second list,
	extract the head of the first list,
	extract the tail of the second list,
	test that the two extracted elements are different
*/
Test(boolListTest,headTailDiffTest) {
	pboollist b1=boollist_create();
	pboollist b2=boollist_create();
	boollist_add_to_tail(b1,1);
	boollist_add_to_tail(b2,0);
	bool v1=boollist_head(b1);
	bool v2=boollist_tail(b2);
	cr_expect(v1!=v2);
}

/*
	create a boollist,
	try to extract head and tail from it,
	verify that head and tail are null
*/
Test(boolListTest,emptyTest) {
	pboollist b1=boollist_create();
	bool v1=boollist_head(b1);
	bool v2=boollist_tail(b1);
	cr_expect(v1==0);
	cr_expect(v2==0);
}

/*
	create a boollist,
	verify that the size of the list is 0,
	verify that the list is empty
*/
Test(boolListTest,emptyListTest) {
	pboollist b1=boollist_create();
	cr_expect(boollist_size(b1)==0);
	cr_expect(boollist_is_empty(b1)==1);
}

/*
	create two boollists,
	enter some elements in the lists,
	verify that the sizes of the lists are correct
*/
Test(boolListTest,listSizeTest) {
	pboollist b1=boollist_create();
	pboollist b2=boollist_create();
	boollist_add_to_tail(b1,1);
	boollist_add_to_tail(b2,0);
	boollist_add_to_head(b2,1);
	cr_expect(boollist_size(b1)==1);
	cr_expect(boollist_size(b2)==2);
}

/*
	create a boollist,
	enter some elements in the lists,
	create a copy of the list,
	verify that the sizes of the 2 lists are equal,
	verify that the heads ant the tails of the two lists are equal
*/
Test(boolListTest,listCopyTest) {
	pboollist b1=boollist_create();
	boollist_add_to_tail(b1,1);
	boollist_add_to_tail(b1,0);
	boollist_add_to_head(b1,1);
	pboollist b2=boollist_copy(b1);
	cr_expect(boollist_size(b1)==boollist_size(b2));
	bool v1=boollist_head(b1);
	bool v2=boollist_head(b2);
	cr_expect(v1==v2);
	bool v3=boollist_tail(b1);
	bool v4=boollist_tail(b2);
	cr_expect(v3==v4);
}

/*
	create a boollist,
	enter a value in the head of the list,
	verify that this value was saved correctly
*/
Test(boolListTest,headSingleValueTest) {
	pboollist b1=boollist_create();
	boollist_add_to_head(b1,1);
	bool v1=boollist_head(b1);
	cr_expect(v1==1);
}

/*
	create a boollist,
	enter a value in the tail of the list,
	verify that this value was saved correctly
*/
Test(boolListTest,tailSingleValueTest) {
	pboollist b1=boollist_create();
	boollist_add_to_tail(b1,1);
	bool v1=boollist_tail(b1);
	cr_expect(v1==1);
}

/*
	create two boollists,
	enter many values in these lists,
	verify that them values were saved correctly
*/
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
