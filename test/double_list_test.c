//gcc double_list_test.c -o double_list_test -l criterion -I '/home/lorenzo/PIntron/include'


#include "double_list.h"
#include "util.h"
#include <stdlib.h>
#include "log.h"

#include "../src/double_list.c"
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
	create two doublelists,
	add 79.82 to the tail of the first list,
	add 79.82 to the head of the second list,
	extract the head of the first list,
	extract the tail of the second list,
	test that the two extracted elements are equal
*/
Test(doubleListTest,headTailTest) {
	pdoublelist d1=doublelist_create();
	pdoublelist d2=doublelist_create();
	doublelist_add_to_tail(d1,79.82);
	doublelist_add_to_head(d2,79.82);
	double v1=doublelist_head(d1);
	double v2=doublelist_tail(d2);
	cr_expect(v1==v2);
}

/*
	create two doublelists,
	add 79.82 to the tail of the first list,
	add 14.97 to the tail of the second list,
	extract the head of the first list,
	extract the tail of the second list,
	test that the two extracted elements are different
*/
Test(doubleListTest,headTailDiffTest) {
	pdoublelist d1=doublelist_create();
	pdoublelist d2=doublelist_create();
	doublelist_add_to_tail(d1,79.82);
	doublelist_add_to_head(d2,14.97);
	double v1=doublelist_head(d1);
	double v2=doublelist_tail(d2);
	cr_expect(v1!=v2);
}

/*
	create a doublelist,
	try to extract head and tail from it,
	verify that head and tail are null
*/
Test(doubleListTest,emptyTest) {
	pdoublelist d1=doublelist_create();
	double v1=doublelist_head(d1);
	double v2=doublelist_tail(d1);
	cr_expect(v1==0);
	cr_expect(v2==0);
}

/*
	create a doublelist,
	verify that the size of the list is 0,
	verify that the list is empty
*/
Test(doubleListTest,emptyListTest) {
	pdoublelist d1=doublelist_create();
	cr_expect(doublelist_size(d1)==0);
	cr_expect(doublelist_is_empty(d1)==1);
}

/*
	create two doublelists,
	enter some elements in the lists,
	verify that the sizes of the lists are correct
*/
Test(doubleListTest,listSizeTest) {
	pdoublelist d1=doublelist_create();
	pdoublelist d2=doublelist_create();
	doublelist_add_to_tail(d1,79.14);
	doublelist_add_to_tail(d2,13.78);
	doublelist_add_to_head(d2,23.67);
	cr_expect(doublelist_size(d1)==1);
	cr_expect(doublelist_size(d2)==2);
}

/*
	create a double,
	enter some elements in the lists,
	create a copy of the list,
	verify that the sizes of the 2 lists are equal,
	verify that the heads and the tails of the two lists are equal
*/
Test(doubleListTest,listCopyTest) {
	pdoublelist d1=doublelist_create();
	doublelist_add_to_tail(d1,79.16);
	doublelist_add_to_tail(d1,0.27);
	doublelist_add_to_head(d1,19.45);
	pdoublelist d2=doublelist_copy(d1);
	cr_expect(doublelist_size(d1)==doublelist_size(d2));
	double v1=doublelist_head(d1);
	double v2=doublelist_head(d2);
	cr_expect(v1==v2);
	double v3=doublelist_tail(d1);
	double v4=doublelist_tail(d2);
	cr_expect(v3==v4);
}

/*
	create a doublelist,
	enter a value in the head of the list,
	verify that this value was saved correctly
*/
Test(doubleListTest,headSingleValueTest) {
	pdoublelist d1=doublelist_create();
	doublelist_add_to_head(d1,83.87);
	double v1=doublelist_head(d1);
	cr_expect(v1==83.87);
}

/*
	create a doublelist,
	enter a value in the tail of the list,
	verify that this value was saved correctly
*/
Test(doubleListTest,tailSingleValueTest) {
	pdoublelist d1=doublelist_create();
	doublelist_add_to_tail(d1,196.124);
	double v1=doublelist_tail(d1);
	cr_expect(v1==196.124);
}

/*
	create two doublelists,
	enter many values in these lists,
	verify that them values were saved correctly
*/
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
