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

/*
	create two intlists,
	add 79 to the tail of the first list,
	add 79 to the head of the second list,
	extract the head of the first list,
	extract the tail of the second list,
	test that the two extracted elements are equal
*/
Test(intListTest,headTailTest) {
	pintlist i1=intlist_create();
	pintlist i2=intlist_create();
	intlist_add_to_tail(i1,79);
	intlist_add_to_head(i2,79);
	int v1=intlist_head(i1);
	int v2=intlist_tail(i2);
	cr_expect(v1==v2);
}

/*
	create two intlists,
	add 79 to the tail of the first list,
	add 87 to the head of the second list,
	extract the head of the first list,
	extract the tail of the second list,
	test that the two extracted elements are different
*/
Test(intListTest,headTailDiffTest) {
	pintlist i1=intlist_create();
	pintlist i2=intlist_create();
	intlist_add_to_tail(i1,79);
	intlist_add_to_head(i2,87);
	int v1=intlist_head(i1);
	int v2=intlist_tail(i2);
	cr_expect(v1!=v2);
}

/*
	create a intlist,
	try to extract head and tail from it,
	verify that head and tail are null
*/
Test(intListTest,emptyTest) {
	pintlist i1=intlist_create();
	int v1=intlist_head(i1);
	int v2=intlist_tail(i1);
	cr_expect(v1==0);
	cr_expect(v2==0);
}

/*
	create a intlist,
	verify that the size of the list is 0,
	verify that the list is empty
*/
Test(intListTest,emptyListTest) {
	pintlist i1=intlist_create();
	cr_expect(intlist_size(i1)==0);
	cr_expect(intlist_is_empty(i1)==1);
}

/*
	create two intlists,
	enter some elements in the lists,
	verify that the sizes of the lists are correct
*/
Test(intListTest,listSizeTest) {
	pintlist i1=intlist_create();
	pintlist i2=intlist_create();
	intlist_add_to_tail(i1,18);
	intlist_add_to_tail(i2,20);
	intlist_add_to_head(i2,91);
	cr_expect(intlist_size(i1)==1);
	cr_expect(intlist_size(i2)==2);
}

/*
	create a intlist,
	enter some elements in the lists,
	create a copy of the list,
	verify that the sizes of the 2 lists are equal,
	verify that the heads and the tails of the two lists are equal
*/
Test(intListTest,listCopyTest) {
	pintlist i1=intlist_create();
	intlist_add_to_tail(i1,1);
	intlist_add_to_tail(i1,0);
	intlist_add_to_head(i1,1);
	pintlist i2=intlist_copy(i1);
	cr_expect(intlist_size(i1)==intlist_size(i2));
	int v1=intlist_head(i1);
	int v2=intlist_head(i2);
	cr_expect(v1==v2);
	int v3=intlist_tail(i1);
	int v4=intlist_tail(i2);
	cr_expect(v3==v4);
}

/*
	create a intlist,
	enter a value in the head of the list,
	verify that this value was saved correctly
*/
Test(intListTest,headSingleValueTest) {
	pintlist i1=intlist_create();
	intlist_add_to_head(i1,83);
	int v1=intlist_head(i1);
	cr_expect(v1==83);
}

/*
	create a intlist,
	enter a value in the tail of the list,
	verify that this value was saved correctly
*/
Test(intListTest,tailSingleValueTest) {
	pintlist i1=intlist_create();
	intlist_add_to_tail(i1,196);
	int v1=intlist_tail(i1);
	cr_expect(v1==196);
}

/*
	create two intlists,
	enter many values in these lists,
	verify that them values were saved correctly
*/
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
