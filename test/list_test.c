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
	create an empty list,
	verify that the size of the list is 0,
	verify that the list is empty
*/
Test(listTest,isEmptyTest) {
	plist l1=list_create();
	cr_expect(list_size(l1)==0);
	cr_expect(list_is_empty(l1)==1);
}

/*
	create two empty lists,
	merge them into the first one,
	verify that the size of the merged list is 0,
	verify that the merged list is empty
*/
Test(listTest,mergeEmptyTest) {
	plist l1=list_create();
	plist l2=list_create();
	list_merge(l1,l2);
	cr_expect(list_size(l1)==0);
	cr_expect(list_is_empty(l1)==1);
}

/*
	create two empty lists,
	merge them into a new one,
	verify that the size of the new list is 0,
	verify that the new list is empty
*/
Test(listTest,mergeNewEmptyTest) {
	plist l1=list_create();
	plist l2=list_create();
	plist l3=list_merge_new(l1,l2);
	cr_expect(list_size(l3)==0);
	cr_expect(list_is_empty(l3)==1);
}

/*
	create an empty list,
	try to use the head function on it,
	verify that this will return 0
*/
Test(listTest,headEmptyTest) {
	plist l1=list_create();
	cr_expect(list_head(l1)==0);
}

/*
	create an empty list,
	try to use the tail function on it,
	verify that this will return 0
*/
Test(listTest,tailEmptyTest) {
	plist l1=list_create();
	cr_expect(list_tail(l1)==0);
}

/*
	create an empty list,
	try to use the remove_from_head function on it,
	verify that this will cause a failure
*/
Test(listTest,removeHeadEmptyTest,.signal = SIGSEGV) {
	plist l1=list_create();
	cr_expect(list_remove_from_head(l1));
}

/*
	create an empty list,
	try to use the remove_from_tail function on it,
	verify that this will cause a failure
*/
Test(listTest,removeTailEmptyTest,.signal = SIGSEGV) {
	plist l1=list_create();
	cr_expect(list_remove_from_tail(l1));
}
