//gcc bit_vector_test.c -o bit_vector_test -l criterion -I '/home/lorenzo/PIntron/include'

#include "bit_vector.h"
#include <stddef.h>
#include <signal.h>
#include "util.h"
#include <stdio.h>
#include <limits.h>
#include <string.h>

#define _ASSERT_VALID_BV( bv )						\
  my_assert(bv!=NULL);

#define _ASSERT_VALID_POS( bv, i )				\
  my_assert((i)<(bv)->n)
//  my_assert(0<=(i) && (i)<(bv)->n)

#include "../src/bit_vector.c"
#include <criterion/criterion.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int BV_comp(pbit_vect bv, pbit_vect bv2)
{
	
	my_assert(bv->n<10001);
	int flag=1;
	for (unsigned int i= 0; i<bv->n; ++i) {
		if(!(BV_get(bv,i)==BV_get(bv2,i)))
			flag=0;
	}
	return flag;
}

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
	create a vector of 10 elements called v1,
	create a vector v2 cloned from v1,
	test that v1 and v2 are equals
*/
Test(BV_test,BV_cloneTest) {
	pbit_vect v1=BV_create(10);
	pbit_vect v2=BV_clone(v1);
	cr_expect(BV_comp(v1,v2)==1);
}

/*
	create a vector of 10 elements called v1,
	create a vector v2 cloned from v1,
	set an element of v1 with a value,
	test that v1 and v2 are different
*/
Test(BV_test,BV_cloneDiffTest) {
	pbit_vect v1=BV_create(10);
	pbit_vect v2=BV_clone(v1);
	BV_set(v2,1,3);
	cr_expect(BV_comp(v1,v2)==0);
}

/*
	create a vector v1 of 1 element,
	set true value in position 0,
	verify that true value is in position 0
*/
Test(BV_test,BV_singleValueTest) {
	pbit_vect v1=BV_create(1);
	BV_set(v1,0,1);
	cr_expect(BV_get(v1,0)==1);
}

/*
	create a vector v1 of 0 element,
	verify that this will cause a failure
*/
Test(BV_test,BV_errorSizeTest,.signal = SIGSEGV) {
	pbit_vect v1=BV_create(0);
}

/*
	create a vector v1 of 3 elements,
	set false value in position 0,
	set true value in position 1,
	set false value in position 2,
	verify all that values
*/
Test(BV_test,BV_multipleValueTest) {
	pbit_vect v1=BV_create(3);
	BV_set(v1,0,0);
	BV_set(v1,1,1);
	BV_set(v1,2,0);
	cr_expect(BV_get(v1,0)==0);
	cr_expect(BV_get(v1,1)==1);
	cr_expect(BV_get(v1,2)==0);
}

/*
	create a vector v1 of 3 elements,
	set true value in position 0,
	set true value in position 1,
	verify that the vector contains only true values
*/
Test(BV_test,BV_allTrueTest) {
	pbit_vect v1=BV_create(2);
	BV_set(v1,0,1);
	BV_set(v1,1,1);
	cr_expect(BV_all_true(v1)==1);
}
