//gcc min_factorization_test.c -o min_factorization_test -l criterion -I '/home/lorenzo/PIntron/include'

#include "min_factorization.h"
#include "color_matrix.h"
#include "list.h"
#include <assert.h>
#include <stdint.h>

#include "../src/min_factorization.c"
#include "../src/my_time.c"
#include "../src/util.c"
#include "../src/list.c"
#include "../src/ext_array.c"
#include "../src/types.c"
#include <criterion/criterion.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bit_vector.h"

#include "util.h"
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include "../src/bool_list.c"
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
	create a vector of 2 elements,
	insert two true values in it,
	verify that the vector contains two true values
*/
Test(min_factorizationTest,count_true_test1) {
	pbit_vect v1=BV_create(2);
	BV_set(v1,0,1);
	BV_set(v1,1,1);
	cr_expect(count_true(v1)==2);
}

/*
	create a vector of 3 elements,
	insert two true values and one false value in it,
	verify that the vector contains two true values
*/
Test(min_factorizationTest,count_true_test2) {
	pbit_vect v1=BV_create(3);
	BV_set(v1,0,1);
	BV_set(v1,1,1);
	BV_set(v1,2,0);
	cr_expect(count_true(v1)==2);
}

/*
	create a vector of 4 elements,
	insert two true values and one false value in it,
	verify that the vector contains two true values
*/
Test(min_factorizationTest,count_true_test3) {
	pbit_vect v1=BV_create(4);
	BV_set(v1,0,1);
	BV_set(v1,1,1);
	BV_set(v1,2,0);
	cr_expect(count_true(v1)==2);
}

/*
	create a vector of 10 elements,
	verify that the vector contains zero true values
*/
Test(min_factorizationTest,count_true_test4) {
	pbit_vect v1=BV_create(10);
	cr_expect(count_true(v1)==0);
}

/*
	create a vector of 10 elements,
	insert one true value in the last position,
	verify that the vector contains one true values
*/
Test(min_factorizationTest,count_true_test5) {
	pbit_vect v1=BV_create(10);
	BV_set(v1,9,1);
	cr_expect(count_true(v1)==1);
}

/*
	create a vector of 10 elements,
	insert two true values in arbitrary positions,
	verify that the vector contains two true values
*/
Test(min_factorizationTest,count_true_test6) {
	pbit_vect v1=BV_create(10);
	BV_set(v1,7,1);
	BV_set(v1,4,1);
	cr_expect(count_true(v1)==2);
}
