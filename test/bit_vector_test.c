//gcc bit_vector_test.c -o bit_vector_test -l criterion -I '/home/lorenzo/PIntron/include'

#include "bit_vector.h"

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
	//eseguo i test
	struct criterion_test_set *tests = criterion_initialize();
	int result = 0;
	if (criterion_handle_args(argc, argv, true))
	result = !criterion_run_all_tests(tests);
	criterion_finalize(tests);
	return result;
}

/*
	creo un vettore v1 di 10 elementi,
	creo un vettore v2 clonato dal v1
	testo che v1 e v2 siano uguali,
	in questo modo ho testato le funzioni
	BV_create, BV_clone e BV_comp
*/
Test(BV_test,BV_cloneTest) {
	pbit_vect v1=BV_create(10);
	pbit_vect v2=BV_clone(v1);
	//BV_set(v2,1,3);
	cr_expect(BV_comp(v1,v2)==1);
}

/*
	creo un vettore v1 di 1 elementi,
	setto in posizione 0 del vettore il valore true,
	testo che il valore sia stato salvato correttamente,
	in questo modo ho testato le funzioni BV_set e BV_get
*/
Test(BV_test,BV_singleValueTest) {
	pbit_vect v1=BV_create(1);
	BV_set(v1,0,1);
	cr_expect(BV_get(v1,0)==1);
}


Test(BV_test,BV_multipleValueTest) {
	pbit_vect v1=BV_create(3);
	BV_set(v1,0,0);
	BV_set(v1,1,1);
	BV_set(v1,2,0);
	cr_expect(BV_get(v1,0)==0);
	cr_expect(BV_get(v1,1)==1);
	cr_expect(BV_get(v1,2)==0);
}

Test(BV_test,BV_allTrueTest) {
	pbit_vect v1=BV_create(2);
	BV_set(v1,0,1);
	BV_set(v1,1,1);
	cr_expect(BV_all_true(v1)==1);
}
