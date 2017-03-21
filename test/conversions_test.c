//gcc conversions_test.c -o conversions_test -l criterion -I '/home/lorenzo/PIntron/include'

#include "list.h"
#include "log.h"
#include "util.h"
#include <stdlib.h>

#include "../src/conversions.c"
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
	create required variables for the function get_ABS_coord,
	set the variables with expected values,
	verify that the returned value of the function is correct
*/
Test(convTest,singleValueTest) {
	int gas=3;
	int gae=5;
	int str=0;
	int crd=2;
	cr_expect(get_ABS_coord(gas,gae,str,crd)==4);
}

/*
	create required variables for the function get_ABS_coord,
	set the variables with expected values,
	verify that the returned value of the function is correct
*/
Test(convTest,singleValueTest2) {
	int gas=1;
	int gae=598;
	int str=1;
	int crd=29;
	cr_expect(get_ABS_coord(gas,gae,str,crd)==29);
}

/*
	create required variables for the function get_ABS_coord,
	set the variables with unexpected values,
	verify that the returned value of the function is correct
*/
Test(convTest,wrongValueTest) {
	int gas=5;
	int gae=3;
	int str=1;
	int crd=2;
	cr_expect(get_ABS_coord(gas,gae,str,crd)==6);
}

/*
	create required variables for the function get_ABS_coord,
	set the variables with unexpected values (negative values),
	verify that the returned value of the function is correct
*/
Test(convTest,wrongValueTest2) {
	int gas=-5;
	int gae=-3;
	int str=1;
	int crd=-2;
	cr_expect(get_ABS_coord(gas,gae,str,crd)==-8);
}

/*
	create required variables for the function get_ABS_coord,
	set the variables with unexpected values (negative values) and wrong str,
	verify that the returned value of the function is correct
*/
Test(convTest,wrongValueTest3) {
	int gas=-5;
	int gae=-3;
	int str=-7;
	int crd=-2;
	cr_expect(get_ABS_coord(gas,gae,str,crd)==0);
}
