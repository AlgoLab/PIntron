//gcc color_matrix_test.c -o color_matrix_test -l criterion -I '/home/lorenzo/PIntron/include'

#include "color_matrix.h"
#include "list.h"
#include "types.h"
#include "bit_vector.h"
#include "log.h"

#include "color_matrix.c"
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

Test(colorMatrixTest,test1) {
	
	cr_expect(v1==v2);
}
