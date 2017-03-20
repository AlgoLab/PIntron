//gcc aug_suffix_tree_test.c -o aug_suffix_tree_test -l criterion -I '/home/lorenzo/PIntron/include' -I '/home/lorenzo/PIntron/stree_src'

#include "aug_suffix_tree.h"

#include <stdio.h>
#include "ext_array.h"
#include "util.h"
#include "log.h"

#include "../src/aug_suffix_tree.c"
#include "../src/ext_array.c"
#include "../src/bool_list.c"
#include "../src/util.c"
#include <criterion/criterion.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

size_t
lst_edge_get_length(const LST_Edge * const edge)
{
  if (!edge)
	 return 0;

  return *(edge->range.end_index) - edge->range.start_index + 1;
}

int PGen_destroy_test(ppreproc_gen pg) {
  if (pg == NULL)
	 return;
  if (pg->gen != NULL)
	 pfree(pg->gen);
  if (pg->alph != NULL)
	 pfree(pg->alph);
  if (pg->alph_occ != NULL)
	 pfree(pg->alph_occ);
  if (pg->keys != NULL)
	 pfree(pg->keys);
  pfree(pg);
  return 1;
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
	create a suffix tree,
	verify that can be destroyed
*/
Test(augSuffixTreeTest,createDestroy) {
	ppreproc_gen s1=PGen_create();
	cr_expect(PGen_destroy_test(s1)==1);
}
