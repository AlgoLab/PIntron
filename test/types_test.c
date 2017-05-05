//gcc types_test.c -o types_test -l criterion -I '/home/lorenzo/PIntron/include'

#include "types.h"
#include "list.h"
#include "bit_vector.h"
#include "../src/my_time.c"
#include "../src/util.c"
#include "../src/list.c"
#include "../src/ext_array.c"
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
#include <signal.h>
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
#include "../src/types.c"
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
	create an empty ppointer,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,pointerTest) {
	ppointer pp1=pointer_create();
	pointer_destroy(pp1);
}

/*
	create an empty EST,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,estTest) {
	pEST pe1=EST_create();
	EST_destroy(pe1);
}

/*
	create an empty EST,
	try to destroy just factorization of it,
	verify that this will work correctly
*/
Test(typesTest,estFactTest) {
	pEST pe1=EST_create();
	EST_destroy_just_factorizations(pe1);
}

/*
	create an empty EST_info,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,estInfoTest) {
	pEST_info pei1=EST_info_create();
	EST_info_destroy(pei1);
}

/*
	create an empty factorization,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,factTest) {
	pfactorization pf1=factorization_create();
	factorization_destroy(pf1);
}

/*
	create an empty embedding plist,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,plistEmbTest) {
	plist pem1=embedding_create();
	embedding_destroy(pem1);
}

/*
	create an empty factor,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,factorTest) {
	pfactor f1=factor_create();
	factor_destroy(f1);
}

/*
	create an empty genomic intron with correct values,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,genomicIntronTest) {
	pgenomic_intron gi1=genomic_intron_create(1,3);
	genomic_intron_destroy(gi1);
}

/*
	create an empty genomic intron with wrong values,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,genomicIntronWrongTest) {
	pgenomic_intron gi1=genomic_intron_create(93,5);
	genomic_intron_destroy(gi1);
}

/*
	create an empty genomic intron with limit values,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,genomicIntronLimitTest) {
	pgenomic_intron gi1=genomic_intron_create(7,7);
	genomic_intron_destroy(gi1);
}

/*
	create an empty genomic intron with 0 values,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,genomicIntronZeroTest) {
	pgenomic_intron gi1=genomic_intron_create(0,0);
	genomic_intron_destroy(gi1);
}

/*
	create an empty genomic intron with high values,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,genomicIntronHighTest) {
	pgenomic_intron gi1=genomic_intron_create(86721589,1246549534);
	genomic_intron_destroy(gi1);
}

/*
	create an empty genomic intron with negative values,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,genomicIntronNegativeTest) {
	pgenomic_intron gi1=genomic_intron_create(-3,-7);
	genomic_intron_destroy(gi1);
}

/*
	create an empty genomic intron with one negative value,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,genomicIntronFirstNegativeTest) {
	pgenomic_intron gi1=genomic_intron_create(-3,7);
	genomic_intron_destroy(gi1);
}

/*
	create an empty genomic intron with one negative value,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,genomicIntronSecondNegativeTest) {
	pgenomic_intron gi1=genomic_intron_create(3,-7);
	genomic_intron_destroy(gi1);
}

/*
	create two empty genomic intron with correct values,
	verify that they are equal
*/
Test(typesTest,genomicIntronCompareTest) {
	pgenomic_intron gi1=genomic_intron_create(1,3);
	pgenomic_intron gi2=genomic_intron_create(1,3);
	cr_expect(genomic_intron_compare(&gi1,&gi2));
}
/*
	create two empty genomic intron with wrong values,
	verify that they are equal
*/
Test(typesTest,genomicIntronCompareWrongTest) {
	pgenomic_intron gi1=genomic_intron_create(93,5);
	pgenomic_intron gi2=genomic_intron_create(93,5);
	cr_expect(genomic_intron_compare(&gi1,&gi2));
}

/*
	create two empty genomic intron with limit values,
	verify that they are equal
*/
Test(typesTest,genomicIntronCompareLimitTest) {
	pgenomic_intron gi1=genomic_intron_create(7,7);
	pgenomic_intron gi2=genomic_intron_create(7,7);
	cr_expect(genomic_intron_compare(&gi1,&gi2));
}

/*
	create two empty genomic intron with 0 values,
	verify that they are equal
*/
Test(typesTest,genomicIntronCompareZeroTest) {
	pgenomic_intron gi1=genomic_intron_create(0,0);
	pgenomic_intron gi2=genomic_intron_create(0,0);
	cr_expect(genomic_intron_compare(&gi1,&gi2));
}

/*
	create two empty genomic intron with high values,
	verify that they are equal
*/
Test(typesTest,genomicIntronCompareHighTest) {
	pgenomic_intron gi1=genomic_intron_create(86721589,1246549534);
	pgenomic_intron gi2=genomic_intron_create(86721589,1246549534);
	cr_expect(genomic_intron_compare(&gi1,&gi2));
}

/*
	create two empty genomic intron with negative values,
	verify that they are equal
*/
Test(typesTest,genomicIntronCompareNegativeTest) {
	pgenomic_intron gi1=genomic_intron_create(-3,-7);
	pgenomic_intron gi2=genomic_intron_create(-3,-7);
	cr_expect(genomic_intron_compare(&gi1,&gi2));
}

/*
	create two empty genomic intron with partial negative values,
	verify that they are equal
*/
Test(typesTest,genomicIntronCompareFirstNegativeTest) {
	pgenomic_intron gi1=genomic_intron_create(-3,7);
	pgenomic_intron gi2=genomic_intron_create(-3,7);
	cr_expect(genomic_intron_compare(&gi1,&gi2));
}

/*
	create two empty genomic intron with partial negative values,
	verify that they are equal
*/
Test(typesTest,genomicIntronCompareSecondNegativeTest) {
	pgenomic_intron gi1=genomic_intron_create(3,-7);
	pgenomic_intron gi2=genomic_intron_create(3,-7);
	cr_expect(genomic_intron_compare(&gi1,&gi2));
}

/*
	create an empty burset frequency with correct values,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,bursetFrequencyCorrectTest) {
	pburset_frequency bf1=burset_frequency_create(1,3);
	burset_frequency_destroy(bf1);
}

/*
	create an empty burset frequency with wrong values,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,bursetFrequencyWrongTest) {
	pburset_frequency bf1=burset_frequency_create(93,5);
	burset_frequency_destroy(bf1);
}

/*
	create an empty burset frequency with limit values,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,bursetFrequencyLimitTest) {
	pburset_frequency bf1=burset_frequency_create(7,7);
	burset_frequency_destroy(bf1);
}

/*
	create an empty burset frequency with 0 values,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,bursetFrequencyZeroTest) {
	pburset_frequency bf1=burset_frequency_create(0,0);
	burset_frequency_destroy(bf1);
}

/*
	create an empty burset frequency with high values,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,bursetFrequencyHighTest) {
	pburset_frequency bf1=burset_frequency_create(86721589,1246549534);
	burset_frequency_destroy(bf1);
}

/*
	create an empty burset frequency with negative values,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,bursetFrequencyNegativeTest) {
	pburset_frequency bf1=burset_frequency_create(-3,-7);
	burset_frequency_destroy(bf1);
}

/*
	create an empty burset frequency with one negative value,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,bursetFrequencyFirstNegativeTest) {
	pburset_frequency bf1=burset_frequency_create(-3,7);
	burset_frequency_destroy(bf1);
}

/*
	create an empty burset frequency with one negative value,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,bursetFrequencySecondNegativeTest) {
	pburset_frequency bf1=burset_frequency_create(3,-7);
	burset_frequency_destroy(bf1);
}

/*
	create two empty burset frequency with correct values,
	verify that they are equal
*/
Test(typesTest,bursetFrequencyCompareTest) {
	pburset_frequency bf1=burset_frequency_create(1,3);
	pburset_frequency bf2=burset_frequency_create(1,3);
	cr_expect(burset_frequency_compare(&bf1,&bf2));
}
/*
	create two empty burset frequency with wrong values,
	verify that they are equal
*/
Test(typesTest,bursetFrequencyCompareWrongTest) {
	pburset_frequency bf1=burset_frequency_create(93,5);
	pburset_frequency bf2=burset_frequency_create(93,5);
	cr_expect(burset_frequency_compare(&bf1,&bf2));
}

/*
	create two empty burset frequency with limit values,
	verify that they are equal
*/
Test(typesTest,bursetFrequencyCompareLimitTest) {
	pburset_frequency bf1=burset_frequency_create(7,7);
	pburset_frequency bf2=burset_frequency_create(7,7);
	cr_expect(burset_frequency_compare(&bf1,&bf2));
}

/*
	create two empty burset frequency with 0 values,
	verify that they are equal
*/
Test(typesTest,bursetFrequencyCompareZeroTest) {
	pburset_frequency bf1=burset_frequency_create(0,0);
	pburset_frequency bf2=burset_frequency_create(0,0);
	cr_expect(burset_frequency_compare(&bf1,&bf2));
}

/*
	create two empty burset frequency with high values,
	verify that they are equal
*/
Test(typesTest,bursetFrequencyCompareHighTest) {
	pburset_frequency bf1=burset_frequency_create(86721589,1246549534);
	pburset_frequency bf2=burset_frequency_create(86721589,1246549534);
	cr_expect(burset_frequency_compare(&bf1,&bf2));
}

/*
	create two empty burset frequency with negative values,
	verify that they are equal
*/
Test(typesTest,bursetFrequencyCompareNegativeTest) {
	pburset_frequency bf1=burset_frequency_create(-3,-7);
	pburset_frequency bf2=burset_frequency_create(-3,-7);
	cr_expect(burset_frequency_compare(&bf1,&bf2));
}

/*
	create two empty burset frequency with partial negative values,
	verify that they are equal
*/
Test(typesTest,bursetFrequencyCompareFirstNegativeTest) {
	pburset_frequency bf1=burset_frequency_create(-3,7);
	pburset_frequency bf2=burset_frequency_create(-3,7);
	cr_expect(burset_frequency_compare(&bf1,&bf2));
}

/*
	create two empty burset frequency with partial negative values,
	verify that they are equal
*/
Test(typesTest,bursetFrequencyCompareSecondNegativeTest) {
	pburset_frequency bf1=burset_frequency_create(3,-7);
	pburset_frequency bf2=burset_frequency_create(3,-7);
	cr_expect(burset_frequency_compare(&bf1,&bf2));
}

/*
	create an empty pintron,
	try to destroy it with intron_destroy,
	verify that this will work correctly
*/
Test(typesTest,intronTest) {
	pintron pin1=intron_create();
	intron_destroy(pin1);
}

/*
	create an empty pintron,
	try to destroy it with intron_destroy_free,
	verify that this will work correctly
*/
Test(typesTest,intronFreeTest) {
	pintron pin1=intron_create();
	intron_destroy_free(pin1);
}

/*
	create an empty pintron,
	try to destroy it with intron_destroy_free_if_false,
	verify that this will work correctly
*/
Test(typesTest,intronFreeFalseTest) {
	pintron pin1=intron_create();
	intron_destroy_free_if_false(pin1);
}

/*
	create an empty subtree factorizations,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,subtreeFactorizationsTest) {
	psubtree_factorizations stf1=subtree_factorizations_create();
	subtree_factorizations_destroy(stf1);
}

/*
	create an empty subtree embeddings,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,subtreeEmbeddingsTest) {
	psubtree_embeddings ste1=subtree_embeddings_create();
	subtree_embeddings_destroy(ste1);
}

/*
	create an empty pairing,
	try to destroy it with pairing_destroy,
	verify that this will work correctly
*/
Test(typesTest,pairingDestroyTest) {
	ppairing pair1=pairing_create();
	pairing_destroy(pair1);
}

/*
	create two empty pairing,
	verify that they are equal
*/
Test(typesTest,pairingCompareTest) {
	ppairing pair1=pairing_create();
	ppairing pair2=pairing_create();
	cr_expect(!pairing_compare(&pair1,&pair2));
}

/*
	create an empty pairing,
	clone it to another pairing,
	verify that they are equal
*/
Test(typesTest,pairingCloneTest) {
	ppairing pair1=pairing_create();
	ppairing pair2=pairing_simple_copy(pair1);
	cr_expect(!pairing_compare(&pair1,&pair2));
}

/*
	create an empty pairing,
	try to destroy it with pairing_destroy2,
	verify that this will work correctly
*/
Test(typesTest,pairingDestroy2Test) {
	ppairing pair1=pairing_create();
	pairing_destroy_2(pair1);
}

/*
	create an empty GEN ESTS,
	try to destroy it and verify that this will work correctly
*/
Test(typesTest,genEstsTest) {
	pGEN_ESTS gest1=GEN_ESTS_create();
	GEN_ESTS_destroy(gest1);
}

/*
	create an empty GEN ESTS,
	try to destroy it,
	verify that this will cause a failure
*/
Test(typesTest,estMegTest,.signal = SIGSEGV) {
	pEST_MEG em1=EST_MEG_create();
	EST_MEG_destroy(em1);
}

/*
	create an empty list,
	try to destroy it with vi_destroy,
	verify that this will work correctly
*/
Test(typesTest,listDestroyTest) {
	plist l1=list_create();
	vi_destroy(l1);
}
