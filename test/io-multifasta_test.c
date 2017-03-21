//gcc io-multifasta_test.c -o io-multifasta_test -l criterion -I '/home/lorenzo/PIntron/include'

#include "io-multifasta.h"
#include "log.h"
#include "util.h"
#define LEN_BUFFER 10000000
#define LEN_STRAND_ARRAY 10
#define LEN_ABSCOORD_ARRAY 100

#include "../src/io-multifasta.c"
#include "../src/list.c"
#include "../src/bool_list.c"
#include "../src/util.c"
#include "../src/types.c"
#include "../src/bit_vector.c"
#include "../src/ext_array.c"
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
	verify that the function get_Complement returns the input if it can not find the letter
*/
Test(ioMultifastaTest,getComplementErrorTest) {
	cr_expect(get_Complement('z')=='z');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest) {
	cr_expect(get_Complement('A')=='T');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest2) {
	cr_expect(get_Complement('a')=='t');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest3) {
	cr_expect(get_Complement('T')=='A');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest4) {
	cr_expect(get_Complement('t')=='a');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest5) {
	cr_expect(get_Complement('C')=='G');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest6) {
	cr_expect(get_Complement('c')=='g');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest7) {
	cr_expect(get_Complement('G')=='C');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest8) {
	cr_expect(get_Complement('g')=='c');
}

/*	DUPLICATE LINE OF CODE FOUND AT LINE 57 OF io-multifasta.c
	verify that the function get_Complement return the correct complement

Test(ioMultifastaTest,getComplementTestX) {
	cr_expect(get_Complement('A')=='T');
}*/

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest9) {
	cr_expect(get_Complement('R')=='Y');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest10) {
	cr_expect(get_Complement('r')=='y');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest11) {
	cr_expect(get_Complement('Y')=='R');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest12) {
	cr_expect(get_Complement('y')=='r');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest13) {
	cr_expect(get_Complement('M')=='K');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest14) {
	cr_expect(get_Complement('m')=='k');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest15) {
	cr_expect(get_Complement('K')=='M');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest16) {
	cr_expect(get_Complement('k')=='m');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest17) {
	cr_expect(get_Complement('B')=='V');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest18) {
	cr_expect(get_Complement('b')=='v');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest19) {
	cr_expect(get_Complement('V')=='B');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest20) {
	cr_expect(get_Complement('v')=='b');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest21) {
	cr_expect(get_Complement('D')=='H');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest22) {
	cr_expect(get_Complement('d')=='h');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest23) {
	cr_expect(get_Complement('H')=='D');
}

/*
	verify that the function get_Complement return the correct complement
*/
Test(ioMultifastaTest,getComplementTest24) {
	cr_expect(get_Complement('h')=='d');
}
