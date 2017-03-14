//gcc testPIntronOutput.c -o testPIntronOutput -l criterion
//./testPIntronOutput
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <criterion/criterion.h>

int confronta(char f1[],char f2[]){
	int flag=1;
	FILE *fp1,*fp2;
	int l;
	char r1,r2;
	if((fp1=fopen(f1,"r"))==NULL){
		flag=-1;
		//errore apertura file 1
	}
	if((fp2=fopen(f2,"r"))==NULL){
		flag=-1;
		//errore apertura file 2
	}
	l=1;
	while(!feof(fp1)){
		r1=fgetc(fp1); 
		r2=fgetc(fp2);      
		if(r1=='\n'){
			++l;
		}
		if(r1!=r2){                
			flag=0;  
			//file diversi, n e' il numero della prima riga diversa
			break;
		}
	}
	if(flag==1){
		//file uguali
	}
	if((fclose(fp1)==EOF)||(fclose(fp2)==EOF)){
		flag=-2;
		//errore chiusura
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

Test(input,confrontaInput) {
	char f1[]="exampleOutput/intermediate_files/ests.txt";
	char f2[]="executionOutput/intermediate_files/ests.txt";
	cr_expect(confronta(f1,f2)==1);

	char f3[]="exampleOutput/intermediate_files/genomic.txt";
	char f4[]="executionOutput/intermediate_files/genomic.txt";
	cr_expect(confronta(f3,f4)==1);
}

Test(output,confrontaOutput) {
	char f1[]="exampleOutput/pintron-all-isoforms.gtf";
	char f2[]="executionOutput/pintron-all-isoforms.gtf";
	cr_expect(confronta(f1,f2)==1);

	char f3[]="exampleOutput/pintron-full-output.json";
	char f4[]="executionOutput/pintron-full-output.json";
	cr_expect(confronta(f3,f4)==1);
}
