//gcc testPIntronOutput.c -o testPIntronOutput -l criterion

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <criterion/criterion.h>

int compare(char f1[],char f2[]){
	int flag=1;
	FILE *fp1,*fp2;
	int l;
	char r1,r2;
	if((fp1=fopen(f1,"r"))==NULL){
		flag=-1;
		//error opening file 1
	}
	if((fp2=fopen(f2,"r"))==NULL){
		flag=-1;
		//error opening file 2
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
			//files are different, n in the number of first different row
			break;
		}
	}
	if(flag==1){
		//files are equals
	}
	if((fclose(fp1)==EOF)||(fclose(fp2)==EOF)){
		flag=-2;
		//error closing files
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


Test(output,confrontaOutput) {
	char f1[]="exampleOutput/pintron-all-isoforms.gtf";
	char f2[]="executionOutput/pintron-all-isoforms.gtf";
	cr_expect(compare(f1,f2)==1);

	char f3[]="exampleOutput/pintron-full-output.json";
	char f4[]="executionOutput/pintron-full-output.json";
	cr_expect(compare(f3,f4)==1);
}
