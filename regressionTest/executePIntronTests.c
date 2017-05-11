//gcc executePIntronTests.c -o executePIntronTests -l criterion

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <criterion/criterion.h>
#include <criterion/parameterized.h>
#include <dirent.h>
#include <malloc.h>

char* strdup(const char* org){
	if(org==NULL)
		return NULL;
	char* newstr=malloc(strlen(org)+1);
	char* p;
	if(newstr==NULL)
		return NULL;
	p=newstr;
	while(*org)
		*p++=*org++; //copy the string 
	return newstr;
}



#define MAX_ROWS (256)
#define MAX_COLUMNS (256)
static char words1[MAX_ROWS][MAX_COLUMNS]={{'\0','\0'}};
static char words2[MAX_ROWS][MAX_COLUMNS]={{'\0','\0'}};

void bubbleSortWordsArray1(int wordCount)
{
	int c; //outer index through rows
	int d; //inner index through rows
	char swap[MAX_COLUMNS]={'\0'};

	for(c=0;c<(wordCount-1);c++)
	{
		for(d=0;d<(wordCount-c-1);d++)
		{
			if(0>strcmp(words1[d],words1[d+1]))
				{
					//words need to be swapped
					strcpy(swap,words1[d]);
					strcpy(words1[d],words1[d+1]);
					strcpy(words1[d+1],swap);
				} 
		}
	}
}

void bubbleSortWordsArray2(int wordCount)
{
	int c; //outer index through rows
	int d; //inner index through rows
	char swap[MAX_COLUMNS]={'\0'};

	for(c=0;c<(wordCount-1);c++)
	{
		for(d=0;d<(wordCount-c-1);d++)
		{
			if(0>strcmp(words2[d],words2[d+1]))
				{
					//words need to be swapped
					strcpy(swap,words2[d]);
					strcpy(words2[d],words2[d+1]);
					strcpy(words2[d+1],swap);
				}
			}
	}
}

char * my_strcat(char *dest, char *src)
{
	char *new_string = malloc(strlen(dest) + strlen(src) + 1);
	strcpy(new_string, dest);
	strcat(new_string, src);
	return(new_string);
}

int main(int argc, char *argv[]){
	system("cd regressionTest;for i in $(ls -d */); do echo ${i%%/}; done > output.txt");
	system("make dist");

	FILE * fp;
	FILE * fpgg;
	char * line = NULL;
	char * linegg = NULL;
	size_t len = 0;
	size_t lengg = 0;
	ssize_t read;
	ssize_t readgg;
	
	char *PI1="bin/pintron --bin-dir=bin/ --genomic='regressionTest/";
	//line
	char *PI2="/genomic.txt' --EST='regressionTest/";
	//line
	char *PI3="/ests.txt' --organism=human --gene=";
	//gene.txt content extracted from test* subfolder for each test
	//regressionTest/+line+/gene.txt
	char *PI4=" --output='regressionTest/";
	//line
	char *PI5="/executionOutput/full.json' --logfile='regressionTest/";
	//line
	char *PI6="/executionOutput/pintron-pipeline-log.txt' general-logfile='regressionTest/";
	//line
	char *PI7="/executionOutput/pintron-log.txt' --gtf='regressionTest/";
	//line
	char *PI8="/executionOutput/pintron-all-isoforms.gtf'";
	
	char *PI=NULL;
	
	char *FG1="regressionTest/";
	//line
	char *FG2="/gene.txt";
	//complete string:
	char *FG=NULL;

	char *MKd1="cd regressionTest/";
	//line
	char *MKd2="/;mkdir -p executionOutput";
	//complete string for make missing dirs
	char *MKdF=NULL;

	fp = fopen("regressionTest/output.txt", "r");
	while ((read = getline(&line,&len,fp)) != -1) {
		line[strlen(line) - 1] = 0;

		FG=my_strcat(FG1,line);
		FG=my_strcat(FG,FG2);

		fpgg = fopen(FG, "r");
		while ((readgg = getline(&linegg,&lengg,fpgg)) != -1) {
			linegg[strlen(linegg) - 1] = 0;}
		fclose(fpgg);

		MKdF=my_strcat(MKd1,line);
		MKdF=my_strcat(MKdF,MKd2);
		system(MKdF);

		PI=my_strcat(PI1,line);
		PI=my_strcat(PI,PI2);
		PI=my_strcat(PI,line);
		PI=my_strcat(PI,PI3);
		PI=my_strcat(PI,linegg);
		PI=my_strcat(PI,PI4);
		PI=my_strcat(PI,line);
		PI=my_strcat(PI,PI5);
		PI=my_strcat(PI,line);
		PI=my_strcat(PI,PI6);
		PI=my_strcat(PI,line);
		PI=my_strcat(PI,PI7);
		PI=my_strcat(PI,line);
		PI=my_strcat(PI,PI8);
		
		system(PI);
	}
	fclose(fp);
}


