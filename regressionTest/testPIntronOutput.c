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
			//files are different, n is the number of first different row
			break;
		}
	}
	if(flag==1){
		//files are equal
	}
	if((fclose(fp1)==EOF)||(fclose(fp2)==EOF)){
		flag=-2;
		//error closing files
	}
	return flag;
}

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

char* returnInfoNextToPattern(char *f1,char *str,int num){
	FILE *fp1;
	int l=1;
	int find=0;
	char temp[10000];
	char *pch;
	if((fp1=fopen(f1,"r"))==NULL){
		perror("fopen failed");
		exit(EXIT_FAILURE);
		return("E-1");
	}

	while(fgets(temp,10000,fp1)!=NULL){
		if((strstr(temp,str))!=NULL){
			find++;
			pch=strtok(temp,":");
			pch=strtok(NULL,", ");

			if((num==1)&&(find==1))
				return(strdup(pch));
			if((num==2)&&(find==2))
				return(strdup(pch));
			if((num==3)&&(find==3))
				return(strdup(pch));
			if((num==4)&&(find==4))
				return(strdup(pch));
			if((num==5)&&(find==5))
				return(strdup(pch));
			if((num==6)&&(find==6))
				return(strdup(pch));
			if((num==7)&&(find==7))
				return(strdup(pch));
			if((num==8)&&(find==8))
				return(strdup(pch));
			if((num==9)&&(find==9))
				return(strdup(pch));
			if((num==10)&&(find==10))
				return(strdup(pch));
			if((num==11)&&(find==11))
				return(strdup(pch));
			if((num==12)&&(find==12))
				return(strdup(pch));
			if((num==13)&&(find==13))
				return(strdup(pch));
			if((num==14)&&(find==14))
				return(strdup(pch));
			if((num==15)&&(find==15))
				return(strdup(pch));
		}
		++l;
	}
	if(find==0){
		//couldn't find pattern
	}
	if(fclose(fp1)==EOF){
		return("E-2");
		//error closing files
	}
   	return("E1");
}

int compareJson(char f1[],char f2[]){
	int c1=!strcmp(returnInfoNextToPattern(f1,"sequence_id",1),returnInfoNextToPattern(f2,"sequence_id",1));
	int c2=!strcmp(returnInfoNextToPattern(f1,"strand",1),returnInfoNextToPattern(f2,"strand",1));
	int c6=!strcmp(returnInfoNextToPattern(f1,"acceptor_alignment_error",1),returnInfoNextToPattern(f2,"acceptor_alignment_error",1));
	int c14=!strcmp(returnInfoNextToPattern(f1,"pattern",1),returnInfoNextToPattern(f2,"pattern",1));
	int c15=!strcmp(returnInfoNextToPattern(f1,"prefix",1),returnInfoNextToPattern(f2,"prefix",1));
	int c16=!strcmp(returnInfoNextToPattern(f1,"relative_end",1),returnInfoNextToPattern(f2,"relative_end",1));
	int c17=!strcmp(returnInfoNextToPattern(f1,"relative_start",1),returnInfoNextToPattern(f2,"relative_start",1));
	int c18=!strcmp(returnInfoNextToPattern(f1,"repeat_sequence",1),returnInfoNextToPattern(f2,"repeat_sequence",1));
	int c19=!strcmp(returnInfoNextToPattern(f1,"suffix",1),returnInfoNextToPattern(f2,"suffix",1));
	int c20=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_end",1),returnInfoNextToPattern(f2,"acceptor_factor_end",1));
	int c21=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_prefix",1),returnInfoNextToPattern(f2,"acceptor_factor_prefix",1));
	int c22=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_start",1),returnInfoNextToPattern(f2,"acceptor_factor_start",1));
	int c23=!strcmp(returnInfoNextToPattern(f1,"donor_factor_start",1),returnInfoNextToPattern(f2,"donor_factor_start",1));
	int c24=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_end",2),returnInfoNextToPattern(f2,"acceptor_factor_end",2));
	int c25=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_prefix",2),returnInfoNextToPattern(f2,"acceptor_factor_prefix",2));
	int c26=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_start",2),returnInfoNextToPattern(f2,"acceptor_factor_start",2));
	int c27=!strcmp(returnInfoNextToPattern(f1,"donor_factor_start",2),returnInfoNextToPattern(f2,"donor_factor_start",2));
	int c28=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_end",3),returnInfoNextToPattern(f2,"acceptor_factor_end",3));
	int c41=!strcmp(returnInfoNextToPattern(f1,"donor_score",2),returnInfoNextToPattern(f2,"donor_score",2));
	int c42=!strcmp(returnInfoNextToPattern(f1,"length",2),returnInfoNextToPattern(f2,"length",3));//+1 in index of f2
	int c43=!strcmp(returnInfoNextToPattern(f1,"number_of_supporting_transcripts",2),returnInfoNextToPattern(f2,"number_of_supporting_transcripts",2));
	int c44=!strcmp(returnInfoNextToPattern(f1,"pattern",2),returnInfoNextToPattern(f2,"pattern",2));
	int c45=!strcmp(returnInfoNextToPattern(f1,"prefix",2),returnInfoNextToPattern(f2,"prefix",2));
	int c46=!strcmp(returnInfoNextToPattern(f1,"relative_end",2),returnInfoNextToPattern(f2,"relative_end",2));
	int c29=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_prefix",3),returnInfoNextToPattern(f2,"acceptor_factor_prefix",3));
	int c30=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_start",3),returnInfoNextToPattern(f2,"acceptor_factor_start",3));
	int c31=!strcmp(returnInfoNextToPattern(f1,"donor_factor_start",3),returnInfoNextToPattern(f2,"donor_factor_start",3));
	int c32=!strcmp(returnInfoNextToPattern(f1,"BPS_position",1),returnInfoNextToPattern(f2,"BPS_position",1));//-1 in index of f1 and f2
	int c33=!strcmp(returnInfoNextToPattern(f1,"BPS_score",2),returnInfoNextToPattern(f2,"BPS_score",2));
	int c36=!strcmp(returnInfoNextToPattern(f1,"acceptor_alignment_error",2),returnInfoNextToPattern(f2,"acceptor_alignment_error",2));
	int c37=!strcmp(returnInfoNextToPattern(f1,"acceptor_exon_prefix",2),returnInfoNextToPattern(f2,"acceptor_exon_prefix",2));
	int c38=!strcmp(returnInfoNextToPattern(f1,"acceptor_score",2),returnInfoNextToPattern(f2,"acceptor_score",2));
	int c39=!strcmp(returnInfoNextToPattern(f1,"donor_alignment_error",2),returnInfoNextToPattern(f2,"donor_alignment_error",2));
	int c40=!strcmp(returnInfoNextToPattern(f1,"donor_exon_suffix",2),returnInfoNextToPattern(f2,"donor_exon_suffix",2));	
	int c47=!strcmp(returnInfoNextToPattern(f1,"relative_start",2),returnInfoNextToPattern(f2,"relative_start",2));
	int c48=!strcmp(returnInfoNextToPattern(f1,"repeat_sequence",2),returnInfoNextToPattern(f2,"repeat_sequence",2));
	int c49=!strcmp(returnInfoNextToPattern(f1,"suffix",2),returnInfoNextToPattern(f2,"suffix",2));
	int c50=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_end",4),returnInfoNextToPattern(f2,"acceptor_factor_end",4));
	int c51=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_prefix",4),returnInfoNextToPattern(f2,"acceptor_factor_prefix",4));
	int c52=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_start",4),returnInfoNextToPattern(f2,"acceptor_factor_start",4));
	int c53=!strcmp(returnInfoNextToPattern(f1,"donor_factor_start",4),returnInfoNextToPattern(f2,"donor_factor_start",4));
	int c55=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_prefix",5),returnInfoNextToPattern(f2,"acceptor_factor_prefix",5));
	int c56=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_start",5),returnInfoNextToPattern(f2,"acceptor_factor_start",5));
	int c57=!strcmp(returnInfoNextToPattern(f1,"donor_factor_start",5),returnInfoNextToPattern(f2,"donor_factor_start",5));
	int c58=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_end",6),returnInfoNextToPattern(f2,"acceptor_factor_end",6));
	int c59=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_prefix",6),returnInfoNextToPattern(f2,"acceptor_factor_prefix",6));
	int c60=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_start",6),returnInfoNextToPattern(f2,"acceptor_factor_start",6));
	int c61=!strcmp(returnInfoNextToPattern(f1,"donor_factor_start",6),returnInfoNextToPattern(f2,"donor_factor_start",6));
	int c62=!strcmp(returnInfoNextToPattern(f1,"BPS_position",2),returnInfoNextToPattern(f2,"BPS_position",2));//-1 in index of f1 and f2
	int c63=!strcmp(returnInfoNextToPattern(f1,"BPS_score",3),returnInfoNextToPattern(f2,"BPS_score",3));
	int c66=!strcmp(returnInfoNextToPattern(f1,"acceptor_alignment_error",3),returnInfoNextToPattern(f2,"acceptor_alignment_error",3));
	int c67=!strcmp(returnInfoNextToPattern(f1,"acceptor_exon_prefix",3),returnInfoNextToPattern(f2,"acceptor_exon_prefix",3));
	int c68=!strcmp(returnInfoNextToPattern(f1,"acceptor_score",3),returnInfoNextToPattern(f2,"acceptor_score",3));
	int c69=!strcmp(returnInfoNextToPattern(f1,"donor_alignment_error",3),returnInfoNextToPattern(f2,"donor_alignment_error",3));
	int c70=!strcmp(returnInfoNextToPattern(f1,"donor_exon_suffix",3),returnInfoNextToPattern(f2,"donor_exon_suffix",3));
	int c71=!strcmp(returnInfoNextToPattern(f1,"donor_score",3),returnInfoNextToPattern(f2,"donor_score",3));
	int c72=!strcmp(returnInfoNextToPattern(f1,"length",3),returnInfoNextToPattern(f2,"length",4));//+1 in index of f2
	int c73=!strcmp(returnInfoNextToPattern(f1,"number_of_supporting_transcripts",3),returnInfoNextToPattern(f2,"number_of_supporting_transcripts",3));
	int c74=!strcmp(returnInfoNextToPattern(f1,"pattern",3),returnInfoNextToPattern(f2,"pattern",3));
	int c75=!strcmp(returnInfoNextToPattern(f1,"prefix",3),returnInfoNextToPattern(f2,"prefix",3));
	int c76=!strcmp(returnInfoNextToPattern(f1,"relative_end",3),returnInfoNextToPattern(f2,"relative_end",3));
	int c77=!strcmp(returnInfoNextToPattern(f1,"relative_start",3),returnInfoNextToPattern(f2,"relative_start",3));
	int c78=!strcmp(returnInfoNextToPattern(f1,"repeat_sequence",3),returnInfoNextToPattern(f2,"repeat_sequence",3));
	int c80=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_prefix",7),returnInfoNextToPattern(f2,"acceptor_factor_prefix",7));
	int c81=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_start",7),returnInfoNextToPattern(f2,"acceptor_factor_start",7));
	int c82=!strcmp(returnInfoNextToPattern(f1,"donor_factor_start",7),returnInfoNextToPattern(f2,"donor_factor_start",7));
	int c83=!strcmp(returnInfoNextToPattern(f1,"BPS_position",3),returnInfoNextToPattern(f2,"BPS_position",3));//-1 in index of f1 and f2
	int c84=!strcmp(returnInfoNextToPattern(f1,"BPS_score",4),returnInfoNextToPattern(f2,"BPS_score",4));
	int c87=!strcmp(returnInfoNextToPattern(f1,"acceptor_alignment_error",4),returnInfoNextToPattern(f2,"acceptor_alignment_error",4));
	int c88=!strcmp(returnInfoNextToPattern(f1,"acceptor_exon_prefix",4),returnInfoNextToPattern(f2,"acceptor_exon_prefix",4));
	int c89=!strcmp(returnInfoNextToPattern(f1,"acceptor_score",4),returnInfoNextToPattern(f2,"acceptor_score",4));
	int c90=!strcmp(returnInfoNextToPattern(f1,"donor_alignment_error",4),returnInfoNextToPattern(f2,"donor_alignment_error",4));
	int c91=!strcmp(returnInfoNextToPattern(f1,"donor_exon_suffix",4),returnInfoNextToPattern(f2,"donor_exon_suffix",4));
	int c92=!strcmp(returnInfoNextToPattern(f1,"donor_score",4),returnInfoNextToPattern(f2,"donor_score",4));
	int c93=!strcmp(returnInfoNextToPattern(f1,"length",4),returnInfoNextToPattern(f2,"length",5));//+1 in index of f2
	int c94=!strcmp(returnInfoNextToPattern(f1,"number_of_supporting_transcripts",4),returnInfoNextToPattern(f2,"number_of_supporting_transcripts",4));
	int c95=!strcmp(returnInfoNextToPattern(f1,"pattern",4),returnInfoNextToPattern(f2,"pattern",4));
	int c96=!strcmp(returnInfoNextToPattern(f1,"prefix",4),returnInfoNextToPattern(f2,"prefix",4));
	int c97=!strcmp(returnInfoNextToPattern(f1,"relative_end",4),returnInfoNextToPattern(f2,"relative_end",4));
	int c98=!strcmp(returnInfoNextToPattern(f1,"relative_start",4),returnInfoNextToPattern(f2,"relative_start",4));
	int c99=!strcmp(returnInfoNextToPattern(f1,"repeat_sequence",4),returnInfoNextToPattern(f2,"repeat_sequence",4));
	int c110=!strcmp(returnInfoNextToPattern(f1,"annotated CDS?",1),returnInfoNextToPattern(f2,"annotated_CDS?",1));
	int c111=!strcmp(returnInfoNextToPattern(f1,"3utr length",1),returnInfoNextToPattern(f2,"3UTR_length",1));
	int c117=!strcmp(returnInfoNextToPattern(f1,"cumulative genome length",1),returnInfoNextToPattern(f2,"cumulative_length",1));
	int c118=!strcmp(returnInfoNextToPattern(f1,"cumulative transcript length",1),returnInfoNextToPattern(f2,"cumulative_length_on_transcript",1));
	int c119=!strcmp(returnInfoNextToPattern(f1,"transcript length",1),returnInfoNextToPattern(f2,"length_on_transcript",1));
	int c126=!strcmp(returnInfoNextToPattern(f1,"annotated CDS?",2),returnInfoNextToPattern(f2,"annotated_CDS?",2));
	int c127=!strcmp(returnInfoNextToPattern(f1,"3utr length",2),returnInfoNextToPattern(f2,"3UTR_length",2));
	int c130=!strcmp(returnInfoNextToPattern(f1,"cumulative transcript length",2),returnInfoNextToPattern(f2,"cumulative_length_on_transcript",2));	
	int c101=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_prefix",8),returnInfoNextToPattern(f2,"acceptor_factor_prefix",8));
	int c102=!strcmp(returnInfoNextToPattern(f1,"acceptor_factor_start",8),returnInfoNextToPattern(f2,"acceptor_factor_start",8));
	int c103=!strcmp(returnInfoNextToPattern(f1,"donor_factor_start",8),returnInfoNextToPattern(f2,"donor_factor_start",8));

	if(c1&&c2&&c6&&c14&&c15&&c16&&c17&&c18&&c19&&c20&&c21&&c22&&c23&&c24&&c25&&c26&&c27&&c28&&c52&&c53
		&&c32&&c33&&c36&&c37&&c38&&c39&&c40&&c41&&c42&&c43&&c44&&c45&&c46&&c47&&c48&&c49&&c50&&c51&&c55
		&&c62&&c63&&c66&&c67&&c68&&c69&&c70&&c71&&c72&&c73&&c74&&c75&&c76&&c77&&c78&&c80&&c81&&c82&&c61
		&&c83&&c84&&c87&&c88&&c89&&c90&&c91&&c92&&c93&&c94&&c95&&c96&&c97&&c98&&c99&&c101&&c102&&c103
		&&c110&&c111&&c117&&c118&&c119&&c126&&c127&&c29&&c30&&c29&&c130&&c57&&c58&&c56&&c31&&c60&&c59
		)
		return 1;
	else
		return 0;
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

int compareGtfCr(char f1[],char f2[]){
	FILE *fp1 = NULL;
	FILE *fp2 = NULL;

	if((fp1=fopen(f1,"r"))==NULL){
		perror("fopen failed");
		exit(EXIT_FAILURE);
	}
	//implied else, fopen successful
	//read each line from file into entry in words array
	int i1=0;
	while(fgets(words1[i1],MAX_COLUMNS,fp1))
	{
		//remove trailing newline from string
		words1[i1][strlen(words1[i1])-1]='\0';
		i1++;
	}
	//'i1' contains number of valid entries in words1[][]
	//sort the array of strings
	bubbleSortWordsArray1(i1);
	
	if((fp2=fopen(f2,"r"))==NULL)
	{
		perror("fopen failed");
		exit(EXIT_FAILURE);
	}

	//implied else, fopen successful
	//read each line from file into entry in words array
	int i2=0;
	while(fgets(words2[i2],MAX_COLUMNS,fp2))
	{
		//remove trailing newline from string
		words2[i2][strlen(words2[i2])-1]='\0';
		i2++;
	}
	//'i2' contains number of valid entries in words2[][]
	//sort the array of strings
	bubbleSortWordsArray2(i2);

	fclose(fp1);
	fclose(fp2);

	if((strcmp(words1[0],words2[0])!=0)&&(strcmp(words1[2],words2[2])!=0)&&(strcmp(words1[3],words2[3])!=0)&&
		(strcmp(words1[4],words2[4])!=0)&&(strcmp(words1[5],words2[5])!=0)&&(strcmp(words1[6],words2[6])!=0)&&
		(strcmp(words1[7],words2[7])!=0)&&(strcmp(words1[8],words2[8])!=0)&&(strcmp(words1[9],words2[9])!=0)&&
		(strcmp(words1[10],words2[10])!=0)&&(strcmp(words1[11],words2[11])!=0)&&(strcmp(words1[12],words2[12])!=0)&&
		(strcmp(words1[13],words2[13])!=0)&&(strcmp(words1[14],words2[14])!=0)&&(strcmp(words1[15],words2[15])!=0)&&
		(strcmp(words1[16],words2[16])!=0)&&(strcmp(words1[17],words2[17])!=0)&&(strcmp(words1[18],words2[18])!=0)&&
		(strcmp(words1[19],words2[19])!=0)&&(strcmp(words1[20],words2[20])!=0)&&(strcmp(words1[21],words2[21])!=0)&&
		(strcmp(words1[22],words2[22])!=0)&&(strcmp(words1[23],words2[23])!=0)&&(strcmp(words1[24],words2[24])!=0)&&
		(strcmp(words1[25],words2[25])!=0)&&(strcmp(words1[26],words2[26])!=0)&&(strcmp(words1[27],words2[27])!=0)&&
		(strcmp(words1[28],words2[28])!=0)&&(strcmp(words1[29],words2[29])!=0)&&(strcmp(words1[30],words2[30])!=0)&&
		(strcmp(words1[31],words2[31])!=0)&&(strcmp(words1[32],words2[32])!=0)&&(strcmp(words1[33],words2[33])!=0)&&
		(strcmp(words1[34],words2[34])!=0)&&(strcmp(words1[35],words2[35])!=0)&&(strcmp(words1[36],words2[36])!=0)&&
		(strcmp(words1[37],words2[37])!=0)&&(strcmp(words1[38],words2[38])!=0)&&(strcmp(words1[39],words2[39])!=0)&&
		(strcmp(words1[40],words2[40])!=0)&&(strcmp(words1[41],words2[41])!=0)&&(strcmp(words1[42],words2[42])!=0)&&
		(strcmp(words1[43],words2[43])!=0)&&(strcmp(words1[44],words2[44])!=0)&&(strcmp(words1[45],words2[45])!=0)&&
		(strcmp(words1[46],words2[46])!=0)&&(strcmp(words1[47],words2[47])!=0)&&(strcmp(words1[48],words2[48])!=0)&&
		(strcmp(words1[49],words2[49])!=0)&&(strcmp(words1[50],words2[50])!=0)&&(strcmp(words1[51],words2[51])!=0)&&
		(strcmp(words1[52],words2[52])!=0)&&(strcmp(words1[53],words2[53])!=0)&&(strcmp(words1[54],words2[54])!=0)&&
		(strcmp(words1[55],words2[55])!=0)&&(strcmp(words1[56],words2[56])!=0)&&(strcmp(words1[57],words2[57])!=0)&&
		(strcmp(words1[58],words2[58])!=0)&&(strcmp(words1[59],words2[59])!=0)&&(strcmp(words1[60],words2[60])!=0)&&
		(strcmp(words1[61],words2[61])!=0)&&(strcmp(words1[62],words2[62])!=0)&&(strcmp(words1[63],words2[63])!=0)&&
		(strcmp(words1[64],words2[64])!=0)&&(strcmp(words1[65],words2[65])!=0)&&(strcmp(words1[66],words2[66])!=0)&&
		(strcmp(words1[67],words2[67])!=0)&&(strcmp(words1[68],words2[68])!=0)&&(strcmp(words1[69],words2[69])!=0)&&
		(strcmp(words1[70],words2[70])!=0)&&(strcmp(words1[71],words2[71])!=0)&&(strcmp(words1[72],words2[72])!=0)&&
		(strcmp(words1[73],words2[73])!=0)&&(strcmp(words1[74],words2[74])!=0)&&(strcmp(words1[75],words2[75])!=0)&&
		(strcmp(words1[76],words2[76])!=0)&&(strcmp(words1[77],words2[77])!=0)&&(strcmp(words1[78],words2[78])!=0)&&
		(strcmp(words1[79],words2[79])!=0)&&(strcmp(words1[80],words2[80])!=0)&&(strcmp(words1[81],words2[81])!=0)&&
		(strcmp(words1[82],words2[82])!=0)&&(strcmp(words1[83],words2[83])!=0)&&(strcmp(words1[84],words2[84])!=0)&&
		(strcmp(words1[85],words2[85])!=0)&&(strcmp(words1[86],words2[86])!=0)&&(strcmp(words1[87],words2[87])!=0)&&
		(strcmp(words1[88],words2[88])!=0)&&(strcmp(words1[89],words2[89])!=0)&&(strcmp(words1[90],words2[90])!=0)&&
		(strcmp(words1[91],words2[91])!=0)&&(strcmp(words1[92],words2[92])!=0)&&(strcmp(words1[93],words2[93])!=0)&&
		(strcmp(words1[94],words2[94])!=0)&&(strcmp(words1[95],words2[95])!=0)&&(strcmp(words1[96],words2[96])!=0)&&
		(strcmp(words1[97],words2[97])!=0)&&(strcmp(words1[98],words2[98])!=0)&&(strcmp(words1[99],words2[99])!=0)&&
		(strcmp(words1[100],words2[100])!=0)&&(strcmp(words1[101],words2[101])!=0)&&(strcmp(words1[102],words2[102])!=0)&&
		(strcmp(words1[103],words2[103])!=0)&&(strcmp(words1[104],words2[104])!=0)&&(strcmp(words1[105],words2[105])!=0)&&
		(strcmp(words1[106],words2[106])!=0)&&(strcmp(words1[107],words2[107])!=0)&&(strcmp(words1[108],words2[108])!=0)&&
		(strcmp(words1[109],words2[109])!=0)&&(strcmp(words1[110],words2[110])!=0)&&(strcmp(words1[111],words2[111])!=0)&&
		(strcmp(words1[112],words2[112])!=0)&&(strcmp(words1[113],words2[113])!=0)&&(strcmp(words1[114],words2[114])!=0)&&
		(strcmp(words1[119],words2[119])!=0)&&(strcmp(words1[120],words2[120])!=0)&&(strcmp(words1[121],words2[121])!=0)&&
		(strcmp(words1[122],words2[122])!=0)&&(strcmp(words1[123],words2[123])!=0)&&(strcmp(words1[124],words2[124])!=0)&&
		(strcmp(words1[125],words2[125])!=0)&&(strcmp(words1[126],words2[126])!=0)&&(strcmp(words1[127],words2[127])!=0)&&
		(strcmp(words1[128],words2[128])!=0)&&(strcmp(words1[129],words2[129])!=0)&&(strcmp(words1[130],words2[130])!=0)&&
		(strcmp(words1[131],words2[131])!=0)&&(strcmp(words1[132],words2[132])!=0)&&(strcmp(words1[133],words2[133])!=0))
		return 0;
	return 1;
}

int compareGtf(char f1[],char f2[]){
	FILE *fp1 = NULL;
	FILE *fp2 = NULL;

	if((fp1=fopen(f1,"r"))==NULL)
	{
		perror("fopen failed");
		exit(EXIT_FAILURE);
	}

	//implied else, fopen successful
	//read each line from file into entry in words array
	int i1=0;
	while(fgets(words1[i1],MAX_COLUMNS,fp1))
	{
		//remove trailing newline from string
		words1[i1][strlen(words1[i1])-1]='\0';
		i1++;
	}
	//'i1' contains number of valid entries in words1[][]
	//sort the array of strings
	bubbleSortWordsArray1(i1);
	
	if((fp2=fopen(f2,"r"))==NULL)
	{
		perror("fopen failed");
		exit(EXIT_FAILURE);
	}

	//implied else, fopen successful
	//read each line from file into entry in words array
	int i2=0;
	while(fgets(words2[i2],MAX_COLUMNS,fp2))
	{
		//remove trailing newline from string
		words2[i2][strlen(words2[i2])-1]='\0';
		i2++;
	}
	//'i2' contains number of valid entries in words2[][]
	//sort the array of strings
	bubbleSortWordsArray2(i2);

	fclose(fp1);
	fclose(fp2);
	
	int i;
	if(strcmp(words1[0],words2[0])!=0)
			return 0;
	for(i=2;i<255;i++)
	{
		if(strcmp(words1[i],words2[i])!=0)
			return 0;
	}
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


Test(output,findPatternTest_Test788) {
	char f1[]="regressionTest/test-788/referenceOutput/pintron-all-isoforms.gtf";
	char f2[]="regressionTest/test-788/executionOutput/pintron-all-isoforms.gtf";
	cr_expect(compareGtf(f1,f2)==1);

	char f3[]="regressionTest/test-788/referenceOutput/full.json";
	char f4[]="regressionTest/test-788/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}

Test(output,compareOutput_TestAMBN) {
	char f1[]="regressionTest/test-AMBN/referenceOutput/pintron-all-isoforms.gtf";
	char f2[]="regressionTest/test-AMBN/executionOutput/pintron-all-isoforms.gtf";
	cr_expect(compareGtfCr(f1,f2)==1);

	char f3[]="regressionTest/test-AMBN/referenceOutput/full.json";
	char f4[]="regressionTest/test-AMBN/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}

Test(output,compareOutput_TestCPB2) {
	char f1[]="regressionTest/test-CPB2/referenceOutput/pintron-all-isoforms.gtf";
	char f2[]="regressionTest/test-CPB2/executionOutput/pintron-all-isoforms.gtf";
	cr_expect(compareGtf(f1,f2)==1);

	char f3[]="regressionTest/test-CPB2/referenceOutput/full.json";
	char f4[]="regressionTest/test-CPB2/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}

Test(output,compareOutput_TestFeatureJsonWithESTalignments) {
	char f3[]="regressionTest/test-feature-json-with-EST-alignments/referenceOutput/full.json";
	char f4[]="regressionTest/test-feature-json-with-EST-alignments/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}

Test(output,compareOutput_TestGtf1) {
	char f3[]="regressionTest/test_gtf1/referenceOutput/full.json";
	char f4[]="regressionTest/test_gtf1/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}

Test(output,compareOutput_TestGtf2) {
	char f3[]="regressionTest/test_gtf2/referenceOutput/full.json";
	char f4[]="regressionTest/test_gtf2/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}

Test(output,compareOutput_TestGtf3) {
	char f1[]="regressionTest/test_gtf3/referenceOutput/pintron-all-isoforms.gtf";
	char f2[]="regressionTest/test_gtf3/executionOutput/pintron-all-isoforms.gtf";
	cr_expect(compareGtfCr(f1,f2)==1);

	char f3[]="regressionTest/test_gtf3/referenceOutput/full.json";
	char f4[]="regressionTest/test_gtf3/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}

Test(output,compareOutput_TestGtf4) {
	char f1[]="regressionTest/test_gtf4/referenceOutput/pintron-all-isoforms.gtf";
	char f2[]="regressionTest/test_gtf4/executionOutput/pintron-all-isoforms.gtf";
	cr_expect(compareGtf(f1,f2)==1);

	char f3[]="regressionTest/test_gtf4/referenceOutput/full.json";
	char f4[]="regressionTest/test_gtf4/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}

Test(output,compareOutput_TestGtf5) {
	char f3[]="regressionTest/test_gtf5/referenceOutput/full.json";
	char f4[]="regressionTest/test_gtf5/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}

Test(output,compareOutput_TestGtf6) {
	char f3[]="regressionTest/test_gtf6/referenceOutput/full.json";
	char f4[]="regressionTest/test_gtf6/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}

Test(output,compareOutput_TestGtf7) {
	char f3[]="regressionTest/test_gtf7/referenceOutput/full.json";
	char f4[]="regressionTest/test_gtf7/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}

Test(output,compareOutput_TestGtf8) {
	char f3[]="regressionTest/test_gtf8/referenceOutput/full.json";
	char f4[]="regressionTest/test_gtf8/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}

Test(output,compareOutput_TestIssue2) {
	char f1[]="regressionTest/test-issue-2/referenceOutput/pintron-all-isoforms.gtf";
	char f2[]="regressionTest/test-issue-2/executionOutput/pintron-all-isoforms.gtf";
	cr_expect(compareGtfCr(f1,f2)==1);

	char f3[]="regressionTest/test-issue-2/referenceOutput/full.json";
	char f4[]="regressionTest/test-issue-2/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}

Test(output,compareOutput_TestIssue13) {
	char f3[]="regressionTest/test-issue-13/referenceOutput/full.json";
	char f4[]="regressionTest/test-issue-13/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}

Test(output,compareOutput_TestIssue31) {
	char f3[]="regressionTest/test-issue-31/referenceOutput/full.json";
	char f4[]="regressionTest/test-issue-31/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}

Test(output,compareOutput_TestIssue39) {
	char f3[]="regressionTest/test-issue-39/referenceOutput/full.json";
	char f4[]="regressionTest/test-issue-39/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}

Test(output,compareOutput_TestMattia1) {
	char f1[]="regressionTest/test-mattia1/referenceOutput/pintron-all-isoforms.gtf";
	char f2[]="regressionTest/test-mattia1/executionOutput/pintron-all-isoforms.gtf";
	cr_expect(compareGtf(f1,f2)==1);

	char f3[]="regressionTest/test-mattia1/referenceOutput/full.json";
	char f4[]="regressionTest/test-mattia1/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}

Test(output,compareOutput_TestMattia2) {
	char f3[]="regressionTest/test-mattia2/referenceOutput/full.json";
	char f4[]="regressionTest/test-mattia2/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}

Test(output,compareOutput_TestMattia3) {
	char f1[]="regressionTest/test-mattia3/referenceOutput/pintron-all-isoforms.gtf";
	char f2[]="regressionTest/test-mattia3/executionOutput/pintron-all-isoforms.gtf";
	cr_expect(compareGtf(f1,f2)==1);

	char f3[]="regressionTest/test-mattia3/referenceOutput/full.json";
	char f4[]="regressionTest/test-mattia3/executionOutput/full.json";
	cr_expect(compareJson(f3,f4)==1);
}
