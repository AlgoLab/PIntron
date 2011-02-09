/**
 *
 *
 *                              PIntron
 *
 * A novel pipeline for computational gene-structure prediction based on
 * spliced alignment of expressed sequences (ESTs and mRNAs).
 *
 * Copyright (C) 2010  Yuri Pirola
 *
 * Distributed under the terms of the GNU Affero General Public License (AGPL)
 *
 *
 * This file is part of PIntron.
 *
 * PIntron is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PIntron is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with PIntron.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "util.h"

#include <stddef.h>
#include <stdbool.h>


void read_rows(FILE *fp,char* my_string)
{
  int n=0;
  char s1[2000],s2[2000],s3[2000],s4[2000],s5[2000],s6[20000],s7[200000],s8[300000],s9[300000];
  size_t bytes_read;
  size_t n_bytes=300000;

  FILE *dest1;
  FILE *dest2;
  bool sxdel=false,ciclo=true;
  int str_conv;

  int coord1,coord2;

  dest1= stdout;
  dest2= stderr;

  while(n!=9){
	 bytes_read=my_getline(&my_string,&n_bytes,fp);
	 n= sscanf(my_string,"%s %s %s %s %s %s %s %s %s",s1,s2,s3,s4,s5,s6,s7,s8,s9);
  }

  while(ciclo){
	 sxdel=false;

		if(n==9){

		  fprintf(dest1,">%s\n",s2);
		  fprintf(dest2,">%s\n",s2);

		  if(strcmp(s9,"X")!=0){
			 fprintf(dest1,"%s\n",s9);
			 fprintf(dest2,"%s\t %s\t %s\t %s\n",s4,s5,s6,s7);
		  }else{
			 sxdel=true;
			 str_conv=atoi(s5);
		  }
		  bytes_read= custom_getline(&my_string,&n_bytes,fp);
		  n=sscanf(my_string,"%s %s %s %s %s %s %s %s %s",s1,s2,s3,s4,s5,s6,s7,s8,s9);
		  my_string[strlen(my_string)-1]='\0';

		  while((n!=9)&&(strcmp(s1,"Distinct")!=0)){

			 bytes_read= custom_getline(&my_string,&n_bytes,fp);
			 my_string[strlen(my_string)-1]='\0';
			 n=sscanf(my_string,"%s %s %s %s %s %s %s %s %s",s1,s2,s3,s4,s5,s6,s7,s8,s9);

			 if((n==7)&&(strcmp(s7,"X")!=0)){
				fprintf(dest1,"%s\n",s7);

				if(sxdel==false)fprintf(dest2,"%s\t %s\t %s\t %s\n",s2,s3,s4,s5);

				else{
				  coord1=atoi(s2);
				  coord2=atoi(s3);

				  coord1=coord1-str_conv;
				  coord2=coord2-str_conv;

				  fprintf(dest2,"%d\t %d\t %s\t %s\n",coord1,coord2,s4,s5);
				}
			 }
		  }
		}else{
		  ciclo=false;
		}
  }
}

void test_create(FILE *fp)
{
 size_t bytes_read;
 size_t n_bytes=10000;
 char* my_string;

 my_string=c_palloc(n_bytes+1);

 while(!(feof(fp))){
	bytes_read= custom_getline(&my_string,&n_bytes,fp);
	my_string[strlen(my_string)-1]='\0';
	if(strcmp(my_string,"EST factorizations:")==0) read_rows(fp,my_string);
 }
}

int main()
{
  test_create(stdin);
  return 0;
}

