/**
 *
 *
 *                              PIntron
 *
 * A novel pipeline for computational gene-structure prediction based on
 * spliced alignment of expressed sequences (ESTs and mRNAs).
 *
 * Copyright (C) 2010  Yuri Pirola, Raffaella Rizzi
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
#include "util.h"
#include "log.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>
//#include <sys/types.h>
//#include <sys/stat.h>


#ifdef __APPLE__

ssize_t custom_getline(char **lineptr, size_t *n, FILE *stream) {
  char *p;
  ssize_t size;
  int c;

  if (lineptr == NULL) {
	 return -1;
  }
  if (stream == NULL) {
	 return -1;
  }
  if (n == NULL) {
	 return -1;
  }
  p= *lineptr;
  size= *n;

  c = fgetc(stream);
  if (c == EOF) {
	 return -1;
  }
  ssize_t readchar= 0;
  while(c != EOF) {
	 FINETRACE("Reading %u-th character %c", readchar, c);
	 if (readchar >= size) {
		size= size + 128;
		TRACE("Extending to %lu bytes.", size);
		p= realloc(p, size);
		if (p == NULL) {
		  return -1;
		}
	 }
	 p[readchar] = c;
	 ++readchar;
	 if (c == '\n') {
		break;
	 }
	 c = fgetc(stream);
  }

  if (readchar<size) {
	 p[readchar]= '\0';
  } else {
	 size= size+1;
	 p= realloc(p, size);
	 p[readchar]= '\0';
  }
  *lineptr= p;
  *n= size;

  return readchar;
}

#endif

char* c_palloc(size_t size) {
  return (char*)palloc(size*sizeof(char));
}

void pfree(const void* const p) {
  if (p==NULL) {
	 FATAL("Cannot free a NULL pointer.");
	 fail();
  }
  free((void *)p);
  FINETRACE("Pointer %p deallocated", p);
}

char* alloc_and_copy(const char* const source) {
  size_t len= strlen(source);
  char* ris= c_palloc(len+1);
  strncpy(ris, source, len+1);
  return ris;
}

void noop_free(void* el) {
  (void)el;
}

char* substring(const int index, const char* const string){
  my_assert(index>=0);
  my_assert(string != NULL);
  FINETRACE("index value is: %d", index);
  FINETRACE("string value is: |%s|", string);
  size_t slen= strlen(string);
  char* ris = c_palloc(slen-index+1);
  memcpy(ris, string+index, (slen-index)*sizeof(char));
  ris[slen-index] = '\0';
  FINETRACE("the resulting substring is: |%s|", ris);
  return ris;
}//end subString

ssize_t my_getline(char **lineptr, size_t *n, FILE *stream) {
  ssize_t ris= custom_getline(lineptr, n, stream);
  while (ris>0 && (*lineptr)[ris-1]<' ') {
	 --ris;
	 (*lineptr)[ris]= '\0';
  }
  return ris;
}

void print_repetitions(FILE* f, const char c, int rep) {
  char* s= c_palloc(rep+1);
  for (int i= 0; i<rep; ++i)
	 s[i]= c;
  s[rep]='\0';
  fprintf(f, "%s", s);
  pfree(s);
}

void resource_usage_log(void) {
  struct rusage *rusage= PALLOC(struct rusage);
  if (getrusage(RUSAGE_SELF, rusage)==0) {
	 INFO("User time:   %10lds %7ldmicrosec.", rusage->ru_utime.tv_sec, rusage->ru_utime.tv_usec);
	 INFO("System time: %10lds %7ldmicrosec.", rusage->ru_stime.tv_sec, rusage->ru_stime.tv_usec);
#ifdef __APPLE__
	 INFO("Mem. used:         NaN KB");
#else
	 char buf[1000];
	 snprintf(buf, 1000, "/proc/%u/statm", (unsigned)getpid());
	 FILE* pf = fopen(buf, "r");
	 if (pf) {
		unsigned int size; //       total program size
		int ret_val= fscanf(pf, "%u", &size);
		if (ret_val==1) {
		  INFO("Mem. used:   %10uKB", size);
		} else {
		  INFO("Mem. used:         NaN KB");
		}
	 }
	 fclose(pf);
#endif
  }
  pfree(rusage);
}

FILE* open_statm_file(void) {
  static char buf[1000];
  snprintf(buf, 1000, "/proc/%u/statm", (unsigned)getpid());
  FILE* pf = fopen(buf, "r");
  if (pf==NULL) {
	 FATAL("Impossible to open the statm file of the current process %s.", buf);
	 fail();
  }
  return pf;
}

void log_info(FILE* const logfile, char* description) {
  my_assert(logfile!=NULL);
#ifdef __APPLE__
  static struct timeval tv;
  if (description==NULL)
	 description= "not-specified";
  gettimeofday(&tv, NULL);
  INFO("log-information %s\t%lu\tNaN", description, (unsigned long)tv.tv_sec);
  fprintf(logfile, "%s\t%lu\tNaN\n", description, (unsigned long)tv.tv_sec);
#else
  static FILE* statmfile= NULL;
  static size_t bsize= 100;
  static char* buff= NULL;
  static struct timeval tv;
  if (buff==NULL) {
	 buff= c_palloc(bsize);
  } else {
	 buff[0]= '\0';
  }
  if (statmfile==NULL) {
	 statmfile= open_statm_file();
  } else {
	 rewind(statmfile);
  }
  if (description==NULL)
	 description= "not-specified";
  my_getline(&buff, &bsize, statmfile);
  gettimeofday(&tv, NULL);
  INFO("log-information %s\t%lu\t%s", description, (unsigned long)tv.tv_sec, buff);
  fprintf(logfile, "%s\t%lu\t%s\n", description, (unsigned long)tv.tv_sec, buff);
#endif
}
