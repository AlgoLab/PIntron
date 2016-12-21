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
#include "simpl_info.h"
#include "util.h"
#include "log.h"
#include "bit_vector.h"
#include <stdio.h>

void psimpl_destroy(psimpl p)
{
  if(p!=NULL){
	 BV_destroy(p->factors_used);
	 BV_destroy(p->factors_not_used);
	 BV_destroy(p->ests_ok);
	 pfree(p);
  }
}

psimpl psimpl_create(void)
{
  psimpl p=(psimpl)palloc(sizeof(struct _simpl_info));
  return p;
}

int countTrue(pbit_vect bv)
{
  int cont=0;
  for(unsigned int i=0;i<bv->n;i++){
	 if(BV_get(bv,i)==true){
		cont=cont+1;
	 }
  }
  return cont;
}

#if defined (LOG_MSG) && (LOG_LEVEL_INFO <= LOG_THRESHOLD)

void psimpl_print(psimpl p)
{
  INFO("Simplified factors: %d", countTrue(p->factors_used));
  INFO("Remaining factors:  %d", countTrue(p->factors_not_used));
  INFO("No. of sequences factorized during simplification: %d", countTrue(p->ests_ok));
}

#endif
