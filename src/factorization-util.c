/**
 *
 *
 *                              PIntron
 *
 * A novel pipeline for computational gene-structure prediction based on
 * spliced alignment of expressed sequences (ESTs and mRNAs).
 *
 * Copyright (C) 2012  Yuri Pirola
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
#include "factorization-util.h"

#include "log.h"
#include "types.h"

void
_impl_print_factorization_on_log_full(const int log_level,
												  plist factorization, const char* const gen_seq){
  plistit plist_it_f;

  plist_it_f=list_first(factorization);
  size_t exon_n= 0;
  while(listit_has_next(plist_it_f)){
	 pfactor pfact=(pfactor)listit_next(plist_it_f);
	 ++exon_n;
	 ALWAYS_LOG(log_level, "   %3zu) E_START=%8d  E_END=%8d  --  G_START=%10d  G_END=%10d  %.2s...%.2s",
					exon_n, pfact->EST_start, pfact->EST_end, pfact->GEN_start, pfact->GEN_end,
					(exon_n==1)||(gen_seq==NULL)?"  ":gen_seq+pfact->GEN_start-2,
					!listit_has_next(plist_it_f)||(gen_seq==NULL)?"  ":gen_seq+pfact->GEN_end+1);
  }
  listit_destroy(plist_it_f);
}

void
_impl_print_factorizations_on_log_full(const int log_level,
													plist factorization_list, const char* const gen_seq){
  plistit plist_it_id= list_first(factorization_list);
  size_t factorization_no= 1;
  while(listit_has_next(plist_it_id)){
	 plist factorization=(plist)listit_next(plist_it_id);
	 ALWAYS_LOG(log_level, "  Factorization %-3zu [%2zu exon(s)]",
					factorization_no, list_size(factorization));
	 _impl_print_factorization_on_log_full(log_level, factorization, gen_seq);
	 ++factorization_no;
  }
  listit_destroy(plist_it_id);
}


void
_impl_print_factorization_on_log(const int log_level,
											plist factorization) {
  _impl_print_factorization_on_log_full(log_level, factorization, NULL);
}

void
_impl_print_factorizations_on_log(const int log_level,
											 plist factorization_list) {
  _impl_print_factorizations_on_log_full(log_level, factorization_list, NULL);
}

