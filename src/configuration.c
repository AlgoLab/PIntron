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
/*
** configuration.c
*/

#include "configuration.h"
#include "util.h"
#include "log.h"
#include "options.h"
#include <unistd.h>


#define __SAVE_CONFIG_FILE__ "./config-dump.ini"


static pconfiguration
check_and_copy(struct gengetopt_args_info* args) {
  DEBUG("Creating the struct for the configuration parameters.");
  pconfiguration config= PALLOC(struct _configuration);

  fail_if(args->min_factor_length_arg<=0);
  config->min_factor_len= args->min_factor_length_arg;
  INFO("CONFIG: The minimum factor lenght is %d", config->min_factor_len);

  fail_if(args->min_intron_length_arg<0);
  config->min_intron_length=args->min_intron_length_arg;
  INFO("CONFIG: The minimum length for an intron is %d", config->min_intron_length);

  fail_if(args->max_intron_length_arg<0);
  config->max_intron_length=args->max_intron_length_arg;
  INFO("CONFIG: The maximum length for an intron is %d", config->max_intron_length);

  fail_if(args->min_string_depth_rate_arg<0.0 ||
			 args->min_string_depth_rate_arg>1.0);
  config->min_string_depth_rate= args->min_string_depth_rate_arg;
  INFO("CONFIG: The minimum string depth rate is %f", config->min_string_depth_rate);

  fail_if(args->max_prefix_discarded_rate_arg<0.0 ||
			 args->max_prefix_discarded_rate_arg>1.0);
  config->max_prefix_discarded_rate= args->max_prefix_discarded_rate_arg;
  INFO("CONFIG: The maximum prefix that can be discarded from an EST is %f (%%) "
		 "(to be used in source links)", config->max_prefix_discarded_rate);

  fail_if(args->max_suffix_discarded_rate_arg<0.0 ||
			 args->max_suffix_discarded_rate_arg>1.0);
  config->max_suffix_discarded_rate= args->max_suffix_discarded_rate_arg;
  INFO("CONFIG: The maximum suffix that can be discarded from an EST is %f (%%) "
		 "(to be used in sink links)", config->max_suffix_discarded_rate);

  fail_if(args->max_prefix_discarded_arg<0);
  config->max_prefix_discarded= args->max_prefix_discarded_arg;
  INFO("CONFIG: The maximum prefix that can be discarded from an EST is %d "
		 "(to be used during factorization construction)",
		 config->max_prefix_discarded);

  fail_if(args->max_suffix_discarded_arg<0);
  config->max_suffix_discarded= args->max_suffix_discarded_arg;
  INFO("CONFIG: The maximum suffix that can be discarded from an EST is %d "
		 "(to be used during factorization construction)",
		 config->max_suffix_discarded);

  fail_if(args->min_distance_of_splice_sites_arg<0);
  config->max_site_difference= args->min_distance_of_splice_sites_arg;
  INFO("CONFIG: The maximum difference (nt) in order to consider the same ss "
		 "two splicing sites: %d", config->max_site_difference);

  fail_if(args->max_no_of_factorizations_arg<0);
  config->max_number_of_factorizations= args->max_no_of_factorizations_arg;
  INFO("CONFIG: The maximum number of factorizations allowed for an EST that is "
		 "not an artifact (if set to 0, checking not performed): %d",
		 config->max_number_of_factorizations);

  fail_if(args->max_difference_of_coverage_arg<0.0 ||
			 args->max_difference_of_coverage_arg>1.0);
  config->max_coverage_diff= args->max_difference_of_coverage_arg;
  INFO("CONFIG: The maximum coverage difference for accepting a factorization: %f",
		 config->max_coverage_diff);

  fail_if(args->max_difference_of_no_of_exons_arg<-1);
  config->max_exonNUM_diff= args->max_difference_of_no_of_exons_arg;
  INFO("CONFIG: The minimum exon number difference for accepting a factorization: %d",
		 config->max_exonNUM_diff);

  fail_if(args->max_difference_of_gap_length_arg<-1);
  config->max_gapLength_diff= args->max_difference_of_gap_length_arg;
  INFO("CONFIG: The minimum gap length difference for accepting a factorization: %d",
		 config->max_gapLength_diff);

  config->retain_externals=
	 (args->retain_externals_arg == retain_externals_arg_true) ?
	 1 : 0;
  INFO("CONFIG: The external factors are retained: %s",
		 (config->retain_externals == 1)?("yes"):("no"));

  //The maximum number of pairings that can compose the vertex
  //set of a MEG
  // AND
  //The maximum frequency of the shortest pairings in a vertex
  //set of a MEG
  //Rationale: if more than max_freq*tot pairings have minimum
  //length and we have more than max pairings,
  //than it is plausible we are allowing too short pairings.
  fail_if(args->max_pairings_in_CMEG_arg<0);
  fail_if(args->max_shortest_pairing_frequence_arg < 0.0 ||
			 args->max_shortest_pairing_frequence_arg>1.0);
  config->max_pairings_in_MEG= args->max_pairings_in_CMEG_arg;
  config->max_freq_shortest_pairing= args->max_shortest_pairing_frequence_arg;
  INFO("CONFIG: We are allowing MEG with at most %d vertices OR "
		 "MEGs that have at most %f%% of minimum length vertices.",
		 config->max_pairings_in_MEG,
		 config->max_freq_shortest_pairing*100);

  fail_if(args->suff_pref_length_intron_arg<=0);
  config->suffpref_length_for_intron=args->suff_pref_length_intron_arg;
  INFO("CONFIG: The suffix/prefix intron length for intron gap alignment: %d.",
		 config->suffpref_length_for_intron);

  fail_if(args->suff_pref_length_est_arg<=0);
  config->suffpref_length_on_est=args->suff_pref_length_est_arg;
  INFO("CONFIG: The suffix/prefix est factor length for intron gap alignment: %d.",
		 config->suffpref_length_on_est);

  fail_if(args->suff_pref_length_genomic_arg<=0);
  config->suffpref_length_on_gen=args->suff_pref_length_genomic_arg;
  INFO("CONFIG: The suffix/prefix exon length for intron gap alignment: %d.",
		 config->suffpref_length_on_gen);

  config->trans_red= !args->no_transitive_reduction_flag;
  INFO("CONFIG: Perform the transitive reduction? %s.",
		 config->trans_red?"yes":"no");

  config->short_edge_comp= !args->no_short_edge_compaction_flag;
  INFO("CONFIG: Perform the short-edge compaction? %s.",
		 config->short_edge_comp?"yes":"no");

  fail_if(args->max_single_factorization_time_arg<0);
  config->max_single_factorization_time= args->max_single_factorization_time_arg;
  INFO("CONFIG: Maximum time for computing a factorization of a single transcript: %u.",
		 config->max_single_factorization_time);

  return config;
}

#define FIELD(opt, suff) args_info.opt##_##suff

#define COPY_char_VALUE(opt)										\
  {																		\
	 FIELD(opt,orig)= alloc_and_copy(FIELD(opt,arg));		\
	 FIELD(opt,given)= 1;											\
  }																		\

#define COPY_int_VALUE(opt)								\
  {																\
	 char buff[25];											\
	 snprintf(buff, 24, "%d", FIELD(opt,arg));		\
	 FIELD(opt,orig)= alloc_and_copy(buff);			\
	 FIELD(opt,given)= 1;									\
  }

#define COPY_double_VALUE(opt)								\
  {																	\
	 char buff[15];												\
	 snprintf(buff, 14, "%.10f", FIELD(opt,arg));		\
	 size_t p= strlen(buff);									\
	 while (p>1 && buff[p-1]=='0' && buff[p-2]!='.') {	\
		--p;															\
		buff[p]='\0';												\
	 }																	\
	 FIELD(opt,orig)= alloc_and_copy(buff);				\
	 FIELD(opt,given)= 1;										\
  }

pconfiguration
config_clone(pconfiguration src) {

  pconfiguration config= PALLOC(struct _configuration);

  config->min_factor_len= src->min_factor_len;
  config->min_intron_length= src->min_intron_length;
  config->max_intron_length= src->max_intron_length;
  config->min_string_depth_rate= src->min_string_depth_rate;
  config->max_prefix_discarded_rate= src->max_prefix_discarded_rate;
  config->max_suffix_discarded_rate= src->max_suffix_discarded_rate;
  config-> max_prefix_discarded= src-> max_prefix_discarded;
  config-> max_suffix_discarded= src-> max_suffix_discarded;
  config->max_site_difference= src->max_site_difference;
  config->max_number_of_factorizations= src->max_number_of_factorizations;
  config->max_coverage_diff= src->max_coverage_diff;
  config->max_exonNUM_diff= src->max_exonNUM_diff;
  config->max_gapLength_diff= src->max_gapLength_diff;
  config->retain_externals= src->retain_externals;
  config->max_pairings_in_MEG= src->max_pairings_in_MEG;
  config->max_freq_shortest_pairing= src->max_freq_shortest_pairing;
  config->suffpref_length_on_est= src->suffpref_length_on_est;
  config->suffpref_length_for_intron= src->suffpref_length_for_intron;
  config->suffpref_length_on_gen= src->suffpref_length_on_gen;
  config->trans_red= src->trans_red;
  config->short_edge_comp= src->short_edge_comp;
  config->max_single_factorization_time= src->max_single_factorization_time;

  return config;
}


void config_destroy(pconfiguration config) {
  DEBUG("Destroying the struct for the configuration parameters.");
  my_assert(config!=NULL);
  pfree(config);
}

pconfiguration config_create(int argc, char** argv) {
  struct gengetopt_args_info args_info;

// Set parameters to the command line values
  struct cmdline_parser_params *params;
  params = cmdline_parser_params_create();
  params->check_required= 0;
  params->initialize= 1;
  params->override= 1;
  if (cmdline_parser_ext(argc, argv, &args_info, params) != 0) {
	 FATAL("Impossible to initialize parameter reader. Terminating.");
	 fail();
  }

// Set parameters to configuration file values (without ovverriding)
  if (access(args_info.config_file_arg, R_OK)==0) {
	 INFO("The configuration file %s is present. "
			"Parameters from the configuration file will "
			"not override the command line!", args_info.config_file_arg);
	 params->initialize= 0;
	 params->override= 0;
	 cmdline_parser_config_file(args_info.config_file_arg, &args_info, params);
  } else {
	 INFO("Configuration file %s not found.", args_info.config_file_arg);
  }

  cmdline_parser_required(&args_info, argv[0]);

// Save parameters
// XXX: Add a line for each parameter !!
  COPY_char_VALUE(config_file);
//  COPY_char_VALUE(org);
//  COPY_char_VALUE(gene);
//  COPY_int_VALUE(strand);
//  COPY_char_VALUE(chromosome);
//  COPY_char_VALUE(est_cluster);
  COPY_int_VALUE(min_factor_length);
  COPY_int_VALUE(min_intron_length);
  COPY_int_VALUE(max_intron_length);
  COPY_double_VALUE(min_string_depth_rate);
  COPY_double_VALUE(max_prefix_discarded_rate);
  COPY_double_VALUE(max_suffix_discarded_rate);
  COPY_int_VALUE(max_prefix_discarded);
  COPY_int_VALUE(max_suffix_discarded);
  COPY_int_VALUE(min_distance_of_splice_sites);
  COPY_int_VALUE(max_no_of_factorizations);
  COPY_double_VALUE(max_difference_of_coverage);
  COPY_int_VALUE(max_difference_of_no_of_exons);
  COPY_int_VALUE(max_difference_of_gap_length);
  COPY_int_VALUE(max_pairings_in_CMEG);
  COPY_double_VALUE(max_shortest_pairing_frequence);
  COPY_int_VALUE(suff_pref_length_genomic);
  COPY_int_VALUE(suff_pref_length_est);
  COPY_int_VALUE(suff_pref_length_intron);
  COPY_int_VALUE(max_single_factorization_time);
//  COPY_int_VALUE(max_seq_in_gst);

  args_info.retain_externals_orig=
	 alloc_and_copy((args_info.retain_externals_arg==retain_externals_arg_true) ?
						 "true": "false");
  args_info.retain_externals_given= 1;

  pconfiguration config= check_and_copy(&args_info);

  if (cmdline_parser_file_save(__SAVE_CONFIG_FILE__, &args_info)!=0) {
	 WARN("Configuration NOT saved in file " __SAVE_CONFIG_FILE__ ".");
  } else {
	 INFO("Configuration saved in file " __SAVE_CONFIG_FILE__ ".");
  }

  cmdline_parser_free(&args_info);
  free(params);

  return config;
}
