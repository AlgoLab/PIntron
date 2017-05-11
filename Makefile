####
#
#
#                              PIntron
#
# A novel pipeline for computational gene-structure prediction based on
# spliced alignment of expressed sequences (ESTs and mRNAs).
#
# Copyright (C) 2010  Gianluca Della Vedova, Yuri Pirola, Raffaella Rizzi
#
# Distributed under the terms of the GNU Affero General Public License (AGPL)
#
#
# This file is part of PIntron.
#
# PIntron is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PIntron is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with PIntron.  If not, see <http://www.gnu.org/licenses/>.
#
####
.PHONY: all reall main remain build clean stree est-fact reest-fact min-factorization remin-factorization test-data-create retest-data-create dist prepare-dist install cds-annotation recds-annotation max-transcr remaxtranscr # build-transcr rebuild-transcr

DEFAULT_STATUS=production
DEFAULT_PROF=no
DEFAULT_CHECKS_debug=yes
DEFAULT_CHECKS_production=no
DEFAULT_LOG_GRAPHS=no
DEFAULT_FORCE_32BIT=no
DEFAULT_USE_COLOR=no
DEFAULT_USE_PAR=no

PROG_DIR=/usr/local/bin

SRC_DIR= src
STREE_DIR= stree_src
INCLUDE_DIR= include

DIST_DOC_DIR= dist-docs
DIST_SCRIPTS_DIR= dist-scripts
DIST_SCRIPTS:= $(wildcard $(DIST_SCRIPTS_DIR)/*)
DIST_DIR= dist

ifneq ($(STATUS), debug)
ifneq ($(STATUS), production)
override STATUS := ${DEFAULT_STATUS}
$(warning Status not specified or invalid! Using default status "${STATUS}")
endif
endif

ifneq ($(PROF), yes)
ifneq ($(PROF), no)
override PROF := ${DEFAULT_PROF}
$(warning Profiling status not specified or invalid! Using default profiling status "${PROF}")
endif
endif

ifneq ($(CHECKS), yes)
ifneq ($(CHECKS), no)
override CHECKS := $(DEFAULT_CHECKS_$(STATUS))
$(warning Checks status not specified or invalid! Using default checks status "${CHECKS}")
endif
endif

ifneq ($(LOG_GRAPHS), yes)
ifneq ($(LOG_GRAPHS), no)
override LOG_GRAPHS := $(DEFAULT_LOG_GRAPHS)
$(warning Logging graph status not specified or invalid! Using default logging graph status "${LOG_GRAPHS}")
endif
endif

ifneq ($(USE_COLOR), yes)
ifneq ($(USE_COLOR), no)
override USE_COLOR:= $(DEFAULT_USE_COLOR)
endif
endif

#####################
## 64 BIT CHECK CODE
tmp_is_64_bit:=/$(shell uname -m)/
ifeq ($(tmp_is_64_bit), /x86_64/)
is_64_bit:=yes
else
ifeq ($(tmp_is_64_bit), /ia64/)
is_64_bit:=yes
else
is_64_bit:=no
endif
endif

ifeq ($(is_64_bit), yes)
ifneq ($(FORCE_32BIT), yes)
ifneq ($(FORCE_32BIT), no)
override FORCE_32BIT:=$(DEFAULT_FORCE_32BIT)
$(warning The machine is running 64 bit kernel but the pointer size has not been specified! Using 32bit? $(FORCE_32BIT))
endif
endif
ifeq ($(FORCE_32BIT), yes)
override WLBIT:=-m32
override compbit:=32
else
override WLBIT:=-m64
override compbit:=64
endif
else
override WLBIT:=
override compbit:=32
endif

## END 64 BIT CHECK CODE
#####################

ifeq ($(USE_COLOR), yes)
PHF:=\033[01;32m
PF:=\033[32m
SF:=\033[00m
else
PHF:=
PF:=
SF:=
endif

use_gcc=no
ifeq ($(CC), gcc)
use_gcc=yes
endif
ifeq ($(CC), cc)
use_gcc=yes
endif

ifeq ($(use_gcc), yes)
ifeq ($(USE_COLOR), yes)
have_colorgcc:=/$(shell which colorgcc)/
ifeq ($(have_colorgcc), //)
#we do not have colorgcc
else
override CC=colorgcc
endif
endif
endif

#####################
# MacOS X detection and configuration
#
uname:=/$(shell uname)/
ifeq ($(uname), /Darwin/)
current_sys:=MacOS
full_current_sys:=MacOS
else
current_sys:=Linux
full_current_sys:=Linux-$(compbit)bit
endif
#####################

#####################
# Architecture-dependent flags
#
ifeq ($(current_sys), MacOS)
CFLAGS_ARCH:=-isysroot /Developer/SDKs/MacOSX10.4u.sdk -arch ppc -arch ppc64 -arch i386 -arch x86_64
LDFLAGS_ARCH:=-isysroot /Developer/SDKs/MacOSX10.4u.sdk -arch ppc -arch ppc64 -arch i386 -arch x86_64
else
ifeq ($(current_sys), Linux)
CFLAGS_ARCH:=
LDFLAGS_ARCH:=
# Customize compiler options for production code
ifeq ($(STATUS), production)
OPTP+=-march=native
endif
else
CFLAGS_ARCH:=
LDFLAGS_ARCH:=
endif
endif
#####################


#####################
# TCMalloc detection
# TCMalloc is a (generally faster) drop-in replacement for the standard malloc
# See: http://code.google.com/p/gperftools/
#
HAS_TCMALLOC:=/$(shell echo "void main() {}" | $(CC) -x c -o /tmp/test - -ltcmalloc_minimal 2> /dev/null && echo yes || echo no)/
#####################



#####################
# Common compiler options
#
ifeq ($(CC), icc)
OPTB=-g -std=c99 -fp-model source -Wcheck -diag-enable warn -pipe
else
OPTB=-g -std=c99 -Wall -Wextra -pedantic -pipe
endif
#####################

#####################
# Compiler options for production code
#
ifeq ($(STATUS), production)
ifeq ($(CC), icc)
OPTP=-O2 -xHost -ip -ipo -vec-report1
else
OPTP=-O2
endif
LOG_LEVEL=INFO
else
OPTP=
endif
#####################

#####################
# Compiler options for debugging code
#
ifeq ($(STATUS), debug)
ifeq ($(CC), icc)
OPTD=-O0
else
OPTD=-O0
endif
LOG_LEVEL=DEBUG
else
OPTD=
endif
#####################

#####################
# Compiler options for profiling code
#
#
ifeq ($(PROF), yes)
ifeq ($(CC), icc)
OPTPROF=-pg
else
OPTPROF=-pg -fprofile-arcs
endif
else
OPTPROF=
ifeq ($(STATUS), production)
OPTP+=-fomit-frame-pointer
endif
endif
#####################

#####################
# Compiler options for checking code
#
#
ifeq ($(CHECKS), yes)
DCHECK=-UNDEBUG
else
DCHECK=-DNDEBUG
endif
#####################

#####################
# Compiler options for logging graphs
#
#
ifeq ($(LOG_GRAPHS), yes)
DLOGGRAPHS=-DLOG_GRAPHS
else
DLOGGRAPHS=-ULOG_GRAPHS
endif
#####################


COMPFLAGS=$(ARCH_DEP) $(WLBIT) $(OPTB) $(OPTP) $(OPTD) $(OPTPROF)

INCLUDE=-I. -I$(INCLUDE_DIR)/ -I$(STREE_DIR)/

LIBS=-lm #-lgsl -lgslcblas #-lefence

ifeq ($(STATUS), production)
ifeq ($(HAS_TCMALLOC), /yes/)
$(warning Using TCMalloc)
LIBS+=-ltcmalloc_minimal
endif
endif

override ADD_DFLAGS:=$(patsubst %,-%,$(subst :, ,$(PERS_DEFINE)))

ifeq (XX$(PERS_DEFINE)XX, XXXX)
PERS_DEFINE=NONE
endif

#####################
# Simboli del preprocessore
#
# NDEBUG Disabilita vari controlli
# LOG_MSG Abilita il log dei messaggi
# LOG_THRESHOLD Livello di visualizzazione dei messaggi di log.
# Puo' assumere il valore LOG_LEVEL_XXX con XXX uguale a
# FATAL, ERROR, WARN, INFO, DEBUG, TRACE
# Assume il valore di LOG_LEVEL_INFO se non definito.
#
DFLAGS=$(DCHECK) $(DLOGGRAPHS) -D_GNU_SOURCE -DLOG_MSG -DLOG_THRESHOLD=LOG_LEVEL_$(LOG_LEVEL) $(ADD_DFLAGS)
#####################

ADD_CFLAGS= $(COMPFLAGS) $(DFLAGS) $(INCLUDE) $(CFLAGS)

COMPILING_DESC=$(CC)-$(current_sys)-$(compbit)bit-status_$(STATUS)-prof_$(PROF)-log_$(LOG_LEVEL)-loggraphs_$(LOG_GRAPHS)-checks_$(CHECKS)-persdefine_$(subst :,_next_,$(PERS_DEFINE))

BASE_BIN_DIR= bin
BIN_DIR= $(BASE_BIN_DIR)/_tmp/$(COMPILING_DESC)
BASE_OBJ_DIR= obj
OBJ_DIR= $(BASE_OBJ_DIR)/$(COMPILING_DESC)

ALL_DIR= $(BIN_DIR) $(SRC_DIR) $(OBJ_DIR)



##
# New files must be inserted just after the options.c row
##
base_SOURCE= \
	$(SRC_DIR)/options.c \
	$(SRC_DIR)/log.c \
	$(SRC_DIR)/double_list.c \
	$(SRC_DIR)/int_list.c \
	$(SRC_DIR)/bool_list.c \
	$(SRC_DIR)/my_time.c \
	$(SRC_DIR)/ext_array.c \
	$(SRC_DIR)/configuration.c \
	$(SRC_DIR)/io-factorizations.c \
	$(SRC_DIR)/list.c \
	$(SRC_DIR)/io-multifasta.c \
	$(SRC_DIR)/types.c \
	$(SRC_DIR)/tree.c \
	$(SRC_DIR)/util.c \
	$(SRC_DIR)/bit_vector.c\
	$(SRC_DIR)/io-meg.c\
	$(SRC_DIR)/io-gen-ests.c\
	$(SRC_DIR)/refine-intron.c \
	$(SRC_DIR)/refine.c \
	$(SRC_DIR)/compute-alignments.c \
	$(SRC_DIR)/exon-complexity.c \
	$(SRC_DIR)/conversions.c \
	$(SRC_DIR)/detect-polya.c \


##
# New files must be inserted just after the options.o row
##
base_OBJ= \
	$(OBJ_DIR)/options.o \
	$(OBJ_DIR)/log.o \
	$(OBJ_DIR)/double_list.o \
	$(OBJ_DIR)/int_list.o \
	$(OBJ_DIR)/bool_list.o \
	$(OBJ_DIR)/my_time.o \
	$(OBJ_DIR)/ext_array.o \
	$(OBJ_DIR)/configuration.o \
	$(OBJ_DIR)/io-factorizations.o \
	$(OBJ_DIR)/list.o \
	$(OBJ_DIR)/io-multifasta.o \
	$(OBJ_DIR)/types.o \
	$(OBJ_DIR)/tree.o \
	$(OBJ_DIR)/util.o \
	$(OBJ_DIR)/bit_vector.o\
	$(OBJ_DIR)/io-meg.o\
	$(OBJ_DIR)/io-gen-ests.o\
	$(OBJ_DIR)/refine-intron.o \
	$(OBJ_DIR)/refine.o \
	$(OBJ_DIR)/compute-alignments.o \
	$(OBJ_DIR)/exon-complexity.o \
	$(OBJ_DIR)/conversions.o \
	$(OBJ_DIR)/detect-polya.o \


stree_SOURCE= \
	$(STREE_DIR)/lst_string.c \
	$(STREE_DIR)/lst_stree.c \

stree_OBJ= \
	$(OBJ_DIR)/lst_string.o \
	$(OBJ_DIR)/lst_stree.o \


est_fact_SOURCE= \
	$(SRC_DIR)/classify-intron.c \
	$(SRC_DIR)/est-factorizations.c \
	$(SRC_DIR)/factorization-util.c \
	$(SRC_DIR)/factorization-refinement.c \
	$(SRC_DIR)/max-emb-graph.c \
	$(SRC_DIR)/aug_suffix_tree.c \
	$(SRC_DIR)/meg-simplification.c \
	$(SRC_DIR)/compute-est-fact.c \
	$(SRC_DIR)/main-est-fact.c

est_fact_OBJ= \
	$(OBJ_DIR)/classify-intron.o \
	$(OBJ_DIR)/est-factorizations.o \
	$(OBJ_DIR)/factorization-util.o \
	$(OBJ_DIR)/factorization-refinement.o \
	$(OBJ_DIR)/max-emb-graph.o \
	$(OBJ_DIR)/aug_suffix_tree.o \
	$(OBJ_DIR)/meg-simplification.o \
	$(OBJ_DIR)/compute-est-fact.o \
	$(OBJ_DIR)/main-est-fact.o

est_fact_PROG=$(BIN_DIR)/est-fact


test_data_create_SOURCE= \
	$(SRC_DIR)/test-data-create.c

test_data_create_OBJ= \
	$(OBJ_DIR)/test-data-create.o

test_data_create_PROG= \
	$(BIN_DIR)/test-data-create

min_factorization_SOURCE= \
	$(SRC_DIR)/min_factorization.c\
	$(SRC_DIR)/color_matrix.c\
	$(SRC_DIR)/simplify_matrix.c\
	$(SRC_DIR)/simpl_info.c \
	$(SRC_DIR)/main-min-factorization.c

min_factorization_OBJ= \
	$(OBJ_DIR)/min_factorization.o\
	$(OBJ_DIR)/color_matrix.o\
	$(OBJ_DIR)/simplify_matrix.o\
	$(OBJ_DIR)/simpl_info.o \
	$(OBJ_DIR)/main-min-factorization.o

min_factorization_PROG=\
	$(BIN_DIR)/min-factorization

intron_agreement_SOURCE= \
	$(SRC_DIR)/classify-intron.c \
	$(SRC_DIR)/agree-introns.c \
	$(SRC_DIR)/main-intron-agreement.c

intron_agreement_OBJ= \
	$(OBJ_DIR)/classify-intron.o \
	$(OBJ_DIR)/agree-introns.o \
	$(OBJ_DIR)/main-intron-agreement.o

intron_agreement_PROG=\
	$(BIN_DIR)/intron-agreement

max_transcr_SOURCE= \
	$(SRC_DIR)/log.c \
	$(SRC_DIR)/util.c \
	$(SRC_DIR)/my_time.c \
	$(SRC_DIR)/MaximalTranscripts.c

max_transcr_OBJ= \
	$(OBJ_DIR)/log.o \
	$(OBJ_DIR)/util.o \
	$(OBJ_DIR)/my_time.o \
	$(OBJ_DIR)/MaximalTranscripts.o

max_transcr_PROG= \
	$(BIN_DIR)/maximal-transcripts

# build_transcr_SOURCE= \
# 	$(SRC_DIR)/util.c \
# 	$(SRC_DIR)/my_time.c \
# 	$(SRC_DIR)/BuildTranscripts.c

# build_transcr_OBJ= \
# 	$(OBJ_DIR)/util.o \
# 	$(OBJ_DIR)/my_time.o \
# 	$(OBJ_DIR)/BuildTranscripts.o

# build_transcr_PROG= \
# 	$(BIN_DIR)/build-transcripts


cds_annotation_SOURCE= \
	$(SRC_DIR)/log.c \
	$(SRC_DIR)/util.c \
	$(SRC_DIR)/my_time.c \
	$(SRC_DIR)/CCDS.c

cds_annotation_OBJ= \
	$(OBJ_DIR)/log.o \
	$(OBJ_DIR)/util.o \
	$(OBJ_DIR)/my_time.o \
	$(OBJ_DIR)/CCDS.o

cds_annotation_PROG= \
	$(BIN_DIR)/cds-annotation


all_header_files:=$(wildcard ${INCLUDE_DIR}/*.h)

stree_header_files:=$(wildcard ${STREE_DIR}/*.h)

#####################
# System information
__DATETIME:=$(shell LANG=C date)
__HOST:=$(shell hostname)
have_git:=/$(shell which git)/
is_git_repository:=/$(wildcard .git)/
have_VERSION:=/$(wildcard VERSION)/
ifeq ($(is_git_repository), //)      # We are NOT in the original git repository

ifeq (/$(wildcard VERSION)/,//)          # We do NOT have the VERSION file
___SRC_DESC:=version-unknown
$(warning Impossible to get the source-code version.)
else                                     # We HAVE the VERSION file
___SRC_DESC:=$(shell cat VERSION)
endif # VERSION

else                                 # We ARE in the original git repository

ifeq ($(have_git), //)                   # We do NOT have git
___GIT_REF:=.git/$(shell cat .git/HEAD | cut -d ' ' -f 2)
ifeq ($(___GIT_REF),.git/)
___GIT_REF:= .git/refs/heads/master
endif
ifeq (/$(wildcard $(___GIT_REF))/,//)        # We CANNOT find the current head
ifeq (/$(wildcard VERSION)/,//)                  # We do NOT have the VERSION file
___SRC_DESC:=version-unknown
$(warning Impossible to get the source-code version.)
else                                             # We HAVE the VERSION file
___SRC_DESC:=$(shell cat VERSION)
endif # VERSION
else                                         # We FOUND the current head
___SRC_DESC:=commit-$(shell tail -n1 $(___GIT_REF))
endif # ___GIT_REF
else                                     # We HAVE git
___SRC_DESC:=$(shell git describe)
endif # have_git

endif # is_git_repository

# Sanitize
__C_SRC_DESC:=$(subst ",\\\",$(___SRC_DESC))
__COMPILER_VER:=$(shell $(CC) --version | head -n 1)
DSYSINFO=-D__BUILD_DATETIME="\"$(__DATETIME)\"" -D__BUILD_HOST="\"$(__HOST)\""
DSYSINFO+=-D__SRC_DESC="\"$(__C_SRC_DESC)\""
DSYSINFO+=-D__BUILD_DESC="\"compiler=$(CC) status=$(STATUS) cflags='$(ADD_CFLAGS)'\""
DSYSINFO+=-D__COMPILER_VER="\"$(__COMPILER_VER)\""
#####################

INNER_DIST_DIR= pintron-$(___SRC_DESC)-$(full_current_sys)
FULL_DIST_DIR= $(DIST_DIR)/$(INNER_DIST_DIR)
FULL_SRC_DIST_DIR= $(DIST_DIR)/pintron-$(___SRC_DESC)-src

all	: build

script_copy_to_bin = cp $(script) $(BASE_BIN_DIR)/$(notdir $(basename $(script))) &&

all	:
	@echo "${PHF}All compiled!${SF}"

stree   : .make $(stree_OBJ)

reall	: clean all
	@echo "${PHF}Cleaned and rebuilt all!${SF}"

clean 	:
	@echo "${PHF}Cleaning objects and programs${SF}" ; \
	rm -rf $(BASE_OBJ_DIR)/* $(BASE_BIN_DIR)/* \
		$(SRC_DIR)/options.c $(INCLUDE_DIR)/options.h \
		$(DIST_DIR)/*

build	: .make est-fact min-factorization intron-agreement max-transcr cds-annotation # build-transcr
	@$(foreach script,$(DIST_SCRIPTS),$(script_copy_to_bin)) \
	echo "Configuration: ${COMPFLAGS}"; \
	echo "Compiler:      ${CC}"; \
	echo "Symbols:       ${DFLAGS}"; \
	echo '${PHF}All built!${SF}'


$(SRC_DIR)/options.c: $(SRC_DIR)/options.ggo
	@echo "${PHF} * Generating the configuration parser...${SF}"; \
	gengetopt --output-dir=`dirname $@` \
		--file-name=`basename $@ .c` < $<; \
	mv $(SRC_DIR)/`basename $@ .c`.h $(INCLUDE_DIR)/

.make   :  Makefile
	@echo "${PHF}Makefile modified! Cleaning enforced.${SF}"; \
	rm -rf $(BASE_OBJ_DIR)/* $(BASE_BIN_DIR)/* \
	   $(SRC_DIR)/options.c $(INCLUDE_DIR)/options.h; \
        touch .make

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(all_header_files) $(stree_header_files)
	@echo '${PF} * Compiling ${SF}$<'; \
	mkdir -pv $(dir $@) ; \
	$(CC) $(ADD_CFLAGS)  $(DSYSINFO) $(CFLAGS_ARCH) -o $@ -c $<

#################################################################
$(OBJ_DIR)/%.o: $(STREE_DIR)/%.c $(stree_header_files)
	@echo '${PF} * Compiling ${SF}$<'; \
	mkdir -pv $(dir $@) ; \
	$(CC) $(ADD_CFLAGS) $(CFLAGS_ARCH) -o $@ -c $<
#################################################################



est-fact	: $(est_fact_PROG)
	@ln -f $(est_fact_PROG) $(BASE_BIN_DIR) ; \

reest-fact	: clean est-fact
	@echo '${PHF}Cleaned and rebuilt!${SF}'

$(est_fact_OBJ)	: $(stree_OBJ) $(base_OBJ) $(est_fact_SOURCE)

$(est_fact_PROG)	: $(stree_OBJ) $(base_OBJ) $(est_fact_OBJ)
	@echo '${PHF} * Linking${SF} $(notdir $@)'; \
	mkdir -pv $(BIN_DIR) ; \
	$(CC) -o $(est_fact_PROG) $(ADD_CFLAGS) $(LDFLAGS_ARCH) $(DSYSINFO) $^ $(LIBS); \
	echo '   ${PHF}...done.${SF}'; \



test-data-create	: $(test_data_create_PROG)
	@ln -f $(test_data_create_PROG) $(BASE_BIN_DIR)

retest-data-create	: clean test-data-create
	@echo '${PHF}Cleaned and rebuilt!${SF}'

$(test_data_create_OBJ)	: $(stree_OBJ) $(base_OBJ) $(test_data_create_SOURCE)

$(test_data_create_PROG)	: $(stree_OBJ) $(base_OBJ) $(test_data_create_OBJ)
	@echo '${PHF} * Linking${SF} $(notdir $@)'; \
	mkdir -pv $(BIN_DIR) ; \
	$(CC) -o $(test_data_create_PROG) $(ADD_CFLAGS) $(LDFLAGS_ARCH) $^ $(LIBS) ; \
	echo '   ${PHF}...done.${SF}'; \



min-factorization	: $(min_factorization_PROG)
	@ln -f $(min_factorization_PROG) $(BASE_BIN_DIR)

remin-factorization : clean min-factorization
	@echo '${PHF}Cleaned and rebuilt!${SF}'

$(min_factorization_OBJ)	: $(base_OBJ) $(min_factorization_SOURCE)

$(min_factorization_PROG)	: $(base_OBJ) $(min_factorization_OBJ)
	@echo '${PHF} * Linking${SF} $(notdir $@)'; \
	mkdir -pv $(BIN_DIR) ; \
	$(CC) -o $(min_factorization_PROG) $(ADD_CFLAGS) $(LDFLAGS_ARCH) $^ $(LIBS) ; \
	echo '   ${PHF}...done.${SF}'; \



intron-agreement	: $(intron_agreement_PROG)
	@ln -f $(intron_agreement_PROG) $(BASE_BIN_DIR)

reintron-agreement : clean intron-agreement
	@echo '${PHF}Cleaned and rebuilt!${SF}'

$(intron_agreement_OBJ)	: $(base_OBJ) $(intron_agreement_SOURCE)

$(intron_agreement_PROG)	: $(base_OBJ) $(intron_agreement_OBJ)
	@echo '${PHF} * Linking${SF} $(notdir $@)'; \
	mkdir -pv $(BIN_DIR) ; \
	$(CC) -o $(intron_agreement_PROG) $(ADD_CFLAGS) $(LDFLAGS_ARCH) $^ $(LIBS) ; \
	echo '   ${PHF}...done.${SF}'; \


max-transcr	: $(max_transcr_PROG)
	@ln -f $(max_transcr_PROG) $(BASE_BIN_DIR)

remax-transcr : clean max-transcr
	@echo '${PHF}Cleaned and rebuilt!${SF}'

$(max_transcr_OBJ)	: $(max_transcr_SOURCE)

$(max_transcr_PROG)	: $(max_transcr_OBJ)
	@echo '${PHF} * Linking${SF} $(notdir $@)'; \
	mkdir -pv $(BIN_DIR) ; \
	$(CC) -o $(max_transcr_PROG) $(ADD_CFLAGS) $(LDFLAGS_ARCH) $^ $(LIBS) ; \
	echo '   ${PHF}...done.${SF}'; \


# build-transcr	: $(build_transcr_PROG)
# 	@ln -f $(build_transcr_PROG) $(BASE_BIN_DIR)

# rebuild-transcr : clean build-transcr
# 	@echo '${PHF}Cleaned and rebuilt!${SF}'

# $(build_transcr_OBJ)	: $(build_transcr_SOURCE)

# $(build_transcr_PROG)	: $(build_transcr_OBJ)
# 	@echo '${PHF} * Linking${SF} $(notdir $@)'; \
# 	mkdir -pv $(BIN_DIR) ; \
# 	$(CC) -o $(build_transcr_PROG) $(ADD_CFLAGS) $(LDFLAGS_ARCH) $^ $(LIBS) ; \
# 	echo '   ${PHF}...done.${SF}'; \


cds-annotation	: $(cds_annotation_PROG)
	@ln -f $(cds_annotation_PROG) $(BASE_BIN_DIR)

recds-annotation : clean cds-annotation
	@echo '${PHF}Cleaned and rebuilt!${SF}'

$(cds_annotation_OBJ)	: $(cds_annotation_SOURCE)

$(cds_annotation_PROG)	: $(cds_annotation_OBJ)
	@echo '${PHF} * Linking${SF} $(notdir $@)'; \
	mkdir -pv $(BIN_DIR) ; \
	$(CC) -o $(cds_annotation_PROG) $(ADD_CFLAGS) $(LDFLAGS_ARCH) $^ $(LIBS) ; \
	echo '   ${PHF}...done.${SF}'; \


raw_script_copy = cp $(script) $(FULL_DIST_DIR)/bin/$(notdir $(basename $(script))) && sed -i -e "s/_____%PINTRON_VERSION%_____/ ${___SRC_DESC}/" $(FULL_DIST_DIR)/bin/$(notdir $(basename $(script))) &&

##
#
# Use PAR::Packer to deploy 'self-contained' perl scripts if possible
#
##

ifneq ($(USE_PAR),yes)
ifneq ($(USE_PAR),no)
override USE_PAR := ${DEFAULT_USE_PAR}
$(warning Using default value for PAR::Packer "${USE_PAR}")
endif
endif

ifeq ($(USE_PAR),no)

$(message Not using PAR for deploying Perl scripts.)
perl_script_copy = cp $(script) $(FULL_DIST_DIR)/bin/$(notdir $(basename $(script))) &&

else

# Check if there exist 'pp' script
# Do not ovveride: it can be specified from the command-line
PP_SCRIPT:=$(shell which pp)

ifeq (/$(shell test -x "${PP_SCRIPT}" && echo "OK")/,/OK/)

perl_script_copy = echo "Packaging Perl script $(script) with '${PP_SCRIPT}'..." && ${PP_SCRIPT} -c -P -o $(FULL_DIST_DIR)/bin/$(notdir $(basename $(script))) $(script) &&

else

$(error Script 'pp' for packaging Perl scripts not found. \
Please specify USE_PAR=no or provide the full pathname of 'pp' in PP_SCRIPT (e.g. PP_SCRIPT=/usr/local/bin/pp))

endif

endif # USE_PAR






script_copy = $(if $(findstring .pl,$(suffix $(script))),$(perl_script_copy),$(raw_script_copy))

prepare-dist	: all
	@echo '${PHF} * Preparing distribution...${SF}'; \
	rm -rf $(DIST_DIR)/pintron-*-$(full_current_sys)/ && \
	mkdir -p $(FULL_DIST_DIR) && \
	mkdir -p $(FULL_DIST_DIR)/bin && \
	mkdir -p $(FULL_DIST_DIR)/doc && \
	cp $(est_fact_PROG) $(FULL_DIST_DIR)/bin && \
	cp $(min_factorization_PROG) $(FULL_DIST_DIR)/bin && \
	cp $(intron_agreement_PROG) $(FULL_DIST_DIR)/bin && \
	cp $(max_transcr_PROG) $(FULL_DIST_DIR)/bin && \
	cp $(cds_annotation_PROG) $(FULL_DIST_DIR)/bin && \
	$(foreach script,$(DIST_SCRIPTS),$(script_copy)) \
	cp -r $(DIST_DOC_DIR)/* $(FULL_DIST_DIR)/doc && \
	( cd $(DIST_DIR) && rm -f pintron-latest && ln -s $(INNER_DIST_DIR) pintron-latest; ) && \
	echo '${PHF}   ...preparation done.${SF}' && \
	echo '${PHF}** The software package has been prepared for $(full_current_sys) in directory $(FULL_DIST_DIR)${SF}';

dist	: prepare-dist
	@echo '${PHF} * Compressing...${SF}'; \
	rm -rf $(DIST_DIR)/pintron-*-$(full_current_sys).tar.gz && \
	tar czf $(FULL_DIST_DIR).tar.gz -C $(DIST_DIR) `basename $(FULL_DIST_DIR)` && \
	echo '${PHF}   ...done.${SF}'; \
	echo '${PHF}** The software package has been prepared for $(full_current_sys) in file $(FULL_DIST_DIR).tar.gz${SF}';

install : prepare-dist


doc_symb_link = ln -s $(doc) $(FULL_SRC_DIST_DIR)/ &&

.PHONY: base-prepare-src-dist
base-prepare-src-dist	: clean
	@echo '${PHF} * Preparing source distribution...${SF}'; \
	rm -rf $(DIST_DIR)/pintron-*-src/ && \
	mkdir -p $(FULL_SRC_DIST_DIR) && \
	cp Makefile $(FULL_SRC_DIST_DIR) && \
	cp -r $(SRC_DIR) $(FULL_SRC_DIST_DIR) && \
	cp -r $(STREE_DIR) $(FULL_SRC_DIST_DIR) && \
	cp -r $(INCLUDE_DIR) $(FULL_SRC_DIST_DIR) && \
	cp -r $(DIST_SCRIPTS_DIR) $(FULL_SRC_DIST_DIR) && \
	cp -r $(DIST_DOC_DIR) $(FULL_SRC_DIST_DIR) && \
	$(foreach doc,$(wildcard $(DIST_DOC_DIR)/*),$(doc_symb_link)) \
	echo '${PHF}   ...preparation done.${SF}' && \
	echo '${PHF}** The source of the software package has been prepared in directory $(FULL_SRC_DIST_DIR)${SF}';

.PHONY: prepare-src-dist
prepare-src-dist	: base-prepare-src-dist
	@echo '${PHF} * Determining source version...${SF}'; \
	echo '$(___SRC_DESC)' > $(FULL_SRC_DIST_DIR)/VERSION ; \
	echo '${PHF}   ...done.${SF}'; \

.PHONY: src-dist
src-dist	: prepare-src-dist
	@echo '${PHF} * Compressing...${SF}'; \
	rm -rf $(FULL_SRC_DIST_DIR).tar.gz && \
	tar czf $(FULL_SRC_DIST_DIR).tar.gz -C $(DIST_DIR) `basename $(FULL_SRC_DIST_DIR)` && \
	echo '${PHF}   ...done.${SF}'; \
	echo '${PHF}** The source of software package has been prepared in file $(FULL_SRC_DIST_DIR).tar.gz${SF}';

check-syntax:
	    gcc -o /dev/null -D$(DEFINE) $(CFLAGS) -S ${CHK_SOURCES}



#Makefile part about tests
test_SOURCE= \
	$(CURDIR)/test/aug_suffix_tree_test.c\
	$(CURDIR)/test/bit_vector_test.c\
	$(CURDIR)/test/bool_list_test.c\
	$(CURDIR)/test/BuildTranscripts_test.c\
	$(CURDIR)/test/conversions_test.c\
	$(CURDIR)/test/double_list_test.c\
	$(CURDIR)/test/exon-complexity_test.c\
	$(CURDIR)/test/ext_array_test.c\
	$(CURDIR)/test/int_list_test.c\
	$(CURDIR)/test/io-multifasta_test.c\
	$(CURDIR)/test/list_test.c\
	$(CURDIR)/test/min_factorization_test.c\
	$(CURDIR)/test/refine-intron_test.c\
	$(CURDIR)/test/simpl_info_test.c\
	$(CURDIR)/test/types_test.c\
	$(CURDIR)/test/util_test.c\
	$(CURDIR)/regressionTest/testPIntronOutput.c\

test_EXEC= \
	$(CURDIR)/test/aug_suffix_tree_test\
	$(CURDIR)/test/bit_vector_test\
	$(CURDIR)/test/bool_list_test\
	$(CURDIR)/test/BuildTranscripts_test\
	$(CURDIR)/test/conversions_test\
	$(CURDIR)/test/double_list_test\
	$(CURDIR)/test/exon-complexity_test\
	$(CURDIR)/test/ext_array_test\
	$(CURDIR)/test/int_list_test\
	$(CURDIR)/test/io-multifasta_test\
	$(CURDIR)/test/list_test\
	$(CURDIR)/test/min_factorization_test\
	$(CURDIR)/test/refine-intron_test\
	$(CURDIR)/test/simpl_info_test\
	$(CURDIR)/test/types_test\
	$(CURDIR)/test/util_test\
	$(CURDIR)/regressionTest/executePIntronTests\
	$(CURDIR)/regressionTest/testPIntronOutput\

CFLAGS = -l criterion -I $(INCLUDE_DIR) -I $(STREE_DIR) -I regressionTest

comp-test:$(test_EXEC) $(all_header_files) $(stree_header_files)
	$(CC) $(test_SOURCE) $(CFLAGS)
	
test: $(test_EXEC)
	#preparing files for regressionTest:
	$(CURDIR)/regressionTest/executePIntronTests
	#executing unitTest:
	$(CURDIR)/test/aug_suffix_tree_test
	$(CURDIR)/test/bit_vector_test
	$(CURDIR)/test/bool_list_test
	$(CURDIR)/test/BuildTranscripts_test
	$(CURDIR)/test/conversions_test
	$(CURDIR)/test/double_list_test
	$(CURDIR)/test/exon-complexity_test
	$(CURDIR)/test/ext_array_test
	$(CURDIR)/test/int_list_test
	$(CURDIR)/test/io-multifasta_test
	$(CURDIR)/test/list_test
	$(CURDIR)/test/min_factorization_test
	$(CURDIR)/test/refine-intron_test
	$(CURDIR)/test/simpl_info_test
	$(CURDIR)/test/types_test
	$(CURDIR)/test/util_test
	#executing regressionTest:
	$(CURDIR)/regressionTest/testPIntronOutput
.PHONY: clean-test
clean-test:
	#cleaning unitTest temp files:
	rm $(CURDIR)/test/*_test
	#cleaning PIntron temp files:
	rm $(CURDIR)/ests.txt
	rm $(CURDIR)/genomic.txt
	rm $(CURDIR)/pintron-log.txt
	#cleaning regressionTest temp files:
	rm $(CURDIR)/regressionTest/executePIntronTests
	rm $(CURDIR)/regressionTest/testPIntronOutput
	rm $(CURDIR)/regressionTest/output.txt
	rm -rf $(CURDIR)/regressionTest/test*/executionOutput/
