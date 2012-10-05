#!/usr/bin/env python3
####
#
#
#                              PIntron
#
# A novel pipeline for computational gene-structure prediction based on
# spliced alignment of expressed sequences (ESTs and mRNAs).
#
# Copyright (C) 2010  Gianluca Della Vedova, Yuri Pirola
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


import random
import sys
import re
import os
import subprocess
import time
import logging
import json
import pprint
import traceback
import csv
import hashlib

from optparse import OptionParser


def md5Checksum(filePath):
    fh = open(filePath, 'rb')
    m = hashlib.md5()
    while True:
        data = fh.read(8192)
        if not data:
            break
        m.update(data)
    return m.hexdigest()


class PIntronError(Exception):
    """Base class for exceptions of the PIntron pipeline."""
    pass


class PIntronIOError(PIntronError):
    """Exception raised for errors related to I/O operations.

    Attributes:
        e_file -- file on which the error occurred
        msg    -- explanation of the error
    """

    def __init__(self, e_file, msg):
        self.e_file = e_file
        self.msg = msg

    def __str__(self):
        return '{0} (offending file: "{1}")'.format(self.msg, self.e_file)


def parse_command_line():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    # parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
    #                   default=False, help="print status messages to stdout")
    parser.add_option("-g", "--genomic",
                      dest="genome_filename",
                      default="genomic.txt",
                      help="FILE containing the genomic sequence",
                      metavar="GENOMIC_FILE")
    parser.add_option("-s", "--EST",
                      dest="EST_filename", default="ests.txt",
                      help="FILE containing the ESTs", metavar="ESTs_FILE")
    # parser.add_option("-1", "--no-full-lengths", action="store_true",
    #                   dest="step1", default=False,
    #                   help="do not compute the full-lengths")
    parser.add_option("-o", "--output",
                      dest="output_filename",
                      default="pintron-full-output.json",
                      help="full output file (default = '%default')",
                      metavar="FILE")
    parser.add_option("-z", "--compress", action="store_true",
                      dest="compress", default=False,
                      help="compress output (default = %default)")
    parser.add_option("-a", "--alignments", action="store_true",
                      dest="output_est_alignments", default=False,
                      help="output ESTs/genome alignments (default = %default)")
    parser.add_option("-l", "--logfile",
                      dest="plogfile", default="pintron-pipeline-log.txt",
                      help="log filename of the pipeline steps (default = '%default')",
                      metavar="FILE")
    parser.add_option("--general-logfile",
                      dest="glogfile", default="pintron-log.txt",
                      help="log filename of the pipline orchestration module (default = '%default')",
                      metavar="FILE")
    parser.add_option("-c", "--continue", action="store_true",
                      dest="from_scratch", default=False,
                      help="resume a previosly interrupted computation (default = %default)")
    parser.add_option("-b", "--bin-dir",
                      dest="bindir", default="",
                      help="DIRECTORY containing the programs (default = system PATH)")
    parser.add_option("-n", "--organism",
                      dest="organism", default="unknown",
                      help="Organism originating the ESTs (default = '%default')")
    parser.add_option("-e", "--gene",
                      dest="gene", default="unknown",
                      help="Gene symbol (or ID) of the locus which the ESTs refer to (default = '%default')")
    parser.add_option("-k", "--keep-intermediate-files", action="store_true",
                      dest="no_clean", default=False,
                      help="keep all intermediate or temporary files (default = %default)")
    parser.add_option("-t", "--gtf",
                      dest="gtf_filename",
                      default="pintron-cds-annotated-isoforms.gtf",
                      help="output GTF FILE with the isoforms that have a CDS annotation "
                      "(default = '%default')",
                      metavar="GTF_FILE")
    parser.add_option("--extended-gtf",
                      dest="extended_gtf_filename",
                      default="pintron-all-isoforms.gtf",
                      help="output GTF FILE with all the predicted isoforms "
                      "(default = '%default')",
                      metavar="GTF_FILE")

    # parser.add_option("--strand",
    #                   dest="strand", type="int", default=1,
    #                   help="[Expert use only] print status messages to stdout")
    # parser.add_option("--chromosome",
    #                   dest="chromosome", default=False,
    #                   help="[Expert use only] print status messages to stdout")
    # parser.add_option("--EST-cluster",
    #                   dEST="EST_cluster", default=False,
    #                   help="[Expert use only] print status messages to stdout")
    # parser.add_option("--min-factor-length",
    #                   dest="min_factor_length", type="int", default=15,
    #                   help="[Expert use only] minimum factor length")
    # parser.add_option("--min-intron-length",
    #                   dest="min_intron_length", type="int", default=60,
    #                   help="[Expert use only] minimum intron length")
    # parser.add_option("--max-intron-length",
    #                   dest="max_intron_length", type="int", default=0,
    #                   help="[Expert use only] max intron length")
    # parser.add_option("--min-string-depth-rate",
    #                   dest="min_string_depth_rate", type="float", default=0.2,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--max-prefix-discarded-rate",
    #                   dest="max_prefix_discarded_rate", type="float", default=0.6,
    #                   help="[Expert use only] largest prefix of an EST that might be discarded (expressed as a fraction of the EST length)")
    # parser.add_option("--max-suffix-discarded-rate",
    #                   dest="max_suffix_discarded_rate", type="float", default=0.6,
    #                   help="[Expert use only] largest suffix of an EST that might be discarded (expressed as a fraction of the EST length)")
    # parser.add_option("--max-prefix-discarded",
    #                   dest="max_prefix_discarded", type="int", default=50,
    #                   help="[Expert use only] largest prefix of an EST that might be discarded during factorization construction (expressed as number of bases)")
    # parser.add_option("--max-suffix-discarded",
    #                   dest="max_suffix_discarded", type="int", default=50,
    #                   help="[Expert use only] largest suffix of an EST that might be discarded during factorization construction (expressed as number of bases)")
    # parser.add_option("--min-distance-of-splice-sites",
    #                   dest="min_distance_splice_sites", type="int", default=50,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--max-no-of-factorizations",
    #                   dest="max_factorizations", type="int", default=0,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--max-difference-of-coverage",
    #                   dest="max_difference_coverage", type="float", default=0.05,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--max-difference-of-no-of-exons",
    #                   dest="max_difference_no_exons", type="int", default=5,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--max-difference-of-gap-length",
    #                   dest="max_difference_gap_length", type="int", default=20,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--retain-externals",
    #                   dest="retain_externals", default=True, action="store_true",
    #                   help="[Expert use only] TODO")
    # parser.add_option("--max-pairings-in-CMEG",
    #                   dest="max_pairings_CMEG", type="int", default=80,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--max-shortest-pairing-frequence",
    #                   dest="max_shortest_pairing_frequence", type="float", default=0.4,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--suffix-prefix-length-intron",
    #                   dest="suffix_prefix_length_intron", type="int", default=70,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--suffix-prefix-length-EST",
    #                   dest="suffix_prefix_length_EST", type="int", default=30,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--suffix-prefix-length-genomic",
    #                   dest="suffix_prefix_length_genomic", type="int", default=30,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--no-transitive-reduction",
    #                   dest="no_transitive_reduction", default=False, action="store_true",
    #                   help="[Expert use only] TODO")
    # parser.add_option("--no-short-edge",
    #                   dest="no_short_edge", default=False, action="store_true",
    #                   help="[Expert use only] TODO")
    parser.add_option("--set-max-factorization-time",
                      dest="max_factorization_time", type="int", default=60,
                      help="[Expert use only] Set a time limit (in mins) for the factorization step")
    parser.add_option("--set-max-factorization-memory",
                      dest="max_factorization_memory", type="int", default=3000,
                      help="[Expert use only] Set a limit (in MiB) for the memory used by the factorization step"
                      " (default = 3000 MiB, approx. 3GB)")
    parser.add_option("--set-max-exon-agreement-time",
                      dest="max_exon_agreement_time", type="int", default=15,
                      help="[Expert use only] Set a time limit (in mins) for the exon agreement step")
    parser.add_option("--set-max-intron-agreement-time",
                      dest="max_intron_agreement_time", type="int", default=30,
                      help="[Expert use only] Set a time limit (in mins) for the intron agreement step")
    parser.add_option("--pas-tolerance",
                      dest="pas_tolerance", type="int", default=30,
                      help="[Expert use only] Maximum allowed difference on the exon final coordinate to identify a PAS")

    (options, args) = parser.parse_args()
    if options.bindir:
        options.bindir = os.path.normpath(options.bindir)

    return(options)


# Transform a JSON file into a GTF
def json2gtf(infile, outfile, gene_name, all_isoforms):
    def write_gtf_line(file, seqname, feature, start, end, score, strand, frame, gene, transcript):
        if end < start:
            start, end = end, start
        f.write("\t".join([seqname, "PIntron", feature, str(start), str(end), score, strand, str(frame),
                           "gene_id \"{0}\"; transcript_id \"{0}.{1}\";\n".format(gene, str(transcript))]))

    logging.debug(str(time.localtime()))
    logging.debug(json2gtf)
    with open(infile, 'r', encoding='utf-8') as f:
        entry = json.load(f)

    with open(outfile, 'w', encoding='utf-8') as f:
        for isoform_id, isoform in entry["isoforms"].items():
            for exon in isoform["exons"]:

                logging.debug("Exon data (json): { %s } ", exon)
#                import pdb; pdb.set_trace()
                if all_isoforms or isoform["annotated CDS?"]:
                    write_gtf_line(f, entry['genome']['sequence_id'], "exon",
                                   exon['chromosome start'], exon['chromosome end'],
                                   "0", entry['genome']['strand'], ".", gene_name, isoform_id)
                    if 'genome 5UTR start' in exon:
                        write_gtf_line(f, entry['genome']['sequence_id'], "5UTR",
                                       exon['genome 5UTR start'], exon['genome 5UTR end'],
                                       "0", entry['genome']['strand'], ".", gene_name, isoform_id)
                    if 'genome start codon start' in exon:
                        write_gtf_line(f, entry['genome']['sequence_id'], "start_codon",
                                       exon['genome start codon start'], exon['genome start codon end'],
                                       "0", entry['genome']['strand'], exon['frame start codon'], gene_name, isoform_id)
                    if 'genome CDS start' in exon:
                        write_gtf_line(f, entry['genome']['sequence_id'], "CDS",
                                       exon['genome CDS start'], exon['genome CDS end'],
                                       "0", entry['genome']['strand'], exon['frame CDS'], gene_name, isoform_id)
                    if 'genome stop codon start' in exon:
                        write_gtf_line(f, entry['genome']['sequence_id'], "stop_codon",
                                       exon['genome stop codon start'], exon['genome stop codon end'],
                                       "0", entry['genome']['strand'], exon['frame stop codon'], gene_name, isoform_id)
                    if 'genome 3UTR start' in exon:
                        write_gtf_line(f, entry['genome']['sequence_id'], "3UTR",
                                       exon['genome 3UTR start'], exon['genome 3UTR end'],
                                       "0", entry['genome']['strand'], ".", gene_name, isoform_id)


def compute_json(ccds_file, variant_file, output_file, from_scratch, pas_tolerance, genomic_seq):
    def dump_and_exit(exon, isoform, isoform_id):
        logging.debug("Exon =>")
        logging.debug(exon)
        logging.debug("Isoform (ID " + str(isoform_id) + ")=>")
        logging.debug(isoform)
        raise PIntronError

    # Find the sequence ID
    # It is stored in the first line of the genomic sequence
    with open(genomic_seq, 'r', encoding='utf-8') as f:
        line = f.readline().rstrip("\r\n")
        m = re.search('^>[^\d]*(\d+):.*:([+-]?\d+)$', line)
        sequence_id = line[1:]
        strand = m.group(2)
        if strand == '-1' or strand == '-':
            strand = '-'
        else:
            strand = '+'

    gene = {
        'version': 4,  # Hardcoding version number
        'program_version': options.version,  # Program version
        'isoforms': {},
        'introns': {},
        'factorizations': {},
        'number of processed ESTs': 0,
        'genome': {
            'sequence_id': sequence_id,
            'strand': strand,
        },
    }

    for file in [ccds_file, variant_file]:
        if not os.access(file, os.R_OK):
            # throw exception and die
            logging.exception("*** Fatal error: Could not read " + file + "\n")

    with open('out-after-intron-agree.txt', mode='r', encoding='utf-8') as fd:
        current = ''
        for line in fd:
            l = line.rstrip()
            if l[0] == '>':
                gene['number of processed ESTs'] = gene['number of processed ESTs'] + 1
                new = re.search('\/gb=([A-Z_0-9]+)', l).groups()
                current = new[0]
                gene['factorizations'][current] = {
                    'polyA?': False,
                    'PAS': False,
                    'exons': [],
                    'EST': current,
                }
                if re.search('\/clone_end=([35])', l):
                    new = re.search('\/clone_end=([35])', l).groups()
                    gene['factorizations'][current]['clone end'] = new[0]

            elif re.match('#polya=1', l):
                gene['factorizations'][current]['polyA?'] = True
            elif re.match('#polyad(\S*)=1', l):
                gene['factorizations'][current]['PAS'] = True
            elif re.match('(\d+) (\d+) (\d+) (\d+)( \S+)? \S+$', l):
                new = re.match('(\d+) (\d+) (\d+) (\d+) (\S+) (\S+)$', l).groups()
                # pprint.pprint(l)
                # pprint.pprint(new)
                exon = {
                    'EST start': int(new[0]),
                    'EST end': int(new[1]),
                    'relative start': int(new[2]),
                    'relative end': int(new[3]),
                    'EST sequence': new[4],
                    'genome sequence': new[5],
                }
                gene['factorizations'][current]['exons'].append(exon)
                if gene['factorizations'][current]['PAS']:
                    gene['factorizations'][current]['exon'] = exon

    with open(variant_file, mode='r', encoding='utf-8') as fd:
        for line in fd:
            row = re.split(' /', line.rstrip())
            index = int(re.sub('^.*\#', '', row.pop(0)))
            isoform = {
                'exons': [],
                'polyA?': False,
                'PAS?': False,
                'annotated CDS?': False,
                'Reference frame?': False,
            }
            for t in row:
                (k, v) = re.split('=', t, 2)
                logging.debug("Reading VariantGTF: " + k + "=>" + v + "!")
                if k == "nex":
                    isoform['number exons'] = int(v)
                elif k == "L":
                    isoform["length"] = int(v)
                elif k == "CDS":
                    if v != '..':
                        isoform["annotated CDS?"] = True
                        m = re.match('^(<?)(\d+)\.\.(\d+)(>?)$', v)
                        (a, isoform["CDS start"], isoform["CDS end"], b) = (m.group(1), int(m.group(2)),
                                                                            int(m.group(3)), m.group(4))
                        isoform['CDS length'] = isoform["CDS end"] - isoform["CDS start"] + 1
                        isoform['canonical start codon?'] = False if a == '<' else True
                        isoform['canonical end codon?']   = False if b == '>' else True
                elif k == "RefSeq":
                    m = re.match('^(.*?)(\(?([NY])([NY])\)?)?$', v, flags=re.IGNORECASE)
                    if m:
                        (r, a, b) = (m.group(1), m.group(3), m.group(4))
                        isoform['reference CDS start codon?'] = False if a == 'N' else True
                        isoform['reference CDS end codon?']   = False if b == 'N' else True
                        if r != None and r:
                            isoform['RefSeq'] = r
                elif k == "ProtL":
                    if v != '..' and isoform["annotated CDS?"]:
                        m = re.match('^(>?)(\d+)$', v, flags=re.IGNORECASE)
                        (a, isoform['protein length']) = (m.group(1), int(m.group(2)))
                        isoform['protein missing codon?'] = False if a != '>' else True
                elif k == "Frame":
                    m = re.match('^y', v, flags=re.IGNORECASE)
                    if m != None and isoform["annotated CDS?"]:
                        isoform['Reference frame?'] = True
                elif k == "Type":
                    if v == 'Ref':
                        isoform['Reference frame?'] = True
                        if 'RefSeq' in isoform:
                            isoform['Type'] = isoform['RefSeq'] + " (Reference TR)"
                        else:
                            isoform['Type'] = "(Reference TR)"
                    else:
                        isoform['Type'] = re.sub('\s+$', '', v)
                elif not re.match('^\s*\#', line):
                    raise ValueError("Could not parse GTF file " + variant_file + "(" + k + "=>" + v + ")\n" + line + "\n")
            gene['isoforms'][index] = isoform

    with open(ccds_file, mode='r', encoding='utf-8') as fd:
        gene['number_isoforms'] = int(fd.readline().rstrip())
        gene['length_genomic_sequence'] = int(fd.readline().rstrip())
        isoform_5utr_length = 0
        isoform_3utr_length = 0
        for line in fd:
            l = line.rstrip()
            l = re.sub('\s+', '', l)
            l = re.sub('#.*', '', l)

            if re.match('^>', l):
                # New isoform
                l = l[1:]
                # print(l)
                fields = [int(x) for x in  re.split(':', l)]

                # pprint.pprint(fields)
                index = fields[0]
                if not index in gene['isoforms']:
                    raise ValueError("CCDS file " + ccds_file + "contains isoform with index " + index + " not in variant_file\n")
                if fields[1] > gene['isoforms'][index]['number exons']:
#                    import pdb; pdb.set_trace()
                    raise ValueError("Wrong number of exons: " + str(index) + "\n " + str(fields[1]) + "!= " +
                                     str(isoform['number exons']) + "\n")

                gene['isoforms'][index]['reference?'] = False if fields[2] == 0 else True
                gene['isoforms'][index]['from RefSeq?'] = False if fields[3] == 0 else True
                gene['isoforms'][index]['NMD flag'] = fields[4]
                if gene['isoforms'][index]['annotated CDS?']:
                    isoform_5utr_length = gene['isoforms'][index]['CDS start']
                    isoform_3utr_length = gene['isoforms'][index]['CDS end']
            elif re.match('^(\d+:){5}(-?\d+:)(-?\d+)$', l):
                # Row contains exon metadata
                exon = {}
                (exon["chromosome start"], exon["chromosome end"], exon["relative start"], exon["relative end"],
                 polyA, exon["5utr length"], exon["3utr length"]) = [max(0, int(x)) for x in  re.split(':', l)]
                exon['genome length'] = abs(exon["chromosome end"] - exon["chromosome start"]) + 1
                # if gene['genome']['strand'] == '-':
                #     (exon["chromosome start"], exon["chromosome end"], exon["relative start"], exon["relative end"]) = (exon["chromosome end"], exon["chromosome start"], exon["relative end"], exon["relative end"])
                if (polyA == 1):
                    gene['isoforms'][index]['polyA?'] = True
                # pprint.pprint(exon)
                logging.debug("Reading CCDS_transcripts: Row contains exon metadata { %s }",
                              "; ".join([line.rstrip(),
                                         str(max(exon["relative end"], exon["relative start"])),
                                         str(min(exon["relative end"], exon["relative start"])),
                                         str(exon["5utr length"]),
                                         str(exon["3utr length"]),
                                         str(abs(exon["relative end"] - exon["relative start"]) + 1 - exon["5utr length"] - exon["3utr length"])
                                     ]))
                if int(re.split(':', l)[4]) < 0:
                    del(exon["5utr length"])
                if int(re.split(':', l)[5]) < 0:
                    del(exon["3utr length"])
                gene['isoforms'][index]['exons'].append(exon)

            elif re.match('^[acgtACGT]+$', l):
                last_exon = gene['isoforms'][index]['exons'][-1]['sequence'] = l
                gene['isoforms'][index]['exons'][-1]['transcript length'] = len(l)
            elif not re.match('^\s*\#', line):
                raise ValueError("Could not parse CCDS file " + ccds_file + " at line:\n" + line + "\n")

    # When the strand is negative, the exons are in reverse order
    for isoform in gene['isoforms'].keys():
        gene['isoforms'][isoform]['exons'].reverse()

    with open('predicted-introns.txt', mode='r', encoding='utf-8') as fd:
        index = 1
        for line in fd:
            intron = {"supporting ESTs": {}}
            (intron['relative start'], intron['relative end'],
             intron['chromosome start'], intron['chromosome end'], intron['length'], intron['number supporting EST'], EST_list,
             intron['donor alignment average error'], intron['acceptor alignment average error'], intron['donor score'],
             intron['acceptor score'], intron['BPS score'], intron['BPS position'], intron['type'], intron['pattern'],
             intron['repeat sequence'], intron['donor suffix'], intron['prefix'], intron['suffix'],
             intron['acceptor prefix']) = re.split("\t", line.rstrip())
             # intron['begin donor'] = intron['relative end'] - len(intron['donor suffix']) +1
             # intron['end acceptor'] = intron['relative begin'] - len(intron['acceptor prefix']) -1
            intron["supporting ESTs"] = {i: {} for i in re.split(',', EST_list) if i != ''}
            #pprint.pprint(intron["supporting ESTs"])
            #import pdb; pdb.set_trace()

            for field in ('relative start', 'relative end', 'chromosome start', 'chromosome end', 'length', 'number supporting EST',
                          'BPS position'):
                intron[field] = int(intron[field])
            for field in ('donor alignment average error', 'acceptor alignment average error', 'donor score',
                          'acceptor score', 'BPS score'):
                intron[field] = float(intron[field])

            if intron['BPS position'] < 0:
                del intron['BPS position']

            gene['introns'][index] = intron
            index += 1

    # add introns to each isoform
    for isoform in gene['isoforms'].values():
        isoform['exons'].sort(key=lambda x: x['relative end'])
        isoform['introns'] = []
        pairs = zip(isoform['exons'][1:], isoform['exons'][:-1])
        for pair in pairs:
            list_extremes = [pair[0]['chromosome end'], pair[0]['chromosome start'],
                             pair[1]['chromosome end'], pair[1]['chromosome start']]
            list_extremes.sort()
            left_border = list_extremes[1] + 1
            right_border = list_extremes[2] - 1
            for index in gene['introns'].keys():
                intron = gene['introns'][index]
                if intron['chromosome start'] == left_border and intron['chromosome end'] == right_border or intron['chromosome end'] == left_border and intron['chromosome start'] == right_border:
                    isoform['introns'].append(index)

    # for each intron, add the alignment of the sorrounding exons.
    # Since different factorizations can support the same intron, the first
    # step is to find all pairs of exons supporting an intron
    def supporting_factors(intron):
        pairs = []
        for est in intron["supporting ESTs"].keys():
            factor = gene['factorizations'][est]
            #                    import pdb; pdb.set_trace()
            good_left  = [exon for exon in factor['exons'] if exon['relative end'] == intron['relative start'] - 1]
            good_right = [exon for exon in factor['exons'] if exon['relative start'] == intron['relative end'] + 1]
            if len(good_left) == 1 and len(good_right) == 1:
                pairs.append([est, good_left[0], good_right[0]])
        if len(pairs) != intron['number supporting EST']:
            pprint.pprint(factor)
            print("\n")
            pprint.pprint(intron)
            print("\n")
            pprint.pprint(pairs)
            raise PIntronError
        return(pairs)

    #
    # Each intron has the list of supporting ESTs.
    # For each such EST we provide the suffix/prefix of the prev/next exon
    for index in gene['introns'].keys():
            # donor_exon = [ exon for exon in isoform['exons'] if (exon['relative end'] == intron['relative start'] - 1) ][0]
            # acceptor_exon = [ exon for exon in isoform['exons'] if (exon['relative start'] == intron['relative end'] + 1) ][0]
            # add the alignment to each intron
        for [est, donor_factor, acceptor_factor] in supporting_factors(gene['introns'][index]):
#            import pdb; pdb.set_trace()
            gene['introns'][index]['supporting ESTs'][est] = {
                'EST donor factor suffix': donor_factor['EST sequence'][-len(gene['introns'][index]['donor suffix']):],
                'EST acceptor factor prefix': acceptor_factor['EST sequence'][:len(gene['introns'][index]['acceptor prefix'])],
                'begin EST acceptor factor': acceptor_factor['EST start'],
                'end EST donor factor': donor_factor['EST end'],
                'end EST acceptor factor': acceptor_factor['EST end'],
                'begin EST donor factor': donor_factor['EST start'],
                # 'EST prefix previous exon'  : gene['introns'][index]['acceptor prefix'],
                # 'EST suffix next exon'  : gene['introns'][index]['donor suffix'],
            # 'EST prefix end'   : acceptor_exon['EST end'],
            # 'EST suffix start' : donor_exon['EST start'],
            }

    def same_coordinates(a, b):
        return True if (a['relative start'] == b['relative start'] and
                        30 >= a['relative end'] - b['relative end'] >= -30) else False

    # pprint.pprint(gene)
    for isoform in gene['isoforms'].keys():
        gene['isoforms'][isoform]['transcript sequence'] = ''.join([s['sequence'] for s in gene['isoforms'][isoform]['exons']])
        if gene['isoforms'][isoform]['annotated CDS?']:
            # Check start/stop codon
            allowed_codons = {'first': ["ATG"],
                              'last': ["TGA", "TAG", "TAA"]
                          }
            codon = {
                'first': gene['isoforms'][isoform]['transcript sequence'][gene['isoforms'][isoform]['CDS start'] - 1: gene['isoforms'][isoform]['CDS start'] + 2],
                'last': gene['isoforms'][isoform]['transcript sequence'][gene['isoforms'][isoform]['CDS end'] - 3:gene['isoforms'][isoform]['CDS end']]
            }
            for p in ['first', 'last']:
                logging.debug("Codon  %s in  [%s]", codon[p].upper(), "/".join(allowed_codons[p]))
                if not (codon[p].upper() in allowed_codons[p]):
                    print("Warning JSON: wrong delimiter. Found " + codon[p].upper() + " instead of " +
                          "/".join(allowed_codons[p]) + " as " + p + " codon")
                    print("Isoform:")
                    pprint.pprint(isoform)
        # Check if we have to add PAS
        if not gene['isoforms'][isoform]['polyA?']:
            continue
        exon = gene['isoforms'][isoform]['exons'][-1]
        # If PAS_factorizations has an exon with the same coordinates,
        # we have a PAS
        if any(x for x in gene['factorizations'].values() if x['PAS'] and same_coordinates(x['exon'], exon)):
            gene['isoforms'][isoform]['PAS?'] = True

    # Enrich the JSON file with information that can be used to compute the GTF file
    def check_codon(codon_type, codon_string):
        possible_values = ["ATG"] if codon_type == 'start' else ["TGA", "TAG", "TAA"]
        if not codon_string.upper() in possible_values:
            logging.debug("Warning: wrong " + codon_type + " delimiter. Found " + codon_string + " instead of " + "/".join(possible_values))
            return True
        else:
            return False
    strand = gene['genome']['strand']
    sequence_id = re.sub(':.*', '', gene['genome']['sequence_id'])
    # The genomic sequence length stored in the JSON file
    # cannot be trusted.
    # seq_record=next(SeqIO.parse(genomic_seq, "fasta"))
    # gene['length_genomic_sequence']=len(seq_record)
    # print(gene['length_genomic_sequence'])

    data_strand = {'first': {"label": "5UTR",
                             "codons": ["ATG"],
                         },
                   'last': {"label": "3UTR",
                            "codons": ["TGA", "TAG", "TAA"],
                        }
               }
    for isoform_id, isoform in gene["isoforms"].items():
        if not isoform["annotated CDS?"]:
            continue
        cumulative_genome_length = 0
        cumulative_transcript_length = 0
        read_start_codon_seq = ''
        read_stop_codon_seq = ''
        ordered_codons = ["start", "stop"] if strand == '+' else  ["stop", "start"]
        for exon in isoform["exons"]:
            cumulative_genome_length_old = cumulative_genome_length
            cumulative_transcript_length_old = cumulative_transcript_length
            cumulative_genome_length += exon['genome length']
            exon['cumulative genome length'] = cumulative_genome_length
            cumulative_transcript_length += exon['transcript length']
            exon['cumulative transcript length'] = cumulative_transcript_length
            if cumulative_transcript_length < isoform['CDS start'] - 1:
                # exon is contained in 5UTR
                if strand == '+':
                    exon['genome 5UTR start'], exon['genome 5UTR end'] = exon['chromosome start'], exon['chromosome end']
                else:
                    exon['genome 5UTR start'], exon['genome 5UTR end'] = exon['chromosome end'], exon['chromosome start']
                continue
            if cumulative_transcript_length_old > isoform['CDS end'] + 1:
                # exon is contained in 3UTR
                if strand == '+':
                    exon['genome 3UTR start'], exon['genome 3UTR end'] = exon['chromosome start'], exon['chromosome end']
                else:
                    exon['genome 3UTR start'], exon['genome 3UTR end'] = exon['chromosome end'], exon['chromosome start']
                continue
            if cumulative_transcript_length_old + 1 <= isoform['CDS start'] - 1 <= cumulative_transcript_length:
                # exon contains a 5UTR portion
                if strand == '+':
                    exon['genome 5UTR start'] = exon['chromosome start']
                    exon['genome 5UTR end'] = exon['chromosome start'] + (exon['5utr length'] - 1)
                else:
                    exon['genome 5UTR start'] = exon['chromosome end']
                    exon['genome 5UTR end'] = exon['chromosome end'] - (exon['5utr length'] - 1)
            if cumulative_transcript_length_old + 1 <= isoform['CDS end'] + 1 <= cumulative_transcript_length:
                # exon contains a 3UTR portion
                if strand == '+':
                    exon['genome 3UTR start'] = exon['chromosome end'] - (exon['3utr length'] - 1)
                    exon['genome 3UTR end'] = exon['chromosome end']
                else:
                    exon['genome 3UTR start'] = exon['chromosome start']
                    exon['genome 3UTR end'] = exon['chromosome start'] + (exon['3utr length'] - 1)

            read_codon_len = 0
            if cumulative_transcript_length_old < isoform['CDS start'] <= cumulative_transcript_length:
                # exon contains at least part of the first codon, including the first character
                read_codon_len = min(3, cumulative_transcript_length - isoform['CDS start'] + 1)
                pos = isoform['CDS start'] - cumulative_transcript_length_old - 1
                read_start_codon_seq += exon['sequence'][pos:pos + read_codon_len]
                if len(read_start_codon_seq) == 3 and check_codon("start", read_start_codon_seq):
                    logging.debug("pos = " + str(pos))
                    logging.debug(read_start_codon_seq)
                    logging.debug("read_codon_len = " + str(read_codon_len))
                    logging.debug("type = " + ordered_codons[0])
                    logging.debug(pformat(exon))
                    logging.debug(pformat(isoform))
                    logging.debug("isoform_id = " + isoform_id)
            elif cumulative_transcript_length_old < (isoform['CDS start'] + 1) <= cumulative_transcript_length or cumulative_transcript_length_old < (isoform['CDS start'] + 2) <= cumulative_transcript_length:
                # exon contains at least part of the start codon, but not the first character
                # Note: the first character of the exon is in the first codon, hence 5utr length=0
                read_codon_len = min(isoform['CDS start'] + 2 - cumulative_transcript_length_old,
                                       cumulative_transcript_length - cumulative_transcript_length_old)
                read_start_codon_seq += exon['sequence'][:read_codon_len]
                if len(read_start_codon_seq) == 3 and check_codon("start", read_start_codon_seq):
                    dump_and_exit(exon, isoform, isoform_id)
            if read_codon_len > 0:
                if strand == '+':
                    exon['genome start codon start'] = exon['chromosome start'] + exon['5utr length']
                    exon['genome start codon end'] = exon['chromosome start'] + exon['5utr length'] + read_codon_len - 1
                else:
                    exon['genome start codon start'] = exon['chromosome end'] - exon['5utr length'] - read_codon_len + 1
                    exon['genome start codon end'] = exon['chromosome end'] - exon['5utr length']

            read_codon_len = 0
            if cumulative_transcript_length_old < isoform['CDS end'] <= cumulative_transcript_length:
                # exon contains at least part of the stop codon, including the last character
                read_codon_len = 3 - len(read_stop_codon_seq)
                final_pos = isoform['CDS end'] - cumulative_transcript_length_old
                read_stop_codon_seq += exon['sequence'][final_pos - read_codon_len:final_pos]
                #                import pdb; pdb.set_trace()
                if check_codon("stop", read_stop_codon_seq):
                    logging.debug(read_stop_codon_seq)
                    logging.debug(read_codon_len)
                    logging.debug(final_pos)
                    logging.debug(cumulative_transcript_length_old)
                    logging.debug(isoform['CDS end'])
                    logging.debug(cumulative_transcript_length)
                    dump_and_exit(exon, isoform, isoform_id)
            elif cumulative_transcript_length_old < (isoform['CDS end'] - 2) <= cumulative_transcript_length:
                # exon contains the first character, but not the last of the stop codon
                read_codon_len = cumulative_transcript_length - (isoform['CDS end'] - 3)
                read_stop_codon_seq += exon['sequence'][-read_codon_len:]
            elif cumulative_transcript_length_old < (isoform['CDS end'] - 1) <= cumulative_transcript_length:
                # exon contains only the second character of the stop codon
                read_codon_len = 1
                read_stop_codon_seq += exon['sequence'][0]
            if read_codon_len > 0:
                if strand == '+':
                    exon['genome stop codon start'] = exon['chromosome end'] - exon['3utr length'] - read_codon_len + 1
                    exon['genome stop codon end'] = exon['chromosome end'] - exon['3utr length']
                else:
                    exon['genome stop codon start'] = exon['chromosome start'] + exon['3utr length']
                    exon['genome stop codon end'] = exon['chromosome start'] + exon['3utr length'] + read_codon_len - 1

            if cumulative_transcript_length >= isoform['CDS start'] and cumulative_transcript_length_old < isoform['CDS end'] - 3:
                # exon contains at least a portion of the CDS
                exon['genome CDS start'] = exon['chromosome start'] + exon['5utr length'] if strand == '+' else exon['chromosome end'] - exon['5utr length']
                if 'genome stop codon start' in exon:
                    exon['genome CDS end'] = exon['genome stop codon start'] - 1 if strand == '+' else exon['genome stop codon end'] + 1
                else:
                    exon['genome CDS end'] = exon['chromosome end'] if strand == '+' else exon['chromosome start']

    #Now we can determine the frames
    for isoform_id, isoform in gene["isoforms"].items():
        if not isoform["annotated CDS?"]:
            continue
        if strand == '+':
            exons = range(isoform["number exons"])
        else:
            exons = range(isoform["number exons"] - 1, 0, -1)
        cumulative_transcript_length = 0
        #        import pdb; pdb.set_trace()
        # logging.debug(exons)
        for exon_id in range(isoform["number exons"]):
            frame = (3 - ((cumulative_transcript_length) % 3)) % 3
            exon = isoform['exons'][exon_id]
            if 'genome start codon end' in exon:
                exon['frame start codon'] = frame
            if 'genome CDS end' in exon:
                exon['frame CDS'] = frame
                cumulative_transcript_length += abs(exon['genome CDS end'] - exon['genome CDS start']) + 1
            if 'genome stop codon end' in exon:
                exon['frame stop codon'] = frame
                cumulative_transcript_length += abs(exon['genome stop codon end'] - exon['genome stop codon start']) + 1

    # import pdb; pdb.set_trace()
    # Clean up the data structure and write the json file
    del gene['factorizations']
    with open(output_file, mode='w', encoding='utf-8') as fd:
        fd.write(json.dumps(gene, sort_keys=True, indent=4))


def exec_system_command(command, error_comment, logfile, output_file="",
                        from_scratch=True):
    if not from_scratch or (not output_file == "" and not os.access(output_file, os.R_OK)):
        logging.debug(str(time.localtime()))
        logging.debug(command)

        try:
            retcode = subprocess.call(command + " 2>> " + logfile, shell=True)
            if retcode < 0:
                print(error_comment, -retcode, file=sys.stderr)
        except OSError as e:
            print("Execution failed:", e, file=sys.stderr)
            raise PIntronError


def check_executables(bindir, exes):
    """Check if the executables are in the path or in the specified directory.
    """

    full_exes = {}
    if bindir:
        if bindir[0] == '~':
            bindir = os.environ["HOME"] + bindir[1:]
        paths = [bindir] + os.environ["PATH"].split(os.pathsep)
    else:
        paths = os.environ["PATH"].split(os.pathsep)
    for exe in exes:
        full_exes[exe] = None
        for path in paths:
            if os.access(os.path.join(path, exe), os.X_OK):
                real_path = os.path.realpath(os.path.abspath(os.path.join(path, exe)))
                md5hex = md5Checksum(real_path)
                logging.debug("Using program '{}' in dir '{}' (md5: {})".format(exe, real_path, md5hex))
                full_exes[exe] = real_path
                break
        if full_exes[exe] == None:
            raise PIntronIOError(exe, "Could not find program '{}'!\n"
                                 "Search path: {}".format(exe,
                                                          ":".join(paths)))
    return full_exes


def pintron_pipeline(options):
    """Executes the whole pipeline, using the input options.
    """

    logging.info("PIntron%s", pintron_version)
    logging.info("Copyright (C) 2010,2011  Paola Bonizzoni, Gianluca Della Vedova, Yuri Pirola, Raffaella Rizzi.")
    logging.info("This program is distributed under the terms of the GNU Affero General Public License (AGPL), either version 3 of the License, or (at your option) any later version.")
    logging.info("This program comes with ABSOLUTELY NO WARRANTY. See the GNU Affero General Public License for more details.")
    logging.info("This is free software, and you are welcome to redistribute it under the conditions specified by the license.")

    logging.info("Running: " + " ".join(sys.argv))

    # Check and copy input data
    logging.info("STEP  1:  Checking executables and preparing input data...")

    logging.debug("Using main program 'pintron' in dir '{}' (md5: {})".format(os.path.realpath(os.path.abspath(sys.argv[0])),
                                                                              md5Checksum(sys.argv[0])))
    exes = check_executables(options.bindir, ["est-fact",
                                             "min-factorization",
                                             "intron-agreement",
                                             "gene-structure",
                                             "compact-compositions",
                                             "maximal-transcripts",
                                             "ests2sam",
                                             "cds-annotation"
                                             ])

    if not os.path.isfile(options.genome_filename) or not os.access(options.genome_filename, os.R_OK):
        raise PIntronIOError(options.genome_filename,
                             'Could not read file "' + options.genome_filename + '"!')
    if not os.path.isfile(options.EST_filename) or not os.access(options.EST_filename, os.R_OK):
        raise PIntronIOError(options.EST_filename,
                             'Could not read file "' + options.EST_filename + '"!')
    if (os.access('genomic.txt', os.F_OK) and
         not os.path.samefile('genomic.txt', options.genome_filename) and
         not os.access('genomic.txt', os.W_OK) ):
        raise PIntronIOError('genomic.txt',
                             'Could not write file "genomic.txt"!')
    if (os.access('ests.txt', os.F_OK) and
         not os.path.samefile('ests.txt', options.EST_filename) and
         not os.access('ests.txt', os.W_OK) ):
        raise PIntronIOError('ests.txt',
                             'Could not write file "ests.txt"!')

    if os.path.isfile('genomic.txt') and os.path.samefile('genomic.txt', options.genome_filename):
        logging.debug('Files "%s" and "genomic.txt" refer to the same file: skip copy.',
                      options.genome_filename)
    else:
        exec_system_command(
            command="cp " + options.genome_filename + " genomic.txt ",
            error_comment="Could not prepare genomic input file",
            logfile=options.plogfile,
            output_file='raw-multifasta-out.txt',
            from_scratch=options.from_scratch)

    if os.path.isfile('ests.txt') and os.path.samefile('ests.txt', options.EST_filename):
        logging.debug('Files "%s" and "ests.txt" refer to the same file: skip copy.',
                      options.EST_filename)
    else:
        exec_system_command(
            command="cp " + options.EST_filename + " ests.txt ",
            error_comment="Could not prepare ESTs input file",
            logfile=options.plogfile,
            output_file='raw-multifasta-out.txt',
            from_scratch=options.from_scratch)

    # Compute factorizations
    logging.info("STEP  2:  Pre-aligning transcript data...")

    exec_system_command(
        command="ulimit -t " + str(options.max_factorization_time * 60) + " && ulimit -v " +
        str(options.max_factorization_memory * 1024) + " && " + exes["est-fact"],
        error_comment="Could not compute the factorizations",
        logfile=options.plogfile,
        output_file='raw-multifasta-out.txt',
        from_scratch=options.from_scratch)

    # Min factorization agreement
    logging.info("STEP  3:  Computing a raw consensus gene structure...")

    exec_system_command(
        command="ulimit -t " + str(options.max_exon_agreement_time * 60) + " && " +
        exes["min-factorization"] + " < raw-multifasta-out.txt >out-agree.txt",
        error_comment="Could not minimize the factorizations",
        logfile=options.plogfile,
        output_file='out-agree.txt',
        from_scratch=options.from_scratch)

    # Intron prediction
    logging.info("STEP  4:  Predicting introns...")

    exec_system_command(
        command="ulimit -t " + str(options.max_intron_agreement_time * 60) + " && " +
        exes["intron-agreement"],
        error_comment="Could not compute the factorizations",
        logfile=options.plogfile,
        output_file='out-after-intron-agree.txt',
        from_scratch=options.from_scratch)

    # Compute the Gene structure
    logging.info("STEP  5:  Computing the final gene structure...")

    exec_system_command(
        command=exes["gene-structure"] + " ./",
        error_comment="Could not compute maximal transcripts",
        logfile=options.plogfile,
        output_file='out-after-intron-agree.txt',
        from_scratch=options.from_scratch)

    # TODO: Gene structure browser

    # The computation of the full-length isoforms should not be avoided
    # if options.step1:
    #     sys.exit(0)

    # Transform compositions into exons
    logging.info("STEP  6:  Computing the final transcript alignments...")

    exec_system_command(
        command=exes["compact-compositions"] + " < out-after-intron-agree.txt > build-ests.txt",
        error_comment="Could not transform factorizations into exons",
        logfile=options.plogfile,
        output_file='build-ests.txt',
        from_scratch=options.from_scratch)

    # Compute maximal transcripts
    logging.info("STEP  7:  Computing the final full-length isoforms...")

    exec_system_command(
        command=exes["maximal-transcripts"] + " < build-ests.txt",
        error_comment="Could not compute maximal transcripts",
        logfile=options.plogfile,
        output_file='CCDS_transcripts.txt',
        from_scratch=options.from_scratch)
    exec_system_command(
        command="cp -f TRANSCRIPTS1_1.txt isoforms.txt",
        error_comment="Could not link isoforms",
        logfile=options.plogfile,
        output_file='CCDS_transcripts.txt',
        from_scratch=options.from_scratch)

    # Annotate CDS
    logging.info("STEP  8:  Annotating CDS...")

    exec_system_command(
        command=exes["cds-annotation"] + " ./ ./ " + options.gene + " " + options.organism,
        error_comment="Could not annotate the CDSs",
        logfile=options.plogfile,
        output_file='CCDS_transcripts.txt',
        from_scratch=options.from_scratch)

    # TODO: Transcripts browser
    # Output the desired file
    logging.info("STEP  9:  Saving outputs...")

    json_output = compute_json(ccds_file="CCDS_transcripts.txt",
                             variant_file="VariantGTF.txt",
                             output_file=options.output_filename,
                             from_scratch=options.from_scratch,
                             pas_tolerance=options.pas_tolerance,
                             genomic_seq=options.genome_filename)

    if options.gtf_filename:
        json2gtf(options.output_filename, options.gtf_filename, options.gene, False)
    if options.extended_gtf_filename:
        logging.debug("""WARNING: you are creating a file that is not consistent with the GTF specifications.  See http://mblab.wustl.edu/GTF22.html""")
        json2gtf(options.output_filename, options.extended_gtf_filename, options.gene, True)
    if options.output_est_alignments:
        exec_system_command(
        command=exes["ests2sam"] + " --directory=. --genome=" + options.genome_filename,
        error_comment="Could not create ESTs-genome alignment file",
        logfile=options.plogfile,
        output_file='est-alignments.sam',
        from_scratch=options.from_scratch)

    # Clean mess
    logging.info("STEP 10:  Finalizing...")

    if options.compress:
        exec_system_command("gzip -q9 " + " ".join([options.output_filename,
                                                    options.plogfile,
                                                    options.glogfile]),
                            error_comment="Could not compress final files",
                            logfile="/dev/null",
        output_file=options.output_filename + '.gz')

    if not options.no_clean:
        tempfiles = ("TEMP_COMPOSITION_TRANS1_1.txt", "TEMP_COMPOSITION_TRANS1_2.txt",
                   "TEMP_COMPOSITION_TRANS1_3.txt", "TEMP_COMPOSITION_TRANS1_4.txt",
                   "TRANSCRIPTS1_1.txt", "TRANSCRIPTS1_2.txt", "TRANSCRIPTS1_3.txt", "TRANSCRIPTS1_4.txt",
                   "VariantGTF.txt", "build-ests.txt", "CCDS_transcripts.txt", "config-dump.ini",
                   "genomic-exonforCCDS.txt", "info-pid-*.log", "isoforms.txt", "gene-struct.txt",
                   "meg-edges.txt", "megs.txt", "out-after-intron-agree.txt", "out-agree.txt", "out-fatt.txt",
                   "predicted-introns.txt", "processed-ests.txt", "processed-megs-info.txt",
                   "processed-megs.txt", "raw-multifasta-out.txt", "time-limits")
        subprocess.call("rm -f " + " ".join(tempfiles), shell=True)


def prepare_loggers(options):
    """Prepare loggers.

    Save DEBUG and higher messages to options.glogfile, and INFO and higher messages to stdout
    Code adapted from
    http://docs.python.org/py3k/library/logging.html?highlight=logging#logging-to-multiple-destinations
    """
    logging.basicConfig(filename=options.glogfile,
                        filemode='w',
                        format='%(levelname)s:%(name)s:%(asctime)s%(msecs)d:%(message)s',
                        datefmt='%Y%m%d-%H%M%S',
                        level=logging.DEBUG)
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('[%(levelname)-8s] %(asctime)s - %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)


if __name__ == '__main__':

    try:
        options = parse_command_line()
## Program version
        pintron_version = '_____%PINTRON_VERSION%_____'
        if (pintron_version[0] == '_' and
            pintron_version[1:] == '____%PINTRON_VERSION%_____'):
            options.version = ''
        else:
            options.version = pintron_version
        prepare_loggers(options)
        pintron_pipeline(options)
    except PIntronError as err:
        logging.exception("*** Fatal error caught during the execution of the pipeline! ***\n"
                          "%s", err)
