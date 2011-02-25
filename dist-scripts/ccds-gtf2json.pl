#!/usr/bin/env perl
####
#
#
#                              PIntron
#
# A novel pipeline for computational gene-structure prediction based on
# spliced alignment of expressed sequences (ESTs and mRNAs).
#
# Copyright (C) 2010  Gianluca Della Vedova
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

use strict;
use warnings;
use JSON;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(max min);
use Switch;

our $debug = 0;
my $fileCCDS;
my $fileGTF;
my $fileOUT;
our $version=1; # File format version

GetOptions (
            'debug' => \$debug,
            'ccds=s' => \$fileCCDS,
            'gtf=s' => \$fileGTF,
#            'pivot=i' => \$indiceRigaBase,
            'output=s' => \$fileOUT,
           );
my $usage="Usage: perl ccds-gtf2json.pl [options]
 --ccds= <file>       (mandatory) name of the file containing the CCDSs
 --gtf= <file>        (mandatory) gtf is the  name of the GTF-encoded file
 --output= <file>     name of the output JSON file
 --debug              emits debugging information
";


die $usage if (not defined $fileCCDS or $fileCCDS eq '' or
               not defined $fileGTF or $fileGTF eq '');
for my $f (($fileCCDS, $fileGTF)) {
    die "File $f does not exists\n" unless (-f $f);
}



my $gene={"version" => $version};

open FA, "<", $fileCCDS or die "Could not read $fileCCDS: $!\n";
my $a=<FA>;
my $b=<FA>;
chomp ($a,$b);
# Must nomify integers, for otherwise they are stored as strings
($gene->{number_isoforms},$gene->{length_genomic_sequence})=($a+0,$b+0);


$gene->{isoforms}={};
my $index;
while (my $l=<FA>) {
    chomp $l;
    $l=~s/\s//g;
    $l=~s/#.*//;
    next if $l=~/^\t*$/;
    if ($l=~/^>/) {
        # New isoform
        my $isoform={};
        $l=~/^>(\d+):(\d+):(\d+):(\d+):(-?\d+)$/;
        $index=$1;

        $isoform->{'number exons'}=$2+0;
        $isoform->{'reference?'}=($3 == 0) ? JSON::false : JSON::true;
        $isoform->{'from RefSeq?'}=($4 == 0) ? JSON::false : JSON::true;
        $isoform->{'NMD flag'}=$5+0;
        $isoform->{exons}=[];
        $isoform->{"coding length"}=0;
        $isoform->{'polyA?'} = JSON::true,
        $isoform->{'annotated CDS?'} = JSON::true,
        $gene->{isoforms}->{$index}=$isoform;
        next;
    }
    if ($l=~/^(\d+:){5}(-?\d+:)(-?\d+)$/) {
        # Row contains exon metadata
        my @line=split(':', $l);
        my $exon={
                  "chromosome start" => max($line[0]+0,0),
                  "chromosome end"   => max($line[1]+0,0),
                  "relative start"     => $line[2]+0,
                  "relative end"       => $line[3]+0,
                  "5utr length"     => max($line[5]+0,0),
                  "3utr length"  => max($line[6]+0,0),
                 };
        $gene->{isoforms}->{$index}->{'polyA?'} = JSON::false if ($line[4]+0 == 0);
        $gene->{isoforms}->{$index}->{'annotated CDS?'} = JSON::false if ($line[5]+0 < 0);
        $gene->{isoforms}->{$index}->{"coding length"} +=
            max($exon->{"relative end"}, $exon->{"relative start"}) -
                min($exon->{"relative end"}, $exon->{"relative start"}) + 1 -
                $exon->{"5utr length"} - $exon->{"3utr length"};
        for my $k ("5utr length",  "3utr length") {
            delete $exon->{$k} if ($exon->{$k} < 0);
        }
        push(@{$gene->{isoforms}->{$index}->{exons}}, $exon);
        next;
    }
    # Row contains only the sequence
    if ($l=~/^[acgt]+$/i) {
        my $last_exon=$gene->{isoforms}->{$index}->{exons}->[-1];
        $last_exon->{sequence}=$l;
        next;
    }

    # If we arrive here, then the row is not correct
    $debug=1;
    print "Unparsable line follows:\n";
    print $l."\n";
    die ("Wrong line in CCDS file $fileCCDS\n", $gene);
}


close FA;
debug_print ("closed $fileCCDS\n", $gene);




open FB, "<", $fileGTF or die "Could not read $fileGTF: $!\n";
while (my $l=<FB>) {
    chomp $l;
    my @fields=split(' /', $l);
    my $index=shift @fields;
    $index=~s/^.*\#//;

    die "Cannot find isoform with index $index\n"
        unless (defined $gene->{isoforms}->{$index});
    my $isoform=$gene->{isoforms}->{$index};

    for my $t (@fields) {
        my ($k,$v)=split('=', $t, 2);
        if ($k eq "nex") {
            die "Wrong number of exons: $index\n $v != ".$isoform->{'number exons'}."\n"
                unless ($v+0 == $isoform->{'number exons'});
        } elsif ($k eq "L") {
            $isoform->{"CDS length"}=$v+0;
        } elsif ($k eq "CDS") {
            next if ($v eq '..');
            $v=~/^(<?)(\d+)\.\.(\d+)(>?)$/;
            my ($a,$b);
            ($a,$isoform->{"CDS start"},$isoform->{"CDS end"},$b)=($1,$2+0,$3+0,$4);
            $isoform->{'canonical start codon?'} = ($a eq '<') ? JSON::false : JSON::true;
            $isoform->{'canonical end codon?'}   = ($b eq '>') ? JSON::false : JSON::true;
        } elsif ($k eq "RefSeq") {
            next if ($v eq '..' or $v eq '' or not defined $isoform->{"CDS"});
            $v=~/^(.*?)(\(?([NY])([NY])\)?)?$/i;
            my ($a,$b);
            ($isoform->{RefSeq},$a,$b)=($1,$3,$4);
            $isoform->{'conserved start codon?'} = ($a eq 'N') ? JSON::false : JSON::true;
            $isoform->{'conserved end codon?'}   = ($b eq 'N') ? JSON::false : JSON::true;
            delete $isoform->{RefSeq} if ($isoform->{RefSeq} eq '');
        } elsif ($k eq "ProtL") {
            next if ($v eq '..' or not defined $isoform->{"CDS"});
            $v=~/^(>?)(\d+)$/i;
            my $a;
            ($a,$isoform->{'protein length'})=($1,$2);
            $isoform->{'protein missing codon?'} = ($a eq '>') ? JSON::true : JSON::false;
        } elsif ($k eq "Frame") {
            next if ($v eq '..' or not defined $isoform->{"CDS"});
            $isoform->{Frame}= ($v=~/^y/i) ? JSON::true : JSON::false;
        } elsif ($k eq "Type") {
            my $ref=($v eq 'Ref') ? JSON::true : JSON::false;
            die "Wrong reference for isoform n. $index\n" unless
                ($ref == $isoform->{'reference?'});
            next if ($isoform->{'reference?'});
            $isoform->{Type}=$v;
        } else {
            print STDERR "Cannot parse line in GTF file: $k => $v\n";
        }
    }
}

close FB;
debug_print ("closed $fileGTF\n", $gene);

my $fh;
if (defined $fileOUT) {
    open $fh, ">", $fileOUT or die "Could not write $fileOUT: $!\n";
} else {
    $fh=*STDOUT;
}

my $json_obj = JSON->new->utf8;
$json_obj = $json_obj->pretty(1); # Enable indentation

print $fh $json_obj->encode($gene);

close $fh;


sub debug_print {
    my $comment=shift;
    my $content=shift;
    if ($debug) {
        print "$comment";
        print system("date");
        print "\n";
        print Dumper($content);
    }
}
