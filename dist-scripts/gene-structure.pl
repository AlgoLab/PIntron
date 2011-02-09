#!/usr/bin/env perl
####
#
#
#                              PIntron
#
# A novel pipeline for computational gene-structure prediction based on
# spliced alignment of expressed sequences (ESTs and mRNAs).
#
# Copyright (C) 2010  Raffaella Rizzi
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

#Threshold for deleting terminal edges
my $term_threshold=20;

foreach my $dir (@ARGV) {
    if($dir eq "./"){
        $dir=".";
    }

    if ((!-e $dir."/predicted-introns.txt")) {
        print STDERR "SKIPPING DIRECTORY ${dir}\n";
        next;
    }
    print STDERR "Analyzing directory $dir\n";

    open GEN, "<".$dir."/genomic.txt" or die "Can't open file. $!";
    open IN, "<".$dir."/predicted-introns.txt" or die "Can't open file. $!";
    open INC, "<".$dir."/out-after-intron-agree.txt" or die "Can't open file. $!";
    open OUT, ">".$dir."/gene-struct.txt" or die "Can't create file. $!";

    my $gen_header=<GEN>;
    close GEN;
    chomp $gen_header;
    $gen_header =~ m/^>chr([xXyY\d]+):(\d+):(\d+):([-+]{0,1}1)/i or die "Header in genomic.txt uncorrect!\n";
 
    my $chr_start=($2 < $3)?($2):($3);
    my $chr_end=($2 < $3)?($3):($2);
    my $strand=$4;

    my @pred_intron_list=();

    my %index_hash=();
    my @index_array=();
 
    while (<IN>) {
        chomp;
        if($_ !~ m/\s*(\d+)\s+(\d+)\s+[\d\.]+\s+[\d\.]+\s+\d+\s+\d+\s+([\w\,]+)/){
        die "Wrong format!\n";
      }
      my $est_ids=$3;
      my $left=$1;
      $left-=1;
      my $right=$2;

      my @ests=split /\,/, $est_ids;

      my %est_hash=();
      foreach my $est(@ests){
        $est_hash{$est}=1;
      }

      my @intron=();

      push @intron, $left;
      push @intron, $right;
      push @intron, \%est_hash;

      $index_hash{$left}=-1;
      $index_hash{$right}=-1;

      push @pred_intron_list, \@intron;
    }

    close IN;

    my @intron_edge=();
    my @exon_edge=();
    my @structure_edge=();

    my @keys = sort {$a <=> $b} keys %index_hash;

    #200:($max-$min)=x:regione ==> regione*200/($max-$min)
    my $visual_length=10000;
    my $min=$keys[0];
    my $max=$keys[$#keys];

    my $index=0;
    
    $index_hash{0}=$index;
    $index_array[$index]=0;
    $intron_edge[$index]=0;
    $exon_edge[$index]=0;
    $structure_edge[$index]=0;
    $index++;    
    
    foreach my $key(@keys){
      $index_hash{$key}=$index;
      $index_array[$index]=$key;
      $intron_edge[$index]=0;
      $exon_edge[$index]=0;
      $structure_edge[$index]=0;
      $index++;
    }
    
    my $upper_bound=$chr_end-$chr_start+1;
    $index_hash{$upper_bound}=$index;
    $index_array[$index]=$upper_bound;
    $intron_edge[$index]=0;
    $exon_edge[$index]=0;
    $structure_edge[$index]=0;
    $index++;    

#    my @keys = sort {$a <=> $b} keys %index_hash;
#        foreach my $key(@keys ){
#        print $key, " ", $index_hash{$key}, "\n";
#    }

#     foreach my $ss(@index_array){
#       print $ss, "\n";
#     }
   
   
   my %out_additional_intron_edges=();
   my %in_additional_intron_edges=();

    #Build intron graph
    foreach my $intron_ref(@pred_intron_list){
        my @intron=@{$intron_ref};

        my $index_left=$index_hash{$intron[0]};
        my $index_right=$index_hash{$intron[1]};

        if($index_right <= $index_left){
            die "Failure 1!\n";
        }
        
        my @edges_out;
        if(exists $out_additional_intron_edges{$index_left}){
            @edges_out=@{$out_additional_intron_edges{$index_left}};
          }
        else{
            @edges_out=();
        }
        push @edges_out, $index_right;
        $out_additional_intron_edges{$index_left}=\@edges_out;
        
        
        my @edges_in;
        if(exists $in_additional_intron_edges{$index_right}){
            @edges_in=@{$in_additional_intron_edges{$index_right}};
          }
        else{
            @edges_in=();
        }
        push @edges_in, $index_left;
        $in_additional_intron_edges{$index_right}=\@edges_in;
        
        foreach my $in($index_left..($index_right-1)){
            $intron_edge[$in]=1;
        }
    }
    
#     for(my $i=0; $i<=$#intron_edge; $i++){
#       print $index_array[$i], " ", $intron_edge[$i], "\n";
#     }

    my %out_additional_exon_edges=();
    my %in_additional_exon_edges=();

    #Build exon graph
    for(my $i=0; $i<=$#pred_intron_list; $i++){
        for(my $j=$i+1; $j<=$#pred_intron_list; $j++){
        my @intron1=@{$pred_intron_list[$i]};
        my @intron2=@{$pred_intron_list[$j]};

        my %ests1=%{$intron1[2]};
        my %ests2=%{$intron2[2]};

        my @keys1=keys %ests1;
        my $found=0;
        my $k=0;
        while($k <= $#keys1){
            if($ests1{$keys1[$k]} == 1 && exists $ests2{$keys1[$k]}){
                $found=1;
                $ests1{$keys1[$k]}=0;
            }
            $k++;
        }   

        $intron1[2]=\%ests1;
        $pred_intron_list[$i]=\@intron1;

        if($found == 1){
            my $index_left=-1;
            my $index_right=-1;
            $index_left=$index_hash{$intron1[1]};
            $index_right=$index_hash{$intron2[0]};
            #print $intron1[1], " ", $intron2[0], "\n";
            
             if($index_right <= $index_left){
                die "Failure 2!\n";
            }
            
            my @edges_out;
            if(exists $out_additional_exon_edges{$index_left}){
                @edges_out=@{$out_additional_exon_edges{$index_left}};
            }
            else{
               @edges_out=();
            }
            push @edges_out, $index_right;
            $out_additional_exon_edges{$index_left}=\@edges_out;

           my @edges_in;
           if(exists $in_additional_exon_edges{$index_right}){
                @edges_in=@{$in_additional_exon_edges{$index_right}};
            }
            else{
               @edges_in=();
            }
            push @edges_in, $index_left;
            $in_additional_exon_edges{$index_right}=\@edges_in;
            
            foreach my $in($index_left..($index_right-1)){
                $exon_edge[$in]=1;
            }
        }   
     }
    }
    
#     for(my $i=0; $i<=$#exon_edge; $i++){
#       print $index_array[$i], " ", $exon_edge[$i], "\n";
#     }

    my %left_terminal_exon_edges=();    
    my %right_terminal_exon_edges=();    
    
    #Spliced alignments composed of one exon are not considered
    my $last_exon="";
    my $first_exon=0;
    my @coords;
    my @first_exon_coord;
    while(<INC>){
           if(/^#/){
            next;
           }
   
           if(/^>/){
                if($last_exon ne ""){
                    @coords=split /\s/, $last_exon;
                    my $left=$coords[2];
                    my $right=$coords[3];
                    if($left != $first_exon_coord[0]){
                       if(!exists $index_hash{$first_exon_coord[1]}){
                            die "Failure 3!\n";
                       }
                       if(!exists $index_hash{$left-1}){
                            die "Failure 4!\n";
                        }
                        if(exists $left_terminal_exon_edges{$index_hash{$first_exon_coord[1]}}){
                            my $left_t=$left_terminal_exon_edges{$index_hash{$first_exon_coord[1]}};
                            if($left_t > $first_exon_coord[0]-1){
                                $left_terminal_exon_edges{$index_hash{$first_exon_coord[1]}}=$first_exon_coord[0]-1;
                            }
                        }
                        else{
                            $left_terminal_exon_edges{$index_hash{$first_exon_coord[1]}}=$first_exon_coord[0]-1;
                        }

                        if(exists $right_terminal_exon_edges{$index_hash{$left-1}}){
                            my $right_t=$right_terminal_exon_edges{$index_hash{$left-1}};
                            if($right_t < $right){
                                $right_terminal_exon_edges{$index_hash{$left-1}}=$right;
                            }
                        }
                        else{
                            $right_terminal_exon_edges{$index_hash{$left-1}}=$right;
                        }
                    }
                }
                $first_exon=1;
            }
            else{
                $last_exon=$_;
                if($first_exon == 1){
                     @coords=split /\s/, $last_exon;
                     @first_exon_coord=($coords[2], $coords[3]);
                }
                $first_exon=0;
            }
    }

    if($last_exon ne ""){
         @coords=split /\s/, $last_exon;
         my $left=$coords[2];
         my $right=$coords[3];
         if($left != $first_exon_coord[0]){
            if(!exists $index_hash{$first_exon_coord[1]}){
                 die "Failure 5!\n";
            }
            if(!exists $index_hash{$left-1}){
                 die "Failure 6!\n";
            }
            if(exists $left_terminal_exon_edges{$index_hash{$first_exon_coord[1]}}){
                my $left_t=$left_terminal_exon_edges{$index_hash{$first_exon_coord[1]}};
                if($left_t > $first_exon_coord[0]-1){
                     $left_terminal_exon_edges{$index_hash{$first_exon_coord[1]}}=$first_exon_coord[0]-1;
                }
             }
             else{
                 $left_terminal_exon_edges{$index_hash{$first_exon_coord[1]}}=$first_exon_coord[0]-1;
             }

             if(exists $right_terminal_exon_edges{$index_hash{$left-1}}){
                 my $right_t=$right_terminal_exon_edges{$index_hash{$left-1}};
                 if($right_t < $right){
                   $right_terminal_exon_edges{$index_hash{$left-1}}=$right;
                 }
             }
             else{
                 $right_terminal_exon_edges{$index_hash{$left-1}}=$right;
             }
      }
    }

    close INC;
    
    foreach my $site_index(keys %left_terminal_exon_edges){
        my $exists_in_exon_edge=0;
        if(exists $in_additional_exon_edges{$site_index}){
            my @first_vertex=@{$in_additional_exon_edges{$site_index}};
            foreach my $in(@first_vertex){
                my $left=$index_array[$in];
                if($left < $left_terminal_exon_edges{$site_index}){
                    $exists_in_exon_edge=1;
                }
            }
        }
        
        if($exists_in_exon_edge == 1){
            $left_terminal_exon_edges{$site_index}=-1;
        }
    }
 
   foreach my $site_index(keys %right_terminal_exon_edges){
        my $exists_in_exon_edge=0;
        
        if(exists $out_additional_exon_edges{$site_index}){
            my @last_vertex=@{$out_additional_exon_edges{$site_index}};
            foreach my $out(@last_vertex){
                my $right=$index_array[$out];
                if($right > $right_terminal_exon_edges{$site_index}){
                    $exists_in_exon_edge=1;
                }
            }
        }
        
        if($exists_in_exon_edge == 1){
            $right_terminal_exon_edges{$site_index}=-1;
        }
    }
    
 
    #Hash terminal sites
    my %terminal_index_hash=();
     
    foreach my $site_index(keys %left_terminal_exon_edges){
        if($left_terminal_exon_edges{$site_index} != -1){
             $terminal_index_hash{$left_terminal_exon_edges{$site_index}}=-1;
          }
    }

    foreach my $site_index(keys %right_terminal_exon_edges){
        if($right_terminal_exon_edges{$site_index} != -1){
             $terminal_index_hash{$right_terminal_exon_edges{$site_index}}=-1;
          }
    }

    my @terminal_index_array=();
    my $offset=scalar(@index_array);
    my @t_keys = sort {$a <=> $b} keys %terminal_index_hash;
    $index=$offset;
    foreach my $key(@t_keys){
      $terminal_index_hash{$key}=$index;
      $terminal_index_array[$index-$offset]=$key;
      $index++;
    }
    
     foreach my $l_site_index(keys %left_terminal_exon_edges){
        if($left_terminal_exon_edges{$l_site_index} != -1){
             $left_terminal_exon_edges{$l_site_index}=$terminal_index_hash{$left_terminal_exon_edges{$l_site_index}};
          }
    }

     foreach my $r_site_index(keys %right_terminal_exon_edges){
        if($right_terminal_exon_edges{$r_site_index} != -1){
             $right_terminal_exon_edges{$r_site_index}=$terminal_index_hash{$right_terminal_exon_edges{$r_site_index}};
        }
    }
  
     for(my $i=0; $i<=$#structure_edge; $i++){
        #Intron region
        if($intron_edge[$i] == 1 && $exon_edge[$i] == 0){
            $structure_edge[$i]=1;
        }
        #Exon region
        if($intron_edge[$i] == 0 && $exon_edge[$i] == 1){
            $structure_edge[$i]=2;
        }
        #Alternative region
        if($intron_edge[$i] == 1 && $exon_edge[$i] == 1){
            $structure_edge[$i]=3;
         }
      }
    
    my $site_counter=1;
    for(my $i=1; $i < $#index_array; $i++, $site_counter++){
        my $abs_coord=($strand eq "+1")?($chr_start+$index_array[$i]):($chr_end-$index_array[$i]);
        print OUT "site\#", $site_counter, " ", $index_array[$i], " $abs_coord", " s";
        if($i < $#index_array-1){
            if($structure_edge[$i] != 0){
                my $col="";
                if($structure_edge[$i] == 1){
                    $col="i";
                }elsif($structure_edge[$i] == 2){
                    $col="e";
                }
                else{
                    $col="a";
                }
                print OUT " $col";
            }
            else{
                print OUT " b";
            }
         }
         print OUT "\n";
    }
    
    my @filter_terminal_sites=();
    foreach my $term_site_index(@terminal_index_array){
        push @filter_terminal_sites, -1;
    }
    my $filter_index=0;
    
    #Filtering left terminal sites
    foreach my $splice(sort {$a <=> $b} keys %left_terminal_exon_edges){
        if($left_terminal_exon_edges{$splice} != -1){
            my $edge=$left_terminal_exon_edges{$splice};
            my $term_coord=$terminal_index_array[$edge-$offset];
 
            my $k=0;
 
            while($k <= $#index_array && $index_array[$k] <= $term_coord){
                $k++;    
            }
            
            $k--;
            
            if($k < 0 || $k+1 > $#index_array){
                die "Failure 7!\n";
            }
            
            if($k >= $splice){
                die "Failure 8!\n";
            }
            
            my $intron_coverage=0;
            my $q=$k;
            my $cover_start=$term_coord;      
            while($q < $splice){
                #Intron or blank
                if($structure_edge[$q] <= 1){
                    $intron_coverage+=$index_array[$q+1]-$cover_start;
                }
                $q++;
                $cover_start=$index_array[$q];
            }
            
            if($intron_coverage < $term_threshold){
                $left_terminal_exon_edges{$splice}=-1;
            }
            else{
                if($filter_terminal_sites[$edge-$offset] == -1){
                    $filter_terminal_sites[$edge-$offset]=$filter_index;
                    $filter_index++;
                }
                $left_terminal_exon_edges{$splice}=$filter_terminal_sites[$edge-$offset]+$offset;
            }   
        }
     }
 
    #Filtering right terminal sites
    foreach my $splice(sort {$a <=> $b} keys %right_terminal_exon_edges){
       if($right_terminal_exon_edges{$splice} != -1){
            my $edge=$right_terminal_exon_edges{$splice};
            my $term_coord=$terminal_index_array[$edge-$offset];
 
            my $k=$#index_array;
            while($k >= 0 && $index_array[$k] >= $term_coord){
                $k--;    
            }

            $k++;
            
            if($k-1 < 0 || $k > $#index_array){
                die "Failure 9!\n";
            }
            
            if($k <= $splice){
                die "Failure 10!\n";
            }
            
            my $intron_coverage=0;
            my $q=$k;
            my $cover_end=$term_coord;      
            while($q > $splice){
                #Intron or blank
                if($structure_edge[$q-1] <= 1){
                    $intron_coverage+=$cover_end-$index_array[$q-1];
                }
                $q--;
                $cover_end=$index_array[$q];
            }     
 
            #print "\t coverage $intron_coverage\n";

            if($intron_coverage < $term_threshold){
                $right_terminal_exon_edges{$splice}=-1;
            }    
            else{
                if($filter_terminal_sites[$edge-$offset] == -1){
                    $filter_terminal_sites[$edge-$offset]=$filter_index;
                    $filter_index++;
                }
                $right_terminal_exon_edges{$splice}=$filter_terminal_sites[$edge-$offset]+$offset;
            }   
 
         }
   }
  
   #for(my $i=0; $i <= $#terminal_index_array; $i++, $site_counter++){
   #     my $abs_coord=($strand eq "+1")?($chr_start+$terminal_index_array[$i]):($chr_end-$terminal_index_array[$i]);
   #     print OUT "site\#", $site_counter, " ", $terminal_index_array[$i], " $abs_coord", " t", "\n";
   #}

   for(my $i=0; $i <= $#filter_terminal_sites; $i++){
        if($filter_terminal_sites[$i] != -1){
            my $abs_coord=($strand eq "+1")?($chr_start+$terminal_index_array[$i]):($chr_end-$terminal_index_array[$i]);
            print OUT "site\#", ($filter_terminal_sites[$i]+$offset-1), " ", $terminal_index_array[$i], " $abs_coord", " t", "\n";
        }
   }
    
    my $i_edge_counter=1;
    foreach my $splice(sort {$a <=> $b} keys %out_additional_intron_edges){
        my @edges=@{$out_additional_intron_edges{$splice}};
        my %unique_hash=();
        my @unique_edges=();
        foreach my $edge(@edges){
            if(!exists $unique_hash{$edge}){
                push @unique_edges, $edge;
                $unique_hash{$edge}=1;
            }
        }    
        foreach my $edge(@unique_edges){
            print OUT "edge\#", $i_edge_counter, " ", $splice, " ", $edge, " i", "\n";
            $i_edge_counter++;
        }
    }
 
    my $e_edge_counter=0;
    foreach my $splice(sort {$a <=> $b} keys %out_additional_exon_edges){
        my @edges=@{$out_additional_exon_edges{$splice}};
        my %unique_hash=();
        my @unique_edges=();
        foreach my $edge(@edges){
            if(!exists $unique_hash{$edge}){
                push @unique_edges, $edge;
                $unique_hash{$edge}=1;
            }
        }
        foreach my $edge(@unique_edges){
            print OUT "edge\#", $e_edge_counter+$i_edge_counter, " ", $splice, " ", $edge, " e", "\n";
            $e_edge_counter++;
        }
    }
 
    my $t_edge_counter=0;
    foreach my $splice(sort {$a <=> $b} keys %left_terminal_exon_edges){
        if($left_terminal_exon_edges{$splice} != -1){
            my $edge=$left_terminal_exon_edges{$splice};
            print OUT "edge\#", $t_edge_counter+$e_edge_counter+$i_edge_counter, " ", $edge-1, " ", $splice, " t", "\n";
            $t_edge_counter++;
        }
     }

    foreach my $splice(sort {$a <=> $b} keys %right_terminal_exon_edges){
       if($right_terminal_exon_edges{$splice} != -1){
            my $edge=$right_terminal_exon_edges{$splice};
            print OUT "edge\#", $t_edge_counter+$e_edge_counter+$i_edge_counter, " ", $splice, " ", $edge-1, " t", "\n";
            $t_edge_counter++;
        }
     }
    
    close OUT;
}

exit;
