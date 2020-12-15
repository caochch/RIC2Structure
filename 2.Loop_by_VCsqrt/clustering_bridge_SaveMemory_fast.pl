#!/usr/bin/perl
# file_name:clustering_bridge.pl
# Clustering procedure for long-read ChIA-PET
# Utilizes pairToPair function in Bedtools

use Graph::Undirected;
use List::Util qw(min max sum);
use List::MoreUtils qw/ uniq/;

# Usage:
die "perl $0 pairTopair.results out.cluster\n" if(@ARGV != 3);
my $pair_to_pair_result=shift;
my $output_cluster_name=shift;
my $output_readID_list=shift;

my %read_info;
my %link_target;
open(PAIR, $pair_to_pair_result) or die "Cannot open temp.pairToPair file!\n";
while($line = <PAIR>){
	chomp $line;
	@elements = split(/\t/, $line);
	$uid1 = $elements[6];
	$uid2 = $elements[16];
	$read_info{$uid1}{"head"}=join(':', $elements[0], $elements[1], $elements[2]);
	$read_info{$uid1}{"tail"}=join(':', $elements[3], $elements[4], $elements[5]);
	$read_info{$uid2}{"head"}=join(':', $elements[10], $elements[11], $elements[12]);
	$read_info{$uid2}{"tail"}=join(':', $elements[13], $elements[14], $elements[15]);
	if(exists $link_target{$uid1}){
		my @tmp=@{$link_target{$uid1}};
		push (@tmp,$uid2);
		@{$link_target{$uid1}}=@tmp;
	}
	else{
		@{$link_target{$uid1}}=($uid2);
	}
	if(exists $link_target{$uid2}){
		my @tmp=@{$link_target{$uid2}};
		push (@tmp,$uid1);
		@{$link_target{$uid2}}=@tmp;
	}
	else{
		@{$link_target{$uid2}}=($uid1);
	}
}
close PAIR;

#simplify and cluster
my %used;
my %clusters;
foreach my $id1 (keys %link_target){
	if($used{$id1}){
		next;
	}
	else{
		$cluster_id++;
		my @cluster_reads;
		my @tmp_id=($id1);
		while(1){
			if(@tmp_id){
				my %next_round;
				foreach my $i (@tmp_id){
					if($used{$i}){
						next;
					}
					else{
						$used{$i}=1;
						push (@cluster_reads,$i);
						my @target_tmp=@{$link_target{$i}};
						foreach my $j (@target_tmp){
							$next_round{$j}=1;
						}
					}
				}
				@tmp_id=keys %next_round;
			}
			else{
				last;
			}
		}
		@{$clusters{$cluster_id}}=@cluster_reads;
	}
}

open(OUTPUT, ">$output_cluster_name") or die "Cannot open output file\n!";
open(OUTR,">$output_readID_list") or die "Cannot open output file\n!";
foreach my $id (keys %clusters){
	my @members=@{$clusters{$id}};
        $cluster_size = $#members + 1;
        @cluster_head_chr = ();
        @cluster_head_start = ();
        @cluster_head_end = ();
        @cluster_tail_chr = ();
        @cluster_tail_start = ();
        @cluster_tail_end = ();
        foreach $m (@members){
                ($head_chr, $head_start, $head_end) = split(':', $read_info{$m}{'head'});
                ($tail_chr, $tail_start, $tail_end) = split(':', $read_info{$m}{'tail'});
                push @cluster_head_chr, $head_chr;
                push @cluster_head_start, $head_start;
                push @cluster_head_end, $head_end;
                push @cluster_tail_chr, $tail_chr;
                push @cluster_tail_start, $tail_start;
                push @cluster_tail_end, $tail_end;
                print OUTR $m,"\n";
        }
        @cluster_head_chr_uniq = uniq(@cluster_head_chr);
        $cluster_head_start_site = min(@cluster_head_start);
        $cluster_head_end_site = max(@cluster_head_end);
        @cluster_tail_chr_uniq = uniq(@cluster_tail_chr);
        $cluster_tail_start_site = min(@cluster_tail_start);
        $cluster_tail_end_site = max(@cluster_tail_end);
        print OUTPUT "@cluster_head_chr_uniq\t$cluster_head_start_site\t$cluster_head_end_site\t@cluster_tail_chr_uniq\t$cluster_tail_start_site\t$cluster_tail_end_site\t$cluster_size\n";
}


















