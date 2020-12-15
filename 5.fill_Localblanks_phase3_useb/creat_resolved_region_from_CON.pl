#!/usr/bin/perl
die "perl $0 input.CON > resolved.region.bed\n" if(@ARGV != 1);
my $phase1_final_constrain_file=shift;

#public
my $output_sequence_id="NC_045512.2";


#initialize dynamic solved
my %dynamic_constrains;
my %dynamic_solved;
my %dynamic_basepairs;

open(PFCF,$phase1_final_constrain_file) || die;
while(my $line=<PFCF>){
        chomp $line;
        my @sub=split/\s+/,$line;
        if($line=~/-1/){
                next;
        }
        elsif($line=~/(\d+)\s+(\d+)/){
                my @pair=($1,$2);
                @pair=sort {$a<=>$b} @pair;
                $dynamic_constrains{$pair[0]."\t".$pair[1]}=$pair[0];
                $dynamic_basepairs{$pair[0]}=$pair[1];  #1-based
                $dynamic_basepairs{$pair[1]}=$pair[0];  #1-based
		foreach ($pair[0]..$pair[1]){
			$dynamic_solved{$_}=1;
		}
        }
}
#finished initialize dynamic solved

my $resolved_bases_blocks=print_resolved_region();
print $resolved_bases_blocks;


sub print_resolved_region{
        my @all_bases=sort {$a<=>$b} keys %dynamic_solved;
        my $bedfiles;
        foreach my $i (0..$#all_bases){
                if($i == 0){
                        my $zerobased=$all_bases[$i]-1;
                        $bedfiles="$output_sequence_id\t$zerobased\t";
                }
                elsif($i == $#all_bases){
                        my $tmp_loci=$all_bases[$i];
                        $bedfiles.="$tmp_loci\n";
                }
                else{
                        if($all_bases[$i] == $all_bases[$i-1]+1){
                                next;
                        }
                        else{
                                my $tmp_loci=$all_bases[$i-1];
                                $bedfiles.="$tmp_loci\n";
                                my $zerobased=$all_bases[$i]-1;
                                $bedfiles.="$output_sequence_id\t$zerobased\t";
                        }
                }
        }
        return $bedfiles;
}

