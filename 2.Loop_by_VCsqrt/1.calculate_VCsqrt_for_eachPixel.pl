#!/usr/bin/perl
die "perl $0 in.sort.list virus.add_readLoci.withID.list resolution\n" if(@ARGV != 3);
my $RICseq_sort_list=shift;
my $arms_loci_list=shift;
my $resolution=shift;

my %arms_loci;
open(ALL,$arms_loci_list) || die;
while(my $line=<ALL>){
	chomp $line;
	my @sub=split/\s+/,$line;
	$arms_loci{$sub[0]}=$line;
}

my $num_of_unique;
my %unique;
my %pairwise_unique;
my %pairwise;
my %coverage;
open(RSL,$RICseq_sort_list) || die;
while(my $line=<RSL>){
	chomp $line;
	my @sub=split/\s+/,$line;
	my $win_a=int($sub[3]/$resolution)*$resolution;
	my $win_b=int($sub[7]/$resolution)*$resolution;
	if($win_a eq $win_b){
		$coverage{$win_a}++;
		$pairwise{$win_a}{$win_b}++;
	}
	else{
		$coverage{$win_a}++;
		$coverage{$win_b}++;
		$pairwise{$win_a}{$win_b}++;
		$pairwise{$win_b}{$win_a}++;
	}
	
	my @arm_info=split/\s+/,$arms_loci{$sub[0]};

	#check first
	if($arm_info[2] ne $sub[3] or $arm_info[5] ne $sub[7]){
		warn $line;
		die;
	}

	#count uniq reads
	my $arms=join"\t",@arm_info[8..16];
	if($unique{$arms}){
		next;
	}
	else{
		$unique{$arms}=1;
		$num_of_unique++;
		if($win_a eq $win_b){
			$pairwise_unique{$win_a}{$win_b}++;
		}
		else{
			$pairwise_unique{$win_a}{$win_b}++;
			$pairwise_unique{$win_b}{$win_a}++;
		}
	}
}

warn $num_of_unique," unique pair tags\n";

my $start_win_index=0;
my $end_win_index=5980;

#print "Win\t";
#foreach ($start_win_index..$end_win_index){
#	print $_*$resolution,"\t";
#}
#print "\n";

foreach my $i ($start_win_index..$end_win_index){
	my $win_a=$i*$resolution;
	foreach my $j ($start_win_index..$end_win_index){
		if($j < $i){
			next;
		}
		my $win_b=$j*$resolution;
		if(!$pairwise{$win_a}{$win_b}){
			next;
		}

		print $win_a,"\t",$win_b,"\t";
		print $pairwise{$win_a}{$win_b},"\t";
		print $pairwise_unique{$win_a}{$win_b},"\t";
		my $depth_a=$coverage{$win_a};
		my $depth_b=$coverage{$win_b};
		my $sqrt_depth=sqrt($depth_a*$depth_b);
		my $score=$pairwise{$win_a}{$win_b}/$sqrt_depth;
		print $score,"\n";	
	}
}

