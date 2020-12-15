#!/usr/bin/perl
die "perl $0 all_pixels.list loop.list resolution scoreCutoff\n" if(@ARGV != 4);
my $all_pixels_list=shift;
my $all_loop_bedpe=shift;
my $resolution=shift;
my $score_cutoff=shift;#0.01

my %pixels_score;
open(APL,$all_pixels_list) || die;
while(my $line=<APL>){
	chomp $line;
	my @sub=split/\s+/,$line;
	$pixels_score{$sub[0]}{$sub[1]}=$sub[4];
}

my %loop_loci_info;
my %loops_score_info;
open(ALB,$all_loop_bedpe) || die;
my $head=<ALB>;
print $head;
while(my $line=<ALB>){
	chomp $line;
	my @sub=split/\s+/,$line;
	my $loop_id=$sub[6];
	if($loop_id =~ /EnrichPixel/){
		next;
	}
	$loop_loci_info{$loop_id}=$line;
	my $left_start_index=int($sub[1]/$resolution);
	my $left_end_index=int($sub[2]/$resolution)-1;
	my $right_start_index=int($sub[4]/$resolution);
	my $right_end_index=int($sub[5]/$resolution)-1;
	my $pixel_count;
	my $pixel_significant_count;
	my $pixel_sum;
	my $pixel_significant_sum;
	foreach my $i ($left_start_index..$left_end_index){
		foreach my $j ($right_start_index..$right_end_index){
			my $win_i=$resolution*$i;
			my $win_j=$resolution*$j;
			if($j < $i+4){
				next;
			}
			$pixel_count++;
			$pixel_sum+=$pixels_score{$win_i}{$win_j};
			if($pixels_score{$win_i}{$win_j} >= $score_cutoff){
				$pixel_significant_count++;
				$pixel_significant_sum+=$pixels_score{$win_i}{$win_j};
			}
		}
	}

	#if($loop_id eq "ManualLoops_21"){
	#	print $pixel_count,"\t",$pixel_sum,"\t",$pixel_significant_count,"\t",$pixel_significant_sum,"\n";
	#}

	$loops_score_info{$loop_id}{"Count"}=$pixel_count;
	$loops_score_info{$loop_id}{"Sum"}=$pixel_sum;
	$loops_score_info{$loop_id}{"sigCount"}=$pixel_significant_count;
	$loops_score_info{$loop_id}{"sigSum"}=$pixel_significant_sum;
}

my %father_son_trees;

foreach my $f_id (keys %loops_score_info){
	foreach my $s_id (keys %loops_score_info){
		if(is_son($loop_loci_info{$f_id},$loop_loci_info{$s_id}) eq "son"){
			$father_son_trees{$f_id}{$s_id}=1;
		}
	}
}

my @test_loops=("ManualLoops_196");
my %already_sorted;
foreach my $id (sort {$loops_score_info{$b}{"sigCount"} <=> $loops_score_info{$a}{"sigCount"}} keys %loops_score_info){
#foreach my $id (@test_loops){
	if($already_sorted{$id}){
		next;
	}
	my $sorted_string=get_sons($id);
	#print "|$sorted_string|\n";
	$sorted_string=~s/^\s+//;
	my @final_sorted_sons=split/\s+/,$sorted_string;
		
	foreach (@final_sorted_sons){
		print $loop_loci_info{$_},"\t";
		print $loops_score_info{$_}{"sigCount"},"\t";
		print $_."\tisSonOf\t".$id,"\n";
	}
}	


sub get_sons{
	my $f_id=shift;
	#print $f_id,"\tfather\n";
	if(exists $father_son_trees{$f_id}){
		$already_sorted{$f_id}=1;
		my $out;
		foreach my $s (sort {$loops_score_info{$b}{"sigCount"} <=> $loops_score_info{$a}{"sigCount"}} keys %{$father_son_trees{$f_id}}){
			if($already_sorted{$s}){
				#print $s,"\t",$loops_score_info{$s}{"sigCount"},"\tinput already assigned\n";
				next;
			}
			else{
				#print $s,"\t",$loops_score_info{$s}{"sigCount"},"\tinput\n";
			}
			$out.=get_sons($s)."\t";
		}
		return $out."\t".$f_id;
	}
	else{
		if($already_sorted{$f_id}){
			#print $f_id,"\tleaf but already used\n";
			return "";
		}
		else{
			#print $f_id,"\tleaf\n";
			$already_sorted{$f_id}=1;
			return $f_id;
		}
	}
}

sub is_son{
	my $f_loci=shift;
	my $s_loci=shift;
	my @f_loci_info=split/\s+/,$f_loci;
	my @s_loci_info=split/\s+/,$s_loci;
	
	my @f_four=($f_loci_info[1],$f_loci_info[2],$f_loci_info[4],$f_loci_info[5]);
	my @s_four=($s_loci_info[1],$s_loci_info[2],$s_loci_info[4],$s_loci_info[5]);

	my $f_region_len=$f_four[3]-$f_four[0];

	if($s_four[0] > $f_four[3]){
		return "notSon";
	}
	if($s_four[3] < $f_four[0]){
		return "notSon";
	}
	if($s_four[0] >= $f_four[1] and $s_four[3] <= $f_four[2]){
		return "son";
	}

	if($f_region_len > 600){	#predict short & strong loops at boundary first
		if($s_four[0] >= $f_four[0]-0.1*$f_region_len and $s_four[3] <= $f_four[1]+0.1*$f_region_len and $s_loci_info[6] !~ /EnrichPixel/){
			return "son"
		}
		if($s_four[0] >= $f_four[2]-0.1*$f_region_len and $s_four[3] <= $f_four[3]+0.1*$f_region_len and $s_loci_info[6] !~ /EnrichPixel/){
			return "son";
		}
	}
	return "notSon";
}


