#!/usr/bin/perl
die "perl $0 all_pixels.list loop.list resolution scoreCutoff\n" if(@ARGV != 4);
my $all_pixels_list=shift;
my $all_loop_bedpe=shift;
my $resolution=shift;
my $score_cutoff=shift;#0.01

my @loop_length_grades=(3000,100000);

my %already_used;
foreach my $i (0..$#loop_length_grades){
	$len=$loop_length_grades[$i];
	open(TMP,">tmp.loop_short_than_$len.bedpe") || die;
	open(ALB,$all_loop_bedpe) || die;
	my $head=<ALB>;
	print TMP $head;
	while(my $line=<ALB>){
	        chomp $line;
		if($already_used{$line}){
			next;
		}
	        my @sub=split/\s+/,$line;
		my @four=($sub[1],$sub[2],$sub[4],$sub[5]);
		@four=sort {$a<=>$b} @four;
		if($four[3]-$four[1] <= $len){
			print TMP $line,"\n";
			$already_used{$line}=1;
		}
	}
	close TMP;
	close ALB;

	my $sorted_loop=`perl 6.score_each_loops.v4.pl $all_pixels_list tmp.loop_short_than_$len.bedpe $resolution $score_cutoff`;
	chomp $sorted_loop;
	if($i == 0){
		print $sorted_loop,"\n";
	}
	else{
		my @all_lines=split/\n/,$sorted_loop;
		shift @all_lines;
		my $new_all_lines=join"\n",@all_lines;
		print $new_all_lines,"\n";
	}
	`rm -rf tmp.loop_short_than_$len.bedpe`;
}
		

