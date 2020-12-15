#!/usr/bin/perl
die "perl $0 out1.all_pixels.score.list score_cutoff resolution\n" if(@ARGV != 3);
my $all_pixels_list=shift;
my $score_cutoff=shift;
my $resolution=shift;

my $minimal_loop_len=4*$resolution;

print "#chr1\tx1\tx2\tchr2\ty1\ty2\tname\tscore\tstrand1\tstrand2\n";
my $pixel_id;
open(APL,$all_pixels_list) || die;
while(my $line=<APL>){
	chomp $line;
	my @sub=split/\s+/,$line;
	if($sub[1] - $sub[0] < 0){
		die;
	}
	if($sub[1] - $sub[0] < $minimal_loop_len){	#minmal loop len
		next;
	}
	if($sub[3] < 2){	#unique reads
		next;
	}
	if($sub[4] < $score_cutoff){	#score
		next;
	}

	$pixel_id++;
	print "hs1\t$sub[0]\t",$sub[0]+$resolution,"\t";
	print "hs1\t$sub[1]\t",$sub[1]+$resolution,"\t";
	print "EnrichPixel_$pixel_id\t255\t+\t+\n";
}
