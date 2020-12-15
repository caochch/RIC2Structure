#!/usr/bin/perl
die "perl $0 merged_loops.bedpe\n" if(@ARGV != 2);
my $hiccups_bedpe=shift;
my $resolution=shift;

my $window=5000;

my $factor=$window/$resolution;

open(OUT,">loops.resolution$resolution.bedpe") || die;
print OUT "#chr1\tx1\tx2\tchr2\ty1\ty2\tname\tscore\tstrand1\tstrand2\n";

my $id;
open(HB,$hiccups_bedpe) || die;
while(my $line=<HB>){
	if($line=~/^#/){
		print $line;
		next;
	}
	chomp $line;
	my @sub=split/\s+/,$line;
	$sub[0]="hs1";
	$sub[1]=int($sub[1]/$factor);
	$sub[2]=int($sub[2]/$factor);
	$sub[3]="hs1";
	$sub[4]=int($sub[4]/$factor);
	$sub[5]=int($sub[5]/$factor);
	my $new=join"\t",@sub;
	print $new,"\n";

	$id++;
	my $numCollapsed=$sub[20];
	my $centroid1=$sub[21];
	my $centroid2=$sub[22];
	my $radius=$sub[23];
	$radius=$radius+0.501*$window;

	my $x_left=int(0.5+($centroid1-$radius)/$factor);
	my $x_right=int(0.5+($centroid1+$radius)/$factor);
	my $y_left=int(0.5+($centroid2-$radius)/$factor);
	my $y_right=int(0.5+($centroid2+$radius)/$factor);

	print OUT "hs1\t$x_left\t$x_right\ths1\t$y_left\t$y_right\tloop_$id\t255\t+\t+\n";

}
