#!/usr/bin/perl
die "perl $0 windowLen slidingLen\n" if(@ARGV != 2);
my $window=shift;
my $sliding=shift;
my $whole_len=29903;

print "#chr1\tx1\tx2\tchr2\ty1\ty2\tname\tscore\tstrand1\tstrand2\n";

foreach (0..int(($whole_len-$window)/$sliding)+1){
	my $start=$_*$sliding;
	my $end=$start+$window;
	$end = $end > $whole_len ? $whole_len : $end;
	print "hs1\t$start\t$end\ths1\t$start\t$end\tSw_$_\t255\t+\t+\n"
}
