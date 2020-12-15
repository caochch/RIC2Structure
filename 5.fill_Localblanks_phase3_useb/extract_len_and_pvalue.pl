#!/usr/bin/perl
die "perl $0 ../debug.txt\n" if(@ARGV != 1);
my $debug_txt=shift;

$/="running cycle";
open(DT,$debug_txt) || die;
<DT>;
while(my $block=<DT>){
	chomp $block;
	my @lines=split/\n/,$block;
	shift @lines;
	my $seq=$lines[1];
	print length($seq),"\t";
	foreach my $l (@lines){
		if($l=~/selected$/){
			print $l,"\n";
			last;
		}
	}
}

