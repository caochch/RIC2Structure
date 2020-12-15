#!/usr/bin/perl
#my @resolution=(2,5,10,25,50);
my @resolution=(5);

#my @score_series=(0.05,0.03,0.01);
my @score_series=(0.02);


foreach my $r (@resolution){
	#`perl 1.calculate_VCsqrt_for_eachPixel.pl ../0.source/virus.onlyPart.sort.list ../0.source/virus.onlyPart.add_readLoci.withID.list $r > out1.onlyPart.all_pixels.resolution$r.score.list`;
	#`perl 2.compare_with_juicebox.pl out1.onlyPart.all_pixels.resolution$r.score.list ../1.Loop_by_HiCCUPS/virus.onlyPart.resolution$r.vcrt.matrix $r > out2.onlyPart.resolution$r.compare.xls`;
	foreach my $score_cutoff (@score_series){
		`perl 3.select_enriched_pixels.pl out1.onlyPart.all_pixels.resolution$r.score.list $score_cutoff $r > out3.onlyPart.enriched_pixels.resolution$r.score$score_cutoff.bedpe`;
		`perl 4.call_loops_by_cluster_pixels.pl out3.onlyPart.enriched_pixels.resolution$r.score$score_cutoff.bedpe $r > out4.onlyPart.loops_of_enriched_pixels.resolution$r.score$score_cutoff.bedpe`;
		`perl 5.format_as_region_loops.pl out4.onlyPart.loops_of_enriched_pixels.resolution$r.score$score_cutoff.bedpe > out5.onlyPart.loops_of_enriched_pixels.resolution$r.score$score_cutoff.loop_format.bedpe`;
		`perl meta_run_step6.pl out1.onlyPart.all_pixels.resolution$r.score.list out4.onlyPart.loops_of_enriched_pixels.resolution$r.score$score_cutoff.bedpe $r $score_cutoff > out6.onlyPart.loops_or_enriched_pixels.resolution$r.score$score_cutoff.sorted.bedpe`;
		`perl 7.format_as_region_beds.pl out6.onlyPart.loops_or_enriched_pixels.resolution$r.score$score_cutoff.sorted.bedpe > out7.onlyPart.loops_or_enriched_pixels.resolution$r.score$score_cutoff.sorted.bed`;
	}
}
