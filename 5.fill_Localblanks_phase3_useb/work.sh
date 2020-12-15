ln -s ../3b.predict_local_phase1_Score0.03/virus.phase1b_final.CON ./
perl creat_resolved_region_from_CON.pl virus.phase1b_final.CON > log.resolved_bases.phase1b_final.real.bed
perl 0.creat_bedpe_from_region.pl log.resolved_bases.phase1b_final.real.bed > out1.loops_for_solved.bedpe
perl 1.fill_local_blanks.useRNAstructure.v1.pl virus.phase1b_final.CON out1.loops_for_solved.bedpe > debug.txt
