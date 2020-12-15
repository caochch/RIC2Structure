ln -s ../3a.predict_local_phase1_Score0.05/virus.phase1a_final.CON ./
rm -rf tmp.* log.* debug.txt
perl 1.Heuristics_predict.useRNAstructure.contain_input_Constraints.pl ../2.Loop_by_VCsqrt/out3.onlyPart.enriched_pixels.resolution5.score0.03.bedpe ../2.Loop_by_VCsqrt/out6.onlyPart.loops_or_enriched_pixels.resolution5.score0.03.sorted.bedpe virus.phase1a_final.CON > debug.txt

