# Commands used on 2026-10-21 for knowledge/topics/rrna_annotation_methods.md (all run from one scratch working dir).
# Inputs: /mnt/n/yzz-分类工作流/projects/mito_benchmark/{gb,bare,eval/dev} and benchmark/sets/mito_bench_2026-10.tsv (dev split, 93 records);
# Infernal 1.1.5 from the mitos env; models: MitoZ v1.1_12snew.cm / v1.1_16snew.cm and MITOS2 refseq89m rrnS.cm / rrnL.cm.
# 1. python3 -I truth_flank.py            -> truth_flank.pkl (NCBI rRNA ends and flank gaps, dev)
# 2. for each dev genome FASTA (cm/fa/<acc>.fa, written from bare/<acc>/<acc>.gbf):
#      cmsearch -g  --cpu 1 --tblout cm/out/<acc>.<mz12|mz16|mt12|mt16>.g.tbl -o /dev/null <model.cm> cm/fa/<acc>.fa
#      cmsearch     --cpu 1 --tblout cm/out/<acc>.<model>.l.tbl            ...   (local mode)
# 3. python3 -I cm_eval.py ; cm_hist.py ; gap_vs_cm.py ; end_vs_nb.py ; rule_sim.py ; rule_sim_truthnb.py ; rr_err.py
# 4. python3 -I fetch_panel.py ; panel_an.py ; panel_sum.py ; panel_kmer.py   (RefSeq panel, 16 random per clade, benchmark accessions excluded)
# 5. cmalign -g --cpu 5 -o pal/<mt|mz>_<rrnS|rrnL>.sto --outformat pfam <model.cm> pal/<rrnS|rrnL>.fa ; python3 -I pal_cols.py ; tal_an.py (dev truth)
# 6. python3 -I clade_build.py rrnS pal/mt_rrnS.sto (and rrnL) ; copy the ECM* lines of the MITOS2 CM into each clade CM (no cmcalibrate: predicted 2 h per model);
#    cmsearch -g -T 10 with the clade CM ; python3 -I clade_eval.py
# 7. cmsearch -g -A aln/<acc>.<S|L>.sto --tblout aln/<acc>.<S|L>.tbl <MITOS2 model> ; python3 -I clip_eval.py
