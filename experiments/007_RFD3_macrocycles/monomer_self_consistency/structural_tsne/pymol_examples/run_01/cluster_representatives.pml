reinitialize
bg_color white
set antialias, 2
set ray_shadows, off
set specular, 0
set ambient, 0.45
set direct, 0.55
set reflect, 0.15

load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_14/length_14_macrocycle_14_0_model_88.cif.gz, cluster_01_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_14/length_14_macrocycle_14_0_model_88/sequence_05/prediction.pdb, cluster_01_afcyc
pair_fit cluster_01_afcyc and name N+CA+C+O, cluster_01_design and name N+CA+C+O
alter cluster_01_design, chain='A'
alter cluster_01_afcyc, chain='B'
create cluster_01, (cluster_01_design or cluster_01_afcyc)
delete cluster_01_design
delete cluster_01_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_14/length_14_macrocycle_14_0_model_141.cif.gz, cluster_02_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_14/length_14_macrocycle_14_0_model_141/sequence_01/prediction.pdb, cluster_02_afcyc
pair_fit cluster_02_afcyc and name N+CA+C+O, cluster_02_design and name N+CA+C+O
alter cluster_02_design, chain='A'
alter cluster_02_afcyc, chain='B'
create cluster_02, (cluster_02_design or cluster_02_afcyc)
delete cluster_02_design
delete cluster_02_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_16/length_16_macrocycle_16_0_model_13.cif.gz, cluster_03_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_16/length_16_macrocycle_16_0_model_13/sequence_07/prediction.pdb, cluster_03_afcyc
pair_fit cluster_03_afcyc and name N+CA+C+O, cluster_03_design and name N+CA+C+O
alter cluster_03_design, chain='A'
alter cluster_03_afcyc, chain='B'
create cluster_03, (cluster_03_design or cluster_03_afcyc)
delete cluster_03_design
delete cluster_03_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_14/length_14_macrocycle_14_0_model_135.cif.gz, cluster_04_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_14/length_14_macrocycle_14_0_model_135/sequence_07/prediction.pdb, cluster_04_afcyc
pair_fit cluster_04_afcyc and name N+CA+C+O, cluster_04_design and name N+CA+C+O
alter cluster_04_design, chain='A'
alter cluster_04_afcyc, chain='B'
create cluster_04, (cluster_04_design or cluster_04_afcyc)
delete cluster_04_design
delete cluster_04_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_16/length_16_macrocycle_16_0_model_154.cif.gz, cluster_05_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_16/length_16_macrocycle_16_0_model_154/sequence_02/prediction.pdb, cluster_05_afcyc
pair_fit cluster_05_afcyc and name N+CA+C+O, cluster_05_design and name N+CA+C+O
alter cluster_05_design, chain='A'
alter cluster_05_afcyc, chain='B'
create cluster_05, (cluster_05_design or cluster_05_afcyc)
delete cluster_05_design
delete cluster_05_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_10/length_10_macrocycle_10_0_model_13.cif.gz, cluster_06_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_10/length_10_macrocycle_10_0_model_13/sequence_05/prediction.pdb, cluster_06_afcyc
pair_fit cluster_06_afcyc and name N+CA+C+O, cluster_06_design and name N+CA+C+O
alter cluster_06_design, chain='A'
alter cluster_06_afcyc, chain='B'
create cluster_06, (cluster_06_design or cluster_06_afcyc)
delete cluster_06_design
delete cluster_06_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_12/length_12_macrocycle_12_0_model_87.cif.gz, cluster_07_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_12/length_12_macrocycle_12_0_model_87/sequence_08/prediction.pdb, cluster_07_afcyc
pair_fit cluster_07_afcyc and name N+CA+C+O, cluster_07_design and name N+CA+C+O
alter cluster_07_design, chain='A'
alter cluster_07_afcyc, chain='B'
create cluster_07, (cluster_07_design or cluster_07_afcyc)
delete cluster_07_design
delete cluster_07_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_16/length_16_macrocycle_16_0_model_133.cif.gz, cluster_08_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_16/length_16_macrocycle_16_0_model_133/sequence_01/prediction.pdb, cluster_08_afcyc
pair_fit cluster_08_afcyc and name N+CA+C+O, cluster_08_design and name N+CA+C+O
alter cluster_08_design, chain='A'
alter cluster_08_afcyc, chain='B'
create cluster_08, (cluster_08_design or cluster_08_afcyc)
delete cluster_08_design
delete cluster_08_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_10/length_10_macrocycle_10_0_model_198.cif.gz, cluster_09_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_10/length_10_macrocycle_10_0_model_198/sequence_04/prediction.pdb, cluster_09_afcyc
pair_fit cluster_09_afcyc and name N+CA+C+O, cluster_09_design and name N+CA+C+O
alter cluster_09_design, chain='A'
alter cluster_09_afcyc, chain='B'
create cluster_09, (cluster_09_design or cluster_09_afcyc)
delete cluster_09_design
delete cluster_09_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_08/length_08_macrocycle_08_0_model_165.cif.gz, cluster_10_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_08/length_08_macrocycle_08_0_model_165/sequence_06/prediction.pdb, cluster_10_afcyc
pair_fit cluster_10_afcyc and name N+CA+C+O, cluster_10_design and name N+CA+C+O
alter cluster_10_design, chain='A'
alter cluster_10_afcyc, chain='B'
create cluster_10, (cluster_10_design or cluster_10_afcyc)
delete cluster_10_design
delete cluster_10_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_18/length_18_macrocycle_18_0_model_41.cif.gz, cluster_11_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_18/length_18_macrocycle_18_0_model_41/sequence_05/prediction.pdb, cluster_11_afcyc
pair_fit cluster_11_afcyc and name N+CA+C+O, cluster_11_design and name N+CA+C+O
alter cluster_11_design, chain='A'
alter cluster_11_afcyc, chain='B'
create cluster_11, (cluster_11_design or cluster_11_afcyc)
delete cluster_11_design
delete cluster_11_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_12/length_12_macrocycle_12_0_model_136.cif.gz, cluster_12_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_12/length_12_macrocycle_12_0_model_136/sequence_03/prediction.pdb, cluster_12_afcyc
pair_fit cluster_12_afcyc and name N+CA+C+O, cluster_12_design and name N+CA+C+O
alter cluster_12_design, chain='A'
alter cluster_12_afcyc, chain='B'
create cluster_12, (cluster_12_design or cluster_12_afcyc)
delete cluster_12_design
delete cluster_12_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_08/length_08_macrocycle_08_0_model_62.cif.gz, cluster_13_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_08/length_08_macrocycle_08_0_model_62/sequence_02/prediction.pdb, cluster_13_afcyc
pair_fit cluster_13_afcyc and name N+CA+C+O, cluster_13_design and name N+CA+C+O
alter cluster_13_design, chain='A'
alter cluster_13_afcyc, chain='B'
create cluster_13, (cluster_13_design or cluster_13_afcyc)
delete cluster_13_design
delete cluster_13_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_16/length_16_macrocycle_16_0_model_108.cif.gz, cluster_14_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_16/length_16_macrocycle_16_0_model_108/sequence_06/prediction.pdb, cluster_14_afcyc
pair_fit cluster_14_afcyc and name N+CA+C+O, cluster_14_design and name N+CA+C+O
alter cluster_14_design, chain='A'
alter cluster_14_afcyc, chain='B'
create cluster_14, (cluster_14_design or cluster_14_afcyc)
delete cluster_14_design
delete cluster_14_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_18/length_18_macrocycle_18_0_model_185.cif.gz, cluster_15_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_18/length_18_macrocycle_18_0_model_185/sequence_02/prediction.pdb, cluster_15_afcyc
pair_fit cluster_15_afcyc and name N+CA+C+O, cluster_15_design and name N+CA+C+O
alter cluster_15_design, chain='A'
alter cluster_15_afcyc, chain='B'
create cluster_15, (cluster_15_design or cluster_15_afcyc)
delete cluster_15_design
delete cluster_15_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_16/length_16_macrocycle_16_0_model_40.cif.gz, cluster_16_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_16/length_16_macrocycle_16_0_model_40/sequence_07/prediction.pdb, cluster_16_afcyc
pair_fit cluster_16_afcyc and name N+CA+C+O, cluster_16_design and name N+CA+C+O
alter cluster_16_design, chain='A'
alter cluster_16_afcyc, chain='B'
create cluster_16, (cluster_16_design or cluster_16_afcyc)
delete cluster_16_design
delete cluster_16_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_12/length_12_macrocycle_12_0_model_128.cif.gz, cluster_17_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_12/length_12_macrocycle_12_0_model_128/sequence_08/prediction.pdb, cluster_17_afcyc
pair_fit cluster_17_afcyc and name N+CA+C+O, cluster_17_design and name N+CA+C+O
alter cluster_17_design, chain='A'
alter cluster_17_afcyc, chain='B'
create cluster_17, (cluster_17_design or cluster_17_afcyc)
delete cluster_17_design
delete cluster_17_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_16/length_16_macrocycle_16_0_model_111.cif.gz, cluster_18_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_16/length_16_macrocycle_16_0_model_111/sequence_03/prediction.pdb, cluster_18_afcyc
pair_fit cluster_18_afcyc and name N+CA+C+O, cluster_18_design and name N+CA+C+O
alter cluster_18_design, chain='A'
alter cluster_18_afcyc, chain='B'
create cluster_18, (cluster_18_design or cluster_18_afcyc)
delete cluster_18_design
delete cluster_18_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_16/length_16_macrocycle_16_0_model_10.cif.gz, cluster_19_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_16/length_16_macrocycle_16_0_model_10/sequence_08/prediction.pdb, cluster_19_afcyc
pair_fit cluster_19_afcyc and name N+CA+C+O, cluster_19_design and name N+CA+C+O
alter cluster_19_design, chain='A'
alter cluster_19_afcyc, chain='B'
create cluster_19, (cluster_19_design or cluster_19_afcyc)
delete cluster_19_design
delete cluster_19_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_10/length_10_macrocycle_10_0_model_98.cif.gz, cluster_20_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_10/length_10_macrocycle_10_0_model_98/sequence_05/prediction.pdb, cluster_20_afcyc
pair_fit cluster_20_afcyc and name N+CA+C+O, cluster_20_design and name N+CA+C+O
alter cluster_20_design, chain='A'
alter cluster_20_afcyc, chain='B'
create cluster_20, (cluster_20_design or cluster_20_afcyc)
delete cluster_20_design
delete cluster_20_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_18/length_18_macrocycle_18_0_model_157.cif.gz, cluster_21_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_18/length_18_macrocycle_18_0_model_157/sequence_06/prediction.pdb, cluster_21_afcyc
pair_fit cluster_21_afcyc and name N+CA+C+O, cluster_21_design and name N+CA+C+O
alter cluster_21_design, chain='A'
alter cluster_21_afcyc, chain='B'
create cluster_21, (cluster_21_design or cluster_21_afcyc)
delete cluster_21_design
delete cluster_21_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_10/length_10_macrocycle_10_0_model_19.cif.gz, cluster_22_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_10/length_10_macrocycle_10_0_model_19/sequence_02/prediction.pdb, cluster_22_afcyc
pair_fit cluster_22_afcyc and name N+CA+C+O, cluster_22_design and name N+CA+C+O
alter cluster_22_design, chain='A'
alter cluster_22_afcyc, chain='B'
create cluster_22, (cluster_22_design or cluster_22_afcyc)
delete cluster_22_design
delete cluster_22_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_14/length_14_macrocycle_14_0_model_86.cif.gz, cluster_23_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_14/length_14_macrocycle_14_0_model_86/sequence_07/prediction.pdb, cluster_23_afcyc
pair_fit cluster_23_afcyc and name N+CA+C+O, cluster_23_design and name N+CA+C+O
alter cluster_23_design, chain='A'
alter cluster_23_afcyc, chain='B'
create cluster_23, (cluster_23_design or cluster_23_afcyc)
delete cluster_23_design
delete cluster_23_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_16/length_16_macrocycle_16_0_model_183.cif.gz, cluster_24_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_16/length_16_macrocycle_16_0_model_183/sequence_06/prediction.pdb, cluster_24_afcyc
pair_fit cluster_24_afcyc and name N+CA+C+O, cluster_24_design and name N+CA+C+O
alter cluster_24_design, chain='A'
alter cluster_24_afcyc, chain='B'
create cluster_24, (cluster_24_design or cluster_24_afcyc)
delete cluster_24_design
delete cluster_24_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_14/length_14_macrocycle_14_0_model_50.cif.gz, cluster_25_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_14/length_14_macrocycle_14_0_model_50/sequence_08/prediction.pdb, cluster_25_afcyc
pair_fit cluster_25_afcyc and name N+CA+C+O, cluster_25_design and name N+CA+C+O
alter cluster_25_design, chain='A'
alter cluster_25_afcyc, chain='B'
create cluster_25, (cluster_25_design or cluster_25_afcyc)
delete cluster_25_design
delete cluster_25_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_18/length_18_macrocycle_18_0_model_152.cif.gz, cluster_26_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_18/length_18_macrocycle_18_0_model_152/sequence_01/prediction.pdb, cluster_26_afcyc
pair_fit cluster_26_afcyc and name N+CA+C+O, cluster_26_design and name N+CA+C+O
alter cluster_26_design, chain='A'
alter cluster_26_afcyc, chain='B'
create cluster_26, (cluster_26_design or cluster_26_afcyc)
delete cluster_26_design
delete cluster_26_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_14/length_14_macrocycle_14_0_model_130.cif.gz, cluster_27_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_14/length_14_macrocycle_14_0_model_130/sequence_07/prediction.pdb, cluster_27_afcyc
pair_fit cluster_27_afcyc and name N+CA+C+O, cluster_27_design and name N+CA+C+O
alter cluster_27_design, chain='A'
alter cluster_27_afcyc, chain='B'
create cluster_27, (cluster_27_design or cluster_27_afcyc)
delete cluster_27_design
delete cluster_27_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_10/length_10_macrocycle_10_0_model_101.cif.gz, cluster_28_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_10/length_10_macrocycle_10_0_model_101/sequence_01/prediction.pdb, cluster_28_afcyc
pair_fit cluster_28_afcyc and name N+CA+C+O, cluster_28_design and name N+CA+C+O
alter cluster_28_design, chain='A'
alter cluster_28_afcyc, chain='B'
create cluster_28, (cluster_28_design or cluster_28_afcyc)
delete cluster_28_design
delete cluster_28_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_16/length_16_macrocycle_16_0_model_179.cif.gz, cluster_29_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_16/length_16_macrocycle_16_0_model_179/sequence_02/prediction.pdb, cluster_29_afcyc
pair_fit cluster_29_afcyc and name N+CA+C+O, cluster_29_design and name N+CA+C+O
alter cluster_29_design, chain='A'
alter cluster_29_afcyc, chain='B'
create cluster_29, (cluster_29_design or cluster_29_afcyc)
delete cluster_29_design
delete cluster_29_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_12/length_12_macrocycle_12_0_model_171.cif.gz, cluster_30_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_12/length_12_macrocycle_12_0_model_171/sequence_07/prediction.pdb, cluster_30_afcyc
pair_fit cluster_30_afcyc and name N+CA+C+O, cluster_30_design and name N+CA+C+O
alter cluster_30_design, chain='A'
alter cluster_30_afcyc, chain='B'
create cluster_30, (cluster_30_design or cluster_30_afcyc)
delete cluster_30_design
delete cluster_30_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_10/length_10_macrocycle_10_0_model_131.cif.gz, cluster_31_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_10/length_10_macrocycle_10_0_model_131/sequence_07/prediction.pdb, cluster_31_afcyc
pair_fit cluster_31_afcyc and name N+CA+C+O, cluster_31_design and name N+CA+C+O
alter cluster_31_design, chain='A'
alter cluster_31_afcyc, chain='B'
create cluster_31, (cluster_31_design or cluster_31_afcyc)
delete cluster_31_design
delete cluster_31_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_18/length_18_macrocycle_18_0_model_196.cif.gz, cluster_32_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_18/length_18_macrocycle_18_0_model_196/sequence_05/prediction.pdb, cluster_32_afcyc
pair_fit cluster_32_afcyc and name N+CA+C+O, cluster_32_design and name N+CA+C+O
alter cluster_32_design, chain='A'
alter cluster_32_afcyc, chain='B'
create cluster_32, (cluster_32_design or cluster_32_afcyc)
delete cluster_32_design
delete cluster_32_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_18/length_18_macrocycle_18_0_model_45.cif.gz, cluster_33_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_18/length_18_macrocycle_18_0_model_45/sequence_05/prediction.pdb, cluster_33_afcyc
pair_fit cluster_33_afcyc and name N+CA+C+O, cluster_33_design and name N+CA+C+O
alter cluster_33_design, chain='A'
alter cluster_33_afcyc, chain='B'
create cluster_33, (cluster_33_design or cluster_33_afcyc)
delete cluster_33_design
delete cluster_33_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_12/length_12_macrocycle_12_0_model_1.cif.gz, cluster_34_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_12/length_12_macrocycle_12_0_model_1/sequence_02/prediction.pdb, cluster_34_afcyc
pair_fit cluster_34_afcyc and name N+CA+C+O, cluster_34_design and name N+CA+C+O
alter cluster_34_design, chain='A'
alter cluster_34_afcyc, chain='B'
create cluster_34, (cluster_34_design or cluster_34_afcyc)
delete cluster_34_design
delete cluster_34_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_10/length_10_macrocycle_10_0_model_177.cif.gz, cluster_35_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_10/length_10_macrocycle_10_0_model_177/sequence_02/prediction.pdb, cluster_35_afcyc
pair_fit cluster_35_afcyc and name N+CA+C+O, cluster_35_design and name N+CA+C+O
alter cluster_35_design, chain='A'
alter cluster_35_afcyc, chain='B'
create cluster_35, (cluster_35_design or cluster_35_afcyc)
delete cluster_35_design
delete cluster_35_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_08/length_08_macrocycle_08_0_model_160.cif.gz, cluster_36_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_08/length_08_macrocycle_08_0_model_160/sequence_01/prediction.pdb, cluster_36_afcyc
pair_fit cluster_36_afcyc and name N+CA+C+O, cluster_36_design and name N+CA+C+O
alter cluster_36_design, chain='A'
alter cluster_36_afcyc, chain='B'
create cluster_36, (cluster_36_design or cluster_36_afcyc)
delete cluster_36_design
delete cluster_36_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_08/length_08_macrocycle_08_0_model_23.cif.gz, cluster_37_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_08/length_08_macrocycle_08_0_model_23/sequence_04/prediction.pdb, cluster_37_afcyc
pair_fit cluster_37_afcyc and name N+CA+C+O, cluster_37_design and name N+CA+C+O
alter cluster_37_design, chain='A'
alter cluster_37_afcyc, chain='B'
create cluster_37, (cluster_37_design or cluster_37_afcyc)
delete cluster_37_design
delete cluster_37_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_12/length_12_macrocycle_12_0_model_51.cif.gz, cluster_38_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_12/length_12_macrocycle_12_0_model_51/sequence_02/prediction.pdb, cluster_38_afcyc
pair_fit cluster_38_afcyc and name N+CA+C+O, cluster_38_design and name N+CA+C+O
alter cluster_38_design, chain='A'
alter cluster_38_afcyc, chain='B'
create cluster_38, (cluster_38_design or cluster_38_afcyc)
delete cluster_38_design
delete cluster_38_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_12/length_12_macrocycle_12_0_model_196.cif.gz, cluster_39_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_12/length_12_macrocycle_12_0_model_196/sequence_05/prediction.pdb, cluster_39_afcyc
pair_fit cluster_39_afcyc and name N+CA+C+O, cluster_39_design and name N+CA+C+O
alter cluster_39_design, chain='A'
alter cluster_39_afcyc, chain='B'
create cluster_39, (cluster_39_design or cluster_39_afcyc)
delete cluster_39_design
delete cluster_39_afcyc
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/rfd3_backbones/length_10/length_10_macrocycle_10_0_model_8.cif.gz, cluster_40_design
load /Users/lkv206/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/monomer_self_consistency/afcyc_predictions/run_01/length_10/length_10_macrocycle_10_0_model_8/sequence_05/prediction.pdb, cluster_40_afcyc
pair_fit cluster_40_afcyc and name N+CA+C+O, cluster_40_design and name N+CA+C+O
alter cluster_40_design, chain='A'
alter cluster_40_afcyc, chain='B'
create cluster_40, (cluster_40_design or cluster_40_afcyc)
delete cluster_40_design
delete cluster_40_afcyc

disable

color mocha_blue, cluster_* and elem C and chain A
color mocha_peach, cluster_* and elem C and chain B