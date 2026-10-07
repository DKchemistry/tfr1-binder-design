reinitialize
set auto_zoom, off
load ligandmpnn_sequences/run_01/length_08/backbones/length_08_macrocycle_08_0_model_0_1.pdb, target_0
load afcyc_predictions/run_01/length_08/length_08_macrocycle_08_0_model_0/sequence_01/prediction.pdb, pred_0
align pred_0 and backbone, target_0 and backbone
load ligandmpnn_sequences/run_01/length_08/backbones/length_08_macrocycle_08_0_model_1_1.pdb, target_1
load afcyc_predictions/run_01/length_08/length_08_macrocycle_08_0_model_1/sequence_01/prediction.pdb, pred_1
align pred_1 and backbone, target_1 and backbone
load ligandmpnn_sequences/run_01/length_08/backbones/length_08_macrocycle_08_0_model_2_1.pdb, target_2
load afcyc_predictions/run_01/length_08/length_08_macrocycle_08_0_model_2/sequence_01/prediction.pdb, pred_2
align pred_2 and backbone, target_2 and backbone
load ligandmpnn_sequences/run_01/length_08/backbones/length_08_macrocycle_08_0_model_3_1.pdb, target_3
load afcyc_predictions/run_01/length_08/length_08_macrocycle_08_0_model_3/sequence_01/prediction.pdb, pred_3
align pred_3 and backbone, target_3 and backbone
load ligandmpnn_sequences/run_01/length_08/backbones/length_08_macrocycle_08_0_model_4_1.pdb, target_4
load afcyc_predictions/run_01/length_08/length_08_macrocycle_08_0_model_4/sequence_01/prediction.pdb, pred_4
align pred_4 and backbone, target_4 and backbone
load ligandmpnn_sequences/run_01/length_08/backbones/length_08_macrocycle_08_0_model_5_1.pdb, target_5
load afcyc_predictions/run_01/length_08/length_08_macrocycle_08_0_model_5/sequence_01/prediction.pdb, pred_5
align pred_5 and backbone, target_5 and backbone
load ligandmpnn_sequences/run_01/length_08/backbones/length_08_macrocycle_08_0_model_6_1.pdb, target_6
load afcyc_predictions/run_01/length_08/length_08_macrocycle_08_0_model_6/sequence_01/prediction.pdb, pred_6
align pred_6 and backbone, target_6 and backbone
load ligandmpnn_sequences/run_01/length_08/backbones/length_08_macrocycle_08_0_model_7_1.pdb, target_7
load afcyc_predictions/run_01/length_08/length_08_macrocycle_08_0_model_7/sequence_01/prediction.pdb, pred_7
align pred_7 and backbone, target_7 and backbone
load ligandmpnn_sequences/run_01/length_08/backbones/length_08_macrocycle_08_0_model_8_1.pdb, target_8
load afcyc_predictions/run_01/length_08/length_08_macrocycle_08_0_model_8/sequence_01/prediction.pdb, pred_8
align pred_8 and backbone, target_8 and backbone
load ligandmpnn_sequences/run_01/length_08/backbones/length_08_macrocycle_08_0_model_9_1.pdb, target_9
load afcyc_predictions/run_01/length_08/length_08_macrocycle_08_0_model_9/sequence_01/prediction.pdb, pred_9
align pred_9 and backbone, target_9 and backbone
hide everything
show cartoon
color gray70, target_*
color cyan, pred_*
zoom
