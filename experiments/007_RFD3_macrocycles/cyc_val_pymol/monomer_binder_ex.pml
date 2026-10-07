set auto_zoom, off
reinitialize

# style
bg_color white
set antialias, 2
set ray_shadows, off
set specular, 0
set ambient, 0.45
set direct, 0.55
set reflect, 0.15
set ray_trace_mode, 1

# pastel color scheme
set_color pastel_blue, [0.45, 0.68, 0.86]
set_color pastel_orange, [0.90, 0.62, 0.38]
set_color pastel_green, [0.48, 0.76, 0.56]
set_color pastel_purple, [0.70, 0.55, 0.82]
set_color pastel_gold, [0.88, 0.78, 0.38]

# Alpha helix, they seem rare
load ../test_2/macrocycle_monomer_10K/macrocycle_monomer_10K_macrocycle_12_0_model_131.cif.gz

# Beta Sheet
load ../test_2/macrocycle_monomer_10K/macrocycle_monomer_10K_macrocycle_10_0_model_3.cif.gz

# Loop 
load ../test_2/macrocycle_monomer_10K/macrocycle_monomer_10K_macrocycle_10_0_model_2.cif.gz


color pastel_blue,   (macrocycle_monomer_10K_macrocycle_12_0_model_131 and chain A and elem C)
color pastel_orange, (macrocycle_monomer_10K_macrocycle_10_0_model_3 and chain A and elem C)
# color pastel_green,  (tfr1_macrocycle_209-212_2 and chain B and elem C)
# color pastel_purple, (tfr1_macrocycle_209-212_3 and chain B and elem C)
color pastel_gold,    (macrocycle_monomer_10K_macrocycle_10_0_model_2 and chain A and elem C)

disable 
enable macrocycle_monomer_10K_macrocycle_12_0_model_131
set_view (\
     0.925132692,    0.309741616,    0.219526529,\
     0.370566338,   -0.611012399,   -0.699532866,\
    -0.082540944,    0.728508532,   -0.680047035,\
     0.000000000,    0.000000000,  -41.698249817,\
     3.817532539,   -4.354623795,    6.197949886,\
  -5378.455078125, 5461.851074219,  -20.000000000 )
viewport 1000,1000
# ray 1000,1000
# png monomer_alpha_helix.png, 1000, 1000, ray=1

disable macrocycle_monomer_10K_macrocycle_12_0_model_131
enable macrocycle_monomer_10K_macrocycle_10_0_model_3 
set_view (\
    -0.681706011,   -0.722375214,    0.115987375,\
     0.661128998,   -0.676126540,   -0.325215518,\
     0.313350290,   -0.145021260,    0.938498676,\
     0.000014905,   -0.000001421,  -38.639804840,\
    -6.739200115,    1.816856384,    3.085157394,\
    29.183401108,   48.095882416,  -20.000000000 )


png monomer_beta_strand.png, 1000, 1000, ray=1

disable 
enable macrocycle_monomer_10K_macrocycle_10_0_model_2
set_view (\
    -0.834045112,   -0.550502896,   -0.036285043,\
     0.351182759,   -0.580488563,    0.734648645,\
    -0.425487757,    0.599983335,    0.677480340,\
     0.000000000,   -0.000000000,  -47.153404236,\
    -0.150466919,    8.547145844,   -4.510023594,\
    37.176101685,   57.130706787,  -20.000000000 )

#png monomer_loop.png, 1000, 1000, ray=1