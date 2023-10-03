t06_L5_theta40-50_twobeam_mesh4-50_mu1_lambda3_E100_nu03_TopFree_Disp.txt
Displacement field of stable state S2 from ABAQUS. The first column is nodal number, the second columm is displacement
along x direction, the third column is displacement along y direction.

t06_L5_theta40-50_twobeam_mesh4-50_indent5_n200_Quasi_mu1_lambda3_E100_nu03.inp
Input file generated from ABAQUS. This file is used to input nodal and element information in finite element code.

t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_stressfree_topfree_Re1e-5.txt
Nodal coordinates of stable state S1 (stress-free). The format is [x1;y1;x2;y2;...]. The sub-index is nodal index.

t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_shape_topfree_Re1e-5.txt
Nodal coordinates of stable state S2 (self-stressed). The format is [x1;y1;x2;y2;...]. The sub-index is nodal index.

t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_NEB_N=11_EngDes_BITSS_kk=4e-3_Re1e-5.txt
NEB results with 11 images (excluding two end states). 11 columns in total. Each column represents nodal coordinates of one image.
The format follows [x1;y1;x2;y2;...]. The sub-index is nodal index.

t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_NEB_N=5_EngDes_BITSS_kk=1e-3_Re1e-5.txt
NEB results with 5 images (excluding two end states). 5 columns in total. Each column represents nodal coordinates of one image.
The format follows [x1;y1;x2;y2;...]. The sub-index is nodal index.

t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_NEB_N=5_EngDes_BITSS_kk=1e-3_Interm_Re1e-5.txt
Iterative NEB results with 5 images (excluding two end states). Each row represents one image's nodal coordinate. Every 5 rows represents one iterative
step. The format follows [x1,y1,x2,y2,...]. The sub-index is nodal index.

t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_EnergyLandscape_N=11_Re1e-5.txt
This file summarizes information of y-coordinate of top edge's center (first column), tilting angle (second column) and total energy (third column). 
For top to bottom, each row represents S1, S2, first image in BITSS, second image in BITSS, saddle point (average of two images in BITSS) and 11 
images in NEB (excluding two end images).

t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_EnergyDescent_dist005-005_BITSS1_Re1e-5.txt
Iterative results of energy descent method that is initiated at the first image from BITSS method. The data is combined in a single column.
For each iterative step, it records the nodal coordinates whose format follows [x1;y1;x2;y2;,,,]. The sub-index is nodal index.

t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS2_alpha10_beta01_dist005-005_iter3_Re1e-5.txt
Nodal coordinates of the second image from BITSS method. The format is [x1;y1;x2;y2;...]. The sub-index is nodal index.

t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS1_alpha10_beta01_dist005-005_iter3_Re1e-5.txt
Nodal coordinates of the first image from BITSS method. The format is [x1;y1;x2;y2;...]. The sub-index is nodal index.

t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS_pt2_alpha10_beta01_dist005-005_iter3_Interm_Re1e-5.txt
Iterative results of the second image's nodal coordinates during BITSS algorithm. All iterative nodal coordinates are combined in a single column.
For top to bottom, the iterative index increases. For each iterative result, the format follows [x1;y1;x2;y2;...]. The sub-index is nodal index.

t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS_pt2_alpha10_beta01_dist005-005_iter3_Interm_Energy_Re1e-5.txt
The iterative results of the second image's y-coordinate of the top edge's center (first column), the tilting angle (second column) and total energy 
(third column) during BITSS method. Each row represents one iterative step. 

t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS_pt1_alpha10_beta01_dist005-005_iter3_Interm_Re1e-5.txt
Iterative results of the first image's nodal coordinates during BITSS algorithm. All iterative nodal coordinates are combined in a single column.
For top to bottom, the iterative index increases. For each iterative result, the format follows [x1;y1;x2;y2;...]. The sub-index is nodal index.

t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS_pt1_alpha10_beta01_dist005-005_iter3_Interm_Energy_Re1e-5.txt
The iterative results of the first image's y-coordinate of the top edge's center (first column), the tilting angle (second column) and total energy 
(third column) during BITSS method. Each row represents one iterative step. 

t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS_alpha10_beta01_dist005-005_iter3_springstiff_Re1e-5.txt
Iterative results of the constrained distance between two images (first column), the actual distance between two images (second column), 
spring stiffness ke (third column), spring stiffness kd (fourth column), total energy of the first image (fifth column) and the total energy of the
second image (sixth column) during BITSS method. Each row represents one iterative step. 







