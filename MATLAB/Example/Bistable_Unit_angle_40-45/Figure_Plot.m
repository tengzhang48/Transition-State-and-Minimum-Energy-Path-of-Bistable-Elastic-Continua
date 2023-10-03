%%% Figure 8a and 8b
% clear;                                                 
% abaqusfile = 'Data_Files\t06_L5_theta40-50_twobeam_mesh4-50_indent5_n200_Quasi_mu1_lambda3_E100_nu03.inp';
% [coord,connect] = inp2mat(abaqusfile);
% nnode = size(coord,1); % total node number
% nelem = size(connect,1); % total element number
% type1 = [1:400;861:1260];
% type2 = 401:860;
% 
% % Input Files
% % S1
% xy_1 = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_stressfree_topfree_Re1e-5.txt');
% % S2 
% xy_2 = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_shape_topfree_Re1e-5.txt');
% % first image in BITSS
% xy_BITSS1 = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS1_alpha10_beta01_dist005-005_iter3_Re1e-5.txt');
% % second image in BITSS
% xy_BITSS2 = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS2_alpha10_beta01_dist005-005_iter3_Re1e-5.txt');
% 
% coord_def_1(:,1) = xy_1(1:2:end,1);
% coord_def_1(:,2) = xy_1(2:2:end,1);
% 
% coord_def_2(:,1) = xy_2(1:2:end,1);
% coord_def_2(:,2) = xy_2(2:2:end,1);
% 
% coord_def_BITSS1(:,1) = xy_BITSS1(1:2:end,:);
% coord_def_BITSS1(:,2) = xy_BITSS1(2:2:end,:);
% 
% coord_def_BITSS2(:,1) = xy_BITSS2(1:2:end,:);
% coord_def_BITSS2(:,2) = xy_BITSS2(2:2:end,:); 
% 
% coord_def_saddle = (coord_def_BITSS1+coord_def_BITSS2)./2;
% 
% TR = triangulation(connect,coord_def_1);
% px = [TR.Points(TR.ConnectivityList(:,1),1),TR.Points(TR.ConnectivityList(:,2),1),TR.Points(TR.ConnectivityList(:,3),1)]';
% py = [TR.Points(TR.ConnectivityList(:,1),2),TR.Points(TR.ConnectivityList(:,2),2),TR.Points(TR.ConnectivityList(:,3),2)]';
% s1 = patch(px,py,'r'); % plot S1
% s1.FaceAlpha = 1;
% s1.EdgeAlpha = 1;
% hold on;
% 
% TR = triangulation(connect,coord_def_2);
% px = [TR.Points(TR.ConnectivityList(:,1),1),TR.Points(TR.ConnectivityList(:,2),1),TR.Points(TR.ConnectivityList(:,3),1)]';
% py = [TR.Points(TR.ConnectivityList(:,1),2),TR.Points(TR.ConnectivityList(:,2),2),TR.Points(TR.ConnectivityList(:,3),2)]';
% s2 = patch(px,py,'b'); % plot S2
% s2.FaceAlpha = 1;
% s2.EdgeAlpha = 1;
% hold on;
% 
% TR = triangulation(connect,coord_def_saddle);
% px = [TR.Points(TR.ConnectivityList(:,1),1),TR.Points(TR.ConnectivityList(:,2),1),TR.Points(TR.ConnectivityList(:,3),1)]';
% py = [TR.Points(TR.ConnectivityList(:,1),2),TR.Points(TR.ConnectivityList(:,2),2),TR.Points(TR.ConnectivityList(:,3),2)]';
% s2 = patch(px,py,[0,128,0]./255); % plot saddle point
% s2.FaceAlpha = 1;
% s2.EdgeAlpha = 1;
% hold on;
% 
% axis equal;hold on;axis off;
%%%

%%% Figure 8d
% clear;
% xx = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_EnergyLandscape_N=11_Re1e-5.txt'); % s1,s2,BITSS1,BITSS2,BITSS average,NEB10
% 
% color_BITSS1 = [128,0,128]./255;
% color_BITSS2 = [255,165,0]./255;
% 
% EBITSS1 = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS_pt1_alpha10_beta01_dist005-005_iter3_Interm_Energy_Re1e-5.txt');
% EBITSS2 = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS_pt2_alpha10_beta01_dist005-005_iter3_Interm_Energy_Re1e-5.txt');
% nBITSS = size(EBITSS1,1);
% 
% neb = [xx(1,1:2:3);xx(6:16,1:2:3);xx(2,1:2:3)]';
% fnplt(cscvn(neb),'k',4);hold on;
% 
% scatter(EBITSS1(1:3:50,1),EBITSS1(1:3:50,3),400,'MarkerFaceColor',color_BITSS1,'MarkerEdgeColor',color_BITSS1);hold on;
% scatter(EBITSS2(1:3:50,1),EBITSS2(1:3:50,3),400,'MarkerFaceColor',color_BITSS2,'MarkerEdgeColor',color_BITSS2);hold on;
% 
% scatter(xx(1,1),xx(1,3),400,'MarkerFaceColor','r','MarkerEdgeColor','r');hold on;
% scatter(xx(2,1),xx(2,3),400,'MarkerFaceColor','b','MarkerEdgeColor','b');hold on;
% color_BITSS = [0,128,0]./255;
% scatter(xx(5,1),xx(5,3),600,'pentagram','MarkerFaceColor',color_BITSS,'MarkerEdgeColor',color_BITSS);hold on;
% 
% for i = 1:11
%     color_NEB = [255./12.*(12-i),200,255./12.*i]./255;
%     scatter(xx(5+i,1),xx(5+i,3),400,'MarkerFaceColor',color_NEB,'MarkerEdgeColor',color_NEB);hold on;
% end
% 
% xlabel('{\it y_{m}} (mm)');
% ylabel('{\it E} (mJ)');
% set(gca,'FontSize',60,'FontWeight','bold','FontName','Times New Roman','LineWidth',4);
% pbaspect([1,1,1]);
%%%

%%% Figure 8d
% clear;
% xx = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_EnergyLandscape_N=11_Re1e-5.txt'); % s1,s2,BITSS1,BITSS2,BITSS average,NEB10
% % format y coordinate; tilting angle; total energy
% for i = 1:11
%     color_NEB = [255./12.*(12-i),200,255./12.*i]./255;
%     scatter(1+i,xx(5+i,3),400,'MarkerFaceColor',color_NEB,'MarkerEdgeColor',color_NEB);hold on;
% end
% 
% scatter(1,xx(1,3),400,'MarkerFaceColor','r','MarkerEdgeColor','r');hold on;
% scatter(13,xx(2,3),400,'MarkerFaceColor','b','MarkerEdgeColor','b');hold on;
% 
% xlabel('Image label');
% ylabel('{\it E} (mJ)');
% set(gca,'FontSize',60,'xlim',[0 14],'ylim',[0 0.35],'FontWeight','bold','FontName','Times New Roman','LineWidth',4);
% pbaspect([1,1,1]);
%%%

%%% Figure 8c
% clear;
% abaqusfile = 'Data_Files\t06_L5_theta40-50_twobeam_mesh4-50_indent5_n200_Quasi_mu1_lambda3_E100_nu03.inp';
% [coord,connect] = inp2mat(abaqusfile);
% nnode = size(coord,1); % total node number
% nelem = size(connect,1); % total element number
% 
% xy_1 = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_stressfree_topfree_Re1e-5.txt');
% xy_2 = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_shape_topfree_Re1e-5.txt');
% xy_NEB = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_NEB_N=5_EngDes_BITSS_kk=1e-3_Re1e-5.txt');
% 
% coord_1(:,1) = xy_1(1:2:end,1);
% coord_1(:,2) = xy_1(2:2:end,1);
% coord_2(:,1) = xy_2(1:2:end,1);
% coord_2(:,2) = xy_2(2:2:end,1);
% 
% TR = triangulation(connect,coord_1);
% px = [TR.Points(TR.ConnectivityList(:,1),1),TR.Points(TR.ConnectivityList(:,2),1),TR.Points(TR.ConnectivityList(:,3),1)]';
% py = [TR.Points(TR.ConnectivityList(:,1),2),TR.Points(TR.ConnectivityList(:,2),2),TR.Points(TR.ConnectivityList(:,3),2)]';
% s1 = patch(px,py,'r'); % plot S1
% s1.FaceAlpha = 0.5;
% s1.EdgeAlpha = 0;
% hold on;
% 
% TR = triangulation(connect,coord_2);
% px = [TR.Points(TR.ConnectivityList(:,1),1),TR.Points(TR.ConnectivityList(:,2),1),TR.Points(TR.ConnectivityList(:,3),1)]';
% py = [TR.Points(TR.ConnectivityList(:,1),2),TR.Points(TR.ConnectivityList(:,2),2),TR.Points(TR.ConnectivityList(:,3),2)]';
% s2 = patch(px,py,'b'); % plot S2
% s2.FaceAlpha = 0.5;
% s2.EdgeAlpha = 0;
% hold on;
% 
% for i = 1:5
%     coord_NEB(:,1) = xy_NEB(1:2:end,i);
%     coord_NEB(:,2) = xy_NEB(2:2:end,i);
%     TR = triangulation(connect,coord_NEB);
%     px = [TR.Points(TR.ConnectivityList(:,1),1),TR.Points(TR.ConnectivityList(:,2),1),TR.Points(TR.ConnectivityList(:,3),1)]';
%     py = [TR.Points(TR.ConnectivityList(:,1),2),TR.Points(TR.ConnectivityList(:,2),2),TR.Points(TR.ConnectivityList(:,3),2)]';
%     ss = patch(px,py,[255./6.*(6-i),200,255./6.*i]./255);
%     ss.FaceAlpha = 0.8;
%     ss.EdgeAlpha = 0;
%     hold on;
% end
% 
% daspect([1 1 1]);
% axis off;
%%%

function [coord,connect] = inp2mat(file)

% Load inp file
fid = fopen(file);
for i = 1:9
    fgetl(fid);
end
data = (fscanf(fid,'%i, %f, %f \n',[3 inf]))';
coord = data(:,2:3);
fgetl(fid);
data = (fscanf(fid,'%i, %i, %i, %i \n',[4 inf]))';
connect = data(:,2:4);
fclose(fid);
end