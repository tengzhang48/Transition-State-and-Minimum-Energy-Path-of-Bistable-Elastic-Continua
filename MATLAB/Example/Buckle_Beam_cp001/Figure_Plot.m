%%% Figure 3a
% clear;                                                 
% abaqusfile = 'Data_Files\BuckleBeam_y1_Left_mesh4-50.inp';
% [coord,connect] = inp2mat(abaqusfile);
% nnode = size(coord,1); % total node number
% 
% xy_1 = load('Data_Files\BuckleBeam_y1_Left_mesh4-50.txt'); % first stable state
% xy_2 = load('Data_Files\BuckleBeam_y1_Right_mesh4-50.txt'); % second stable state
% xy_BITSS1 = load('Data_Files\BuckleBeam_y1_mesh4-50_BITSS_pt1_alpha10_beta01_dist005-005_re1e-8_iter3_dist0.txt'); % first image in BITSS
% xy_BITSS2 = load('Data_Files\BuckleBeam_y1_mesh4-50_BITSS_pt2_alpha10_beta01_dist005-005_re1e-8_iter3_dist0.txt'); % second image in BITSS
% 
% % iterative results of first image in BITSS
% xy_BITSS1_int = load('Data_Files\BuckleBeam_y1_mesh4-50_BITSS_pt1_alpha10_beta01_dist005-005_re1e-8_iter3_Interm_dist0.txt');
% nBITSS1 = size(xy_BITSS1_int,1)./2./nnode; % number of iterative results for first image
% xy_BITSS1_int = reshape(xy_BITSS1_int,2*nnode,nBITSS1);
% 
% % iterative results of second image in BITSS
% xy_BITSS2_int = load('Data_Files\BuckleBeam_y1_mesh4-50_BITSS_pt2_alpha10_beta01_dist005-005_re1e-8_iter3_Interm_dist0.txt');
% nBITSS2 = size(xy_BITSS2_int,1)./2./nnode;
% xy_BITSS2_int = reshape(xy_BITSS2_int,2*nnode,nBITSS2);
% 
% color_BITSS1 = [128,0,128]./255;
% color_BITSS2 = [255,165,0]./255;
% 
% coord_def_1(:,1) = xy_1(1:2:end,1);
% coord_def_1(:,2) = xy_1(2:2:end,1);
% 
% coord_def_2(:,1) = xy_2(1:2:end,1);
% coord_def_2(:,2) = xy_2(2:2:end,1);
% 
% ii=5; % choose the iterative index during BITSS algorithm
% 
% coord_def_BITSS1_int(:,1) = xy_BITSS1_int(1:2:end,ii);
% coord_def_BITSS1_int(:,2) = xy_BITSS1_int(2:2:end,ii);
% 
% coord_def_BITSS2_int(:,1) = xy_BITSS2_int(1:2:end,ii);
% coord_def_BITSS2_int(:,2) = xy_BITSS2_int(2:2:end,ii);
% 
% TR = triangulation(connect,coord_def_1);
% px = [TR.Points(TR.ConnectivityList(:,1),1),TR.Points(TR.ConnectivityList(:,2),1),TR.Points(TR.ConnectivityList(:,3),1)]';
% py = [TR.Points(TR.ConnectivityList(:,1),2),TR.Points(TR.ConnectivityList(:,2),2),TR.Points(TR.ConnectivityList(:,3),2)]';
% s1 = patch(px,py,'r'); % plot first stable state
% s1.FaceAlpha = 1;
% s1.EdgeAlpha = 0;
% hold on;
% 
% TR = triangulation(connect,coord_def_2);
% px = [TR.Points(TR.ConnectivityList(:,1),1),TR.Points(TR.ConnectivityList(:,2),1),TR.Points(TR.ConnectivityList(:,3),1)]';
% py = [TR.Points(TR.ConnectivityList(:,1),2),TR.Points(TR.ConnectivityList(:,2),2),TR.Points(TR.ConnectivityList(:,3),2)]';
% s2 = patch(px,py,'b'); % plot second stable state
% s2.FaceAlpha = 1;
% s2.EdgeAlpha = 0;
% hold on;
% 
% TR = triangulation(connect,coord_def_BITSS1_int);
% px = [TR.Points(TR.ConnectivityList(:,1),1),TR.Points(TR.ConnectivityList(:,2),1),TR.Points(TR.ConnectivityList(:,3),1)]';
% py = [TR.Points(TR.ConnectivityList(:,1),2),TR.Points(TR.ConnectivityList(:,2),2),TR.Points(TR.ConnectivityList(:,3),2)]';
% patch(px,py,color_BITSS1);hold on; % plot iterative result of first image in BITSS
% 
% TR = triangulation(connect,coord_def_BITSS2_int);
% px = [TR.Points(TR.ConnectivityList(:,1),1),TR.Points(TR.ConnectivityList(:,2),1),TR.Points(TR.ConnectivityList(:,3),1)]';
% py = [TR.Points(TR.ConnectivityList(:,1),2),TR.Points(TR.ConnectivityList(:,2),2),TR.Points(TR.ConnectivityList(:,3),2)]';
% patch(px,py,color_BITSS2);hold on; % plot iterative result of second image in BITSS
% 
% axis equal;hold on;axis off;
%%%%%%

%%% Figure 3b and 3c
% clear;                                                 
% abaqusfile = 'Data_Files\BuckleBeam_y1_Left_mesh4-50.inp';
% [coord,connect] = inp2mat(abaqusfile);
% nnode = size(coord,1); % total node number
% middle_node = 56; % node number of beam center
% 
% color_BITSS1 = [128,0,128]./255;
% color_BITSS2 = [255,165,0]./255;
% 
% % input iterative results of BITSS. Format dist0, actual distance, ke, kd, U1 and U2
% BITSS_stiffness = load('Data_Files\BuckleBeam_y1_mesh4-50_BITSS_alpha10_beta01_dist005-005_re1e-8_iter3_springstiff_dist0.txt');
% disti = BITSS_stiffness(:,1); % constrained distance between two images
% dista = BITSS_stiffness(:,2); % actual distance between two images
% ke = BITSS_stiffness(:,3); % spring stiffness ke
% kd = BITSS_stiffness(:,4); % spring stiffness kd
% E1 = BITSS_stiffness(:,5); % energy of first image
% E2 = BITSS_stiffness(:,6); % energy of second image
% 
% xy_BITSS1_int = load('Data_Files\BuckleBeam_y1_mesh4-50_BITSS_pt1_alpha10_beta01_dist005-005_re1e-8_iter3_Interm_dist0.txt');
% nBITSS1 = size(xy_BITSS1_int,1)./2./nnode; % number of iterations in BITSS
% xy_BITSS1_int = reshape(xy_BITSS1_int,2*nnode,nBITSS1); % iterative results of first image
% xy_BITSS2_int = load('Data_Files\BuckleBeam_y1_mesh4-50_BITSS_pt2_alpha10_beta01_dist005-005_re1e-8_iter3_Interm_dist0.txt');
% nBITSS2 = size(xy_BITSS2_int,1)./2./nnode;
% xy_BITSS2_int = reshape(xy_BITSS2_int,2*nnode,nBITSS2); % iterative results of second image
% 
% % input energy landscape, format s1,s2,BITSS1,BITSS2,BITSS average,NEB15
% xx = load('Data_Files\BuckleBeam_y1_mesh4-50_EnergyLandscape_N=15.txt'); 
% 
% % Figure 3c
% y_mep = -xx(1,1):0.1:xx(1,1); % beam center coordinate along MEP
% e_mep = spline([xx(1,1);xx(6:20,1);xx(2,1)],[xx(1,2);xx(6:20,2);xx(2,2)],y_mep); % spline interpolation of MEP 
% plot(y_mep,e_mep,'k','LineWidth',4);hold on; % plot MEP
% scatter(-xy_BITSS1_int(2*middle_node-1,1:5:end),E1(1:5:end,:),200,'MarkerEdgeColor',color_BITSS1,'MarkerFaceColor',color_BITSS1);hold on;
% scatter(-xy_BITSS2_int(2*middle_node-1,1:5:end),E2(1:5:end,:),200,'MarkerEdgeColor',color_BITSS2,'MarkerFaceColor',color_BITSS2);hold on;
% scatter(xx(5,1),xx(5,2),500,'pentagram','MarkerEdgeColor',[0,128,0]./255,'MarkerFaceColor',[0,128,0]./255);hold on;
% scatter(xx(1,1),xx(1,2),200,'MarkerEdgeColor','r','MarkerFaceColor','r');hold on;
% scatter(xx(2,1),xx(2,2),200,'MarkerEdgeColor','b','MarkerFaceColor','b');hold on;
% xlabel('{\it y_{c}} (mm)');ylabel('{\it E} (mJ)');
% %
% 
% % Figure 3b
% yyaxis left
% yaxl = gca;
% % output BITSS iterative data every 3 iterations
% scatter(1:3:nBITSS1,ke(1:3:end,:),200,'MarkerEdgeColor','r','MarkerFaceColor','r');hold on;
% 
% yyaxis right
% yaxr = gca;
% scatter(1:3:nBITSS1,kd(1:3:end,:),200,'MarkerEdgeColor','b','MarkerFaceColor','b');hold on;
% 
% xlabel('Iteration');
% yyaxis left
% ylabel('{\it k_{e}} (mJ^{-1})');
% yyaxis right
% ylabel('{\it k_{d}} (N/mm)');
% 
% ax = gca;
% ax.YAxis(1).Color = 'r';
% ax.YAxis(1).Scale = 'log';
% ax.YAxis(2).Color = 'b';
% ax.YAxis(2).Scale = 'linear';
% ax.FontSize = 40;
% ax.XLim = [0 65];
% ax.FontWeight = 'bold';
% ax.FontName = 'Times New Roman';
% ax.LineWidth = 4;
% %
% 
% set(gca,'FontSize',40,'FontWeight','bold','FontName','Times New Roman','LineWidth',4);
% pbaspect([1,1,1]);
%%%

%%% Figure 3e
% clear;
% abaqusfile = 'Data_Files\BuckleBeam_y1_Left_mesh4-50.inp';
% [coord,connect] = inp2mat(abaqusfile);
% nnode = size(coord,1); % total node number
% xy_1 = load('Data_Files\BuckleBeam_y1_Left_mesh4-50.txt');
% xy_2 = load('Data_Files\BuckleBeam_y1_Right_mesh4-50.txt');
% xy_NEB = load('Data_Files\BuckleBeam_y1_mesh4-50_NEB_N=7_EngDes_BITSS_kk=1e-5.txt');
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
% s1.FaceAlpha = 1;
% s1.EdgeAlpha = 0;
% hold on;
% 
% TR = triangulation(connect,coord_2);
% px = [TR.Points(TR.ConnectivityList(:,1),1),TR.Points(TR.ConnectivityList(:,2),1),TR.Points(TR.ConnectivityList(:,3),1)]';
% py = [TR.Points(TR.ConnectivityList(:,1),2),TR.Points(TR.ConnectivityList(:,2),2),TR.Points(TR.ConnectivityList(:,3),2)]';
% s2 = patch(px,py,'b'); % plot S2
% s2.FaceAlpha = 1;
% s2.EdgeAlpha = 0;
% hold on;
% 
% for i = 7:7 % plot (i+1)-th image in NEB results
%     coord_NEB(:,1) = xy_NEB(1:2:end,1*i);
%     coord_NEB(:,2) = xy_NEB(2:2:end,1*i);
%     TR = triangulation(connect,coord_NEB);
%     px = [TR.Points(TR.ConnectivityList(:,1),1),TR.Points(TR.ConnectivityList(:,2),1),TR.Points(TR.ConnectivityList(:,3),1)]';
%     py = [TR.Points(TR.ConnectivityList(:,1),2),TR.Points(TR.ConnectivityList(:,2),2),TR.Points(TR.ConnectivityList(:,3),2)]';
%     ss = patch(px,py,[255./8.*(8-i),200,255./8.*i]./255);
%     ss.FaceAlpha = 1;
%     ss.EdgeAlpha = 0;
%     hold on;
% end
% 
% axis equal;axis off;
%%%

%%% Figure 3f
% clear;
% xx = load('Data_Files\BuckleBeam_y1_mesh4-50_EnergyLandscape_N=15.txt'); % s1,s2,BITSS1,BITSS2,BITSS average,NEB10
% Indent = load('Data_Files\BuckleBeam_y1_mesh4-50_IndentMiddle_Right-to-Left.txt');
% Indent_color = [255 165 0]./255;
% y_mep = xx(2,1):0.1:xx(1,1);
% e_mep = spline([xx(1,1);xx(6:20,1);xx(2,1)],[xx(1,2);xx(6:20,2);xx(2,2)],y_mep); % spline interpolation of MEP
% 
% plot(y_mep,e_mep,'r','LineWidth',4);hold on;
% plot(Indent(:,1),Indent(:,2),'Color',Indent_color,'LineWidth',4,'LineStyle','--');hold on;
% scatter(xx(1,1),xx(1,2),200,'MarkerFaceColor','r','MarkerEdgeColor','r');hold on;
% scatter(xx(2,1),xx(2,2),200,'MarkerFaceColor','b','MarkerEdgeColor','b');hold on;
% for i = 1:15
%     color_NEB = [255./16.*(16-i),200,255./16.*i]./255;
%     scatter(xx(5+i,1),xx(5+i,2),200,'MarkerFaceColor',color_NEB,'MarkerEdgeColor',color_NEB);hold on;
% end
% scatter(xx(5,1),xx(5,2),500,'pentagram','MarkerFaceColor',[0,128,0]./255,'MarkerEdgeColor',[0,128,0]./255);hold on;
% 
% xlabel('{\it y_{c}} (mm)');ylabel('{\it E} (mJ)');
% set(gca,'FontSize',40,'xlim',[-6 6],'FontWeight','bold','FontName','Times New Roman','LineWidth',4);
% pbaspect([1,1,1]);
%%%

%%% Figure 3d
clear;
abaqusfile = 'Data_Files\BuckleBeam_y1_Left_mesh4-50.inp';
[coord,connect] = inp2mat(abaqusfile);
nnode = size(coord,1); % total node number
xy_1 = load('Data_Files\BuckleBeam_y1_Left_mesh4-50.txt'); % S1
xy_2 = load('Data_Files\BuckleBeam_y1_Right_mesh4-50.txt'); % S2
xy_NEB = load('Data_Files\BuckleBeam_y1_mesh4-50_NEB_N=3_EngDes_BITSS_kk=1e-5.txt'); % NEB results with 5 images
xy_NEB_Interm = load('Data_Files\BuckleBeam_y1_mesh4-50_NEB_N=3_EngDes_BITSS_kk=1e-5_Interm.txt'); % iterative results during NEB algorithm

coord_1(:,1) = xy_1(1:2:end,1);
coord_1(:,2) = xy_1(2:2:end,1);
coord_2(:,1) = xy_2(1:2:end,1);
coord_2(:,2) = xy_2(2:2:end,1);

TR = triangulation(connect,coord_1);
px = [TR.Points(TR.ConnectivityList(:,1),1),TR.Points(TR.ConnectivityList(:,2),1),TR.Points(TR.ConnectivityList(:,3),1)]';
py = [TR.Points(TR.ConnectivityList(:,1),2),TR.Points(TR.ConnectivityList(:,2),2),TR.Points(TR.ConnectivityList(:,3),2)]';
s1 = patch(px,py,'r'); % plot S1
s1.FaceAlpha = 1;
s1.EdgeAlpha = 0;
hold on;

TR = triangulation(connect,coord_2);
px = [TR.Points(TR.ConnectivityList(:,1),1),TR.Points(TR.ConnectivityList(:,2),1),TR.Points(TR.ConnectivityList(:,3),1)]';
py = [TR.Points(TR.ConnectivityList(:,1),2),TR.Points(TR.ConnectivityList(:,2),2),TR.Points(TR.ConnectivityList(:,3),2)]';
s2 = patch(px,py,'b'); % plot S2
s2.FaceAlpha = 1;
s2.EdgeAlpha = 0;
hold on;

for i = 1:3 % plot NEB results with 5 images
    coord_NEB(:,1) = xy_NEB(1:2:end,i);
    coord_NEB(:,2) = xy_NEB(2:2:end,i);
    TR = triangulation(connect,coord_NEB);
    px = [TR.Points(TR.ConnectivityList(:,1),1),TR.Points(TR.ConnectivityList(:,2),1),TR.Points(TR.ConnectivityList(:,3),1)]';
    py = [TR.Points(TR.ConnectivityList(:,1),2),TR.Points(TR.ConnectivityList(:,2),2),TR.Points(TR.ConnectivityList(:,3),2)]';
    ss = patch(px,py,[255./4.*(4-i),200,255./4.*i]./255);
    ss.FaceAlpha = 1;
    ss.EdgeAlpha = 1;
    hold on;
end

ii = 4; % choose iterative index
ii_total = size(xy_NEB_Interm,1)./3;
for i = 1:3
    coord_NEB_Interm(:,1) = xy_NEB_Interm(i*ii,1:2:end)';
    coord_NEB_Interm(:,2) = xy_NEB_Interm(i*ii,2:2:end)';
    TR = triangulation(connect,coord_NEB_Interm);
    px = [TR.Points(TR.ConnectivityList(:,1),1),TR.Points(TR.ConnectivityList(:,2),1),TR.Points(TR.ConnectivityList(:,3),1)]';
    py = [TR.Points(TR.ConnectivityList(:,1),2),TR.Points(TR.ConnectivityList(:,2),2),TR.Points(TR.ConnectivityList(:,3),2)]';
    ss = patch(px,py,[255./4.*(4-i),200,255./4.*i]./255);
    ss.FaceAlpha = 1;
    ss.EdgeAlpha = 1;
    hold on;
end

daspect([1 1 1]);
axis off;
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