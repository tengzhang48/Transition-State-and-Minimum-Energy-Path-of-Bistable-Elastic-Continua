%%% Post processing results
% organize it in the format of traditional fem
% such as x = [x_1^1,x_2^1,x_1^2,x_2^2,..., x_1^n, x_2^n]
%         f = [f_1^1,f_2^1,f_1^2,f_2^2,..., f_1^n, f_2^n]
% Solve the nonlinear equation kdx = f

clear;
clc

% initial coordinate and element connection
% Read abaqus node and element information

abaqusfile = 'Data_Files\t06_L5_theta40-50_twobeam_mesh4-50_indent5_n200_Quasi_mu1_lambda3_E100_nu03.inp';
[coord,connect] = inp2mat(abaqusfile);
nnode = size(coord,1);   % total node number
nelem = size(connect,1); % total element number
mu_1 = 1; lambda_1 = 3; type1 = [1:400;861:1260];
mu_2 = 38.5; lambda_2 = 57.7; type2 = 401:860;
mu = zeros(nelem,1);     % material properties are assigned to each element
lambda = zeros(nelem,1);
mu(type1) = mu_1; lambda(type1) = lambda_1;
mu(type2) = mu_2; lambda(type2) = lambda_2;

%%% Input results files
% stable states
xy_1 = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_stressfree_topfree_Re1e-5.txt');
xy_2 = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_shape_topfree_Re1e-5.txt');
% two images in BITSS
xy_BITSS1 = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS1_alpha10_beta01_dist005-005_iter3_Re1e-5.txt');
xy_BITSS2 = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS2_alpha10_beta01_dist005-005_iter3_Re1e-5.txt');
% iterative results in BITSS
xy_BITSS1_Interm = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS_pt1_alpha10_beta01_dist005-005_iter3_Interm_Re1e-5.txt');
nBITSS1 = size(xy_BITSS1_Interm,1)./2./nnode; % number of iterations
xy_BITSS1_Interm = reshape(xy_BITSS1_Interm,2*nnode,nBITSS1);
xy_BITSS2_Interm = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS_pt2_alpha10_beta01_dist005-005_iter3_Interm_Re1e-5.txt');
nBITSS2 = size(xy_BITSS2_Interm,1)./2./nnode;
xy_BITSS2_Interm = reshape(xy_BITSS2_Interm,2*nnode,nBITSS2);
% NEB results
xy_NEB = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_NEB_N=11_EngDes_BITSS_kk=4e-3_Re1e-5.txt');

total_node = [1:nnode]';                        % total node

%%% set up boundary conditions
l = abs(coord(:,2)-6.2) < 1e-5; % top edge
top_node = total_node (l);

l = abs(coord(:,2)-0) < 1e-3; % bottom edge
bottom_node = total_node (l);

disp_fix_1 = [bottom_node,ones(length(bottom_node),1),zeros(length(bottom_node),1)]; % fix the bottom edge in x direction
disp_fix_2 = [bottom_node,2*ones(length(bottom_node),1),zeros(length(bottom_node),1)]; % fix the bottom edge in y direction 
disp_fix = [disp_fix_1;disp_fix_2];
%%%

Indent_node1 = 151;
Indent_node2 = 152; % indent nodes are used to measure y coordinate of top edge
Angle_node1 = 7;
Angle_node2 = 8; % two angle nodes are used to measure tilting angle phi of top edge
xy_bb = (xy_BITSS1+xy_BITSS2)./2; % average of two images in BITSS, taken as saddle point

%%% calculate y coordinate, tilting angle and total energy of stable
%%% states, BITSS and NEB results
% S1
y1 = (xy_1(2*Indent_node1) + xy_1(2*Indent_node2))./2;
t1 = atan((xy_1(2*Angle_node2)-xy_1(2*Angle_node1))./(xy_1(2*Angle_node2-1)-xy_1(2*Angle_node1-1)));
E1 = total_energy(xy_1,coord,connect,mu,lambda);
% S2
y2 = (xy_2(2*Indent_node1) + xy_2(2*Indent_node2))./2;
t2 = atan((xy_2(2*Angle_node2)-xy_2(2*Angle_node1))./(xy_2(2*Angle_node2-1)-xy_2(2*Angle_node1-1)));
E2 = total_energy(xy_2,coord,connect,mu,lambda);
% first image in BITSS
yb1 = (xy_BITSS1(2*Indent_node1) + xy_BITSS1(2*Indent_node2))./2;
tb1 = atan((xy_BITSS1(2*Angle_node2)-xy_BITSS1(2*Angle_node1))./(xy_BITSS1(2*Angle_node2-1)-xy_BITSS1(2*Angle_node1-1)));
Eb1 = total_energy(xy_BITSS1,coord,connect,mu,lambda);
% second image in BITSS
yb2 = (xy_BITSS2(2*Indent_node1) + xy_BITSS2(2*Indent_node2))./2;
tb2 = atan((xy_BITSS2(2*Angle_node2)-xy_BITSS2(2*Angle_node1))./(xy_BITSS2(2*Angle_node2-1)-xy_BITSS2(2*Angle_node1-1)));
Eb2 = total_energy(xy_BITSS2,coord,connect,mu,lambda);
% saddle point
ybb = (yb1 + yb2)./2;
tbb = atan((xy_bb(2*Angle_node2)-xy_bb(2*Angle_node1))./(xy_bb(2*Angle_node2-1)-xy_bb(2*Angle_node1-1)));
Ebb = total_energy(xy_bb,coord,connect,mu,lambda);
% NEB
N = 11;
yNEB = zeros(N,1);ENEB = zeros(N,1);tNEB = zeros(N,1);
for ii = 1:N
    yNEB(ii) = (xy_NEB(2*Indent_node1,ii) + xy_NEB(2*Indent_node2,ii))./2;
    tNEB(ii) = atan((xy_NEB(2*Angle_node2,ii)-xy_NEB(2*Angle_node1,ii))./(xy_NEB(2*Angle_node2-1,ii)-xy_NEB(2*Angle_node1-1,ii)));
    ENEB(ii) = total_energy(xy_NEB(:,ii),coord,connect,mu,lambda);
end
% iterative BITSS results
yBITSS1 = zeros(nBITSS1,1); EBITSS1 = zeros(nBITSS1,1); tBITSS1 = zeros(nBITSS1,1);
yBITSS2 = zeros(nBITSS2,1); EBITSS2 = zeros(nBITSS2,1); tBITSS2 = zeros(nBITSS2,1);
for jj = 1:nBITSS1
    yBITSS1(jj) = (xy_BITSS1_Interm(2*Indent_node1,jj)+xy_BITSS1_Interm(2*Indent_node2,jj))./2;
    yBITSS2(jj) = (xy_BITSS2_Interm(2*Indent_node1,jj)+xy_BITSS2_Interm(2*Indent_node2,jj))./2;
    tBITSS1(jj) = atan((xy_BITSS1_Interm(2*Angle_node2,jj)-xy_BITSS1_Interm(2*Angle_node1,jj))./...
        (xy_BITSS1_Interm(2*Angle_node2-1,jj)-xy_BITSS1_Interm(2*Angle_node1-1,jj)));
    tBITSS2(jj) = atan((xy_BITSS2_Interm(2*Angle_node2,jj)-xy_BITSS2_Interm(2*Angle_node1,jj))./...
        (xy_BITSS2_Interm(2*Angle_node2-1,jj)-xy_BITSS2_Interm(2*Angle_node1-1,jj)));
    EBITSS1(jj) = total_energy(xy_BITSS1_Interm(:,jj),coord,connect,mu,lambda);
    EBITSS2(jj) = total_energy(xy_BITSS2_Interm(:,jj),coord,connect,mu,lambda);
end
%%%

%%% File Output
filename = 't06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_EnergyLandscape_N=11_Re1e-5.txt';
fid=fopen(filename,'wt');
fprintf(fid,'%.12f,%.12f,%.12f\n',y1,t1,E1);
fprintf(fid,'%.12f,%.12f,%.12f\n',y2,t2,E2);
fprintf(fid,'%.12f,%.12f,%.12f\n',yb1,tb1,Eb1);
fprintf(fid,'%.12f,%.12f,%.12f\n',yb2,tb2,Eb2);
fprintf(fid,'%.12f,%.12f,%.12f\n',ybb,tbb,Ebb);
for ii = 1:11
    fprintf(fid,'%.12f,%.12f,%.12f\n',yNEB(ii),tNEB(ii),ENEB(ii));
end
fclose(fid);
%%%

%%% BITSS iterative results output
filename = 't06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS_pt1_alpha10_beta01_dist005-005_iter3_Interm_Energy_Re1e-5.txt';
fid = fopen(filename,'wt');
for ii = 1:nBITSS1
    fprintf(fid,'%.12f,%.12f,%.12f\n',yBITSS1(ii),tBITSS1(ii),EBITSS1(ii));
end
fclose(fid);

filename = 't06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS_pt2_alpha10_beta01_dist005-005_iter3_Interm_Energy_Re1e-5.txt';
fid = fopen(filename,'wt');
for ii = 1:nBITSS2
    fprintf(fid,'%.12f,%.12f,%.12f\n',yBITSS2(ii),tBITSS2(ii),EBITSS2(ii));
end
fclose(fid);
%%%

function strain_energy = total_energy(xy,coord,connect,mu,lambda)

force_true=0;
hessian_true=0; % no output of force and hessian matrix

% Define stored variables
Ne= size(connect,1);                            % total element number
nnode = size(coord,1);                          % total node number
strain_energy = 0.0;                            % Total strain energy

% fxy = zeros(2*nnode,1);                         % Force vector
% kxy = zeros(2*nnode,2*nnode);                   % Stiffness matrix
% kxy_s = sparse(2*nnode,2*nnode);                   % Stiffness sparse matrix

% This format is to use parallel simulations in matlab
x_ev = ones(Ne,3);                              % deformed x vector
y_ev = ones(Ne,3);                              % Current y vector
coordx_ev = zeros(Ne,3);                        % Undeformed x vector
coordy_ev = zeros(Ne,3);

for ii = 1:3
    x_ev(:,ii) = xy(2*connect(:,ii)-1,1);
    y_ev(:,ii) = xy(2*connect(:,ii),1);
    coordx_ev(:,ii) = coord(connect(:,ii),1);
    coordy_ev(:,ii) = coord(connect(:,ii),2);
end

parfor ii = 1:Ne
    
    xa = x_ev(ii,:);
    ya = y_ev(ii,:);
    % undeformed coordinate
    coordx = coordx_ev(ii,:);
    coordy = coordy_ev(ii,:);
    % element material properties
    mu_e=mu(ii,1);
    lambda_e=lambda(ii,1);
    % energy and forces
    [s_e,~,~] = element_energy_force_stiffness(xa,ya,coordx,coordy,mu_e,lambda_e,force_true,hessian_true);
    % assemble energy and forces
    strain_energy = strain_energy + s_e;
 
end
end

function fxy = total_force_boundary(xy,coord,connect,mu,lambda,disp_fix,force_applied,force_gradient)

% force=-d(strain_energy)/dx
% gradient=d(strain_energy)/dx
% The default output is force

force_true=1; % output of force
hessian_true=0; % no output of hessian matrix

[~,fxy,~,~] = total_energy_force_stiffness(xy,coord,connect,mu,lambda,force_true,hessian_true);
% Set the force equal to 0 for node with prescribed displacement
% find the force index
force_ind = 2*(disp_fix(:,1)-1)+disp_fix(:,2);
fxy(force_ind,1) = 0;


% Add the external forces
if ~isempty(force_applied)
    force_ind = 2*(force_applied(:,1)-1)+force_applied(:,2);
    fxy(force_ind,1) = fxy(force_ind,1) + force_applied(:,3);
end

if force_gradient<1
    % the output is gradient of the energy
    fxy = -fxy;
end

end

function kxy = total_stiffness_boundary(xy,coord,connect,mu,lambda,disp_fix)

% force=-d(strain_energy)/dx
% gradient=d(strain_energy)/dx
% The default output is force

force_true=1;
hessian_true=1;

[~,fxy,kxy,~] = total_energy_force_stiffness(xy,coord,connect,mu,lambda,force_true,hessian_true);
% Set the force equal to 0 for node with prescribed displacement
% find the force index
force_ind = 2*(disp_fix(:,1)-1)+disp_fix(:,2);

% Update the stiffness matrix
kxy(force_ind,:) = 0;
kxy(:,force_ind) = 0;

for ii = 1:length(force_ind)
    kxy(force_ind(ii),force_ind(ii)) = 1;
end
end

function [fxy,kxy] = total_force_stiffness_boundary(xy,coord,connect,mu,lambda,disp_fix,force_applied,force_gradient)

% force=-d(strain_energy)/dx
% gradient=d(strain_energy)/dx
% The default output is force

force_true=1;
hessian_true=1;

[~,fxy,kxy,~] = total_energy_force_stiffness(xy,coord,connect,mu,lambda,force_true,hessian_true);
% Set the force equal to 0 for node with prescribed displacement
% find the force index
force_ind = 2*(disp_fix(:,1)-1)+disp_fix(:,2);
fxy(force_ind,1) = 0;

% Update the stiffness matrix
kxy(force_ind,:) = 0;
kxy(:,force_ind) = 0;

for ii = 1:length(force_ind)
    kxy(force_ind(ii),force_ind(ii)) = 1;
end

% Add the external forces
if ~isempty(force_applied)
    force_ind = 2*(force_applied(:,1)-1)+force_applied(:,2);
    fxy(force_ind,1) = fxy(force_ind,1) + force_applied(:,3);
end

if force_gradient<1
    % the output is gradient of the energy
    fxy = -fxy;
end

end


function [strain_energy,fxy_s,kxy_s,strain_energy_element] = total_energy_force_stiffness(xy,coord,connect,mu,lambda,force_true,hessian_true)

% Define stored variables
Ne= size(connect,1);                            % total element number
nnode = size(coord,1);                          % total node number
strain_energy = 0.0;                            % Total strain energy
strain_energy_element = zeros(Ne,1);            % Strain energy in each element

% fxy = zeros(2*nnode,1);                         % Force vector
% kxy = zeros(2*nnode,2*nnode);                   % Stiffness matrix
% kxy_s = sparse(2*nnode,2*nnode);                   % Stiffness sparse matrix

% This format is to use parallel simulations in matlab
fxy_ev = zeros(Ne,6);                           % Force x and y vector
x_ev = ones(Ne,3);                              % deformed x vector
y_ev = ones(Ne,3);                              % Current y vector
coordx_ev = zeros(Ne,3);                        % Undeformed x vector
coordy_ev = zeros(Ne,3);
kxy_ev = zeros(Ne,36);                          % element stiffness matrix

for ii = 1:3
    x_ev(:,ii) = xy(2*connect(:,ii)-1,1);
    y_ev(:,ii) = xy(2*connect(:,ii),1);
    coordx_ev(:,ii) = coord(connect(:,ii),1);
    coordy_ev(:,ii) = coord(connect(:,ii),2);
end

parfor ii = 1:Ne
    
    xa = x_ev(ii,:);
    ya = y_ev(ii,:);
    % undeformed coordinate
    coordx = coordx_ev(ii,:);
    coordy = coordy_ev(ii,:);
    % element material properties
    mu_e=mu(ii,1);
    lambda_e=lambda(ii,1);
    % energy and forces
    [s_e,fxy_e,kxy_e] = element_energy_force_stiffness(xa,ya,coordx,coordy,mu_e,lambda_e,force_true,hessian_true);
    % assemble energy and forces
    strain_energy = strain_energy + s_e;
    if (force_true>0)
        fxy_ev(ii,:) = fxy_e;
    end
    % assemble stiffness matrix
    if (hessian_true>0)
        kxy_ev(ii,:) = kxy_e(:);
    end
    % Store the energy in each element
    strain_energy_element(ii) = s_e;
    
end

% assign the element force to global force and global stiffness using
% sparse matrix
fxy_s = sparse(2*nnode,1);
kxy_s = sparse(2*nnode,2*nnode);

node_map_f = zeros(Ne,6);
if (force_true>0&&hessian_true<1)
    for ii = 1:3
        node_ii = connect(:,ii);
        node_map_f(:,2*ii-1:2*ii) = [2*node_ii-1,2*node_ii];
    end
    node_map_f = node_map_f(:);
    fxy_s = sparse(node_map_f,ones(length(node_map_f),1),fxy_ev(:),2*nnode,1);
end
    
if (force_true>0&&hessian_true>0)
    node_map_kc = zeros(Ne,6,6);
    node_map_kr = zeros(Ne,6,6);
    for ii=1:3
        node_ii = connect(:,ii);
        node_map_f(:,2*ii-1:2*ii) = [2*node_ii-1,2*node_ii];
        for jj = 1:3
            node_jj = connect(:,jj);
            
            node_map_kc(:,[2*ii-1:2*ii],2*jj-1) = [2*node_ii-1,2*node_ii];
            node_map_kc(:,[2*ii-1:2*ii],2*jj) = [2*node_ii-1,2*node_ii];
            node_map_kr(:,2*ii-1,[2*jj-1:2*jj]) = [2*node_jj-1,2*node_jj];
            node_map_kr(:,2*ii,[2*jj-1:2*jj]) = [2*node_jj-1,2*node_jj];
            
        end
    end
    node_map_kr = node_map_kr(:,:);
    node_map_kr = node_map_kr(:);
    node_map_kc = node_map_kc(:,:);
    node_map_kc = node_map_kc(:);
    node_map_f = node_map_f(:);
    fxy_s = sparse(node_map_f,ones(length(node_map_f),1),fxy_ev(:),2*nnode,1);
    kxy_s = sparse(node_map_kc,node_map_kr,kxy_ev(:),2*nnode,2*nnode);
end
end

function [strain_energy,fxy,kxy] = element_energy_force_stiffness(xa,ya,coordx,coordy,mu,lambda,force_true,hessian_true)

xc = sum(xa)/3;                 % x (radial) coordinate of the triangle center
coordc = sum(coordx)/3;         % coordc is the x (radial) coordinate of the triangle center
fxy = zeros(6,1);               % f_1^1,f_2^1,f_1,
kxy = zeros(6,6);
% Gaussian interation
Nint = 1;
[xi,w] = Gaussian_point_weight(Nint,1,2);
z1 = xi(1,1);
z2 = xi(2,1);
[dNdx,yita] = shape_function_derivative(z1,z2,coordx,coordy);

% calculate energy
Fxy = [xa;ya]*dNdx; % (F_xy)deformation gradient
Cxy = transpose(Fxy)*Fxy;
I1xy = trace(Cxy);
Jv =det(Fxy);
strain_energy = strain_energy_density(I1xy,Jv,mu,lambda);
strain_energy = strain_energy*yita*w(1);              % consider the volume

if (force_true>0&&hessian_true<1)
    % Using complex step to calculate the force and stiffness matrix numerically
    dh = 1e-5; % Change the value dh if the simulation does not converge
    Fxy_h = zeros(6,4);
    for aa = 1:3
        Fxy_h(2*aa-1,:) = Fxy(:);
        Fxy_h(2*aa,:) = Fxy(:);
        Fxy_h(2*aa-1,1) = Fxy_h(2*aa-1,1) + dNdx(aa,1)*dh*1i;
        Fxy_h(2*aa-1,3) = Fxy_h(2*aa-1,3) + dNdx(aa,2)*dh*1i;
        Fxy_h(2*aa,2) = Fxy_h(2*aa,2) + dNdx(aa,1)*dh*1i;
        Fxy_h(2*aa,4) = Fxy_h(2*aa,4) + dNdx(aa,2)*dh*1i;
    end
    
    for aa = 1:6
        % calculate fx xa(xaaa) = xa(aa)+dh*i
        Fxy_aa = Fxy_h(aa,:);
        Fxy_aa = reshape(Fxy_aa,2,2);
        Cxy_aa = transpose(Fxy_aa)*Fxy_aa;
        I1xy_aa = trace(Cxy_aa);
        Jv_aa = det(Fxy_aa);
        strain_energy_aa = strain_energy_density(I1xy_aa,Jv_aa,mu,lambda);
        strain_energy_aa = strain_energy_aa*yita*w(1);              % consider the volume
        % fx = -du/dx
        fxy(aa) = -imag(strain_energy_aa/dh);
        % kx = d^2u/dx^2
        
    end
end

if(force_true>0&&hessian_true>0)
    % Using complex step to calculate the force and stiffness matrix numerically
    dh = 1e-5; % Change the value dh if the simulation does not converge
    Fxy_h = zeros(6,4);
    for aa = 1:3
        Fxy_h(2*aa-1,:) = Fxy(:);
        Fxy_h(2*aa,:) = Fxy(:);
        Fxy_h(2*aa-1,1) = Fxy_h(2*aa-1,1) + dNdx(aa,1)*dh*1i;
        Fxy_h(2*aa-1,3) = Fxy_h(2*aa-1,3) + dNdx(aa,2)*dh*1i;
        Fxy_h(2*aa,2) = Fxy_h(2*aa,2) + dNdx(aa,1)*dh*1i;
        Fxy_h(2*aa,4) = Fxy_h(2*aa,4) + dNdx(aa,2)*dh*1i;
    end
    
    for aa = 1:6
        % calculate fx xa(xaaa) = xa(aa)+dh*i
        Fxy_aa = Fxy_h(aa,:);
        Fxy_aa = reshape(Fxy_aa,2,2);
        Cxy_aa = transpose(Fxy_aa)*Fxy_aa;
        I1xy_aa = trace(Cxy_aa);
        Jv_aa = det(Fxy_aa);
        strain_energy_aa = strain_energy_density(I1xy_aa,Jv_aa,mu,lambda);
        strain_energy_aa = strain_energy_aa*yita*w(1);              % consider the volume
        % fx = -du/dx
        fxy(aa) = -imag(strain_energy_aa/dh);
        % kx = d^2u/dx^2
        kxy(aa,aa) = 2*(real(strain_energy-strain_energy_aa)/dh)/dh;
    
        for bb = aa+1:6
            
            Fxy_bb = Fxy_h(bb,:);
            Fxy_bb = reshape(Fxy_bb,2,2);
            Cxy_bb = transpose(Fxy_bb)*Fxy_bb;
            I1xy_bb = trace(Cxy_bb);
            Jv_bb = det(Fxy_bb);
            strain_energy_bb = strain_energy_density(I1xy_bb,Jv_bb,mu,lambda);
            strain_energy_bb = strain_energy_bb*yita*w(1);              % consider the volume
            Fxy_ab = Fxy_aa + Fxy_bb - Fxy;
            Cxy_ab = transpose(Fxy_ab)*Fxy_ab;
            I1xy_ab = trace(Cxy_ab);
            Jv_ab = det(Fxy_ab);
            strain_energy_ab = strain_energy_density(I1xy_ab,Jv_ab,mu,lambda);
            strain_energy_ab = strain_energy_ab*yita*w(1);              % consider the volume
            kxy(aa,bb) = (real(strain_energy_aa + strain_energy_bb - strain_energy - strain_energy_ab)/dh)/dh;
            kxy(bb,aa) = kxy(aa,bb);
            
        end
    end
    
end
end


function strain_energy = strain_energy_density(I1xy,Jv,mu,lambda)
strain_energy = 1/2*mu*(I1xy-2);
strain_energy =  strain_energy + (-mu*log(Jv)+1/2*lambda*(log(Jv))^2);
end

function N = shape_function(z1,z2)

% quadratic shape function in reduced space

N(1) = xi(1);
N(2) = xi(2);
N(3) = 1.-xi(1)-xi(2);
N = [N1; N2; N3];

end

function [dNdx,yita] = shape_function_derivative(z1,z2,xa,ya)

% derivative of shape function in reduced space

dNdz = [ 1,   0;
    0,   1;
    -1,  -1;];
%
dxdz = [xa*dNdz;ya*dNdz];
yita = det(dxdz);                                                       % Jacobian of the mapping

% derivative of shape function in physical space
dNdx = dNdz/dxdz;

end


function [xi,w] = Gaussian_point_weight(Nint,shape,dim)
% triangle element

w = [0.5];                                                      % Gaussian weights
xi(1,1) = 1./3.;
xi(2,1) = 1./3.;

end

% Read abaqus inp code
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