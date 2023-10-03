%%% Apply NEB to capture MEP
clear;
clc

abaqusfile = 'Data_Files\t06_L5_theta40-50_twobeam_mesh4-50_indent5_n200_Quasi_mu1_lambda3_E100_nu03.inp';
[coord,connect] = inp2mat(abaqusfile);
nnode = size(coord,1); % total node number
nelem = size(connect,1); % total element number
mu_1 = 1; lambda_1 = 3; type1 = [1:400;861:1260];
mu_2 = 38.5; lambda_2 = 57.7; type2 = 401:860;
mu = zeros(nelem,1);     % material properties are assigned to each element
lambda = zeros(nelem,1);
mu(type1) = mu_1; lambda(type1) = lambda_1;
mu(type2) = mu_2; lambda(type2) = lambda_2;

%%% Input files
% S1
xy_1 = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_stressfree_topfree_Re1e-5.txt');
% S2
xy_2 = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_shape_topfree_Re1e-5.txt');
% first image in BITSS
xy_BITSS1 = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS1_alpha10_beta01_dist005-005_iter3_Re1e-5.txt');
% second image in BITSS
xy_BITSS2 = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS2_alpha10_beta01_dist005-005_iter3_Re1e-5.txt');
% iterative results from energy descent method
xy_EngDes1 = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_EnergyDescent_dist005-005_BITSS1_Re1e-5.txt');
nDes1 = size(xy_EngDes1,1)./2./nnode;
xy_EngDes1 = reshape(xy_EngDes1,2*nnode,nDes1);

total_node = [1:nnode]';                        % total node

l = abs(coord(:,2)-6.2) < 1e-5; % top edge
top_node = total_node (l);

l = abs(coord(:,2)-0) < 1e-3; % bottom edge
bottom_node = total_node (l);

% Displacement boundary conditions
% [node, direction, values]
disp_fix_1 = [bottom_node,ones(length(bottom_node),1),zeros(length(bottom_node),1)]; % fix the bottom edge in x direction
disp_fix_2 = [bottom_node,2*ones(length(bottom_node),1),zeros(length(bottom_node),1)]; % fix the bottom edge in y direction 
disp_fix = [disp_fix_1;disp_fix_2];

disp_ind = 2*(disp_fix(:,1)-1)+disp_fix(:,2); % find the fixed node index
force_applied = [];
force_gradient = 0;

efn = @(xy)total_energy(xy,coord,connect,mu,lambda);
% efn returns total energy
gfn = @(xy)total_force_boundary(xy,coord,connect,mu,lambda,disp_fix,force_applied,force_gradient);
% gfn returns energy gradient 
hfn = @(xy)total_stiffness_boundary(xy,coord,connect,mu,lambda,disp_fix);
% hfn returns hessian matrix

u1 = efn(xy_1); % S1 energy
u2 = efn(xy_2); % S2 energy

%%% start with low number of images 
N = 5;
xy_init = zeros(2*nnode,N);
xy_init(:,1:4) = xy_EngDes1(:,40:-10:10); % Use iterative results in energy descent method as initial guess
xy_init(:,5) = (xy_BITSS1+xy_BITSS2)./2; % Use BITSS results as initial guess
x0 = xy_init(:);
%%% 

%%% Linear interpolation of low number of images as initial guess
% N = 11;
% xy_init = zeros(2*nnode,N);
% xy_init1 = load('Data_Files\t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_NEB_N=5_EngDes_BITSS_kk=1e-3_Re1e-5.txt');
% xy_init(:,1) = (xy_1 + xy_init1(:,1))./2;
% xy_init(:,end) = (xy_2 + xy_init1(:,end))./2;
% xy_init(:,2:2:10) = xy_init1;
% xy_init(:,3:2:9) = (xy_init1(:,1:4) + xy_init1(:,2:5))./2;
% x0 = xy_init(:);
%%%

%%% run
tic;
x1 = NEB_Jacobian_NewTangt(x0,xy_1,u1,xy_2,u2,N,2*nnode,efn,gfn,hfn); % NEB without image climbing
timeElapsed1 = toc;

tic;
x2 = NEB_Jacobian_NewTangt_Saddle(x1,xy_1,u1,xy_2,u2,N,2*nnode,efn,gfn,hfn); % climbing image algorithm
timeElapsed2 = toc;
%%%

x2 = reshape(x2,2*nnode,N);

for i = 1:N
    coord_def(:,1) = x2(1:2:end,i);
    coord_def(:,2) = x2(2:2:end,i);
    TR = triangulation(connect,coord_def);
    triplot(TR,'Color','r');hold on; % plot NEB results
end

coord_def(:,1) = xy_1(1:2:end,1);
coord_def(:,2) = xy_1(2:2:end,1);
TR = triangulation(connect,coord_def);
triplot(TR,'k');hold on; % plot S1

coord_def(:,1) = xy_2(1:2:end,1);
coord_def(:,2) = xy_2(2:2:end,1);
TR = triangulation(connect,coord_def);
triplot(TR,'k');hold on;axis equal;axis off; % plot S2

coord_def(:,1) = xy_BITSS1(1:2:end,1);
coord_def(:,2) = xy_BITSS1(2:2:end,1);
TR = triangulation(connect,coord_def);
triplot(TR,'b');hold on;axis equal;axis off; % plot first image in BITSS

coord_def(:,1) = xy_BITSS2(1:2:end,1);
coord_def(:,2) = xy_BITSS2(2:2:end,1);
TR = triangulation(connect,coord_def);
triplot(TR,'b');hold on;axis equal;axis off; % plot second image in BITSS

filename = 't06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_NEB_N=5_EngDes_BITSS_kk=1e-3_Re1e-5.txt';
fid = fopen(filename,'wt');
for ii = 1:2*nnode
    fprintf(fid,'%.12f,%.12f,%.12f,%.12f,%.12f\n',...
        x2(ii,1),x2(ii,2),x2(ii,3),x2(ii,4),x2(ii,5));
end
fclose(fid);

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