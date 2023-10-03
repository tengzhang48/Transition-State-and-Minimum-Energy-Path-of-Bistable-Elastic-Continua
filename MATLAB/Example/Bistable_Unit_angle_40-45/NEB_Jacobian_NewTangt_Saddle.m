function NEB_Jacobian_NewTangt = NEB_Jacobian_NewTangt_Saddle(x0,x1,u1,x2,u2,N,n,U,DU,DDU) 
% x=[x1;x2;...;xN] N images
% tangent vectors updated estimation based on XI-1, XI and XI+1.
% climbing image algorithm is applied
% N number of images
% n degrees of freedom
% U energy function, input one image, return total energy
% DU energy gradient, input one image, return energy gradient
% DDU hessian matrix, input one image, return hessian matrix 
% x1 x2 two stable states
% u1 u2 energy of x1 x2
% x0 is the initial guess

global kk;
kk = 4e-3; % set up spring stiffness

%%% Set up shared variables with outfun
%%% Output results during iterations
history.x = [];
history.fval = [];
function stop = outfun(x,optimValues,state)
stop = false;
abaqusfile = 'Data_Files\t06_L5_theta40-50_twobeam_mesh4-50_indent5_n200_Quasi_mu1_lambda3_E100_nu03.inp';
[coord,connect] = inp2mat(abaqusfile); % input nodal and element information
switch state
    case 'init'
        hold on
    case 'iter'
        % Concatenate current point and objective function
        % value with history. x must be a row vector.
        history.fval = [history.fval; optimValues.fval];
        history.x = [history.x; x];

        tempx = x;
        [temp_force,~] = FB(tempx,x1,u1,x2,u2,N,n,U,DU,DDU); % calculate force on images
        tempx = reshape(tempx,n,N);
        [temp_tau,~] = tangent(tempx,x1,x2,N,n); % calculate the tangent vector
        temp_tau = temp_tau(:);
        temp_force_perpendicular = temp_force-(temp_tau'*temp_force).*temp_tau; % perpendicular force on images
        temp_force_parallel = (temp_tau'*temp_force).*temp_tau; % parallel force on images
        temp_force_perpendicular_magn = norm(temp_force_perpendicular); % magnitude of perpendicular force
        temp_force_parallel_magn = norm(temp_force_parallel); % magnitude of parallel force
        fprintf('perpendicular %f parallel %f spring constant %f\n',temp_force_perpendicular_magn,temp_force_parallel_magn,kk);
        % output force magnitude and spring constant
        temp_force_parallel = reshape(temp_force_parallel,n,N);
        temp_force_perpendicular = reshape(temp_force_perpendicular,n,N);
        figure(1);
        for i = 1:N
            coord_def(:,1) = tempx(1:2:end-1,i);
            coord_def(:,2) = tempx(2:2:end,i);
            TR = triangulation(connect,coord_def);
            triplot(TR);hold on; % plot each image during iteration
            quiver(tempx(1:2:end-1,i),tempx(2:2:end,i),temp_force_parallel(1:2:end-1,i),temp_force_parallel(2:2:end,i),0,'r');hold on;
            quiver(tempx(1:2:end-1,i),tempx(2:2:end,i),temp_force_perpendicular(1:2:end-1,i),temp_force_perpendicular(2:2:end,i),0,'g');hold on;
            % plot the distribution of parallel and perpendicular force
        end
        coord_def(:,1) = x1(1:2:end-1);
        coord_def(:,2) = x1(2:2:end);
        TR = triangulation(connect,coord_def);
        triplot(TR,'k');hold on; % plot first stable shape
        coord_def(:,1) = x2(1:2:end-1);
        coord_def(:,2) = x2(2:2:end);
        TR = triangulation(connect,coord_def);
        triplot(TR,'k');hold on; % plot second stable shape
        hold off;
    case 'done'
        hold off
    otherwise
end
end 
%%%

%%% run
options = optimoptions(@fsolve,'OutputFcn',@outfun,'Display','iter','SpecifyObjectiveGradient',true,...
    'MaxIterations',100,'MaxFunctionEvaluations',100,'StepTolerance',1e-6,'FunctionTolerance',1e-6);
[xsol,fval,exitflag,output,JAC] = fsolve(@(x)FB(x,x1,u1,x2,u2,N,n,U,DU,DDU),x0,options);
NEB_Jacobian_NewTangt = xsol;
%%%

end

function [FB,FB_Jac] = FB(x,x1,u1,x2,u2,N,n,U,DU,DDU) % force balance with Jacobian matrix
% FB: force from potential energy and elastic band on images, climbing
% image algorithm is used
% FB_Jac: Jacobian matrix of FB regarding x for computing acceleration

global kk;

x = reshape(x,[n,N]); % re-arrange the input into N column vectors
[tau,t] = tangent(x,x1,x2,N,n); % calculate the tangent vector from images

uu = zeros(1,N); % energy of N images
du = zeros(n,N); % energy gradient of N images
ddu = zeros(n,n,N); % laplace operator energy of N images

for i = 1:N
    uu(i) = U(x(:,i));
    du(:,i) = DU(x(:,i));
    ddu(:,:,i) = DDU(x(:,i));
end

[~,index_saddle] = max(uu); % find the image with the highest energy

FN = zeros(n,N); % perpendicular force for N images
FP = zeros(n,N); % parallel force for N images

for i = 1:N  
    FN(:,i) = -du(:,i)+du(:,i)'*tau(:,i).*tau(:,i); % perpendicular force
end

if (1 ~= index_saddle) % for the first image
    FP(:,1) = -kk.*(norm(x(:,2)-x(:,1))-norm(x(:,1)-x1)).*tau(:,1);
else
    FP(:,1) = du(:,1)'*tau(:,1).*tau(:,1);
end

if (N ~= index_saddle) % for the last image
    FP(:,N) = -kk.*(norm(x2-x(:,N))-norm(x(:,N)-x(:,N-1))).*tau(:,N);
else
    FP(:,N) = du(:,N)'*tau(:,N).*tau(:,N);
end

for i=2:N-1
    if (i ~= index_saddle)
        FP(:,i) = -kk.*(norm(x(:,i+1)-x(:,i))-norm(x(:,i)-x(:,i-1))).*tau(:,i);
    else
        FP(:,i) = du(:,i)'*tau(:,i).*tau(:,i); % climbing image 
    end
end

FB = FN+FP; % total force exerted on N images
FB = FB(:); % reshape into a single column vector

FB_Jac = zeros(n*N,n*N); % Jacobian matrix,n*N dimension

%%% calculate Jacobian matrix 
for i = 2:N-1
    dtIdXI = eye(n).*(1./norm(x(:,i)-x(:,i-1))-1./norm(x(:,i+1)-x(:,i)))+(x(:,i+1)-x(:,i)).*(x(:,i+1)-x(:,i))'./...
        norm(x(:,i+1)-x(:,i)).^3-(x(:,i)-x(:,i-1)).*(x(:,i)-x(:,i-1))'./norm(x(:,i)-x(:,i-1)).^3;
    dtauIdXI = dtIdXI./norm(t(:,i))-1./norm(t(:,i)).^3.*t(:,i).*(dtIdXI*t(:,i))';
    dFIdXI = -ddu(:,:,i)+tau(:,i).*(ddu(:,:,i)*tau(:,i)+dtauIdXI*du(:,i))'+(du(:,i)'*tau(:,i)).*dtauIdXI...
    +kk./norm(x(:,i)-x(:,i-1)).*tau(:,i).*(x(:,i)-x(:,i-1))'-kk./norm(x(:,i+1)-x(:,i)).*tau(:,i).*(x(:,i)-x(:,i+1))'...
    +kk.*(norm(x(:,i)-x(:,i-1))-norm(x(:,i+1)-x(:,i))).*dtauIdXI;
    
    dtIdXIn1 = -eye(n)./norm(x(:,i)-x(:,i-1))+1./norm(x(:,i)-x(:,i-1)).^3.*(x(:,i)-x(:,i-1)).*(x(:,i)-x(:,i-1))';
    dtauIdXIn1 = 1./norm(t(:,i)).*dtIdXIn1-1./norm(t(:,i)).^3.*t(:,i).*(dtIdXIn1*t(:,i))';
    dFIdXIn1 = tau(:,i).*(dtauIdXIn1*du(:,i))'+du(:,i)'*tau(:,i).*dtauIdXIn1+kk./norm(x(:,i)-x(:,i-1)).*tau(:,i).*...
        (x(:,i-1)-x(:,i))'+kk.*(norm(x(:,i)-x(:,i-1))-norm(x(:,i+1)-x(:,i))).*dtauIdXIn1;
    
    dtIdXIp1 = eye(n)./norm(x(:,i+1)-x(:,i))-1./norm(x(:,i+1)-x(:,i)).^3.*(x(:,i+1)-x(:,i)).*(x(:,i+1)-x(:,i))';
    dtauIdXIp1 = 1./norm(t(:,i)).*dtIdXIp1-1./norm(t(:,i)).^3.*t(:,i).*(dtIdXIp1*t(:,i))';
    dFIdXIp1 = tau(:,i).*(dtauIdXIp1*du(:,i))'+du(:,i)'*tau(:,i).*dtauIdXIp1-kk./norm(x(:,i+1)-x(:,i)).*tau(:,i).*...
        (x(:,i+1)-x(:,i))'+kk.*(norm(x(:,i)-x(:,i-1))-norm(x(:,i+1)-x(:,i))).*dtauIdXIp1;
    
    FB_Jac(i*n-n+1:i*n,i*n-n+1:i*n) = dFIdXI;
    FB_Jac(i*n-n+1:i*n,i*n-2*n+1:i*n-n) = dFIdXIn1;
    FB_Jac(i*n-n+1:i*n,i*n+1:i*n+n) = dFIdXIp1;
end

dt1dX1 = eye(n).*(1./norm(x(:,1)-x1)-1./norm(x(:,2)-x(:,1)))+(x(:,2)-x(:,1)).*(x(:,2)-x(:,1))'./...
    norm(x(:,2)-x(:,1)).^3-(x(:,1)-x1).*(x(:,1)-x1)'./norm(x(:,1)-x1).^3;
dtau1dX1 = dt1dX1./norm(t(:,1))-1./norm(t(:,1)).^3.*t(:,1).*(dt1dX1*t(:,1))';
dF1dX1 = -ddu(:,:,1)+tau(:,1).*(ddu(:,:,1)*tau(:,1)+dtau1dX1*du(:,1))'+(du(:,1)'*tau(:,1)).*dtau1dX1...
    +kk./norm(x(:,1)-x1).*tau(:,1).*(x(:,1)-x1)'-kk./norm(x(:,2)-x(:,1)).*tau(:,1).*(x(:,1)-x(:,2))'...
    +kk.*(norm(x(:,1)-x1)-norm(x(:,2)-x(:,1))).*dtau1dX1;

dt1dX2 = eye(n)./norm(x(:,2)-x(:,1))-1./norm(x(:,2)-x(:,1)).^3.*(x(:,2)-x(:,1)).*(x(:,2)-x(:,1))';
dtau1dX2 = 1./norm(t(:,1)).*dt1dX2-1./norm(t(:,1)).^3.*t(:,1).*(dt1dX2*t(:,1))';
dF1dX2 = tau(:,1).*(dtau1dX2*du(:,1))'+du(:,1)'*tau(:,1).*dtau1dX2-kk./norm(x(:,2)-x(:,1)).*tau(:,1).*...
    (x(:,2)-x(:,1))'+kk.*(norm(x(:,1)-x1)-norm(x(:,2)-x(:,1))).*dtau1dX2;

FB_Jac(1:n,1:n) = dF1dX1;
FB_Jac(1:n,n+1:2*n) = dF1dX2;

dtNdXN = eye(n).*(1./norm(x(:,N)-x(:,N-1))-1./norm(x2-x(:,N)))+(x2-x(:,N)).*(x2-x(:,N))'./...
    norm(x2-x(:,N)).^3-(x(:,N)-x(:,N-1)).*(x(:,N)-x(:,N-1))'./norm(x(:,N)-x(:,N-1)).^3;
dtauNdXN = dtNdXN./norm(t(:,N))-1./norm(t(:,N)).^3.*t(:,N).*(dtNdXN*t(:,N))';
dFNdXN = -ddu(:,:,N)+tau(:,N).*(ddu(:,:,N)*tau(:,N)+dtauNdXN*du(:,N))'+(du(:,N)'*tau(:,N)).*dtauNdXN...
    +kk./norm(x(:,N)-x(:,N-1)).*tau(:,N).*(x(:,N)-x(:,N-1))'-kk./norm(x2-x(:,N)).*tau(:,N).*(x(:,N)-x2)'...
    +kk.*(norm(x(:,N)-x(:,N-1))-norm(x2-x(:,N))).*dtauNdXN;

dtNdXNn1 = -eye(n)./norm(x(:,N)-x(:,N-1))+1./norm(x(:,N)-x(:,N-1)).^3.*(x(:,N)-x(:,N-1)).*(x(:,N)-x(:,N-1))';
dtauNdXNn1 = 1./norm(t(:,N)).*dtNdXNn1-1./norm(t(:,N)).^3.*t(:,N).*(dtNdXNn1*t(:,N))';
dFNdXNn1 = tau(:,N).*(dtauNdXNn1*du(:,N))'+du(:,N)'*tau(:,N).*dtauNdXNn1+kk./norm(x(:,N)-x(:,N-1)).*tau(:,N).*...
    (x(:,N-1)-x(:,N))'+kk.*(norm(x(:,N)-x(:,N-1))-norm(x2-x(:,N))).*dtauNdXNn1;

FB_Jac(N*n-n+1:N*n,N*n-n+1:N*n) = dFNdXN;
FB_Jac(N*n-n+1:N*n,N*n-2*n+1:N*n-n) = dFNdXNn1;
%%%
end

function [tau,t] = tangent(x,x1,x2,N,n) % calculate the tangent vector for N images
tau = zeros(n,N); % tau vector with unit magnitude
t = zeros(n,N); % t vector,magnitude is not one

for i = 2:N-1   % t = (x(i+1)-x(i))/||x(i+1)-x(i)||+(x(i)-x(i-1))/||x(i)-x(i-1)||
    t(:,i) = (x(:,i+1)-x(:,i))./norm(x(:,i+1)-x(:,i))+(x(:,i)-x(:,i-1))./norm(x(:,i)-x(:,i-1));
    tau(:,i) = t(:,i)./norm(t(:,i));
end

t(:,1) = (x(:,2)-x(:,1))./norm(x(:,2)-x(:,1))+(x(:,1)-x1)./norm(x(:,1)-x1);
tau(:,1) = t(:,1)./norm(t(:,1)); % the first image

t(:,N) = (x2-x(:,N))./norm(x2-x(:,N))+(x(:,N)-x(:,N-1))./norm(x(:,N)-x(:,N-1));
tau(:,N) = t(:,N)./norm(t(:,N)); % the last image
end

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

