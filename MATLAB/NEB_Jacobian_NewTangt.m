function NEB_Jacobian_NewTangt=NEB_Jacobian_NewTangt(x0,x1,u1,x2,u2,N,n,U,DU,DDU) 
% x=[x1;x2;...;xN] N images
% tangent vectors updated estimation based on XI-1, XI and XI+1
% no climbing image algorithm
% N number of images
% n degrees of freedom
% U energy function, input one image, return total energy
% DU energy gradient, input one image, return energy gradient
% DDU hessian matrix, input one image, return hessian matrix
% x1 x2 two stable states
% u1 u2 energy of x1 x2,
% x0 is the initial guess

global kk; % spring constant is a fixed value
kk=1e-5;

%%% Set up shared variables with outfun
%%% Output results during iterations
history.x = [];
history.fval = [];
function stop = outfun(x,optimValues,state)
stop = false;
abaqusfile = 'Data_Files\BuckleBeam_y1_Left_mesh4-50.inp';
[coord,connect] = inp2mat(abaqusfile); % input nodal and element information
switch state
    case 'init'
        hold on
    case 'iter'
        history.fval = [history.fval; optimValues.fval];
        history.x = [history.x; x];
        tempx = x;
        [temp_force,~] = FB(tempx,x1,u1,x2,u2,N,n,U,DU,DDU); % calculate total force on images
        tempx = reshape(tempx,n,N);
        [temp_tau,~] = tangent(tempx,x1,x2,N,n); % calculate tangent vectors of images
        temp_tau = temp_tau(:); % rearranged in a single column vector
        temp_force_perpendicular = temp_force-(temp_tau'*temp_force).*temp_tau; % perpendicular force
        temp_force_parallel = (temp_tau'*temp_force).*temp_tau; % parallel force
        temp_force_perpendicular_magn = norm(temp_force_perpendicular); % magnitude of perpendicular force
        temp_force_parallel_magn = norm(temp_force_parallel); % magnitude of parallel force
        fprintf('perpendicular %f parallel %f spring constant %f\n',temp_force_perpendicular_magn,temp_force_parallel_magn,kk);
        temp_force_parallel = reshape(temp_force_parallel,n,N);
        temp_force_perpendicular = reshape(temp_force_perpendicular,n,N);
        figure(1);
        for i = 1:N

            %%% record iterative results 
            filename = 'BuckleBeam_y1_mesh4-50_NEB_N=3_EngDes_BITSS_kk=1e-5_Interm.txt';
            fid = fopen(filename,'at');
            for ii = 1:n
                fprintf(fid,'%.12f,',tempx(ii,i));
            end
            fprintf(fid,'\n');
            fclose(fid);
            %%%

            coord_def(:,1) = tempx(1:2:end-1,i);
            coord_def(:,2) = tempx(2:2:end,i);
            TR = triangulation(connect,coord_def);
            triplot(TR);hold on; % plot N images
            quiver(tempx(1:2:end-1,i),tempx(2:2:end,i),temp_force_parallel(1:2:end-1,i),temp_force_parallel(2:2:end,i),0,'r');hold on;
            quiver(tempx(1:2:end-1,i),tempx(2:2:end,i),temp_force_perpendicular(1:2:end-1,i),temp_force_perpendicular(2:2:end,i),0,'g');hold on;
            % plot the distribution of perpendicular and parallel force
        end
        coord_def(:,1) = x1(1:2:end-1);
        coord_def(:,2) = x1(2:2:end);
        TR = triangulation(connect,coord_def);
        triplot(TR,'k');hold on; % plot the first stable state
        coord_def(:,1) = x2(1:2:end-1);
        coord_def(:,2) = x2(2:2:end);
        TR = triangulation(connect,coord_def);
        triplot(TR,'k');hold on; % plot the second stable state
        hold off;
    case 'done'
        hold off
    otherwise
end
end
%%%

%%% run
options = optimoptions(@fsolve,'OutputFcn',@outfun,'Display','iter','SpecifyObjectiveGradient',true,...
    'MaxIterations',5000,'MaxFunctionEvaluations',50000,'StepTolerance',1e-6,'FunctionTolerance',1e-6);
[xsol,fval,exitflag,output,JAC] = fsolve(@(x)FB(x,x1,u1,x2,u2,N,n,U,DU,DDU),x0,options);
NEB_Jacobian_NewTangt = xsol;
%%%

end

function [FB,FB_Jac]=FB(x,x1,u1,x2,u2,N,n,U,DU,DDU) % force balance with Jacobian matrix
% FB: force from potential energy and elastic band on images
% FB_Jac: Jacobian matrix of FB regarding x for computing acceleration

global kk;

x = reshape(x,[n,N]); % arrange the input into N column vectors
[tau,t] = tangent(x,x1,x2,N,n); % calculate tangent vectors

uu = zeros(1,N); % energy of N images
du = zeros(n,N); % energy gradient of N images
ddu = zeros(n,n,N); % laplace operator energy of N images

for i = 1:N
    uu(i) = U(x(:,i));
    du(:,i) = DU(x(:,i));
    ddu(:,:,i) = DDU(x(:,i));
end

FN = zeros(n,N); % perpendicular force
FP = zeros(n,N); % parallel force

for i = 1:N  
    FN(:,i) = -du(:,i)+du(:,i)'*tau(:,i).*tau(:,i); % perpendicular force
end

FP(:,1) = -kk.*(norm(x(:,2)-x(:,1))-norm(x(:,1)-x1)).*tau(:,1); % paralle force for the first image
FP(:,N) = -kk.*(norm(x2-x(:,N))-norm(x(:,N)-x(:,N-1))).*tau(:,N); % parallel force for the last image

for i = 2:N-1
    FP(:,i) = -kk.*(norm(x(:,i+1)-x(:,i))-norm(x(:,i)-x(:,i-1))).*tau(:,i); % parallel force
end

FB = FN+FP; % total force exerted on N images
FB = FB(:); % reshape into a single column vector

FB_Jac = zeros(n*N,n*N); % Jacobian matrix,n*N dimension
%%% calculate Jabocian matrix FB_Jac
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

function [tau,t]=tangent(x,x1,x2,N,n) % calculate tangent vectors of N images 
tau = zeros(n,N); % tau vector with unit magnitude
t = zeros(n,N); % t vector,magnitude is not one

for i = 2:N-1   % t=(x(i+1)-x(i))/||x(i+1)-x(i)||+(x(i)-x(i-1))/||x(i)-x(i-1)||
    t(:,i) = (x(:,i+1)-x(:,i))./norm(x(:,i+1)-x(:,i))+(x(:,i)-x(:,i-1))./norm(x(:,i)-x(:,i-1));
    tau(:,i) = t(:,i)./norm(t(:,i));
end

t(:,1) = (x(:,2)-x(:,1))./norm(x(:,2)-x(:,1))+(x(:,1)-x1)./norm(x(:,1)-x1);
tau(:,1) = t(:,1)./norm(t(:,1)); % tangent vector for the first image

t(:,N) = (x2-x(:,N))./norm(x2-x(:,N))+(x(:,N)-x(:,N-1))./norm(x(:,N)-x(:,N-1));
tau(:,N) = t(:,N)./norm(t(:,N)); % tangent vector for the last image
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

