function [x1,x2] = bitss_hessian(x1, x2, energy_fn, gradient_fn, hessian_fn, varargin)
% Perform the BITSS method to find a transition state.
% Returns the final two images on two sides of the transition state.
%
% Arguments:
%   x1, x2: Two starting states for BITSS. Typically minima.
%   energy_fn: The energy function to use. Takes one state and returns energy.
%   gradient_fn: The gradient of the energy. Returns the gradient in a column vector.
%   hessian_fn: Hessian matrix of the energy,used in calculating Total
%   hessian matrix for computing acceleration

% Optional Arguments (Using name-value pairs):
%   dist_step: The fraction change in the distance each BITSS step. (default = 0.5)
%   max_iter: The maximum number of iterations in the BITSS method. (default = 100)
%   convergence_type: Choose the type of the convergence test, 'distance' or 'gradient. (default = 'distance')
%   convergence_value: Criteria for the convergence test. (distance: default = initial_dist/100, 'gradient':1e-5)
%   min_max_iter: Maximum number of iterations for minimisation. (default = 10000)
%   min_convergence_value: Convergence criteria for gradinet size. (default = 1e-5)

  %%% Initialise
  x=[x1;x2];
  x10=x1;x20=x2; % two stable states set as initials
  ndof = length(x1); % degrees of freedom  
  dist0 = norm(x1 - x2); % initial distance between two states
  U1 = energy_fn(x1); % energy of initial states
  U2 = energy_fn(x2);

  global num_iter; % number of calling minimization function
  num_iter = 0; 

  global ke kd; % spring stiffness in BITSS ke and kd
  compute_coef(x1, x2, energy_fn, gradient_fn(x1), gradient_fn(x2),dist0); % Compute ke and kd in the initial step
  kke = full(ke);
  kkd = full(kd);

  %%% Record the set distance and actual distance between states dist, spring stiffness ke and
  %%% kd, energy of two states U1 and U2 in the initial step
  filename = 't06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS_alpha10_beta01_dist005-005_iter3_springstiff_Re1e-5.txt';
  fid = fopen(filename,'at');
  fprintf(fid,'%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n',dist0,dist0,kke,kkd,U1,U2);
  fclose(fid);
  %%%

  %%% Record the first image during iterations in the initial step
  filename = 't06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS_pt1_alpha10_beta01_dist005-005_iter3_Interm_Re1e-5.txt';
  fid = fopen(filename,'at');
  for ii = 1:ndof
      fprintf(fid,'%.12f\n',x1(ii));
  end
  fclose(fid);
  %%%

  %%% Record the second image during iterations in the initial step
  filename = 't06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS_pt2_alpha10_beta01_dist005-005_iter3_Interm_Re1e-5.txt';
  fid = fopen(filename,'at');
  for ii = 1:ndof
      fprintf(fid,'%.12f\n',x2(ii));
  end
  fclose(fid);
  %%%

  %%% Set parameters
  dist_step = 0.05; 
  max_iter = 200;
  convergence_type = 'distance';
  convergence_value = 0.05 * dist0;
  for iarg = 1:2:length(varargin)
    switch varargin{iarg}
      case 'dist_step'
        dist_step = varargin{iarg+1};
      case 'max_iter'
        max_iter = varargin{iarg+1};
      case 'convergence_type'
        convergence_type = varargin{iarg+1};
        if (~ismember(convergence_type, {'distance','gradient'}))
          error("convergence_type should be 'distance' or 'gradient'");
        elseif ((convergence_value == 0) & strcmp(convergence_type, 'gradient'))
          convergence_value = 1e-5;
        end
      case 'convergence_value'
        convergence_value = varargin{iarg+1};
    end
  end
  %%%

  %%% Run
  abaqusfile = 'Data_Files\t06_L5_theta40-50_twobeam_mesh4-50_indent5_n200_Quasi_mu1_lambda3_E100_nu03.inp'; % Input node and element information
  [coord,connect] = inp2mat(abaqusfile);
  for iter = 1:max_iter
      fprintf('iteration number %d\n',iter);
      dist0 = (1 - dist_step) * dist0; % Calculate the constrained distance dist0
      options = optimoptions(@fsolve,'Display','iter','SpecifyObjectiveGradient',true,...
          'MaxIterations',200,'MaxFunctionEvaluations',200,'StepTolerance',1e-8,'FunctionTolerance',1e-8); % Setup parameters for fsolve
      [xsol,fval_sol,exitflag_sol,output_sol,JAC_sol] = fsolve(@(x)total_gradient_Hessian(x,energy_fn,...
          gradient_fn,hessian_fn,dist0,ndof),x,options); % Minimize BITSS energy through fsolve function
      x = xsol;
      x1 = x(1:ndof);
      x2 = x(ndof+1:end);

      %%% Plot the configuration of x1 and x2
      figure(1);
      coord_def(:,1) = x1(1:2:end-1);
      coord_def(:,2) = x1(2:2:end);
      TR = triangulation(connect,coord_def);
      triplot(TR,'r');hold on; % plot first image

      coord_def(:,1) = x2(1:2:end-1);
      coord_def(:,2) = x2(2:2:end);
      TR = triangulation(connect,coord_def);
      triplot(TR,'r');hold on; % plot second image
      %%%

      %%% Plot the configuration of stable states x10 and x20
      coord_def(:,1) = x10(1:2:end-1);
      coord_def(:,2) = x10(2:2:end);
      TR = triangulation(connect,coord_def);
      triplot(TR,'k');hold on;

      coord_def(:,1) = x20(1:2:end-1);
      coord_def(:,2) = x20(2:2:end);
      TR = triangulation(connect,coord_def);
      triplot(TR,'k');hold on;
      hold off;
      %%%

      %%% Record the intermediate images x1 and x2, actual distance between x1 and
      %%% x2 adist, constrained distance dist0, spring stiffness ke and kd,
      %%% and energy of two states U1 and U2 during iterations
      filename = 't06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS_pt1_alpha10_beta01_dist005-005_iter3_Interm_Re1e-5.txt';
      fid = fopen(filename,'at');
      for ii = 1:ndof
          fprintf(fid,'%.12f\n',x1(ii)); % Output state x1
      end
      fclose(fid);

      filename = 't06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS_pt2_alpha10_beta01_dist005-005_iter3_Interm_Re1e-5.txt';
      fid = fopen(filename,'at');
      for ii = 1:ndof
          fprintf(fid,'%.12f\n',x2(ii)); % Output state x2
      end
      fclose(fid);

      compute_coef(x1, x2, energy_fn, gradient_fn(x1), gradient_fn(x2),dist0); % Compute the spring stiffness ke and kd
      U1 = energy_fn(x1);
      U2 = energy_fn(x2);
      adist = norm(x1-x2); % Distance between two states x1 and x2
      kke = full(ke);
      kkd = full(kd);

      filename='t06_L5_theta40-45_twobeam_mesh4-50_mu1_lambda3_E100_nu03_BITSS_alpha10_beta01_dist005-005_iter3_springstiff_Re1e-5.txt';
      fid = fopen(filename,'at');
      fprintf(fid,'%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n',dist0,adist,kke,kkd,U1,U2);
      fclose(fid);

      fprintf('U1 %f U2 %f Distance0 %f Actual Distance %f\n',U1,U2,dist0,adist);
      fprintf('min_iter %f\n',num_iter);
      %%%
    
      if (bitss_check_convergence(x, dist0, gradient_fn, convergence_type, convergence_value))
          break;
      end
  end
end

function compute_coef(x1, x2, energy_fn, g1, g2,dist0) 
%%% Calculate the spring stiffness ke and kd based on two states x1 and x2
%%% and their gradient g1 and g2. dist0 is constrained distance during each 
%%% iteration. dist0 can be replaced with actual distance between two states
%%% adist

  global ke kd;

  %%%% set up parameters alpha and beta
  alpha = 10;  
  beta = 0.1;
  %%%%

  eb = estimate_barrier(x1, x2, energy_fn); %%% Estimate the energy barrier between x1 and x2
  adist=norm(x1-x2); % Calculate the actual distance between two states
  
  ke = alpha / (2 * eb);
  kd1 = sqrt(sum(g1.^2) + sum(g2.^2)) / (2.8284 * beta * dist0); 
  %%% dist0 can be replaced with adist
  kd2 = eb / (beta * dist0^2);
  %%% dist0 can be replaced with adist
  kd = max(kd1, kd2);
  
  function eb = estimate_barrier(x1, x2, energy_fn)
    emin = max(energy_fn(x1), energy_fn(x2));
    emax = -Inf;
    for i = 1:9
      t = i / 10;
      x = (1 - t) * x1 + t * x2;
      emax = max(emax, energy_fn(x));
    end
    eb = emax - emin;
  end
end

function converged = bitss_check_convergence(x, dist0, gradient_fn, convergence_type, convergence_value)
  if strcmp(convergence_type, 'distance')
    converged = (dist0 < convergence_value);
  elseif strcmp(convergence_type, 'gradient')
    ndof = length(x) / 2;
    xmid = (x(1:ndof) + x(ndof+1:end)) / 2;
    gnorm = norm(gradient_fn(xmid));
    converged = (gnorm < convergence_value);
  end
end

function [dEdX,ddEdXX] = total_gradient_Hessian(x,energy_fn,gradient_fn,hessian_fn,dist0,ndof)
%%% Return the gradient dEdX and heissian matrix ddEdXX of BITSS energy. 
%%% Inputs are two states x=[x1;x2], energy function energy_fn, energy gradient
%%% function gradient_fn, hessian function hessian_fn, constrained distance
%%% between two states dist0, and degrees of freedom ndof
global ke kd;
global num_iter;

num_iter = num_iter+1; % Record the number of calling this function

x1 = x(1:ndof); % first state
x2 = x(ndof+1:end); % second state

E1 = energy_fn(x1);E2 = energy_fn(x2); % energy
G1 = gradient_fn(x1);G2 = gradient_fn(x2); % energy gradient
H1 = hessian_fn(x1);H2 = hessian_fn(x2); % energy hessian matrix
adist = norm(x1-x2); % actual distance between x1 and x2

%%% Update the stiffness ke and kd every three calling of minimization
%%% function
if (mod(min_iter-1,3) == 0); compute_coef(x1, x2, energy_fn, G1, G2,dist0); end
%%%

dEdX1 = (1+2.*ke.*(E1-E2)).*G1+2.*kd.*(1-dist0./adist).*(x1-x2);
dEdX2 = (1+2.*ke.*(E2-E1)).*G2+2.*kd.*(1-dist0./adist).*(x2-x1);
dEdX = [dEdX1;dEdX2]; % BITSS energy gradient

ddEdX1dX1 = (1+2.*ke.*(E1-E2)).*H1+2.*ke.*G1*G1'+2.*kd.*(1-dist0./adist).*eye(ndof)+2.*kd.*dist0./adist.^3.*(x1-x2)*(x1-x2)';
ddEdX2dX2 = (1+2.*ke.*(E2-E1)).*H2+2.*ke.*G2*G2'+2.*kd.*(1-dist0./adist).*eye(ndof)+2.*kd.*dist0./adist.^3.*(x2-x1)*(x2-x1)';
ddEdX1dX2 = -2.*ke.*G1*G2'-2.*kd.*(1-dist0./adist).*eye(ndof)-2.*kd.*dist0./adist.^3.*(x1-x2)*(x1-x2)';
ddEdX2dX1 = -2.*ke.*G2*G1'-2.*kd.*(1-dist0./adist).*eye(ndof)-2.*kd.*dist0./adist.^3.*(x1-x2)*(x1-x2)';

ddEdXX = [ddEdX1dX1,ddEdX1dX2;
    ddEdX2dX1,ddEdX2dX2]; % BITSS hessian matrix
end

function [coord,connect] = inp2mat(file)
% Load inp file to output the nodal 
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