% p_robust_MG0adapt_solver runs the variant of the a-posteriori-steered 
% multigrid (APS-MG) solver adaptive number of post-smoothing developed in 
% [Miraci, Papez, Vohralik. SIAM J. Sci. Comput. (2021)]
% https://doi.org/10.1137/19M1275929 
%
% Linear system: A*x_exact = F
%
% Solver features:
% *) initial guess is x_approx = x_0 = 0 and each iteration consists of:
%    *) MG V-cycle with zero pre- and an adaptive number 
%       of post-smoothing steps, P1-coarse solve 
%    *) smoothing is additive Schwarz/block-Jacobi with respect 
%       to patches of elements around vertices; 
%    *) optimal step-sizes (lambda) in the error correction stage
%
% Stopping criterion employed 
%  || F - A*x_approx || / ||F|| <= 10^(-5)
%
% INPUT:
% *) problem_number: a range of different regularity problems 
%    sine: 2; peak: 4; L-shaped domain: 5;
%    skyscraper (with jump O(1) in diffusion coeff): 12; 
%    skyscraper (with jump O(10^7) in diffusion coeff): 13;
%    checkerboard (with jump O(1) in diffusion coeff): 14; 
%    checkerboard (with jump O(10^6) in diffusion coeff): 15;
% *) m: non-decreasing sequence of polynomial degrees 
%    distributed per level 
%    (the coarsest level is always P1)
%    implicitly sets the number of levels J in the hierarchy
% *) mesh_uniformity: bulk-chasing parameter determining the 
%    type of mesh hierarchy
%    uniformly-refined: 1; graded: smaller than 1;
% *) Hmax: initial mesh-size
% *) theta: adaptivity-parameter for the number 
%    of post-smoothing steps,
% *) maxiter: to limit the maximum number of iterations
%
% OUTPUT:
% *) results: structure with fields:
%    *) x_approx: the last approximation computed by the solver
%    *) n_iter: number of iterations to reach the stopping criterion
%    *) rel_res: relative residual per iteration
%    *) adapt_smoothingsteps_level: sum over iterations of smoothing
%       steps per level
%
% Remark: below we limit the maximum number of 
%    post-smoothing steps by setting maxlocsmooths_adapt = 5.
%    Setting maxlocsmooths_adapt = 1, the solver then
%    runs the non-adaptive variant p_robust_MG01_solver.m
%
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG


function [results] = p_robust_MG0adapt_solver(problem_number, m, mesh_uniformity, Hmax, theta, maxiter)

if nargin < 6
    maxiter = 100;
end
if nargin < 5
    theta = 0.2;
end
if nargin < 4
    Hmax = 0.1;
end
if nargin < 3
    mesh_uniformity = 1;    % uniform mesh refinement
end

maxlocsmooths_adapt = 5; %to avoid possible infinite loops


%% initialize

initialize;

%% solver

fprintf('\n************************************************************\n');
fprintf('***\n');
fprintf('***      A-Posteriori-Steered MG(0,adapt) Solver       \n');
fprintf('***  (adaptive number of post-smoothing steps, theta = %.1f',theta);fprintf(')\n');
fprintf('***\n');
fprintf('***                 %s \n', solution.name);
fprintf('***            number of levels %d ',J-1);fprintf('\n');
fprintf('***   polynomial order per level [');fprintf('%g ', m);fprintf(']\n');
fprintf('***\n');
fprintf('************************************************************\n');


results.adapt_smoothingsteps_level = zeros(1, J);
x_approx = zeros(length(DOF2(J).freenodes),1);
res = F{J};
rel_res_init = 1; 

% iterations
for iter = 1:maxiter

    % calling MG solver with adaptive number of smoothing per level
    [rho_alg, ~, smoothiters] = ...
        smoother_bJ_adaptnumber(x_approx,res,theta,maxlocsmooths_adapt);

    %new iterate
    x_approx = x_approx + rho_alg;
    res = F{J} - A{J}*x_approx;
    relres(iter) = norm(res)/norm(F{J});

    fprintf('Iter %d:     relative residual norm = %e \n', iter, relres(iter));
    fprintf('   adaptive number of smoothings per level: \n');
    disp(smoothiters);
    results.adapt_smoothingsteps_level = results.adapt_smoothingsteps_level + smoothiters;

    %stopping criterion
    if (relres(iter) <= 10^(-5)* rel_res_init)
        fprintf('Accuracy attained in iter %d \n', iter);
        break;
    end

end% of iterations

%% summary

results.x_approx = x_approx;
results.n_iter = iter;
results.rel_res = relres;
%results.adapt_smoothingsteps_level;

end
