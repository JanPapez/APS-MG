% p_robust_MG01_solver runs the non-adaptive variant of the
% a-posteriori-steered multigrid (APS-MG) solver developed in 
% [Miraci, Papez, Vohralik. SIAM J. Sci. Comput. (2021)]
% https://doi.org/10.1137/19M1275929 
%
% Linear system: A*x_exact = F
%
% Solver features:
% *) initial guess is x_approx = x_0 = 0 and each iteration consists of:
%    *) MG V-cycle with zero pre- and one post-smoothing step, 
%       P1-coarse solve 
%    *) smoothing is additive Schwarz/block-Jacobi with respect 
%       to patches of elements around vertices; 
%    *) optimal step-sizes (lambda) in the error correction stage
% 
% Stopping criterion employed 
%  || F -’ A*x_approx || / ||F|| <= 10^(-5)
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
%    implicitly sets the numer of levels J in the hierarchy
% *) mesh_uniformity: bulk-chasing parameter determining the 
%    type of mesh hierarchy
%    uniformly-refined: 1; graded: smaller than 1;
% *) Hmax: initial mesh-size
% *) maxiter: to limit the maximum number of iterations
%
% OUTPUT:
% *) results: structure with fields:
%    *) x_approx: the last approximation computed by the solver
%    *) n_iter: number of iterations to reach the stopping criterion
%    *) rel_res: relative residual per iteration
%
%
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG


function [results] = p_robust_MG01_solver(problem_number, m, mesh_uniformity, Hmax, maxiter)

if nargin < 5
    maxiter = 100;
end
if nargin < 4
    Hmax = 0.1;
end
if nargin < 3
    mesh_uniformity = 1;    % uniform mesh refinement
end


%% initialize

initialize;

%% solver

fprintf('\n****************************************************\n');
fprintf('***\n');
fprintf('***      A-Posteriori-Steered MG(0,1) Solver       \n');
fprintf('***               (non-adaptive)   \n');
fprintf('***\n');
fprintf('***                  %s \n', solution.name);
fprintf('***            number of levels %d ',J-1);
fprintf('\n');
fprintf('***   polynomial order per level [');fprintf('%g ', m);fprintf(']\n');
fprintf('***\n');
fprintf('****************************************************\n\n');


x_approx = zeros(length(DOF2(J).freenodes),1);
res = F{J};
rel_res_init = 1;

% iterations
for iter = 1:maxiter

    % calling MG solver
    [rho_alg, ~] = smoother_bJ(x_approx,res);

    x_approx = x_approx + rho_alg;
    res = F{J} - A{J}*x_approx;
    relres(iter) = norm(res)/norm(F{J});
    fprintf('Iter %d:     relative residual norm = %e \n', iter, relres(iter));

    %stopping criterion
    if (relres(iter) <= 10^(-5)* rel_res_init)
        fprintf('Accuracy attained in iter %d \n', iter);
        break;
    end

end% of iterations
    
%% results

results.x_approx = x_approx;
results.n_iter = iter;
results.rel_res = relres;




