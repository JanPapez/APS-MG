% p_robust_MG0adapt_solver runs the local adaptive smoothing adaptive 
% version of a-posteriori-steered multigrid (APS-MG) solver developed 
% in [Miraci, Papez, Vohralik. Comput. Methods Appl. Math. (2021)] 
% https://doi.org/10.1515/cmam-2020-0024 
%
% Linear system: A*x_exact = F
%
% Solver features:
% *) initial guess is x_approx = x_0 = 0 and each iteration consists of:
%    *) A first MG V-cycle (full-smoothing substep) with zero pre- 
%       and one post-smoothing step, P1-coarse solve 
%    *) A second MG V-cycle (adaptive-smoothing substep) with zero pre- 
%       and one post-smoothing step ONLY locally on adaptively chosen 
%       patches, P1-coarse solve        
%    *) smoothing is adaptively chosen between additive Schwarz (AS)
%       and weighted restricted additive Schwarz (wRAS) with respect 
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
% *) gamma: parameter between 0 and 1 in the theoretical 
%    conditional statement determining whether the 
%    adaptive-smoothing substep takes place or not
% *) theta: adaptivity-parameter determining the portion
%    of patches with high-estimated algebraic error to be selected 
%    for local smoothing
% *) maxiter: to limit the maximum number of iterations
% *) plot: whether the error and estimated error distribution 
%    per patches should be plotted; first entry 1 if plot should 
%    occur, second entry on which iteration the plot is done
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


function [results] = p_robust_MG01_locallyadaptive_solver(problem_number, m, mesh_uniformity, Hmax, gamma, theta, maxiter, plot)
if nargin < 8
    plot = 0;
end
if nargin < 7
    maxiter = 100;
end
if nargin < 6
    theta = 0.95; % unsquared
end
if nargin < 5
    gamma = 0.7; % unsquared
end
if nargin < 4
    Hmax = 0.1;
end
if nargin < 3
    mesh_uniformity = 1;    % uniform mesh refinement
end

threshold_theta = theta^2;
threshold_gamma = gamma^2;
threshold_R = 5;    % for lambda_j \le R condition


% if plotting of the marked patches in some iteration is required
plot_patches = plot(1);
if plot_patches
   plot_algerror_initer = plot(2);
end


%% initialize
initialize;
percentage_selected = zeros(maxiter,J);

% only this part differs from initialization for other solvers
for j = 2:J
    extended(j).local_sizes_cumsum = [1, cumsum(extended(j).local_sizes)];
end


if plot_patches
addpath('scripts_plotting');
end

%% solver

fprintf('\n************************************************************\n');
fprintf('***\n');
fprintf('***      A-Posteriori-Steered MG(0,1) Solver       \n');
fprintf('***  (adaptive local smoothing, theta = %.2f, gamma = %.1f',theta,gamma);fprintf(')\n');
fprintf('***\n');
fprintf('***                 %s \n', solution.name);
fprintf('***            number of levels %d ',J-1);fprintf('\n');
fprintf('***   polynomial order per level [');fprintf('%g ', m);fprintf(']\n');
fprintf('***\n');
fprintf('************************************************************\n');


x_approx = zeros(length(DOF2(J).freenodes),1);
res = F{J};
rel_res_init = 1;

% iterations
for iter = 1:maxiter
    
    % FULL SOLVE
    [rho_alg, local_rhonorms, lambda_j, rho_j, rho_j_loc] = ...
        W_full_smoother(x_approx,res);
    
    %for plotting marked patches
    if plot_patches == 1
        x_approx_prev = x_approx;
        res_prev = res;
    end
    
    %new iterate
    x_approx = x_approx + rho_alg;
    res = F{J} - A{J}*x_approx;
    relres = norm(res)/norm(F{J});
    
    if threshold_theta > 0 % adaptive step
        
        fprintf('Iter %d+1/2: relative residual norm = %e \n', iter-1, relres);
        
        % MARKING
        for lev = 1:J
            lam_local_rhonorms{lev} = local_rhonorms{lev}.*lambda_j(lev);
        end
        
        tmp = cat(1,lam_local_rhonorms{:});
        
        if threshold_theta == 1
            threshold_value = 0;
        else
            [tmp,~] = sort(tmp,'descend');
            sumeta = cumsum(tmp);
            ell = find(sumeta >= sumeta(end)*threshold_theta,1);
            threshold_value = (tmp(ell)+tmp(ell+1))/2;
        end
        
        % here comes the plotting of the marked patches if needed
        if plot_patches == 1
            plot_marked_patches_for_W;
        end
        
        % TEST IF THE ADAPTIVE STEP SHOULD BE PERFORMED
        if threshold_theta == 1
            Test2_rhs = sum(tmp);
        else
            Test2_rhs = sumeta(ell);
        end
        
        if (local_rhonorms{1} >= threshold_value)
            rho_selected = rho_j{1};
            Test2_lhs = full(rho_j{1}'*(A{1}*rho_selected));
        else
            rho_selected = zeros(size(F{1}));
            Test2_lhs = 0;
        end
        
        for j = 2:J
            rho_selected = prolong_matrices(j).p*(prolong_matrices(j).h*rho_selected);
            I = find(lam_local_rhonorms{j} >= threshold_value);
            
            if ~isempty(I)
                loc_index = [];
                for jjj = 1:length(I)
                    loc_index = [loc_index, ...
                        extended(j).local_sizes_cumsum(I(jjj))+1:extended(j).local_sizes_cumsum(I(jjj))+extended(j).local_sizes(I(jjj))];
                end
                rho_selected_j = accumarray(extended(j).ext_indices(loc_index), rho_j_loc{j}(loc_index), [length(rho_selected),1]);
                rho_selected = rho_selected + lambda_j(j)*rho_selected_j;
                Test2_lhs = Test2_lhs + lambda_j(j)*(rho_j{j}'*(A{j}*rho_selected) );
            end
        end
        
        Test2_lambda = min(lambda_j <= threshold_R);
        
        % ADAPTIVE SOLVE
        if ((Test2_lhs <= threshold_gamma * Test2_rhs) && Test2_lambda)
            fprintf('            cond.statement PASSED: alpha = %e \n', Test2_lhs/Test2_rhs);
            [rho_alg, percentage_selected(iter,:)] = ...
                W_adapt_smoother(x_approx,res,lam_local_rhonorms,threshold_value);
            
            %new iterate
            x_approx = x_approx + rho_alg;
            res = F{J} - A{J}*x_approx;
            relres = norm(res)/norm(F{J});
        else
            fprintf('            cond.statement NOT PASSED: alpha = %e \n', Test2_lhs/Test2_rhs);
        end        
    end
    
    fprintf('Iter %d:     relative residual norm = %e \n', iter, relres);
    
    %stopping criterion
    if (relres <= 10^(-5)* rel_res_init)
        fprintf('Accuracy attained in iter %d \n', iter);
        percentage_selected = percentage_selected(1:iter,:);
        break;
    end
       
end% of iterations

results.x_approx = x_approx;
results.n_iter = iter;
results.rel_res = relres;




