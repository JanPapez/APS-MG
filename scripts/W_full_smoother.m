function [rho, local_rhonorms, lambda_j, rho_j, rho_j_loc] = W_full_smoother(xapprox,R)
% block Jacobi smoother with no presmoothing and one postsmoothing steps
% Includes the test if wRAS or AS smoother should be employed and saves
% local quantities to be used for marking patches for adaptive local
% refinement. On each level, the optimal step size is computed.
% Details are given in [Miraci, Papez, Vohralik. Comput. Methods Appl. Math. (2021)]

% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG

global J A F;
global extended prolong_matrices;

local_rhonorms = struct([]);
lambda_j = zeros(J,1);
rho_j = struct([]);
rho_j_loc = struct([]);

% algebraic residual
if nargin < 2
    R = F{J}-A{J}*xapprox;
end
rR = struct([]);
rR{J} = R;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RESTRICTION TO THE COARSEST GRID, no pre-smoothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
for j = J:-1:2
    rRj = rR{j};
    rR{j-1} = ((rRj'* prolong_matrices(j).p) * prolong_matrices(j).h)';
end
 
%%%%%%%%%%%%%%%%%%%%%%% 
%COARSE GRID CORRECTION
%%%%%%%%%%%%%%%%%%%%%%%

R_0 = rR{1};
rho = A{1}\R_0;

local_rhonorms{1} = rho'*R_0;
lambda_j(1) = 1;
rho_j{1} = rho;

%%%%%%%%%%%%%%%
%POST-SMOOTHING
%%%%%%%%%%%%%%%

for j = 2:J
    Aj = A{j};
    rRj = rR{j};
    
    rho = prolong_matrices(j).p*(prolong_matrices(j).h*rho);
    
    % extend the vectors
    Arho = Aj*rho;
    RHS = rRj - Arho;
    RHS_ext = RHS(extended(j).ext_indices);
    
    % local solve
    localcoef_ext = (extended(j).L_chol_ext')\RHS_ext;
    localcoef_ext = extended(j).L_chol_ext\localcoef_ext;
    
    % compute the norms || \nabla rho_j,a ||_{omega_a}
    tmp = cumsum(RHS_ext.*localcoef_ext);
    tmp2 = cumsum(extended(j).local_sizes);
    initzeros = sum(tmp2 == 0);
    tmp2 = tmp2(tmp2>0);
    tmp3 = tmp(tmp2);
    local_rhonorms{j} = diff([zeros(initzeros+1,1); tmp3]);
    
    % compute the norms || \nabla I_pj ( hat_j,a rho_j,a) ||_{omega_a}
    wRAS_tmp = cumsum(RHS_ext.*(localcoef_ext.*extended(j).psiaweight_local_ext));
    wRAS_tmp2 = cumsum(extended(j).local_sizes);
    wRAS_initzeros = sum(wRAS_tmp2 == 0);
    wRAS_tmp2 = wRAS_tmp2(wRAS_tmp2>0);
    wRAS_tmp3 = wRAS_tmp(wRAS_tmp2);
    local_wRAS_rhonorms{j} = diff([zeros(wRAS_initzeros+1,1); wRAS_tmp3]);
    
    % test if wRAS or AS smoothing should be used
    lhs = sqrt(sum(local_rhonorms{j})/3);
    
    rho_j_loc{j} = localcoef_ext.*extended(j).psiaweight_local_ext;
    rho_wRAS_j = accumarray(extended(j).ext_indices, rho_j_loc{j});
    Arho_wRAS_j = Aj*rho_wRAS_j;
    rhs = ( (rRj'* rho_wRAS_j) - (rho' * Arho_wRAS_j) ) / sqrt( rho_wRAS_j' * Arho_wRAS_j);
    
    wRAS_bool_1 = (lhs <= rhs);
    
    sum_rhoa_sq = sum(local_rhonorms{j});
    sum_wRAS_rhoa_sq = sum(local_wRAS_rhonorms{j});
    
    wRAS_bool_2 = (sum_wRAS_rhoa_sq <= sum_rhoa_sq);
    
    % use either wRAS or AS smoother
    if (wRAS_bool_1 && wRAS_bool_2)
        rhoj = rho_wRAS_j;
        Arhoj = Arho_wRAS_j;
    else
        disp('AS is used instead of wRAS in full V-step')
        rho_j_loc{j} = localcoef_ext;
        rhoj = accumarray(extended(j).ext_indices, localcoef_ext);
        Arhoj = Aj*rhoj;
    end
    
    % optimal stepsize
    rhojenergynormsquared = rhoj'* Arhoj;
    lambda_j(j) = ( (rRj'* rhoj) - (rho' * Arhoj) ) / rhojenergynormsquared;
    rho_j{j} = rhoj;
    
    % updating the solution
    rho = rho + lambda_j(j)*rhoj;
end





