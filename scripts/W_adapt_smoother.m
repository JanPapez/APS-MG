function [rho, percentage_selected_iter] = W_adapt_smoother(xapprox,R,lam_local_rhonorms,threshold_value)
% block Jacobi smoother with no presmoothing and one local postsmoothing 
% steps on marked patches to be run after 'W_full_smoother'
% Includes the test if wRAS or AS smoother should be employed. On each 
% level, the optimal step size is computed.
% Details are given in [Miraci, Papez, Vohralik. Comput. Methods Appl. Math. (2021)]

% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG

global J A F;
global prolong_matrices patches meshdata;

percentage_selected_iter = zeros(J,1);

% algebraic residual
if nargin < 4
    R = F{J}-A{J}*xapprox;
end
rR = struct([]);
rR{J} = R;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RESTRICTION TO THE COARSEST GRID, no pre-smoothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
for j = J:-1:3
    rR{j-1} = ((rR{j}'* prolong_matrices(j).p) * prolong_matrices(j).h)';
end

%%%%%%%%%%%%%%%%%%%%%%% 
%COARSE GRID CORRECTION
%%%%%%%%%%%%%%%%%%%%%%%

j = 2;
if (lam_local_rhonorms{1} >= threshold_value)
    res = rR{j};
    R_0 = ((res'* prolong_matrices(j).p)* prolong_matrices(j).h)';
    temp_0 = A{1}\R_0;
    rho = prolong_matrices(j).p * (prolong_matrices(j).h * temp_0);
    percentage_selected_iter(1) = 100;
else
    rho = zeros(size(rR{2}));
end

%%%%%%%%%%%%%%%
%POST-SMOOTHING
%%%%%%%%%%%%%%%

for j = 2:J
    
    Aj = A{j};
    rRj = rR{j};
    
    if j ~= 2
        rho = prolong_matrices(j).p*(prolong_matrices(j).h*rho);
    end
    
    Arho = Aj*rho;
    rho_wRAS_j = zeros(size(rho));
    rho_dAS_j  = zeros(size(rho));
    
    
    I = find(lam_local_rhonorms{j} >= threshold_value);
    percentage_selected_iter(j) = length(I)/meshdata(j).nc*100;
    
    lhs = 0;
    wRAS_lhs = 0;
    
    if ~isempty(I)
        % smoothing on the marked patches only
        for index = I(:)'
            % indices of elements corresponding to the patch
            patchFEM0indices = patches(j).FEM0indices{index};
            weights = patches(j).hatfunction{index};

            localstiffmatrix = Aj(patchFEM0indices,patchFEM0indices);
            
            rhs_loc = (rRj(patchFEM0indices) - Arho(patchFEM0indices));
            localcoef = localstiffmatrix\rhs_loc;
            
            lhs = lhs + localcoef'*localstiffmatrix*localcoef;
            wRAS_lhs = wRAS_lhs + (localcoef.*weights)'*localstiffmatrix*(localcoef.*weights);
            
            rho_wRAS_j(patchFEM0indices) = rho_wRAS_j(patchFEM0indices) + localcoef.*weights;
            rho_dAS_j(patchFEM0indices)  = rho_dAS_j(patchFEM0indices)  + localcoef;
        end
        
        % test if wRAS or AS smoothing should be used (simpler in comparison to full step)
        lhs = sqrt( lhs/3 );
        
        Arho_wRAS_j = Aj*rho_wRAS_j;
        rhs = ( (rRj'* rho_wRAS_j) - (rho' * Arho_wRAS_j) ) / sqrt( rho_wRAS_j' * Arho_wRAS_j);
        
        wRAS_bool_1 = (lhs <= rhs);
        
        % use either wRAS or AS smoother
        if (wRAS_bool_1)
            rhoj = rho_wRAS_j;
            Arhoj = Arho_wRAS_j;
        else
            disp('AS is used instead of wRAS in adaptive V-step')
            rhoj = rho_dAS_j;
            Arhoj = Aj*rhoj;
        end
        
        % optimal stepsize
        rhojenergynormsquared = rhoj'* Arhoj;
        lambdaj = ( (rRj'* rhoj) - (rho' * Arhoj) ) / rhojenergynormsquared;
        
        % updating the solution
        rho = rho + lambdaj*rhoj;
    else
        %do not update rho <- no patches marked on this level
    end 
    
end
