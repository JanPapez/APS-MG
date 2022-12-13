function [rhoprev,local_rhonorms] = compute_rhotilde(xapprox,R)
% function for computing the local decomposition of the error over the 
% patches and levels
% used to be compared with the local distribution of error indicators
%
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG

global J A F;
global extended prolong_matrices;

local_rhonorms = struct([]);

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
    
    if j ~= 2
        rR{j-1} = ((rRj'* prolong_matrices(j).p) * prolong_matrices(j).h)';%@@@
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%% 
%COARSE GRID CORRECTION
%%%%%%%%%%%%%%%%%%%%%%%

res = rR{2};
R_0 = ((res'* prolong_matrices(j).p)* prolong_matrices(j).h)';%@@@
temp_0 = A{1}\R_0;
rhoprev = prolong_matrices(j).p * (prolong_matrices(j).h * temp_0);

local_rhonorms{1} = temp_0'*R_0;

%%%%%%%%%%%%%%%
%POST-SMOOTHING
%%%%%%%%%%%%%%%

for j = 2:J
    
    Aj = A{j};
    rRj = rR{j};
       
    if j ~= 2
        rhoprev = prolong_matrices(j).p*(prolong_matrices(j).h*rhoprev); %@@@
    end
    
    RHS = rRj - Aj*rhoprev;
    rhotildej = Aj\RHS;
    
    RHS_ext = RHS(extended(j).ext_indices);
    localcoef_ext = rhotildej(extended(j).ext_indices);
    
    %compute the norms || \nabla rho_j,a ||_{omega_a}
    tmp = cumsum(RHS_ext.*localcoef_ext);
    tmp2 = cumsum(extended(j).local_sizes);
    initzeros = sum(tmp2 == 0);
    tmp2 = tmp2(tmp2>0);
    tmp3 = tmp(tmp2);
    local_rhonorms{j} = diff([zeros(initzeros+1,1);tmp3]);
    
    % updating the rho
    rhoprev = rhoprev + rhotildej;
    
end



