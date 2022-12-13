function [rhoprev,levelwisenorms,smoothiters] = smoother_bJ_adaptnumber(xapprox,R,theta,maxlocsmooths)
% block Jacobi smoother with no presmoothing and adaptively chosen number 
% of postsmoothing steps;
% The information about the blocks and the block diagonal matrix are given
% in the global variable extended. The blocks overlap because of 
% overlapping of patches. 
% On each level, the optimal step size is computed.
%
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG


global J A F;
global extended prolong_matrices;

if nargin < 4
    maxlocsmooths = 5;
end

% algebraic residual
if isempty(R), R = F{J}-A{J}*xapprox; end

rR = struct([]);
rR{J} = R;

levelwisenorms = zeros(1,J);
sum_levelwisenorms = 0;
smoothiters = ones(1,J);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RESTRICTION TO THE COARSEST GRID, no pre-smoothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
for j = J:-1:2
   
    rRj = rR{j};
    
    if j ~= 2
        rR{j-1} = ((rRj'* prolong_matrices(j).p) * prolong_matrices(j).h)';
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%% 
%COARSE GRID CORRECTION
%%%%%%%%%%%%%%%%%%%%%%%

res = rR{2};
R_0 = ((res'* prolong_matrices(j).p)* prolong_matrices(j).h)';
temp_0 = A{1}\R_0;
rhoprev = prolong_matrices(j).p * (prolong_matrices(j).h * temp_0);

levelwisenorms(1) = temp_0'*R_0;

sum_levelwisenorms = levelwisenorms(1);

%%%%%%%%%%%%%%%
%POST-SMOOTHING one step
%%%%%%%%%%%%%%%

for j = 2:J
    
    Aj = A{j};
    rRj = rR{j};
    
    if j ~= 2
        rhoprev = prolong_matrices(j).p*(prolong_matrices(j).h*rhoprev);
    end
    
    for smoothonlevel = 1:maxlocsmooths
        smoothiters(j) = smoothonlevel;
        Arhoprev = Aj*rhoprev;

        % extend the vectors
        RHS = rRj - Arhoprev;
        RHS_ext = RHS(extended(j).ext_indices);

        localcoef_ext = extended(j).L_chol_ext\(extended(j).L_chol_ext'\RHS_ext);
        
        % optimal stepsize
        rhoj =  accumarray(extended(j).ext_indices, localcoef_ext);
        Arhoj = Aj * rhoj;

        rhojenergynormsquared = rhoj'* Arhoj;
        lambdaj = ( (rRj'* rhoj) - (rhoprev' * Arhoj) ) / rhojenergynormsquared;

        tmp = lambdaj^2*rhojenergynormsquared;
        levelwisenorms(j) = levelwisenorms(j) + tmp;

        % updating the solution
        rhoprev = rhoprev + lambdaj*rhoj;
        
        %check if we smooth again
        if tmp < (theta^2)*sum_levelwisenorms
            sum_levelwisenorms = sum_levelwisenorms + tmp;
            break; 
        end
        sum_levelwisenorms = sum_levelwisenorms + tmp;
              
    end
end


