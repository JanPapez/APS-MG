function [discrerror_K, toterror_K, algerror_K] = evaluate_error_on_elements(x,xapprox,j)
% function for evaluating the norm of the discretization, total and
%   algebraic errors on mesh elements
% 
% INPUT: "exact" solution of the algebraic system
%        approximate vector
%        knowledge of true infinite-dimensional solution is assumed to
%           evaluate the discretization error
% 
% OUPUT: vector of || \nabla( error ) ||_K
% 
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG


global meshdata DOF2 ref_element solution

discrerror_K = zeros(meshdata(j).ne,1);
algerror_K   = zeros(meshdata(j).ne,1);
toterror_K   = zeros(meshdata(j).ne,1);

ref_elem = ref_element{j};

% quadrature nodes and weights
X = ref_elem.quadrature.X;
Y = ref_elem.quadrature.Y;
W = ref_elem.quadrature.Weights;
s = ref_elem.quadrature.s;

nFEM = DOF2(j).nFEM;

for ie = 1:meshdata(j).ne
    
    locvertices = meshdata(j).elements(1:3,ie);
    locFEMindices = DOF2(j).elementFEM(:,ie);
    
    coordinates = meshdata(j).coord(:,locvertices);
    Jloc = [coordinates(:,2)-coordinates(:,1), coordinates(:,3)-coordinates(:,1)];
    dJloc = abs(det(Jloc));
    Wloc = W*dJloc;
    
    quad_nodes = repmat(coordinates(:,1),1,s) + Jloc*([X,Y])';
    
    locFEMgrad = ref_elem.FEMgrad/Jloc;
    
    uhgrad_inqnodes  = [reshape(locFEMgrad(:,1),s,nFEM)*x(locFEMindices), ...
        reshape(locFEMgrad(:,2),s,nFEM)*x(locFEMindices)];
    
    uhigrad_inqnodes = [reshape(locFEMgrad(:,1),s,nFEM)*xapprox(locFEMindices), ...
        reshape(locFEMgrad(:,2),s,nFEM)*xapprox(locFEMindices)];
    
    ugrad_inqnodes = solution.ugrad(quad_nodes(1,:)',quad_nodes(2,:)');
    
    discr_error = uhgrad_inqnodes - ugrad_inqnodes;
    alg_error = uhigrad_inqnodes - uhgrad_inqnodes;
    tot_error =  uhigrad_inqnodes - ugrad_inqnodes;
    
    dotproduct = sum(discr_error.*discr_error,2);
    discrerror_K(ie) = sqrt(Wloc'*dotproduct);
    
    dotproduct = sum(alg_error.*alg_error,2);
    algerror_K(ie) = sqrt(Wloc'*dotproduct);
    
    dotproduct = sum(tot_error.*tot_error,2);
    toterror_K(ie) = sqrt(Wloc'*dotproduct);
    
end