function [A,F,localstiffmatrices,dJloc_all] = globalstiffmatrix_varp(j)
% function for evaluating the local stiffness matrices and assembling the
%   global one
% 
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG

global meshdata DOF2 ref_element solution;


ref_elem = ref_element(j).level;
fh_coef = solution.fh_coef{j};

nFEM = DOF2(j).nFEM;

locI = 1:nFEM;
locJ = reshape(locI(ones(nFEM,1),:),nFEM*nFEM,1);
locI = reshape(locI(ones(nFEM,1),:)',nFEM*nFEM,1);

I = reshape(DOF2(j).elementFEM(locI,:),nFEM*nFEM*meshdata(j).ne,1);
J = reshape(DOF2(j).elementFEM(locJ,:),nFEM*nFEM*meshdata(j).ne,1);
Avalues = zeros(size(J));

IF = reshape(DOF2(j).elementFEM(1:nFEM,:),nFEM*meshdata(j).ne,1);
Fvalues = zeros(size(IF));

localstiffmatrices = zeros(nFEM,nFEM,meshdata(j).ne);
dJloc_all = zeros(1, meshdata(j).ne);

% quadrature weights
W = ref_elem.quadrature.Weights;
s = ref_elem.quadrature.s;

FEMgrad = ref_elem.FEMgrad;
FEMmassmatrix = ref_elem.FEMmassmatrix;

for ie = 1:meshdata(j).ne
    
    tenzorsqrt = solution.tenzorsqrt{meshdata(j).elements(4,ie)};
    
    coordinates = meshdata(j).coord(:,meshdata(j).elements(1:3,ie));
    Jloc = [coordinates(:,2)-coordinates(:,1), coordinates(:,3)-coordinates(:,1)];
    dJloc = abs(det(Jloc));
    Wloc = W*dJloc;
    
    locFEMgrad = (FEMgrad/Jloc)*tenzorsqrt';
    
    locFEMgrad_reshaped = zeros(s,2*nFEM);
    locFEMgrad_reshaped(:,1:2:end-1) = reshape(locFEMgrad(:,1),s,nFEM);
    locFEMgrad_reshaped(:,2:2:end) = reshape(locFEMgrad(:,2),s,nFEM);
    
    Ltemp = locFEMgrad_reshaped'*(diag(Wloc)*locFEMgrad_reshaped);
    L = Ltemp(1:2:end,1:2:end) + Ltemp(2:2:end,2:2:end);
    
    Floc = FEMmassmatrix*dJloc*(fh_coef(:,ie));
    
    Avalues((nFEM*nFEM*(ie-1)+1):(nFEM*nFEM*ie)) = L(:);
    Fvalues((nFEM*(ie-1)+1):(nFEM*ie)) = Floc;
     
    localstiffmatrices(:,:,ie) = L;
    dJloc_all(ie) = dJloc;
end

A = sparse(I,J,Avalues);
F = sparse(IF,ones(size(IF)),Fvalues);

% incorporating the Dirichlet boundary condition
xD = DOF2(j).dirichletvalues;

F = F-A*xD;

end
