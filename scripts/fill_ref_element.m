function [ref_element] = fill_ref_element(DOF)
% function setting, for the reference simplex, the quadrature nodes and 
%   weights, and evaluating the values and gradients of Dubiner FEM basis functions 
%   in quadrature nodes 
%
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG


% computing quadrature nodes and weights
[X,Y,Weights] = simplexquad(2*DOF.m+2,2);
s = length(Weights);

% values and gradients of lagrangean FEM basis functions in quadrature nodes
[FEMvalue_temp, FEMgrad_temp] = phikl_all(X,Y,DOF.m);
FEMvalue = FEMvalue_temp(:);
FEMgrad = [reshape(FEMgrad_temp(:,1:2:end), s*DOF.nFEM, 1), reshape(FEMgrad_temp(:,2:2:end), s*DOF.nFEM, 1)];

ref_element.quadrature.X = X;
ref_element.quadrature.Y = Y;
ref_element.quadrature.Weights = Weights;
ref_element.quadrature.s = s;

ref_element.FEMvalue = FEMvalue;
ref_element.FEMgrad = FEMgrad;

% construction of reference FEM mass matrix (psi_k, psi_ell)_hatK
FEMvalue_reshaped = reshape(ref_element.FEMvalue, s, DOF.nFEM);
ref_element.FEMmassmatrix = (FEMvalue_reshaped'*diag(Weights))*FEMvalue_reshaped;

end