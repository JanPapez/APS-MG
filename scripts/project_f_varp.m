function [fh_coef, osc] = project_f_varp(j)
% function for interpolating the rhs function f as piece-wise polynomial f_h
% assuring
%       (f, v_h) = (f_h, v_h)   for all v_h in V_h
% and computing the oscillation terms
%       osc_K := (f - f_h, f - f_h)_K    K in T_h
%
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG


global meshdata DOF2 ref_element solution;

ref_elem = ref_element(j).level;

fh_coef = zeros(DOF2(j).nFEM,meshdata(j).ne);
osc = zeros(meshdata(j).ne,1);

% quadrature nodes and weights
X = ref_elem.quadrature.X;
Y = ref_elem.quadrature.Y;
W = ref_elem.quadrature.Weights;
s = ref_elem.quadrature.s;

FEMvalue = reshape(ref_elem.FEMvalue, s, DOF2(j).nFEM);

invFEMmassmatrix = inv(ref_elem.FEMmassmatrix);

for ie = 1:meshdata(j).ne
    
    coordinates = meshdata(j).coord(:,meshdata(j).elements(1:3,ie));    
    Jloc = [coordinates(:,2)-coordinates(:,1), coordinates(:,3)-coordinates(:,1)];
    dJloc = abs(det(Jloc));
    
    quad_nodes = repmat(coordinates(:,1),1,s) + Jloc*([X,Y])';
    floc = solution.f(quad_nodes(1,:),quad_nodes(2,:));
    
    temp = invFEMmassmatrix* ((floc*diag(W))*FEMvalue)';
    fh_coef(:,ie) = temp;

    osc(ie) = (W*dJloc)'*((floc' - FEMvalue*temp).^2);
end