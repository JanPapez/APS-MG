function [Iout,refine] = interpolation_matrix_varp_NVB(j, refined)
% function for constructing the interpolation matrix between two
% consecutive levels, with refinement corresponding to refinemesh_mine_nvb
%
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG


global DOF1 DOF2 meshdata

m = DOF1(j).m;
nFEM = (m+1)*(m+2)/2;

% nonlagrangean grid
coordinates_temp = nonlagrangeannodes(m)';

%% lagrange grid on the element refined into 4 triangles
J_1 = 0.5;
coordinates_1 = J_1*coordinates_temp;
J_2 = [-0.5, -0.5; 1, 0.5];
coordinates_2 = repmat([0.5;0],1,nFEM) + J_2*coordinates_temp;
J_3 = [0.5, 0.5; -1 -0.5];
coordinates_3 = repmat([0;1],1,nFEM) + J_3*coordinates_temp;
J_4 = 0.5;
coordinates_4 = repmat([0.5;0],1,nFEM) + J_4*coordinates_temp;
coordinates_all = [coordinates_1, coordinates_2, coordinates_3, coordinates_4];

[I1_loc, ~] = phikl_all(coordinates_all(1,:)',coordinates_all(2,:)',m);

%% lagrange grid on the element refined into 3 triangles (edges 2+3 refined)
coordinates_1 = coordinates_temp*0.5;
coordinates_2 = repmat([1; 0],1,nFEM) + [-1;1]*coordinates_temp(1,:) + [-.5;0]*coordinates_temp(2,:);
coordinates_3 = repmat([.5;0],1,nFEM) + [-.5;1]*coordinates_temp(1,:) + [-.5;.5]*coordinates_temp(2,:);
coordinates_all = [coordinates_1, coordinates_2, coordinates_3];

[I2_loc, ~] = phikl_all(coordinates_all(1,:)',coordinates_all(2,:)',m);

%% lagrange grid on the element refined into 3 triangles (edges 1+3 refined)
coordinates_1 = repmat([.5;0],1,nFEM) + coordinates_temp*0.5;
coordinates_2 = repmat([0;1],1,nFEM) + [0;-1]*coordinates_temp(1,:) + [.5; -1]*coordinates_temp(2,:);
coordinates_3 = repmat([0;1],1,nFEM) + [.5;-1]*coordinates_temp(1,:) + [.5;-.5]*coordinates_temp(2,:);
coordinates_all = [coordinates_1, coordinates_2, coordinates_3];

[I3_loc, ~] = phikl_all(coordinates_all(1,:)',coordinates_all(2,:)',m);

%% lagrange grid on the element refined into 2 triangles (edge 3 refined)
coordinates_1 = repmat([0;1],1,nFEM) + [0;-1]*coordinates_temp(1,:) + [.5;-1]*coordinates_temp(2,:);
coordinates_2 = repmat([1;0],1,nFEM) + [-1;1]*coordinates_temp(1,:) + [-.5;0]*coordinates_temp(2,:);
coordinates_all = [coordinates_1, coordinates_2];

[I4_loc, ~] = phikl_all(coordinates_all(1,:)',coordinates_all(2,:)',m);

%% assembling the global interpolation matrix

ne = meshdata(j-1).ne;
refine = zeros(ne, 5);

interp.row = zeros(ne*nFEM*nFEM*4,1);
interp.column = interp.row;
interp.values = interp.row;
interp.values2 = interp.row;
interp.currindex = 0;

index = 0;
% not refined elements
elem = refined{1}; nelem0 = length(elem);
for ie = 1:nelem0
    ie_prev = elem(ie);
    interp = addblocktosparse(DOF1(j).elementFEM(:,ie),DOF2(j-1).elementFEM(:,ie_prev),eye(nFEM),interp);
end
refine(elem,1) = 0;
refine(elem,2) = 1:nelem0;

index = index + nelem0;

%elements refined into 4 triangles
elem = refined{2}; nelem1 = length(elem);
for ie = 1:nelem1
    ie_prev = elem(ie);
    temp = [DOF1(j).elementFEM(:,index+ie); DOF1(j).elementFEM(:,index+ie+nelem1); DOF1(j).elementFEM(:,index+ie+2*nelem1); DOF1(j).elementFEM(:,index+ie+3*nelem1)];
    interp = addblocktosparse(temp,DOF2(j-1).elementFEM(:,ie_prev),I1_loc,interp);
end
refine(elem,1) = 1;
refine(elem,2:5) = reshape((index+1):(index+4*nelem1),nelem1,4);

index = index + 4*nelem1;

%elements refined into 3 triangles (edges 2+3 refined)
elem = refined{3}; nelem2 = length(elem);
for ie = 1:nelem2
    ie_prev = elem(ie);
    temp = [DOF1(j).elementFEM(:,index+ie); DOF1(j).elementFEM(:,index+ie+nelem2); DOF1(j).elementFEM(:,index+ie+2*nelem2)];
    interp = addblocktosparse(temp,DOF2(j-1).elementFEM(:,ie_prev),I2_loc,interp);
end
refine(elem,1) = 2;
refine(elem,2:4) = reshape(index+1:(index+3*nelem2),nelem2,3);

index = index + 3*nelem2;

%elements refined into 3 triangles (edges 1+3 refined)
elem = refined{4}; nelem3 = length(elem);
for ie = 1:nelem3
    ie_prev = elem(ie);
    temp = [DOF1(j).elementFEM(:,index+ie); DOF1(j).elementFEM(:,index+ie+nelem3); DOF1(j).elementFEM(:,index+ie+2*nelem3)];
    interp = addblocktosparse(temp,DOF2(j-1).elementFEM(:,ie_prev),I3_loc,interp);
end
refine(elem,1) = 3;
refine(elem,2:4) = reshape(index+1:(index+3*nelem3),nelem3,3);

index = index + 3*nelem3;

%elements refined into 2 triangles (edge 3 refined)
elem = refined{5}; nelem4 = length(elem);
for ie = 1:nelem4
    ie_prev = elem(ie);
    temp = [DOF1(j).elementFEM(:,index+ie); DOF1(j).elementFEM(:,index+ie+nelem4)];
    interp = addblocktosparse(temp,DOF2(j-1).elementFEM(:,ie_prev),I4_loc,interp);
end
refine(elem,1) = 4;
refine(elem,2:3) = reshape(index+1:(index+2*nelem4),nelem4,2);

I = sparse(interp.row(1:interp.currindex), interp.column(1:interp.currindex), interp.values(1:interp.currindex));
nI = sparse(interp.row(1:interp.currindex), interp.column(1:interp.currindex), interp.values2(1:interp.currindex));

Iout = dividesparsematrices(I,nI);

end

function [spmatrix] = addblocktosparse(rowindex, columnindex, block, spmatrix)
    rowindex = rowindex(:);
    nrows = length(rowindex);
    columnindex = columnindex(:);
    ncolumns = length(columnindex);
    block = full(block);
    rowindex = repmat(rowindex,1,ncolumns);
    columnindex = repmat(columnindex',nrows,1);
    spmatrix.row(spmatrix.currindex + (1:(nrows*ncolumns))) = rowindex(:);
    spmatrix.column(spmatrix.currindex + (1:(nrows*ncolumns))) = columnindex(:);
    spmatrix.values(spmatrix.currindex + (1:(nrows*ncolumns))) = block(:);
    spmatrix.values2(spmatrix.currindex + (1:(nrows*ncolumns))) = all(block(:),2);
    spmatrix.currindex = spmatrix.currindex + (nrows*ncolumns);
end

function [C] = dividesparsematrices(A,B)
    C = A.*spfun(@(x) 1./x, B);
end
