function [Iout] = p_interpolation(j)
% polynomial interpolation matrix from degree m_(j-1) to m_j
%
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG

global DOF1 DOF2 meshdata

m_j = DOF2(j).m;
m_prevj = DOF2(j-1).m;
nFEM = (m_j+1)*(m_j+2)/2;

% non-lagrangean grid
coordinates_temp = nonlagrangeannodes(m_j)';

%local interpolation matrix
[I1_loc, ~] = phikl_all(coordinates_temp(1,:)',coordinates_temp(2,:)',m_prevj);

% assembling the global interpolation matrix

ne = meshdata(j).ne;

interp.row = zeros(ne*3*nFEM,1);
interp.column = interp.row;
interp.values = interp.row;
interp.values2 = interp.row;
interp.currindex = 0;

% go over the elements
for ie = 1:ne
    ie_prev = ie;
    interp = addblocktosparse(DOF2(j).elementFEM(:,ie),DOF1(j).elementFEM(:,ie_prev),I1_loc,interp);
end

I = sparse(interp.row, interp.column, interp.values);
nI = sparse(interp.row, interp.column, interp.values2);

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
