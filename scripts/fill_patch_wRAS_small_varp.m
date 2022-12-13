function [loc_mat, patches, hats] = fill_patch_wRAS_small_varp
% function for gathering the local information associated with patches;
% namely: value of hat function in DOFs associated with the patch ('hats'),
% global index of DOFs on the patch ('patches'), and local stiffness
% matrices ('loc_mat')
%
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG

global meshdata DOF2 J A ;


for j = 2:J
    
    m = DOF2(j).m;
    nFEM = DOF2(j).nFEM;
    
    % FEM element basis indices not on the edge
    not_edge = zeros(3,nFEM-(m+1+m));
    not_edge(1,:) = setdiff(1:nFEM, unique([cumsum(m+1:-1:1), 1:m+1]));
    not_edge(2,:) = setdiff(1:nFEM, unique([1:m+1, (nFEM+1) - cumsum(1:m+1)]));
    not_edge(3,:) = setdiff(1:nFEM, unique([(nFEM+1) - cumsum(1:m+1), cumsum(m+1:-1:1)]));
    
    % evaluating the hat function in points associated with DOFs
    nodes = nonlagrangeannodes(m);
    X = nodes(:,1);
    Y = nodes(:,2);
    
    [phia_finepatch, ~] = phikl_all(X,Y,1);
    phia_finepatch = phia_finepatch(:,[1 3 2]);
    
    
    nc_j = meshdata(j).nc;
    Aj = A{j};
    
    psiaweight = zeros(DOF2(j).lastFEMindex, 1);
    
    elements_of_patch_all = struct('element',cell(1,meshdata(j).nc));
    for ell = 1:size(meshdata(j).elements,2)
        for ell2 = 1:3
            elements_of_patch_all(meshdata(j).elements(ell2,ell)).element(end+1) = ell;
        end
    end
    
    loc = struct([]);
    hat = struct([]);
    patch = struct([]);
    
    
    for index = 1:nc_j
        % indices of elements corresponding to the patch
        elements_of_patch = elements_of_patch_all(index).element;
        patchFEM0indices = index;
        
        for el = 1:length(elements_of_patch)
            
            ie = elements_of_patch(el);
            which_vertex = find(meshdata(j).elements(1:3,ie) == index);
            
            DOFslocal = DOF2(j).elementFEM(:, ie);
            psiaweight(DOFslocal(:)) = phia_finepatch(:, which_vertex);
            
            elementFEM0indices = DOF2(j).elementFEM(not_edge(which_vertex,:),ie);
            
            patchFEM0indices = [patchFEM0indices; elementFEM0indices];
            
        end
        
        psiaweight_local = psiaweight(patchFEM0indices);
        patchFEM0indices = DOF2(j).freenodesnumbering_toglobal(patchFEM0indices);
        psiaweight_local = psiaweight_local(patchFEM0indices>eps);
        patchFEM0indices = patchFEM0indices(patchFEM0indices>eps);
        
        
        hat(index).vertex = psiaweight_local;
        
        patch(index).vertex =  patchFEM0indices;
        
        localstiffmatrix = Aj(patchFEM0indices,patchFEM0indices);
        [T,p] = chol(localstiffmatrix);
        if p ~= 0 %checking if Hermitian
            disp('error');
        end
        
        loc(index).chol = T;
        
    end
    hats(j).level = hat;
    patches(j).level = patch;
    loc_mat(j).level = loc;
end

end