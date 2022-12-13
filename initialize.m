% script intializing the sequence of meshes, assembling the systems, and 
%  other data structures for p-robust MG solvers
%
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG


global solution meshdata DOF1 DOF2 ref_element element A F x;
global extended patches hats loc_mat prolong_matrices J;

J = length(m);

addpath('scripts');

% choice of the problem
solution = solution_def(problem_number);

% initial mesh
fill_meshdata(Hmax,1);
meshdata(1).I = [];

DOF = struct([]); ref_element = struct([]); osc = struct([]); A = struct([]); F = struct([]); x = struct([]);
element = struct([]);

fill_DOF1(meshdata,m(1),1); 

%% build the mesh hierarchy
for j = 1:J
    if j > 1
        if mesh_uniformity == 1
            marked = 1:meshdata(j-1).ne; % uniform refinement

        else
            [discr_error_K, ~, ~] = evaluate_error_on_elements(x{j-1},x{j-1},j-1);

            [discr_error_K,idx] = sort(discr_error_K,'descend');
            sumeta = cumsum(discr_error_K.^2);
            ell = find(sumeta >= sumeta(end)*mesh_uniformity,1);
            marked = idx(1:ell);
        end

        if mesh_uniformity < 1
            [meshdata(j).coord,meshdata(j).edges,meshdata(j).elements,refined] = ...
                refinemesh_mine_nvb(solution.geometry, meshdata(j-1).coord,meshdata(j-1).edges,meshdata(j-1).elements, marked);
            meshdata(j).nc = length(meshdata(j).coord(1,:));
            meshdata(j).ne = length(meshdata(j).elements(1,:));
        else
            [meshdata(j).coord,meshdata(j).edges,meshdata(j).elements,refined] = ...
                refinemesh_mine(solution.geometry, meshdata(j-1).coord,meshdata(j-1).edges,meshdata(j-1).elements, marked);
            meshdata(j).nc = length(meshdata(j).coord(1,:));
            meshdata(j).ne = length(meshdata(j).elements(1,:));
        end


    end

    % indexing the degrees of freedom
    if j>1
        fill_DOF1(meshdata,m(j-1),j); 
    end

    fill_DOF2(meshdata,m(j),j); 
    if j > 1
        if m(j-1) == m(j)
            p_int = 1;
        else
            p_int = p_interpolation(j); %polynomial interpolation matrix  degree
        end
        DOF2(j).Ip = p_int;
    end

    % data on reference triangle for each level (DOFs change since polynomial degree varies)
    ref_element(j).level = fill_ref_element(DOF2(j));

    if j > 1
        if mesh_uniformity < 1
            [meshdata(j).I,meshdata(j-1).refine] = interpolation_matrix_varp_NVB(j,refined);
            meshdata(j).I1 = interpolation_matrix1_NVB(j,refined); %interpolation matrix for construction of hat functions
        else
            [meshdata(j).I,meshdata(j-1).refine] = interpolation_matrix_varp(j,refined);
            meshdata(j).I1 = interpolation_matrix1(j,refined); %interpolation matrix for construction of hat functions
        end
    end

    if length(A) < J %we have to compute the finest problem for the first time

        % projecting rhs f as piece-wise polynomial and measuring the oscillations
        [tmp, tmp2] = project_f_varp(j);
        solution.fh_coef{j} = tmp;
        osc{j} = tmp2;

        % assembling the stiffness matrix and the right-hand side
        [tempA,tempF,element(j).FEMstiffmatrix, element(j).dJloc] = globalstiffmatrix_varp(j);

        tempx = DOF2(j).dirichletvalues;
        %comment below if the exact algebraic solution is not needed
        tempx(DOF2(j).freenodes) = tempA(DOF2(j).freenodes,DOF2(j).freenodes)\tempF(DOF2(j).freenodes);

        x{j} = tempx; 
        A{j} = tempA(DOF2(j).freenodes,DOF2(j).freenodes);
        F{j} = tempF(DOF2(j).freenodes);

        clear tempx tempA tempF lintempx templinA templinF;
    end
end


%% assembly block-diagonal matrices for solvers (to speed up solution in MATLAB)

[loc_mat, patches, hats] = fill_patch_wRAS_small_varp;

elements_of_patch = struct([]);
elements_of_patch(1).el = [];
for j = 2:J
    elements_of_patch(j).el = find_patch_elements_counterclockwise(j);
end

extended = struct([]);
prolong_matrices = struct([]);
for j = 2:J
    
    prolong_matrices(j).h = meshdata(j).I(DOF1(j).freenodes,DOF2(j-1).freenodes);
    if size(DOF2(j).Ip,2) == 1
        prolong_matrices(j).p = 1;
    else
        prolong_matrices(j).p = DOF2(j).Ip(DOF2(j).freenodes,DOF1(j).freenodes);
    end
    
    temp_1 = patches(j).level;

    extended(j).ext_indices = cat(1,temp_1(:).vertex);
    extended(j).psiaweight_local_ext = cat(1,hats(j).level(:).vertex);

    extended(j).L_chol_ext = blkdiag(loc_mat(j).level(:).chol);
    loc_mat(j).level = [];

    DOF2(j).Ip = [];
    meshdata(j).I = [];

    extended(j).local_sizes = arrayfun(@(x) size(temp_1(x).vertex,1), 1:numel(temp_1));
    temp_1 = [];
end
loc_mat = [];
   
