% run_paper_examples illustrates numerical results obtained in
%
% [1] [Miraci, Papez, Vohralik. SIAM J. Sci. Comput. (2021)]
% https://doi.org/10.1137/20M1349503 
% and
% [2] [Miraci, Papez, Vohralik. Comput. Methods Appl. Math. (2021)] 
% https://doi.org/10.1515/cmam-2020-0024 
%
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG


%Experiment from Table 1 of [1]
problem_number = 5; % L-shaped domain
m = [1,3,3,3]; % polynomial degree per level 
mesh_uniformity = 1; % uniformly-refined meshes
Hmax = 0.145; % initial mesh-size
p_robust_MG01_solver(problem_number, m, mesh_uniformity, Hmax);


%Experiment from Table 2 of [1]
problem_number = 5; % L-shaped domain
m = [1,6,6,6,6,6,6,6,6,6,6]; % polynomial degree per level 
mesh_uniformity = 0.8; % graded meshes
Hmax = inf; % initial mesh-size
p_robust_MG01_solver(problem_number, m, mesh_uniformity, Hmax);


%Experiment from Tables 3 and 4 of [1]
problem_number = 5; % L-shaped domain
m = [1,1,1,1,3]; % polynomial degree per level 
mesh_uniformity = 1; % uniformly-refined meshes
Hmax = 0.145; % initial mesh-size
theta = 0.2; %adaptivity-parameter for the number of post-smoothing steps
p_robust_MG0adapt_solver(problem_number, m, mesh_uniformity, Hmax, theta);


%Experiment from Figure 4 of [2]
problem_number = 5; % L-shaped domain
m = [1,3,3]; % polynomial degree per level 
mesh_uniformity = 1; % uniformly-refined meshes
Hmax = .25; % initial mesh-size
gamma = 0.7;
theta = 0.95; %adaptivity-parameter for the local adaptive smoothing
maxiter = 100; % maximum number of iterations
plot = [1,3]; %plot patches on iteration 3
p_robust_MG01_locallyadaptive_solver(problem_number, m, mesh_uniformity, Hmax, gamma, theta, maxiter, plot);


%Experiment from Figure 5 of [2]
problem_number = 4; % peak problem
m = [1,6,6]; % polynomial degree per level 
mesh_uniformity = 1; % uniformly-refined meshes
Hmax = .145; % initial mesh-size
gamma = 0.7;
theta = 0.95;
maxiter = 100; % maximum number of iterations
plot = [1,4]; %plot patches on iteration 3
p_robust_MG01_locallyadaptive_solver(problem_number, m, mesh_uniformity, Hmax, gamma, theta, maxiter, plot);


%Experiment from Table 1 of [2]
problem_number = 5; % L-shaped domain
m = [1,2,4,6]; % polynomial degree per level 
mesh_uniformity = 1; % uniformly-refined meshes
Hmax = 0.145; % initial mesh-size
gamma = 0.7;
theta = 0.95;
p_robust_MG01_locallyadaptive_solver(problem_number, m, mesh_uniformity, Hmax, gamma, theta);

