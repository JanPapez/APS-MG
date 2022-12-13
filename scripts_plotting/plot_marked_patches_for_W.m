% script for plotting the local distribution of the error and local error
% indicators and showing the patches marked for local adaptive smoothing
% part of 'p_robust_MG0adapt_solver'
%
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG


if iter == plot_algerror_initer
    
    % computing the exact level-wise decomposition and its patch-wise
    % norms
    [rhotilde,local_rhotildenorms] = compute_rhotilde(x_approx_prev,res_prev);
    
    if norm(F{J} - A{J}*(x_approx_prev+rhotilde))/norm(F{J}) > 1e-12
    % approximation to the true error was not computed accurately enough
        disp('Error in computation of rhotilde !')
    end
    
    %PLOT: each subplot has its own scaling
    figure();
    for j = 2:J
        sfh1 = subplot(J-1,2,2*(j-2)+1);
        plot_vertex_values(sqrt(abs(local_rhotildenorms{j})),j);
        str = ['$\| \nabla \tilde{\rho}_',num2str(j-1),'^i \|_{\omega_{',num2str(j-1),'}^{\bf a}}$'];
        tit = title(str);
        set(tit,'interpreter','latex','FontSize', 16)
        c1=colorbar; c1.FontSize = 16;
        sfh1.Position = sfh1.Position + [-0.05 -0.03 0.001 0];
        c1.Position = c1.Position - [0 0 0.02 0];
        
        if problem_number == 5 %remove the corner
            hold on;
            fill([0 0 1 1],[0 -1 -1 0],'white','EdgeColor',[1 1 1]);
        end
        
        Marked = find(lam_local_rhonorms{j} >= threshold_value);
        
        %PLOTTING LOCAL VALUES WITH OPTIMAL STEP SIZES NOW
        
        sfh2 = subplot(J-1,2,2*(j-2)+2);
        plot_vertex_values_marked(sqrt(lam_local_rhonorms{j}),Marked,j);
        str = ['$ (\lambda^i_',num2str(j-1),')^{\frac{1}{2}} \|\nabla \rho_{',num2str(j-1),',{\bf a}}^i \|_{\omega_{',num2str(j-1),'}^{\bf a}}$'];
        tit = title(str);
        set(tit,'interpreter','latex','FontSize',16)
        c2= colorbar; c2.FontSize = 16;
        sfh2.Position = sfh2.Position + [-0.05 -0.03 0.001 0];
        c2.Position = c2.Position - [0 0. 0.02 0];
        
        if problem_number == 5 %remove the corner
            hold on;
            fill([0 0 1 1],[0 -1 -1 0],'white','EdgeColor',[1 1 1]);
        end
    end
    
end