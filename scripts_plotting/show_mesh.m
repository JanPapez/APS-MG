function  show_mesh( meshdata, J )
% function for plotting the discretization mesh
%
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG

pdemesh(meshdata(J).coord,meshdata(J).edges,meshdata(J).elements,...
    'Elementlabels','off','Nodelabels','off')
axis square;
axis off;

end
