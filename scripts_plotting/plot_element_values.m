function plot_element_values(values, j)
% function for plotting a vector of values, each associated to a single 
% element
% 
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG

global meshdata J

if nargin < 2
    j = J;
end

h = pdesurf(meshdata(j).coord, meshdata(j).elements, (values(:))');
set(h, 'edgecolor','none')
colormap('default'), colorbar, view(0,90);

set(gca,'XTick',[]);
set(gca,'YTick',[]);






