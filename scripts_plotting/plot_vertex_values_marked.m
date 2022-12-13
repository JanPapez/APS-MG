function plot_vertex_values_marked(values, marked, j, yncolorbar, whichcolormap)
% function for plotting a vector of values, each associated to a single 
% vertex or patch and marking the set of them in 'marked' by red border
% 
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG

global meshdata solution

[v,c] = voronoin([meshdata(j).coord'; -20 -20; -20 20; 20 20; 20 -20]);
for i = 1:meshdata(j).nc
    patch(v(c{i},1),v(c{i},2),values(i));
end

if ~isempty(marked)
    for i = 1:length(marked) %marked ones
        vertex_marked = marked(i);
        S = patch(v(c{vertex_marked},1),v(c{vertex_marked},2),values(vertex_marked));
        S.EdgeColor = 'red';
        S.LineWidth = 0.7;
    end
end 

switch solution.name
    case 'lshape'
        hold on;
        fill([0 0 1.1 1.1 0], [-1.1 0 0 -1.1 -1.1], 'w');
    case 'crack'
        hold on;
        fill([-1.1 0 -1 -1.1], [-1.1 -1 0 -1.1], 'w');
        fill([-1.1 -1 0 -1.1], [1.1 0 1 1.1], 'w');
        fill([1.1 0 1 1.1], [1.1 1 0 1.1], 'w');
        fill([1.1 1 0 1.1], [-1.1 0 -1 -1.1], 'w');
        box off;
end

if nargin > 3 
    if yncolorbar == 1
        axis square;
        colorbar('southoutside');
    else
        colorbar(yncolorbar);
    end
    if nargin == 4
        colormap(whichcolormap);
    end
end

set(gca,'XTick',[]);
set(gca,'YTick',[]);

xlim([min(meshdata(j).coord(1,:)) max(meshdata(j).coord(1,:))])
ylim([min(meshdata(j).coord(2,:)) max(meshdata(j).coord(2,:))])

hold off