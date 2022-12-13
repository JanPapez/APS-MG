function fill_meshdata(Hmax, j)
% function for initializing the mesh data using Delaunay triangulation
%
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG


global solution meshdata;

[meshdata(j).coord,meshdata(j).edges,meshdata(j).elements] = ...
    initmesh(solution.geometry, 'Hmax', Hmax);

[meshdata(j).coord,meshdata(j).edges,meshdata(j).elements] = ...
    permute_mesh(meshdata(j).coord,meshdata(j).edges,meshdata(j).elements);

meshdata(j).nc = length(meshdata(j).coord(1,:));
meshdata(j).ne = length(meshdata(j).elements(1,:));
meshdata(j).diamOmega = solution.diamOmega;
meshdata(j).CF = 1/(2*pi);

end

function [p,e,t] = permute_mesh(p,e,t)

nt=size(t,2);

% Find longest side of each triangle
ls=3*ones(1,nt);
d1=(p(1,t(1,:))-p(1,t(2,:))).^2+(p(2,t(1,:))-p(2,t(2,:))).^2;
d=(p(1,t(2,:))-p(1,t(3,:))).^2+(p(2,t(2,:))-p(2,t(3,:))).^2;
ii=find(d>d1);
ls(ii)=1*ones(size(ii));
d1=max(d,d1);
d=(p(1,t(3,:))-p(1,t(1,:))).^2+(p(2,t(3,:))-p(2,t(1,:))).^2;
ii=find(d>d1);
ls(ii)=2*ones(size(ii));

% Permute so longest side is 3
ii=find(ls==1);
d=t(1,ii);
t(1,ii)=t(2,ii);
t(2,ii)=t(3,ii);
t(3,ii)=d;
ii=find(ls==2);
d=t(1,ii);
t(1,ii)=t(3,ii);
t(3,ii)=t(2,ii);
t(2,ii)=d;

end