function phi = minimal_angle(coord, elem)
% function for computing the minimal angle of each element
% the minimal angle in the triangulation on level j can be given as
%   phi = minimal_angle(meshdata(j).coord, meshdata(j).elements);
%   min_angle = min(phi);
%
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG

b1=sqrt(sum((coord(:,elem(2,:))-coord(:,elem(1,:))).^2));
b2=sqrt(sum((coord(:,elem(3,:))-coord(:,elem(2,:))).^2));
b3=sqrt(sum((coord(:,elem(1,:))-coord(:,elem(3,:))).^2));

a1=acos((b2.^2+b3.^2-b1.^2)./(2*b2.*b3));
a2=acos((b1.^2+b3.^2-b2.^2)./(2*b1.*b3));
a3=acos((b1.^2+b2.^2-b3.^2)./(2*b1.*b2));
phi = min([a1;a2;a3]);

end