function [out] = nonlagrangeannodes(p)
% function for setting a nonuniform grid for lagrangean degrees of freedom
% for p-order FEM
%
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG

alpha_all = [1 1 1.4152 0.1001 0.2751 0.9808 1.0999 1.2832 1.3648 1.4773 1.4959 1.5743 1.5770 1.6223 1.6258];

if p > 15
    alpha = 1;
else
    alpha = alpha_all(p);
end


% (nodes.m)
% input: p=polynomial order of interpolant
% input: alpha=tunable warping parameter
% output: X,Y vectors of node coordinates in equilateral triangle function [X,Y] = nodes(p, alpha)
% total number of nodes
N = (p+1)*(p+2)/2;
% 1) compute Gauss-Lobatto-Legendre node distribution
gaussX = gll(p);
% 2) create equidistributed nodes on equilateral triangle
L1 = zeros(N,1); 
L2 = zeros(N,1); 
L3 = zeros(N,1); 
sk = 1;

for n=1:p+1
    for m=1:p+2-n
        L1(sk) = (n-1)/p;
        L3(sk) = (m-1)/p;
        sk = sk+1;
    end
end
L2 = 1-L1-L3;
X = -L2 + L3;
Y = (-L2-L3+2*L1)/sqrt(3);
% 3) compute blending function at each node for each edge
blend1 = 4*L2.*L3;
blend2 = 4*L1.*L3;
blend3 = 4*L1.*L2;
% 4) amount of warp for each node, for each edge
warpfactor1 = warpfactor(p, gaussX, L3-L2); warpfactor2 = warpfactor(p, gaussX, L1-L3); warpfactor3 = warpfactor(p, gaussX, L2-L1);
% 5) combine blend & warp
warp1 = blend1.*warpfactor1.*(1 + (alpha*L1).^2); warp2 = blend2.*warpfactor2.*(1 + (alpha*L2).^2); warp3 = blend3.*warpfactor3.*(1 + (alpha*L3).^2);
% 6) accumulate deformations associated with each edge
X = X + 1*warp1 + cos(2*pi/3)*warp2 + cos(4*pi/3)*warp3; Y = Y + 0*warp1 + sin(2*pi/3)*warp2 + sin(4*pi/3)*warp3;

% transformation to the triangle with vertices [0,0], [1,0], [0,1]
coordinates = [X(1), X(p+1), X(end); Y(1), Y(p+1), Y(end)];
Jloc = [coordinates(:,2)-coordinates(:,1), coordinates(:,3)-coordinates(:,1)];
new = Jloc\[X,Y]';

new = new-new(1,1);

X = new(1,:)';
Y = new(2,:)';

out = [Y,X];
return



function warp = warpfactor(p, xnodes, xout) 

warp = zeros(size(xout));
xeq  = linspace(1,-1,p+1)';
for i=1:p+1
    d = (xnodes(i)-xeq(i));
    for j=2:p
        if(i~=j)
            d = d.*(xout-xeq(j))/(xeq(i)-xeq(j));
        end
    end
    % deflate end roots
    if(i~=1)
        d = -d/(xeq(i)-xeq(1));
    end
    if(i~=(p+1))
        d = d/(xeq(i)-xeq(p+1));
    end
    warp = warp+d;
end
return



function x=gll(p)
% polynomial order  + 1
N=p+1;
% Use the Chebyshev-Gauss-Lobatto nodes as the first guess
x=cos(pi*(0:p)/p)';
% The Legendre Vandermonde Matrix
V=zeros(N,N);
xold=2;
% Newton-Raphson iteration
while max(abs(x-xold))>eps
    xold=x;
    V(:,1)=1;
    V(:,2)=x;
    for k=2:p
        V(:,k+1)=( (2*k-1)*x.*V(:,k)-(k-1)*V(:,k-1) )/k;
    end
    x=xold-( x.*V(:,N)-V(:,p) )./( N*V(:,N) );
end
return