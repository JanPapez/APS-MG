function [value, grad] = phikl_all(alpha1,alpha2,m)
% function giving the values and the gradient of Lagrange basis function
% (associated with non-lagrangean grid) in point(s) (in bary.coord.)
% [alpha1,alpha2]
%
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG
    
    n = (m+1)*(m+2)/2;  

    % non-lagrangean grid
    coordinates = nonlagrangeannodes(m);
    
    X = [alpha1(:); coordinates(:,1)];
    Y = [alpha2(:); coordinates(:,2)];
    
    [db,dbx,dby] = dubiner2d(2*X-1,2*Y-1,m);
    
    Vandermonde = db(length(alpha1)+1:end,:);
    Transf = Vandermonde \ eye(n);
    
    outval = db(1:length(alpha1),:)*Transf;
    outdx = 2*dbx(1:length(alpha1),:)*Transf;
    outdy = 2*dby(1:length(alpha1),:)*Transf;
    
    value = outval;
    
    grad = zeros(length(alpha1), 2*size(outdx,2));
    grad(:,1:2:end) = outdx;
    grad(:,2:2:end) = outdy;
    grad(isnan(grad)) = 0;
    
end


function [db,dbx,dby] = dubiner2d(x,y,k)
% function [db,dbx,dby]=dubiner2d(x,y,k)
% input:
% x,y: the coordinates at which we evaluate the Dubiner basis
% k: the degree of the polynomial. Must be a column vector.
% output:
% db: the value of dubiner basis at given coordinates, Nnodes x Nbasis
% dbx: the derivative to x of the Dubiner basis, Nnodes x Nbasis
% dby: the derivative to y of the Dubiner basis, Nnodes x Nbasis
%
% Last modified: June 5 2012, code by John Burkardt
Nnodes=size(x,1);
Nbasis=nchoosek(k+2,2);
db=zeros(Nnodes,Nbasis);
dbx=zeros(Nnodes,Nbasis);
dby=zeros(Nnodes,Nbasis);
eta1(y~=1)= 2*(1+x(y~=1))./(1-y(y~=1))-1;
eta1(y==1)= -1;
eta1=eta1';
eta2=y;
% locate the index with p,q
a=zeros(k+1,k+1);
index=1;
for l=2:k+2
    for i=l-1:-1:1
        a(i,l-i)=index;
        index=index+1;
    end
end
loc=@(p,q) a(p+1,q+1);
%
% Dubiner Basis
JP=JacobiP(k,0,0,eta1);
for p=0:k
    A=JP(:,p+1).*((1-eta2)/2).^p;
    JPQ=JacobiP(k-p,2*p+1,0,eta2);
    for q=0:k-p
        db(:,loc(p,q))=A.*JPQ(:,q+1);
    end
end
%
% Derivative of Dubiner Basis
dbx(:,loc(0,0))=0.*x; 
dby(:,loc(0,0))=0.*y;
if k>0
    dbx(:,loc(1,0))=1+0.*x; 
    dby(:,loc(1,0))=1/2+0.*y;
end
for p=1:k-1
    dbx(:,loc(p+1,0))=(2*p+1)/(p+1)*(1+2.*x+y)/2.*dbx(:,loc(p,0))...
        -p/(p+1)*(1-y).^2/4.*dbx(:,loc(p-1,0))...
        +(2*p+1)/(p+1)*db(:,loc(p,0));
    dby(:,loc(p+1,0))=1/2*(2*p+1)/(p+1)*db(:,loc(p,0))...
        +(2*p+1)/(p+1)*(1+2*x+y)/2.*dby(:,loc(p,0))...
        -p/(p+1)*(1/2* (y-1).*db(:,loc(p-1,0))+(1-y).^2/4.*dby(:,loc(p-1,0)) );
end
for p=0:k-1
    dbx(:,loc(p,1))=dbx(:,loc(p,0)).*(1+2*p+(3+2*p)*y)/2;
    dby(:,loc(p,1))=dby(:,loc(p,0)).*(1+2*p+(3+2*p)*y)/2+(3+2*p)/2*db(:,loc(p,0));
end
for p=0:k-1
    for q=1:k-p-1
        a=(2*q+2*p+2)*(2*q+2*p+3)/2/(q+1)/(q+2*p+2);
        b=(2*q+2*p+2)*(2*p+1)^2/2/(q+1)/(q+2*p+2)/(2*q+2*p+1);
        c=(2*q+2*p+3)*(q+2*p+1)*q/(q+1)/(q+2*p+2)/(2*q+2*p+1);
        dbx(:,loc(p,q+1))=(a*y+b).*dbx(:,loc(p,q))-c*dbx(:,loc(p,q-1));
        dby(:,loc(p,q+1))=(a*y+b).*dby(:,loc(p,q))-c*dby(:,loc(p,q-1))+a*db(:,loc(p,q));
    end
end
end


function v = JacobiP(n,alpha,beta,x)
% Subfunction with evaluation of Jacobi polynomials
% taken from code by John Burkardt,
% distributed under the GNU LGPL license
v=zeros(size(x,1),n+1);
v(:,1) = 1.0;
if ( n == 0 )
    return
end
v(:,2) = ( 1.0 + 0.5 * ( alpha + beta ) ) * x(:) + 0.5 * ( alpha - beta );
for i = 2 : n
    c1 = 2 * i * ( i + alpha + beta ) * ( 2 * i - 2 + alpha + beta );
    c2 = ( 2 * i - 1 + alpha + beta ) * ( 2 * i + alpha + beta ) ...
        * ( 2 * i - 2 + alpha + beta );
    c3 = ( 2 * i - 1 + alpha + beta ) * ( alpha + beta ) * ( alpha - beta );
    c4 = - 2 * ( i - 1 + alpha ) * ( i - 1 + beta ) * ( 2 * i + alpha + beta );
    v(:,i+1) = ( ( c3 + c2 * x(:) ) .* v(:,i) + c4 * v(:,i-1) ) / c1;
end
end



