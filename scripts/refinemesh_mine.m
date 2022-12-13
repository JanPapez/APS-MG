function [p1,e1,t1,refined]=refinemesh_mine(g,p,e,t,u,it,mode)
%REFINEMESH Refine a triangular mesh.
%
%       [P1,E1,T1]=REFINEMESH(G,P,E,T) returns a refined version
%       of the triangular mesh specified by the geometry G, point matrix P,
%       edge matrix E, and triangle matrix T.
%
%       G describes the geometry of the PDE problem. G can
%       either be a Decomposed Geometry Matrix or the name of Geometry
%       MATLAB-file. See either DECSG or PDEGEOM for details.
%
%       The triangular mesh is given by the mesh data P, E, and T.
%       Details can be found under INITMESH.
%
%       [P1,E1,T1,U1]=REFINEMESH(G,P,E,T,U) refines the mesh and also
%       extends the function U to the new mesh by linear interpolation.
%       The number of rows in U should correspond to the number of
%       columns in P, and U1 will have as many rows as there are
%       points in P1. Each column of U is interpolated separately.
%
%       An extra input argument is interpreted as a list of
%       subdomains to refine, if it is a row vector, or a list of
%       triangles to refine, if it is a column vector.
%
%       The default refinement method is regular refinement, where
%       all of the specified triangles are divided into four triangles
%       of the same shape. Longest edge refinement, where
%       the longest edge of each specified triangle is bisected, can
%       be demanded by giving 'longest' as a final parameter. Using
%       'regular' as the final parameter results in regular refinement.
%       Some triangles outside of the specified set may also be refined,
%       in order to preserve the triangulation and its quality.
%
%       See also INITMESH, PDEGEOM

%       A. Nordmark 4-26-94, AN 6-21-94.
%       Copyright 1994-2010 The MathWorks, Inc.

% modified by Jan Papez
%   interpolation of a function is removed
%   permutation of vertices is removed
%   new ordering of refined triangles is considered
%
% Jan Papez, Ani Miraci, December 2022
%       APS-MG MATLAB package https://github.com/JanPapez/APS-MG

refined = cell(1); % the indices of refined elements

np=size(p,2);
ne=size(e,2);
nt=size(t,2);

if nargout==4
  intp=0; %%%
else
  intp=0;
end

if nargin-intp==4
  it=(1:nt)';                           % All triangles
  mode='regular';
end

if (~intp) && nargin==5
  it=u;
  if ischar(it)
    mode=it;
    it=(1:nt)';                         % All triangles
  else
    mode='regular';
  end
end

if (~intp) && nargin==6
  mode=it;
  it=u;
end

if intp && nargin==6
  if ischar(it)
    mode=it;
    it=(1:nt)';                         % All triangles
  else
    mode='regular';
  end
end

if strcmp(mode,'regular')==0 && strcmp(mode,'longest')==0
  error(message('pde:refinemesh:InvalidRefineMode'));
end

if size(it,1)>1                        % Triangles
  it=it';
else                                    % Subdomains
  it=pdesdt(t,it);
end

% Cannot use matrix indices that exceeds the size of a signed int
[comp,maxsize]=computer;
indexproblem=np^2>maxsize;

itt1=ones(1,nt);
itt1(it)=zeros(size(it));
it1=find(itt1);                         % Triangles not yet to be refined
it=find(itt1==0);                       % Triangles whos longest side is to be bisected

% Make a connectivity matrix, with edges to be refined.
% -1 means no point is yet allocated
ip1=t(1,it);
ip2=t(2,it);
if strcmp(mode,'regular')
  ip3=t(3,it);
end
A=sparse(ip1,ip2,-1,np,np);
if strcmp(mode,'regular')
  A=A+sparse(ip2,ip3,-1,np,np);
  A=A+sparse(ip3,ip1,-1,np,np);
end
A=-((A+A.')<0);
newpoints=1;

% loop until no additional hanging nodes are introduced
while newpoints
  newpoints=0;
  n=length(it1);
  ip1=t(1,it1);
  ip2=t(2,it1);
  ip3=t(3,it1);
  m1=zeros(1,n);
  m2=m1;
  m3=m1;
  for i=1:n
    m3(i)=A(ip1(i),ip2(i));
    m1(i)=A(ip2(i),ip3(i));
    m2(i)=A(ip3(i),ip1(i));
  end
  ii=find(m3);
  if ~isempty(ii)
    itt1(it1(ii))=zeros(size(ii));
  end
  ii=find((m1 | m2) & (~m3));
  if ~isempty(ii)
    A=A+sparse(ip1(ii),ip2(ii),-1,np,np);
    A=-((A+A.')<0);
    newpoints=1;
    itt1(it1(ii))=zeros(size(ii));
  end
  it1=find(itt1);                       % Triangles not yet fully refined
  it=find(itt1==0);                     % Triangles fully refined
end

% Find edges to be refined
if ~indexproblem
  ie=full(A(e(1,:)+(e(2,:)-1)*np))==-1;
else
  ie=l_extract(A,e(1,:),e(2,:))==-1;
end

ie1=find(ie==0);                        % Edges not to be refined
ie=find(ie);                            % Edges to be refined

% Get the edge "midpoint" coordinates
x = (p(1,e(1,ie))+p(1,e(2,ie)))/2;
y = (p(2,e(1,ie))+p(2,e(2,ie)))/2;

% Create new points
p1=[p [x;y]];
ip=(np+1):(np+length(ie));
np1=np+length(ie);
% Create new edges
e1=[e(:,ie1) ...
        [e(1,ie);ip;e(3,ie);(e(3,ie)+e(4,ie))/2;e(5:7,ie)] ...
        [ip;e(2,ie);(e(3,ie)+e(4,ie))/2;e(4,ie);e(5:7,ie)]];
% Fill in the new points
if ~indexproblem
  A(e(1,ie)+np*(e(2,ie)-1))=ip;
  A(e(2,ie)+np*(e(1,ie)-1))=ip;
else
  A=l_assign(A,[e(1,ie) e(2,ie)],[e(2,ie) e(1,ie)],[ip ip]);
end

% Generate points on interior edges
[i1,i2]=find(A==-1 & A.'==-1);
i=find(i2>i1);
i1=i1(i);
i2=i2(i);
p1=[p1 ((p(1:2,i1)+p(1:2,i2))/2)];
ip=(np1+1):(np1+length(i));
np1=np1+length(i);
% Fill in the new points
if ~indexproblem
  A(i1+np*(i2-1))=ip;
  A(i2+np*(i1-1))=ip;
else
  A=l_assign(A,[i1 i2],[i2 i1],[ip ip]);
end

% Lastly form the triangles
ip1=t(1,it);
ip2=t(2,it);
ip3=t(3,it);
if ~indexproblem
  mp1=full(A(ip2+np*(ip3-1)));
  mp2=full(A(ip3+np*(ip1-1)));
  mp3=full(A(ip1+np*(ip2-1)));
else
  mp1=l_extract(A,ip2,ip3);
  mp2=l_extract(A,ip3,ip1);
  mp3=l_extract(A,ip1,ip2);
end

% Find out which sides are refined
bm=1*(mp1>0)+2*(mp2>0);
% The number of new triangles
nt1=length(it1)+length(it)+sum(mp1>0)+sum(mp2>0)+sum(mp3>0);
t1=zeros(4,nt1);
t1(:,1:length(it1))=t(:,it1);           % The unrefined triangles
nnt1=length(it1);
refined{1} = it1; %%%

if isempty(bm)
	i = bm;
else
	i=find(bm==3);                          % All sides are refined
end
t1(:,(nnt1+1):(nnt1+length(i)))=[t(1,it(i));mp3(i);mp2(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[mp3(i);t(2,it(i));mp1(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[mp2(i);mp1(i);t(3,it(i));t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[mp1(i);mp2(i);mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
refined{2} = it(i); %%%

if isempty(bm)
	i = bm;
else
	i=find(bm==2);                          % Sides 2 and 3 are refined
end
t1(:,(nnt1+1):(nnt1+length(i)))=[t(1,it(i));mp3(i);mp2(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[t(2,it(i));t(3,it(i));mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[mp3(i);t(3,it(i));mp2(i);t(4,it(i))];
nnt1=nnt1+length(i);
refined{3} = it(i); %%%

if isempty(bm)
	i = bm;
else
	i=find(bm==1);                          % Sides 3 and 1 are refined
end
t1(:,(nnt1+1):(nnt1+length(i)))=[t(2,it(i));mp1(i);mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[t(3,it(i));t(1,it(i));mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[t(3,it(i));mp3(i);mp1(i);t(4,it(i))];
nnt1=nnt1+length(i);
refined{4} = it(i); %%%

if isempty(bm)
	i = bm;
else
	i=find(bm==0);                          % Side 3 is refined
end
t1(:,(nnt1+1):(nnt1+length(i)))=[t(3,it(i));t(1,it(i));mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[t(2,it(i));t(3,it(i));mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
refined{5} = it(i); %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function k=l_extract(A,i,j)

if numel(i)~=numel(j)
  error(message('pde:refinemesh:ijNumel'))
end

k=zeros(size(i));

for l=1:numel(i)
  k(l)=A(i(l),j(l));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A=l_assign(A,i,j,k)

if numel(i)~=numel(j) || numel(i)~=numel(k) 
  error(message('pde:refinemesh:ijkNumel'))
end

for l=1:numel(i)
  A(i(l),j(l))=k(l);
end
