function [X,U,Astar] = generateVarTimeSeriesDynamic(O,H,T,snr,Q,randseed)
% generates VAR(1) synthetic time-series data from Foti et al. 2016 fig 3+4
% INPUTS
%   O        : number of observed variables
%   H        : number of hidden variables
%   T        : number of time steps
%   snr      : snr of process
%   Q        : number of time blocks
%   randseed : seed for random number generator
% OUTPUTS
%   X        : observed time-series
%   U        : hidden time-series
%   Astar    : matrix defining VAR(1) process


% --- parse parameters ---
if nargin < 6, randseed = rand; end
if nargin < 5, Q = 5;           end
if nargin < 4, snr = 2;         end
if nargin < 3, T = 1000;        end
if nargin < 2, H = 5;           end
if nargin < 1, O = 50;          end


% --- construct Astar for VAR(1) process ---
rng(randseed);
A = diag(0.2*ones(O,1));
for i = 1:O
  ind = randperm(O,2);
  while any(ind == i)
    ind = randperm(O,2);
  end
  A(i,ind) = round(rand(1,2))-1/2;
end
B = 2*randn(O,H);
for i = 1:H
  ind = randperm(O,round(0.2*O));
  B(ind,i) = 0;
end
C = zeros(H,O);
D = diag(randn(H,1));
Astar0 = [A B; C D];
Astar0 = Astar0 / max(eig(Astar0));


% --- create dynamic Astar ---
Astar = zeros([size(Astar0) Q]);
Astar(:,:,1) = Astar0;
for i = 2:Q
  A = Astar(1:O,1:O,i-1);
  % find edges to remove
  edge = abs(A) > 0;
  edge_ind = find(edge(:));
  edge_ind = setxor(edge_ind,1:O+1:O^2); % don't pick indices on diagonal
  edge_rc = zeros(2,length(edge_ind));
  for k = 1:length(edge_ind)
    [edge_rc(1,k), edge_rc(2,k)] = ind2sub([O O],edge_ind(k));
  end
  % remove edge 1
  k = randperm(length(edge_ind),1);
  A(edge_rc(1,k),edge_rc(2,k)) = 0;
  A(edge_rc(2,k),edge_rc(1,k)) = 0;
  edge_ind(k) = [];
  edge_rc(:,k) = [];
  % remove edge 2
  k = randperm(length(edge_ind),1);
  A(edge_rc(1,k),edge_rc(2,k)) = 0;
  A(edge_rc(2,k),edge_rc(1,k)) = 0;
  edge_ind(k) = [];
  edge_rc(:,k) = [];
  % find edges to add
  edge = abs(A) > 0;
  noedge_ind = find(~edge(:));
  noedge_rc = zeros(2,length(noedge_ind));
  for k = 1:length(noedge_ind)
    [noedge_rc(1,k), noedge_rc(2,k)] = ind2sub([O O],noedge_ind(k));
  end
  % add edge 1
  k = randperm(length(noedge_ind),1);
  v = round(rand)-1/2;
  A(noedge_rc(1,k),noedge_rc(2,k)) = v;
  A(noedge_rc(2,k),noedge_rc(1,k)) = v;
  noedge_ind(k) = [];
  noedge_rc(:,k) = [];
  % add edge 2
  k = randperm(length(noedge_ind),1);
  v = round(rand)-1/2;
  A(noedge_rc(1,k),noedge_rc(2,k)) = v;
  A(noedge_rc(2,k),noedge_rc(1,k)) = v;
  noedge_ind(k) = [];
  noedge_rc(:,k) = [];
  Astar(:,:,i) = Astar0;
  Astar(1:O,1:O,i) = A;
end
for i = 1:Q
  Astar(:,:,i) = Astar(:,:,i) / max(eig(Astar(:,:,i)));
end


% --- generate data from VAR(1) process ---
Astar_spectralNorms = zeros(1,Q);
for i = 1:Q
  Astar_spectralNorms(i) = norm(Astar(:,:,i),2);
end
sigma = mean(Astar_spectralNorms) / snr;
XU = zeros(O+H,T*Q);
XU(:,1) = randn(O+H,1);
for t = 2:T
  e = sigma*randn(O+H,1);
  XU(:,t) = Astar(:,:,1)*XU(:,t-1) + e;
end
for q = 2:Q
  XU(:,T*(q-1)+1) = randn(O+H,1);
  for t = 2:T
    e = sigma*randn(O+H,1);
    XU(:,T*(q-1)+t) = Astar(:,:,q)*XU(:,T*(q-1)+t-1) + e;
  end
end
X = XU(1:O,:);
U = XU(O+1:end,:);


end

