function [X,U,Astar] = generateVarTimeSeries(O,H,T,snr,randseed)
% generates VAR(1) synthetic time-series data from Foti et al. 2016 fig 3+4
% INPUTS
%   O        : number of observed variables
%   H        : number of hidden variables
%   T        : number of time steps
%   snr      : snr of process
%   randseed : seed for random number generator
% OUTPUTS
%   X     : observed time-series
%   U     : hidden time-series
%   Astar : matrix defining VAR(1) process


  % --- parse parameters ---
  if nargin < 5, randseed = rand; end
  if nargin < 4, snr = 2;    end
  if nargin < 3, T = 1000;   end
  if nargin < 2, H = 5;      end
  if nargin < 1, O = 50;     end


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
  Astar = [A B; C D];
  Astar = Astar / svds(Astar,1); % modified 5/23/2018 -- matt


  % --- generate data from VAR(1) process ---
  sigma = norm(Astar,2) / snr;
  XU = zeros(O+H,T);
  XU(:,1) = randn(O+H,1);
  for t = 2:T
    e = sigma*randn(O+H,1);
    XU(:,t) = Astar*XU(:,t-1) + e;
  end
  X = XU(1:O,:);
  U = XU(O+1:end,:);


end

