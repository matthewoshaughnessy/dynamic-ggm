function [A,D,X,Theta] = generateGraph(p,n,varargin)
% GENERATEGRAPH
%  Generates a synthetic random graph
% INPUTS
%  p        : number of nodes
%  n        : number of samples 
% PARAMETERS
%  type     : type of graph to generate (default: 'ba')
%              - 'ba'     : scale-free graph using Barabasi-Albert 1999
%  randseed : seed for random number generator (default: random prime)
%  nedge    : number of edges (default: p)
% OUTPUTS
%  A        : p x p adjacency matrix with zero diagonal
%  D        : p x p diagonal degree matrix
%  X        : p x n matrix of observations
%  Theta    : p x p precision matrix


  % === parse inputs ===
  parser = inputParser;
  addRequired( parser, 'p',    @(x)isnumeric(x)&&x==round(x));
  addRequired( parser, 'n',    @(x)isnumeric(x)&&x==round(x));
  addParameter(parser, 'type',     'ba',     @(x)any(strcmpi(x,{'ba'})));
  addParameter(parser, 'randseed', [],       @(x)isnumeric(x));
  addParameter(parser, 'nedge',    [],       @(x)isnumeric(x)&&x==round(x));
  parse(parser,p,n,varargin{:});
  params = parser.Results;
  if isempty(params.nedge), params.nedge = params.p; end
  if ~isempty(params.randseed), rng(params.randseed); end
  
  
  % === generate graph ===
  Ans = zeros(params.p);
  switch lower(params.type)
    
    % --- scale-free graph using barabasi-albert (1999) method ---
    case 'ba'
      % create initial 4-node cycle
      ninit = 4;
      iinit = 1:ninit;
      for i = 1:ninit-1
        Ans(iinit(i),iinit(i+1)) = 1;
      end
      Ans(iinit(end),iinit(1)) = 1;
      % generate each additional node
      for i = ninit+1:params.nedge
        d = sum(Ans+Ans');
        pcs = cumsum(d / sum(d));
        r = rand;
        inew = find(r < pcs, 1);
        Ans(i,inew) = 1;
        A = Ans + Ans';
        D = diag(sum(A));
      end
      
    otherwise
      error('invalid graph type');

  end
  
  % === generate observations ===
  L = 1.1*D - A;
  Lambda = diag(diag(inv(L)));
  Theta = sqrt(Lambda)*L*sqrt(Lambda);
  Sigma = inv(Theta);
  X = mvnrnd(zeros(n,p),Sigma).';
  %S = 1/n*X*X';
  

end