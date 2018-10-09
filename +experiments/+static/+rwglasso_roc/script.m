% INPUTS:
%  - lambda_m
%  - lambda_d
%  - lambda_a
%  - pinnov
%  - randseed

p = 10;
n = 10;

rng(randseed);
[A,D,X,Theta] = util.generateGraph(p,n);
S = 1/n*X*X';

Theta_prev = util.addNoise(Theta,pinnov);
f = @(x) x;


%% recover with reweighted l1-penalized maximum likelihood

% parameters
nrw = 3;
Theta_rw = eye(p);

% initialization
Theta_hist = zeros(p,p,nrw);

% --- solve reweighted graphical lasso ---
for irw = 1:nrw
  
  % update weights
  lambda_rw = lambda_m ./ (abs(Theta_rw) + lambda_d*abs(f(Theta_prev)) + lambda_a);
  lambda_rw = lambda_rw - diag(diag(lambda_rw));
  
  % update precision matrix estimate
  cvx_begin quiet
  variable Theta_rw(p,p) symmetric
  minimize ( - log_det(Theta_rw) + trace(S*Theta_rw) + sum(lambda_rw(:).*abs(Theta_rw(:))) )
  subject to
    Theta_rw - 1e-12*eye(p) == semidefinite(p)
  cvx_end
  
  % save data from this iteration
  Theta_hist(:,:,irw) = Theta_rw;
  err(irw)  = util.evaluateGraph(Theta, Theta_rw, 'all', 1e3*cvx_slvtol);
  
end

