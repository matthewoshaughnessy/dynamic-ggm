% INPUTS:
%  - pinnov
%  - lambda0
%  - xi
%  - eta
%  - randseed

% pinnov = 0.5;
% lambda0 = 0.1;
% xi = 1;
% eta = 0.1;
% randseed = 0;

% parameters
p = 15;
n = 65;
nrw = 3;
k = 1;

% generate data
rng(randseed);
[A,D,X,Theta] = util.generateGraph(p,n);
S = 1/n*X*X';
Theta_est = util.addNoise(Theta,pinnov);

%% recover with reweighted l1-penalized maximum likelihood

% initialization
ndiagmask = ones(p)-eye(p);
Theta_prev = eye(p);
Theta_hist = zeros(p,p,nrw);

% --- solve reweighted graphical lasso ---
for irw = 1:nrw
  
  fprintf('Iteration %d/%d...', irw, nrw);
  
  % update weights
  Lambda_rw = xi*(k+1) ./ (lambda0*xi*abs(Theta_prev) + abs(Theta_est) + eta);
  
  % update precision matrix estimate
  diagLambdarwThetarw = diag(diag(Lambda_rw.*abs(Theta_prev)));
  cvx_begin quiet
    variable Theta_rw(p,p) symmetric
    minimize ( - log_det(Theta_rw) + trace(S*Theta_rw) + lambda0*...
      sum( Lambda_rw(:).*abs(Theta_rw(:)).*ndiagmask(:) ) )
    subject to
      Theta_rw - 1e-12*eye(p) == semidefinite(p)
  cvx_end
  
  % save data from this iteration
  Theta_hist(:,:,irw) = Theta_rw;
  Theta_prev = Theta_rw;
  err(irw)  = util.evaluateGraph(Theta, Theta_rw, 'all', 1e3*cvx_slvtol);
  fprintf('done. tpr = %f, fpr = %f.\n', err(irw).tpr, err(irw).fpr);
  
end
