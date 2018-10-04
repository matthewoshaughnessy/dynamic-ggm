% INPUTS:
%  - alpha
%  - randseed

p = 50;
n = 50;

[A,D,X,Theta] = util.generateGraph(p,n,'randseed',randseed);
S = 1/n*X*X';


%% recover with reweighted l1-penalized maximum likelihood

% parameters
nrw = 3;
Theta_rw = eye(p);
e = 0.5*ones(1,p);

% initialization
Theta_hist = zeros(p,p,nrw);
gcost = zeros(nrw,1);
scost = zeros(nrw,1);
rmse = zeros(nrw,1);
theta_ni = @(i) sum(abs(Theta_rw(i,setdiff(1:p,i))));

% --- solve reweighted graphical lasso ---
for irw = 1:nrw
  
  % update weights
  beta = 2*alpha/e(1);
  lambdarw = zeros(p);
  for i = 1:p
    for j = 1:p
      ni = setdiff(1:p,i);
      nj = setdiff(1:p,j);
      lambdarw(i,j) = alpha*(1/(sum(abs(Theta_rw(i,ni)))+e(i)) + ...
        1/(sum(abs(Theta_rw(j,nj)))+e(j)));
    end
  end
  lambdarw = lambdarw - diag(diag(lambdarw));
  
  % update precision matrix estimate
  cvx_begin quiet
  variable Theta_rw(p,p) symmetric
  minimize ( - log_det(Theta_rw) + trace(S*Theta_rw) ...
    + sum(lambdarw(:).*abs(Theta_rw(:))) + beta*sum(abs(diag(Theta_rw))) )
  subject to
  Theta_rw - eps*eye(p) == semidefinite(p)
  cvx_end
  
  cost_surrogate = -log_det(Theta_rw) + trace(S*Theta_rw) + ...
    sum(lambdarw(:).*abs(Theta_rw(:))) + beta*sum(abs(diag(Theta_rw)));
  cost_global = -log_det(Theta_rw) + trace(S*Theta_rw) + ...
    alpha*log(sum(sum(abs(Theta_rw-diag(diag(Theta_rw))))+e)) + ...
    beta*sum(abs(diag(Theta_rw)));
  
  % save data from this iteration
  Theta_hist(:,:,irw) = Theta_rw;
  scost(irw) = cost_surrogate;
  gcost(irw) = cost_global;
  rmse(irw) = norm(Theta(:)-Theta_rw(:)) / norm(Theta(:));
  err(irw)  = util.evaluateGraph(Theta, Theta_rw, 'all', 1e3*cvx_slvtol);
  
end

