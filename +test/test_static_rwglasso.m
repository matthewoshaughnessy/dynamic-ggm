randseed = 0;
p = 30;
n = 30;

[A,D,X,Theta] = util.generateGraph(p,n,'randseed',randseed);
S = 1/n*X*X';


%% recover with reweighted l1-penalized maximum likelihood
nrw = 3;
alpha_sweep = logspace(-2.5,0,6);

% initialization
Theta_rw_hist = zeros(p,p,nrw,length(alpha_sweep));
cost_rw = zeros(nrw,length(alpha_sweep));
rmse_rw = zeros(nrw,length(alpha_sweep));
theta_ni = @(i) sum(abs(Theta_rw(i,setdiff(1:p,i))));

for ia = 1:length(alpha_sweep)
  
  e = 1;
  alpha = alpha_sweep(ia);
  beta = 2*alpha/e;
  Theta_rw = eye(p);
  
  % --- solve reweighted graphical lasso ---
  for irw = 1:nrw
    fprintf('Trial %d/%d: alpha = %.3e, iteration %d/%d...\n', ...
      ia, length(alpha_sweep), alpha, irw, nrw);
    
    % update weights
    e = diag(Theta_rw);
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
    
    % save data from this iteration
    Theta_rw_hist(:,:,irw,ia) = Theta_rw;
    cost_rw(irw,ia) = cvx_optval;
    rmse_rw(irw,ia) = norm(Theta(:)-Theta_rw(:)) / norm(Theta(:));
    err_rw(irw,ia) = util.evaluateGraph(Theta, Theta_rw, 'all', 100*cvx_slvtol);
    
  end
end


%%
cols = cbrewer('seq','Blues',nrw+1);
fpr = reshape([err_rw.fpr],[nrw length(alpha_sweep)]);
tpr = reshape([err_rw.tpr],[nrw length(alpha_sweep)]);
clf;
for i = 1:3
  plot(fpr(i,:),tpr(i,:),'o-'); hold on;
end
%plot(err_ml.fpr,err_ml.tpr,'o','color',parulacols(4));
grid on; xlabel('False positive rate'); ylabel('True positive rate');
legend('Graphical lasso', ...
  'Weighted graphical lasso (1)', ...
  'Weighted graphical lasso (2)');
set(gca,'fontsize',18); grid on;


%%
for i = 1:length(alpha_sweep)
  for j = 1:nrw+1
    subplot(4,length(alpha_sweep),(j-1)*length(alpha_sweep)+i)
    if j <= nrw
      imagesc(abs(Theta_rw_hist(:,:,j,i)));
      xlabel({sprintf('rmse = %.3f', rmse_rw(j,i)), ...
        sprintf('cost = %.4f', cost_rw(j,i))});
    else
      imagesc(abs(Theta));
    end
    axis image; set(gca,'xtick',[],'ytick',[],'fontsize',14); drawnow;
  end
end


%% recover with maximum likelihood
cvx_begin quiet
  variable Theta_ml(p,p) symmetric
  minimize ( - log_det(Theta_ml) + trace(S*Theta_ml) )
  subject to
    Theta_ml - 1e-12*eye(p) == semidefinite(p)
cvx_end
rmse_ml = norm(Theta(:)-Theta_ml(:)) / norm(Theta(:));
err_ml = util.evaluateGraph(Theta, Theta_ml, 'all', 1e-6);

%% recover with l1-penalized maximum likelihood
lambda_gl_sweep = 10^-0.5; %logspace(-2,0,20);
for i = 1:length(lambda_gl_sweep)
  lambda_gl = lambda_gl_sweep(i);
  fprintf('Trial %d/%d: lambda = %.3e...\n', i, length(lambda_gl_sweep), lambda_gl);
  cvx_begin quiet
  variable Theta_gl(p,p) symmetric
    minimize ( - log_det(Theta_gl) + trace(S*Theta_gl) ...
      + lambda_gl*sum(abs(Theta_gl(:))) )
    subject to
      Theta_gl - eps*eye(p) == semidefinite(p)
  cvx_end
  rmse_gl(i) = norm(Theta(:)-Theta_gl(:)) / norm(Theta(:));
  err_gl(i) = util.evaluateGraph(Theta, Theta_gl, 'all', 1e-6);
end
Theta_gl2 = graphicalLasso(S,lambda_gl_sweep,100000,cvx_slvtol);

%%
clf; plot([err_gl.fpr],[err_gl.tpr],'o-'); hold on;
xlabel('False positive rate'); ylabel('True positive rate');
set(gca,'fontsize',18); grid on;
