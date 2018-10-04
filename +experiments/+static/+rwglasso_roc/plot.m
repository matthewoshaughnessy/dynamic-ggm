load +experiments/+static/+rwglasso_roc/results.mat
[nalpha,ntrials] = size(results);
p = length(results{1}.A);
nrw = length(results{1}.err_rw);
fpr = zeros(nrw,nalpha,ntrials);
tpr = zeros(nrw,nalpha,ntrials);
cost = zeros(nrw,nalpha,ntrials);
rmse = zeros(nrw,nalpha,ntrials);
Theta_rw_hist = zeros(p,p,nrw,nalpha,ntrials);
Theta = zeros(p,p,ntrials);
alphas = zeros(1,nalpha);
for it = 1:ntrials
  Theta(:,:,it) = results{1,it}.Theta;
  for ia = 1:nalpha
    alphas(ia) = results{ia,it}.inputs{1};
    for ir = 1:nrw
      fpr(ir,ia,it) = results{ia,it}.err_rw(ir).fpr;
      tpr(ir,ia,it) = results{ia,it}.err_rw(ir).tpr;
      cost(ir,ia,it) = results{ia,it}.cost_rw(ir);
      rmse(ir,ia,it) = results{ia,it}.rmse_rw(ir);
      Theta_rw_hist(:,:,ir,ia,it) = results{ia,it}.Theta_rw_hist(:,:,ir);
    end
  end
end


%%
cols = cbrewer('seq','Blues',nrw+1);
figure(1); clf;
f = @(x) mean(x,3);
for i = 1:nrw
  plot(f(fpr(i,:)),f(tpr(i,:)),'o-'); hold on;
end
%plot(err_ml.fpr,err_ml.tpr,'o','color',parulacols(4));
grid on; xlabel('False positive rate'); ylabel('True positive rate');
legend('Graphical lasso', ...
  'Weighted graphical lasso (1)', ...
  'Weighted graphical lasso (2)', ...
  'location', 'se');
set(gca,'fontsize',24); grid on;


%%
clf;
itrial = 1;
ialpha_disp = round(linspace(1,nalpha,8));
for i = 1:length(ialpha_disp)
  for j = 1:nrw+1
    subplot(4,length(ialpha_disp),(j-1)*length(ialpha_disp)+i)
    if j <= nrw
      imagesc(abs(Theta_rw_hist(:,:,j,ialpha_disp(i),itrial)));
      xlabel({sprintf('rmse = %.3f', rmse(j,ialpha_disp(i))), ...
        sprintf('cost = %.4f', cost(j,ialpha_disp(i)))});
    else
      imagesc(abs(Theta));
    end
    if j == 1
      title(sprintf('\\alpha = %.1e',alphas(ialpha_disp(i))));
    end
    axis image; set(gca,'xtick',[],'ytick',[],'fontsize',14); drawnow;
  end
end

