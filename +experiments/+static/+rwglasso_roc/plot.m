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
for it = 1:ntrials
  Theta(:,:,it) = results{1,it}.Theta;
  for ia = 1:nalpha
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
clf;
for i = 1:nrw
  plot(mean(fpr(i,:),3),mean(tpr(i,:),3),'o-'); hold on;
end
%plot(err_ml.fpr,err_ml.tpr,'o','color',parulacols(4));
grid on; xlabel('False positive rate'); ylabel('True positive rate');
legend('Graphical lasso', ...
  'Weighted graphical lasso (1)', ...
  'Weighted graphical lasso (2)');
set(gca,'fontsize',18); grid on;


%%
clf;
itrial = 1;
ialpha_disp = 1:5:28;
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
    axis image; set(gca,'xtick',[],'ytick',[],'fontsize',14); drawnow;
  end
end

