load +experiments/+static/+rwglasso_roc/results.mat


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

