load +experiments/+static/+rwglasso_roc/results.mat
params = struct( ...
  'lambda_m', logspace(-5,1,100), ...
	'lambda_d', logspace(-2,1,10), ...
	'lambda_a', 1e-4, ...
	'pinnov',   [0 0.02 0.04 0.06 0.08], ...
  'randseed', 1:12);


%%
[nlm,nld,nla,ni,ntrials] = size(results);
nrw = length(results{1}.err);
fpr = zeros(nrw,nlm,nld,nla,ni,ntrials);
tpr = zeros(nrw,nlm,nld,nla,ni,ntrials);
mcc = zeros(nrw,nlm,nld,nla,ni,ntrials);
for it = 1:ntrials
  for ilm = 1:nlm
    for ild = 1:nld
      for ila = 1:nla
        for ii = 1:ni
          for irw = 1:nrw
            r = results{ilm,ild,ila,ii,it};
            fpr(irw,ilm,ild,ila,ii,it) = r.err(irw).fpr;
            tpr(irw,ilm,ild,ila,ii,it) = r.err(irw).tpr;
            mcc(irw,ilm,ild,ila,ii,it) = r.err(irw).mcc;
          end
        end
      end
    end
  end
end
clearvars ii ila ild ilm irw it


%% plot ROC curve
cols = cbrewer('seq','Blues',nrw+1);
figure(1); clf;
ild = 1;
ila = 1;
ii = 1;
f = @(x) mean(x,6);
for i = 1:nrw
  plot(f(fpr(i,:,ild,ila,ii,:)),f(tpr(i,:,ild,ila,ii,:)), ...
    'x-','color',cols(i+1,:),'linew',2); hold on;
end
%plot(err_ml.fpr,err_ml.tpr,'o','color',parulacols(4));
grid on; xlabel('False positive rate'); ylabel('True positive rate');
legend('Graphical lasso', ...
  'Weighted graphical lasso (1)', ...
  'Weighted graphical lasso (2)', ...
  'location', 'se');
set(gca,'fontsize',24); grid on;
%export_fig -transparent +experiments/+static/+rwglasso_roc/plot1.pdf


%% plot MCC
cols = cbrewer('seq','Blues',nrw+1);
figure(1); clf;
ild = 1;
ila = 1;
ii = 1;
f = @(x) mean(x,6,'omitnan');
for i = 1:nrw
  plot(f(mcc(i,:,ild,ila,ii,:)), ...
    'x-','color',cols(i+1,:),'linew',2); hold on;
end
grid on; xlabel('$$\lambda_m$$','interpreter','latex'); ylabel('MCC');
legend('Graphical lasso', ...
  'Weighted graphical lasso (1)', ...
  'Weighted graphical lasso (2)', ...
  'location', 'se');
set(gca,'fontsize',24); grid on;


%%
figure(2); clf;
itrial = 1;
ialpha_disp = round(linspace(1,nalpha,8));
for i = 1:length(ialpha_disp)
  for j = 1:nrw+1
    subplot(4,length(ialpha_disp),(j-1)*length(ialpha_disp)+i)
    if j <= nrw
      imagesc(abs(Theta_hist(:,:,j,ialpha_disp(i),itrial)));
      xlabel({sprintf('rmse = %.3f', rmse(j,ialpha_disp(i))), ...
        sprintf('gcost = %.4f', gcost(j,ialpha_disp(i)))});
    else
      imagesc(abs(Theta));
    end
    if j == 1
      title(sprintf('\\alpha = %.1e',alphas(ialpha_disp(i))));
    end
    if i == 1
      if j <= nrw, ylabel(sprintf('Iteration %d',j));
      else, ylabel('Ground truth'); end
    end
    axis image; set(gca,'xtick',[],'ytick',[],'fontsize',14); drawnow;
  end
end
export_fig -transparent +experiments/+static/+rwglasso_roc/plot2.pdf


%%
figure(3); clf;
cols = cbrewer('seq','Blues',nalpha);
subplot(121);
for i = 1:10:nalpha
  plot(scost(:,i),'o-','color',cols(end-i+1,:)); hold on; grid on;
end
title('Surrogate cost'); xlabel('Iteration'); set(gca,'fontsize',24);
subplot(122);
for i = 1:10:nalpha
  plot(gcost(:,i),'o-','color',cols(end-i+1,:)); hold on; grid on;
end
title('Global cost'); xlabel('Iteration'); set(gca,'fontsize',24);
export_fig -transparent +experiments/+static/+rwglasso_roc/plot3.pdf

