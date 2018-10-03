function plotSL(S,L,A)
% INPUTS
%  S : sparse component (F x O x O)
%  L : low-rank component (F x O x O)
%  A : Observed component of Astar (O x O)
  
  allS = squeeze(all(abs(S) > 1e-3,1));
  clf;
  
  subplot(221);
  imagesc(allS);
  title(sprintf('K_Y (||K_Y||_0 = %d)',sum(allS(:))));
  set(gca,'fontsize',12,'xtick',[],'ytick',[]);
  axis image;
  
  subplot(222);
  spectra = zeros(size(S,2),size(S,1));
  for i = 1:size(L,1)
    spectra(:,i) = svd(squeeze(L(i,:,:)));
  end
  semilogy(spectra,'linew',0.5); grid on; hold on;
  title('Spectra of recovered L');
  set(gca,'fontsize',12,'xtick',0:10:size(S,2));
  xlim([1 size(S,2)]);
  
  subplot(223);
  symA = abs(A)+abs(A)' > 0;
  imagesc(symA);
  title(sprintf('A + A^T (||A+A^T||_0 = %d)',sum(symA(:))));
  set(gca,'fontsize',12,'xtick',[],'ytick',[]);
  axis image;
  
  subplot(224);
  imshowpair(allS,symA);
  title({'purple = FN','green = FP'});
  set(gca,'fontsize',12,'xtick',[],'ytick',[]);
  axis image;

end