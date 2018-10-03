function plotS(S,A)
  
  allS = squeeze(all(abs(S) > 1e-3,1));
  
  subplot(311);
  imagesc(allS);
  title(sprintf('K_Y (||K_Y||_0 = %d)',sum(allS(:))));
  set(gca,'fontsize',16);
  axis image;
  
  subplot(312);
  symA = abs(A)+abs(A)' > 0;
  imagesc(symA);
  title(sprintf('A + A^T (||A+A^T||_0 = %d)',sum(symA(:))));
  set(gca,'fontsize',16);
  axis image;
  
  subplot(313);
  imshowpair(allS,symA);
  axis image;

end