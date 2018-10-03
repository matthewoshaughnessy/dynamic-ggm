function videoSL(S, L, videofilename)

  % --- parse inputs ---
  assert(isequal(size(S),size(L)));
  if nargin < 3
    videofilename = [];
  end

  % --- open video file ---
  if ~isempty(videofilename)
    v = VideoWriter('lvsglasso_out','MPEG-4');
    open(v);
  end

  % --- plot ---
  clim = [min(abs(S(:))) max(abs(S(:)))];
  for i = 1:size(S,1)
    Si = squeeze(S(i,:,:));
    Li = squeeze(L(i,:,:));
    subplot(211);
    imagesc(abs(Si),clim); colorbar;
    title({'Recovered sparse component', ...
      sprintf('(frame %d -- ||S_i||_0 = %d)',i,sum(abs(Si(:))>1e-3))});
    set(gca,'fontsize',16); axis image;
    subplot(212);
    imagesc(abs(Li)); colorbar;
    title({'Recovered low-rank component', sprintf('(frame %d -- rank %d)',i,rank(squeeze(L(i,:,:))))});
    set(gca,'fontsize',16); axis image;
    drawnow;
    if ~isempty(videofilename)
      writeVideo(v,getframe(gcf));
    end
  end

  % --- close video ---
  if ~isempty(videofilename)
    close(v);
  end

end
