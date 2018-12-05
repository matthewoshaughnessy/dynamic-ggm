function [Theta_noisy, A_noisy] = addNoise(Theta,pchange)
% ADDNOISE
%  Adds noise to the precision matrix of a graph
% INPUTS
%  Theta   : p x p precision matrix
%  pchange : probability of a connection switching
% OUTPUTS
%  Theta_noisy : p x p modified precision matrix
  
  
  p = length(Theta);
  mask_triu = triu(ones(p),1);
  ind_triu = find(ones(p) == mask_triu);
  
  Theta_d = diag(diag(Theta));
  Theta_ud = triu(Theta,1);
  
  ind_z  = find(Theta_ud == 0 & mask_triu);
  ind_nz = find(abs(Theta_ud) > 0 & mask_triu);
  
  nchange = round(pchange*length(ind_nz));
  ind_add = ind_z(randperm(length(ind_z),nchange));
  ind_remove = ind_nz(randperm(length(ind_nz),nchange));
  
  Theta_noisy = Theta_ud;
  Theta_noisy(ind_add) = rand(size(ind_add))*(max(Theta(ind_nz))-min(Theta(ind_nz))) + min(Theta(ind_nz));
  Theta_noisy(ind_remove) = 0;
  
  A_noisy = Theta_noisy ~= 0 + Theta_noisy.' ~= 0;
  Theta_noisy = Theta_noisy + Theta_noisy.' + Theta_d;

end

