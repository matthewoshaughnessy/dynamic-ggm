function plotROC(err)
% INPUT
%  err : struct array (generated by +util/evaluateGraph) containing fields
%        tp, tn, fp, and fn

  tpr = [err.tp] ./ ([err.tp] + [err.fn]);
  fpr = [err.fp] ./ ([err.fp] + [err.tn]);
  
  plot(fpr / 

end