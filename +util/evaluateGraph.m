function out = evaluateGraph(A, Ahat, type, tol)
% INPUTS
%  A    : ground truth adjacency matrix
%  A    : inferred adjacency matrix
%  type : type of statistic to produce, one of
%          'mcc'  - matthews correlation coefficient
%          'f1'   - f1 score
%          'prec' - precision
%          'rec'  - recall
%          'tp'   - true positives
%          'fp'   - false positives
%          'fp'   - false positives
%          'fn'   - false negatives
%          'tpr'  - true positive rate
%          'fpr'  - false positive rate
%          'all'  - struct of all statistics
%  tol  : tolerance for existence of edge (if A, Ahat are not 0/1)

  if nargin == 3
    tol = 0;
  end

  A_thresh    = triu(abs(A)    > tol, 1);
  Ahat_thresh = triu(abs(Ahat) > tol, 1);
  
  triu_mask = triu(true(size(A)),1);
  a    = A_thresh(triu_mask);
  ahat = Ahat_thresh(triu_mask);
  
  out.tp = sum( ahat &  a);
  out.tn = sum(~ahat & ~a);
  out.fp = sum( ahat & ~a);
  out.fn = sum(~ahat &  a);
  
  out.tpr = out.tp / (out.tp + out.fn);
  out.fpr = out.fp / (out.fp + out.tn);
  
  out.prec = out.tp / (out.tp + out.fp);
  out.rec  = out.tp / (out.tp + out.fn);
  out.f1 = 2 / (1/out.rec + 1/out.prec);
  out.mcc = (out.tp*out.tn - out.fp*out.fn) / ...
    sqrt((out.tp+out.fp)*(out.tp+out.fn)*(out.tn+out.fp)*(out.tn+out.fn));
  
  if ~strcmpi(type,'all')
    out = out.(type);
  end
  
end