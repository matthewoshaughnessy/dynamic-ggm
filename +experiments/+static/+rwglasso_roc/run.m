% Replicate Figure 2a (ROC curves) from Liu and Ihler, 'Learning scale free
% networks by reweighted l1 regularization,' Proc. AISTATS 2011.

params = struct( ...
  'alpha', logspace(-3,0,28), ...
  'randseed', 0);

varsToStore = {'A','S','Theta','Theta_rw_hist','cost_rw','rmse_rw','err_rw'};

scriptname = 'experiments.static.rwglasso_roc.script';
if strfind(getenv('HOSTNAME'),'.pace.gatech.edu')
  resultsJOBNUM = sweep(scriptname, params, ...
    'varsToStore', varsToStore, ...
    'jobNum', JOBNUM, 'totalJobs', TOTALJOBS);
  save('resultsJOBNUM','resultsJOBNUM');
else
  results = sweep(scriptname, params, ...
    'varsToStore', varsToStore);
  save('+experiments/+static/+rwglasso_roc/results','results');
end

