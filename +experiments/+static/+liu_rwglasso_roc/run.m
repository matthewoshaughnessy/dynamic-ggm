% Replicate Figure 2a (ROC curves) from Liu and Ihler, 'Learning scale free
% networks by reweighted l1 regularization,' Proc. AISTATS 2011.

params = struct( ...
  'alpha', logspace(-1.9,-0.2,140), ...
  'randseed', 1:14);

%varsToStore = {'A','S','Theta','Theta_hist','scost','gcost','rmse','err'};
varsToStore = {'err'};

scriptname = 'experiments.static.liu_rwglasso_roc.script';
if strfind(getenv('HOSTNAME'),'.pace.gatech.edu')
  resultsJOBNUM = sweep(scriptname, params, ...
    'varsToStore', varsToStore, ...
    'jobNum', JOBNUM, 'totalJobs', TOTALJOBS);
  save('resultsJOBNUM','resultsJOBNUM');
else
  results = sweep(scriptname, params, ...
    'varsToStore', varsToStore);
  save('+experiments/+static/+liu_rwglasso_roc/results','results');
end

