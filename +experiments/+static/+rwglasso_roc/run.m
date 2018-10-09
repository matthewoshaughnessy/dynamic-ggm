% Create ROC curves to test LSM-derived reweighted graphical model

params = struct( ...
  'lambda_m', logspace(-5,1,10), ...
	'lambda_d', logspace(-2,1,10), ...
	'lambda_a', logspace(-3,0,4), ...
	'pinnov',   [0 0.02 0.04 0.06 0.08], ...
  'randseed', 1:12);

varsToStore = {'err'};

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

