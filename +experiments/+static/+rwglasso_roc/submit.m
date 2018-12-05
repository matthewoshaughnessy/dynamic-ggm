sweeputil.pacesweep('dggm_static_rwglasso_roc', ...
  '+experiments/+static/+rwglasso_roc/run.m', ...
  'excludedirs',   {'+datasets', '+experiments', '+visualize'}, ...
  'queue',         'rozell', ...
  'walltime',      '5:00:00:00', ...
  'matlabversion', 'r2016b', ...
  'username',      'moshaughnessy6', ...
  'headnode',      'login-d', ...
  'email',         'moshaughnessy6@gatech.edu', ...
  'mem',           '100gb', ...
  'nodes',         1, ...
  'ppn',           28);