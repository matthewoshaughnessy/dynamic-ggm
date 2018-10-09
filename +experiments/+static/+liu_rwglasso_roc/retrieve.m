% --- retrieve files ---
mkdir(fullfile('+experiments/+static/+rwglasso_roc','results'));
system(['scp -r moshaughnessy6@login-d.pace.gatech.edu:' ...
  '~/data/dggm_static_rwglasso_roc/results* +experiments/+static/+rwglasso_roc/results/'], '-echo');


% --- concatenate results ---
resultfiles = dir('+experiments/+static/+rwglasso_roc/results/results*.mat');
resultfiles(strcmp('results.mat',resultfiles)) = [];
if isempty(resultfiles), error('No results returned!'); end
resultfiles = {resultfiles.name};
[~, reindex] = sort(str2double(regexp(resultfiles,'\d+','match','once')));
resultfiles = resultfiles(reindex);
resultfiles(strcmp(resultfiles,'results.mat')) = [];
fprintf('\nConcatenating results: 1/%d...\n', length(resultfiles));
load(['+experiments/+static/+rwglasso_roc/results/' resultfiles{1}]);
results = results1;
for i = 2:length(resultfiles)
  fprintf('Concatenating results: %d/%d.', i, length(resultfiles));
  load(['+experiments/+static/+rwglasso_roc/results/' resultfiles{i}]); fprintf('.');
  assert(isequal(size(results),size(results1)));
  eval(sprintf('ind = ~cellfun(@isempty,results%d);',i)); fprintf('.');
  eval(sprintf('results(ind) = results%d(ind);',i)); fprintf('\n');
end
fprintf('Saving concatenated results...');
save(fullfile('+experiments/+static/+rwglasso_roc','results.mat'),'results','-v7.3');
fprintf('done!\n');


% --- clean up ---
fprintf('Cleaning up...');
system(sprintf('rm -r %s', fullfile('+experiments/+static/+rwglasso_roc','results')));
clear;
fprintf('done!\n');
