function plot_obj_timeseries(modvar)
% modvar is 1 for vis and 2 for aud
% the signed rank test as implemented in Pantazis et al., 2004 
% "Statistical Surface-Based Morphometry Using a Non-Parametric Approach."
% function for plotting boundedline: https://www.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m

% paths
addpath(genpath('/home/iampolina/OR/WORKFLOW/ANALYSIS/rerun_2021/PLOT/plot_boundedline'))
addpath('/home/iampolina/OR/WORKFLOW/ANALYSIS/rerun_2021/stats_pantazis')
filedir = ['/scratch/iampolina/OR/DATA/OBJ2021/3pseudotrials_mvnn/']
files{1} = dir(fullfile(filedir, 'cl_objvis*.mat'));
files{2} = dir(fullfile(filedir, 'cl_obj_aud*.mat'));
names = [{'vis'}, {'aud'}];
% prepare
for mod = modvar
    for v = 1:numel(files{mod})
        load(fullfile(files{mod}(v).folder, files{mod}(v).name))
        DA_end(isnan(DA_end)) = 0;
        objdata{mod}(v,:) = squeeze(DA_end);
        clear DA_end
    end
end
for mod = modvar
    meancurve{mod} = nanmean(objdata{mod}, 1); 
    semcurve{mod} = std(objdata{mod})/(sqrt(size(objdata{mod},1))-1);
    [SV, clusters,clustersize,StatMapPermPV(mod,:,:)] = permutation_cluster_1sample_alld(squeeze(objdata{mod}(:,:))-50, 10000, 0.05, 0.05, 'right');

end

SV(SV==0) = NaN;
CCC = summer;
% plot
[l,p] = boundedline(1:size(meancurve{modvar},2), meancurve{modvar}, semcurve{modvar}, 'r', 'transparency', 0.25)
outlinebounds(l,p);
hold on
xticks([0:20:200])
hold on
xticklabels({'-200','-100','0', '100','200', '300', '400', '500', '600', '700', '800'})
hold on
plot(1:200, SV+59, ['r*'])
hold on
xlim([0 201])
hline(50, 'k', ' ');
print(fullfile(pwd, ['obj_' names{modvar} '_boundedline']), '-dsvg')
end