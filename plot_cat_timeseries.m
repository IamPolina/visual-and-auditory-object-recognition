function plot_cat_timeseries(modvar)
% modvar is 1 for vis and 2 for aud
% the signed rank test as implemented in Pantazis et al., 2004 
% "Statistical Surface-Based Morphometry Using a Non-Parametric Approach."
% function for plotting boundedline: https://www.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m

% paths
addpath('/home/iampolina/OR/WORKFLOW/ANALYSIS/rerun_2021/stats_pantazis')
filedir = ['/scratch/iampolina/OR/DATA/CAT2021/check/']
addpath(genpath('/home/iampolina/OR/WORKFLOW/ANALYSIS/rerun_2021/PLOT/plot_boundedline'))
files{1} = dir(fullfile(filedir, 'cl_cat_vis*.mat'));
files{2} = dir(fullfile(filedir, 'cl_cat_aud*.mat'));
% prepare
for mod = modvar
    for v = 1:numel(files{mod})
        load(fullfile(files{mod}(v).folder, files{mod}(v).name))
        DA_end(isnan(DA_end)) = 0;
        %histogram(nonzeros(DA_end))
        for t = 1:size(DA_end,2)
            catdata{mod}(v, t) = nanmean(squeeze(DA_end(:,t)'));
        end
        clear DA_end
    end
end

for mod = modvar
    meancurve{mod} = nanmean(catdata{mod}, 1); 
    semcurve{mod} = std(catdata{mod})/(sqrt(size(catdata{mod},1))-1);
    [SV, clusters,clustersize,StatMapPermPV(mod,:,:)] = permutation_cluster_1sample_alld(squeeze(catdata{mod}(:,:))-50, 10000, 0.05, 0.05, 'right');
end
SV(SV==0) = NaN;
% plot
CCC = summer;
figure(20)
[l,p] = boundedline(1:size(meancurve{modvar},2), meancurve{modvar}, semcurve{modvar}, 'b', 'transparency', 0.25)
outlinebounds(l,p);
hold on
xticks([0:20:200])
hold on
xticklabels({'-200','-100','0', '100','200', '300', '400', '500', '600', '700', '800'})
hold on
plot(1:200, SV+58, ['b*'])
hold on
xlim([0 201])
hold on
hline(50, 'k', ' ');
print(fullfile(pwd, 'cat_vis_boundedline'), '-dsvg')
end