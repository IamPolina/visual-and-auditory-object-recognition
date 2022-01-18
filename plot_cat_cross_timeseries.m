function plot_cat_cross_timeseries
% paths
addpath('/home/iampolina/OR/WORKFLOW/ANALYSIS/rerun_2021/stats_pantazis')
addpath(genpath('/home/iampolina/OR/WORKFLOW/ANALYSIS/rerun_2021/PLOT/plot_boundedline'))
filedir = ['/scratch/iampolina/OR/DATA/CATCROSS2021/attempt01']
files = dir(fullfile(filedir, 'cl_cat_cross*.mat'));
% prepare

for v = 1:numel(files)
    load(fullfile(files(v).folder, files(v).name))
    %[DA_end,Q] = sub_time_after_dec(DA_end,2);
    DADA(v,:,:) = (squeeze(nanmean(DA_end(:,1,:,:),1)) + squeeze(nanmean(DA_end(:,2,:,:),1))')/2;
    clear DA_end
end
    
meancurve = nanmean(DADA, 1); 
stdcurve = std(DADA);

h = pcolor(squeeze(meancurve));%/sqrt(size(DADA,1))));
h.FaceColor = 'interp'
set(h, 'EdgeColor', 'none'); colorbar
colormap(jet)
caxis([50 58])
%title('standard error')
figure(2)
h1 = pcolor(squeeze(stdcurve))
set(h1, 'EdgeColor', 'none'); colorbar
print(fullfile(pwd, 'cat_dec'), '-dsvg')


%1D plot
% cross decoding at the diagonal

[l,p] = boundedline(1:size(meancurve,2), diag(squeeze(meancurve)), ...
         diag(squeeze(stdcurve))/sqrt(size(DADA,1)), 'g', 'transparency', 0.25)
outlinebounds(l,p);
hold on
xticks([0:20:200])
hold on
xticklabels({'-200','-100','0', '100','200', '300', '400', '500', '600', '700', '800'})
hold on

end