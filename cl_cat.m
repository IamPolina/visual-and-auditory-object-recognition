function cl_cat(subn, modality)
%% this function outputs category-level classification accuracy 
%% input:
% subn - number of the subject
% modality = 'vis' or 'aud'

Numcats = 3;
use_mvnn = 1;
if isequal(modality,'vis'), triggertemp=100; else, triggertemp=200; end

% the indices of the exemplars, sorted by category divisions: 
%[{{moving}, {nonmoving}}; 
% {{big}, {small}; 
%{{made by nature}, {man-made}}]. 
ctg = [{{1:24}, {25:48}}; {{[1:12, 25:36]},{[13:24,37:48]}}; {{[1:6, 13:18, 25:30, 37:42]}, {[7:12, 19:24, 31:36, 43:48]}}];
%% 0) add paths
addpath('/home/iampolina/OR/WORKFLOW/ANALYSIS')
addpath('/home/iampolina/OR/WORKFLOW/ANALYSIS/rerun_2021/mvnn')
addpath('/home/iampolina/OR/WORKFLOW/ANALYSIS/libsvm-3.11')
direeg = '/scratch/iampolina/OR/DATA/EEG';
%% 1). Load  data load conditions X ...
%% X electodes X probe X time
DT=200;
DA_end = NaN(Numcats, DT);
for KK = 1:Numcats % 3 categories
%% 1A). Load data
if subn<11 || subn>24
    subdireeg = dir(fullfile(direeg, ['sub' num2str(subn, '%02d')], 'timelock_EM_excl.mat'));
else
    if isequal(modality, 'vis')
        subdireeg = dir(fullfile(direeg, ['sub' num2str(subn, '%02d') '_2'], 'timelock_EM_excl.mat'));
    else
        subdireeg = dir(fullfile(direeg, ['sub' num2str(subn, '%02d')], 'timelock_EM_excl.mat'));
    end
end
    
fileName = [subdireeg(1).folder, '/', subdireeg(1).name];
load(fileName)
%% 1B). Sort trials by category division 
for cat = 1:2
    counter=0;
    for cond = ctg{KK, cat}{1}+triggertemp                       
        new_ind = find(timelock.trialinfo == cond);
        counter=counter+1;
        CondMatrix(cat, counter,:, :, :)  = squeeze(nanmean(timelock.trial(new_ind, :, :),1));
    end
end

%% 2). Dimensions:
resdir = fullfile('/scratch/iampolina/OR/DATA/CAT2021'); % where to save the results
mkdir(resdir)
DC = size(CondMatrix,1); % dimension condition
DP = size(CondMatrix,2); %dimension probe
DE = size(CondMatrix,3); % dimension electrode
DT = size(CondMatrix,4); % dimension time
%% 3). Set parameters for supertrials
['create supertrials']
K=3;
L=floor(DP/K);
[H,edges] = histcounts([1:DP],K);
edges = ceil(edges);
edges(1) = edges(1) + 1;
%% 4). Prepare output matrices
num_permutations = 100;
DA = NaN(num_permutations,DT);


for perm =1:num_permutations 
    tic
    ['perm: ' num2str(perm)]
    %% 5). Shuffle trials within each category division
    permutedC = NaN(DC,DP,DE,DT); %% CHANGE
    for co = 1:DC
        permutedC(co, :, :, :) = squeeze(CondMatrix(co, randperm(DP),:,:)); %% CHANGE
    end
    %% 6). Use multivariate noise normalization if set
    if use_mvnn
        permutedC = MVNN(permutedC, K);
    end
    %% 7). Sort trials into super-trials
    for co = 1:DC
        for step= 1:K %average by steps
                pseudo_trialD(co,step,:,:)= squeeze(mean(permutedC(co,edges(step):edges(step)+H(step)-1,:,:),2)); %assign pseudo trial to pseudo_trial_D
        end
    end
    
    %% 8). Classification
    for condA=1:DC%M %loop through all conditions %% CONDITIONS
        for condB = condA+1:DC%M  %loop through all conditions >condA+1 (to get all pair-wise combinations
            %['condA:' num2str(condA) ' ' 'condB: ' num2str(condB)]
            for time_point =1:DT % all time points are independent

                    for testpt = 1:K
                        trainpt = setdiff(1:K, testpt);
                        training_data=[squeeze(pseudo_trialD(condA,trainpt,:, time_point)) ; squeeze(pseudo_trialD(condB,trainpt,:,time_point))];
                        labels_train=[ones(1,K-1) 2*ones(1,K-1)];
                        model = libsvmtrain(labels_train', training_data,'-s 0 -t 0 -q');
                        testing_data=[squeeze(pseudo_trialD(condA,testpt,:,time_point))'; squeeze(pseudo_trialD(condB,testpt,:,time_point))'];
                        labels_test= [1 2];
                        [predicted_label, accuracy, decision_values] = libsvmpredict(labels_test', testing_data , model);
                        DAtemp(testpt) = accuracy(1);
                    end
                    DA(perm, time_point) = mean(DAtemp);
                    clear DAtemp
                    
             end % time point
         end % cond A
     end %cond B
toc
end % permutation
DA_end(KK, :)=squeeze(nanmean(DA,[1]));
end
save([resdir, '/', 'cl_cat_' modality '_sub' num2str(subn, '%02d')], 'DA_end');