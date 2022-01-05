function cl_cat(subn)
Numcats = 3;
%modality = 'vis' or 'aud'
modality = 'vis';
use_mvnn = 1;
triggertemp=100;
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

if subn<11 || subn>24
    subdireeg = dir(fullfile(direeg, ['sub' num2str(subn, '%02d')], 'timelock_EM_excl.mat'));
else
    subdireeg = dir(fullfile(direeg, ['sub' num2str(subn, '%02d') '_2'], 'timelock_EM_excl.mat'));
end
    
fileName = [subdireeg(1).folder, '/', subdireeg(1).name];
load(fileName)

for cat = 1:2
    counter=0;
    for cond = ctg{KK, cat}{1}+triggertemp                       
        new_ind = find(timelock.trialinfo == cond);
        counter=counter+1;
        CondMatrix(cat, counter,:, :, :)  = squeeze(nanmean(timelock.trial(new_ind, :, :),1));
    end
end

%% 2). Dimensions:
resdir = fullfile('/scratch/iampolina/OR/DATA/CAT2021');
mkdir(resdir)
DC = size(CondMatrix,1); % dimension condition
DP = size(CondMatrix,2); %dimension probe
DE = size(CondMatrix,3); % dimension electrode
DT = size(CondMatrix,4); % dimension time
%% 3). Create pseudotrials
['create pseudotrials']
K=3;
L=floor(DP/K);
[H,edges] = histcounts([1:DP],K);
edges = ceil(edges);
edges(1) = edges(1) + 1;
%% 4). Prepare output matrices
num_permutations = 100;
DA = NaN(num_permutations,DT);
%% 5). Decoding

for perm =1:num_permutations 
    tic
    ['perm: ' num2str(perm)]
    % Shuffle #1
    permutedC = NaN(DC,DP,DE,DT); %% CHANGE
    for co = 1:DC
        permutedC(co, :, :, :) = squeeze(CondMatrix(co, randperm(DP),:,:)); %% CHANGE
    end
    
    if use_mvnn
        permutedC = MVNN(permutedC, K);
    end

%     % shuffle #2
    %pseudo_trialD=NaN(DC,K, DE,DT); %pre-allocate memory
    for co = 1:DC
        for step= 1:K %average by steps
                pseudo_trialD(co,step,:,:)= squeeze(mean(permutedC(co,edges(step):edges(step)+H(step)-1,:,:),2)); %assign pseudo trial to pseudo_trial_D
        end
    end
    %% 6). Use multivariate noise normalization if set
%     if use_mvnn
%       for t = 1:DT
%           for c = 1:DC
%               for p = 1:K
%                 pseudo_trialD(c,p,:,t) = zscore(pseudo_trialD(c,p,:,t), 1,3);
%               end
%           end
%       end
%         pseudo_trialD = MVNN(pseudo_trialD, K);
    %end
    
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
%plot(1:400, squeeze(nanmean(nanmean(DA_end,2),3)))
%title([num2str(subn)])   
end
save([resdir, '/', 'cl_cat_' modality '_sub' num2str(subn, '%02d')], 'DA_end');

%figure(2); plot(1:200, DA_end); %average across conditionsA&B