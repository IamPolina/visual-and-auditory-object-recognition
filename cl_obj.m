function cl_obj(subn, modality)
%% this function outputs object-level classification accuracy 
%% input:
% subn - number of the subject
% modality = 'vis' or 'aud'
use_nn = 1;
triggertemp=100;
%% 0) add paths
addpath('/home/iampolina/OR/WORKFLOW/ANALYSIS')
addpath('/home/iampolina/OR/WORKFLOW/ANALYSIS/rerun_2021/mvnn')
addpath('/home/iampolina/OR/WORKFLOW/ANALYSIS/libsvm-3.11')
direeg = '/scratch/iampolina/OR/DATA/EEG';
%% 1). Load  data load conditions X ...
%% X electodes X probe X time 
Nobj = 48;
if subn<11 || subn>24 
    subdireeg = dir(fullfile(direeg, ['sub' num2str(subn, '%02d')], 'timelock_EM_excl.mat'));
else
    subdireeg = dir(fullfile(direeg, ['sub' num2str(subn, '%02d') '_2'], 'timelock_EM_excl.mat'));
end
fileName = [subdireeg(1).folder, '/', subdireeg(1).name];
load(fileName)
for cond = (1:Nobj)+triggertemp
    howmany(cond-triggertemp) = sum(timelock.trialinfo==cond);
end
nch = min(howmany);
C = NaN(Nobj,nch,size(timelock.trial,2), size(timelock.trial,3));
for cond = (1:Nobj)+triggertemp
    new_ind = find(timelock.trialinfo == cond);
    CondMatrix(cond-triggertemp,:,:,:) = timelock.trial(new_ind(1:nch), :, :);
end
%% 2). Dimensions:
resdir = fullfile('/scratch/iampolina/OR/DATA/OBJ2021');
mkdir(resdir)
CondMatrix = permute(CondMatrix, [1,3,2,4]);
DC = size(CondMatrix,1); % dimension condition
DE = size(CondMatrix,2); %dimension electrode
DP = size(CondMatrix,3); % dimension probe
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
DA = zeros(num_permutations,DT);
DA_end = zeros(DT);
%% 5). Decoding
for perm =1:num_permutations 
    tic
    ['perm: ' num2str(perm)]
    % Shuffle #1
    permutedC = NaN(DC,DE,DP,DT); %% CHANGE
    for co = 1:DC
        permutedC(co, :, :, :) = squeeze(CondMatrix(co, :, randperm(DP),:)); %% CHANGE
    end
    % shuffle #2
    pseudo_trialD=NaN(DC,DE,K,DT); %pre-allocate memory
    for co = 1:DC
        for step= 1:K %average by steps
                pseudo_trialD(co,:,step,:)= squeeze(mean(permutedC(co,:,edges(step):edges(step)+H(step)-1,:),3)); %assign pseudo trial to pseudo_trial_D
        end
    end
    %% 6). Use noise normalization if set
    if use_nn
      for t = 1:DT
          for c = 1:DC
              for p = 1:K
                pseudo_trialD(c,:,p,t) = zscore(pseudo_trialD(c,:,p,t), 1,2);
              end
          end
      end
    end
    for time_point =1:DT % all time points are independent
            for condA=1:DC%M %loop through all conditions %% CONDITIONS
                for condB = condA:DC%M  %loop through all conditions >condA+1 (to get all pair-wise combinations
            %['condA:' num2str(condA) ' ' 'condB: ' num2str(condB)]
            

                    for testpt = 1:K
                        trainpt = setdiff(1:K, testpt);
                        training_data=[squeeze(pseudo_trialD(condA,:,trainpt, time_point))' ; squeeze(pseudo_trialD(condB,:,trainpt,time_point))'];
                        labels_train=[ones(1,K-1) 2*ones(1,K-1)];
                        model = libsvmtrain(labels_train', training_data,'-s 0 -t 0 -q');
                        testing_data=[squeeze(pseudo_trialD(condA,:,testpt,time_point)); squeeze(pseudo_trialD(condB,:,testpt,time_point))];
                        labels_test= [1 2];
                        [predicted_label, accuracy, decision_values] = libsvmpredict(labels_test', testing_data , model);
                        DAtemp(testpt) = accuracy(1);
                        clear training_data labels_train model testing_data labels_test accuracy
                    end
                    DAtemp2(condA, condB) = mean(DAtemp);
                    clear DAtemp
                    
               end % cond B
            end % cond A
        DA(perm, time_point) = mean(squareform([DAtemp2-50]')+50);
        clear DAtemp2
    end %time point
toc
end % permutation
DA_end=squeeze(nanmean(DA,[1]));  
save([resdir, '/', 'cl_obj_' modality '_sub' num2str(subn, '%02d')], 'DA_end');