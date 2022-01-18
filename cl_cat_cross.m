function cl_cat_cross(subn)
Numcats = 3;
%modality = 'vis' or 'aud'
modality = 'cross';
use_mvnn = 1;
shrink_time = 1;
ctg = [{{1:24}, {25:48}}; {{[1:12, 25:36]},{[13:24,37:48]}}; {{[1:6, 13:18, 25:30, 37:42]}, {[7:12, 19:24, 31:36, 43:48]}}];
DA_end = NaN(Numcats, 2, 50,50);
%% 0) add paths
addpath('/home/iampolina/OR/WORKFLOW/ANALYSIS')
addpath('/home/iampolina/OR/WORKFLOW/ANALYSIS/rerun_2021/mvnn')
addpath('/home/iampolina/OR/WORKFLOW/ANALYSIS/libsvm-3.11')
direeg = '/scratch/iampolina/OR/DATA/EEG';
resdir = fullfile('/scratch/iampolina/OR/DATA/CATCROSS2021/attempt06');
mkdir(resdir)
%% 1). Load  data load conditions X ...
%% X electodes X probe X time
for KK = 1:Numcats % 3 categories
  for mod = 1:2
    triggertemp=mod*100;
  
    if subn<11 || subn>24
        subdireeg = dir(fullfile(direeg, ['sub' num2str(subn, '%02d')], 'timelock_EM_excl.mat'));
    else
        if isequal(triggertemp, 100)
            subdireeg = dir(fullfile(direeg, ['sub' num2str(subn, '%02d') '_2'], 'timelock_EM_excl.mat'));
        else
            subdireeg = dir(fullfile(direeg, ['sub' num2str(subn, '%02d')], 'timelock_EM_excl.mat'));
        end
    end

    fileName = [subdireeg(1).folder, '/', subdireeg(1).name];
    load(fileName)

    for cat = 1:2
        counter=0;
        for cond = ctg{KK, cat}{1}+triggertemp                       
            new_ind = find(timelock.trialinfo == cond);
            counter=counter+1;
            CondMatrix{mod}(cat, counter,:, :, :)  = squeeze(nanmean(timelock.trial(new_ind, :, :),1));
            clear new_ind
        end
    end
    clear timelock subdireeg fileName
    %% 2). Dimensions:
    DC{mod} = size(CondMatrix{mod},1); % dimension condition
    DP{mod} = size(CondMatrix{mod},2); %dimension probe
    DE{mod} = size(CondMatrix{mod},3); % dimension electrode
    DT{mod} = size(CondMatrix{mod},4); % dimension time
    
    
    %% 3). Parameters for supertrials
    ['create pseudotrials']
    K=3;
    L=floor(DP{mod}/K);
    [H{mod},edges{mod}] = histcounts([1:DP{mod}],K);
    edges{mod} = ceil(edges{mod});
    edges{mod}(1) = edges{mod}(1) + 1;
  end
  
%% 4). Prepare output matrices
num_permutations = 100;

%% 5). Decoding

for perm =1:num_permutations 
    tic
    ['perm: ' num2str(perm)]
    for mod = 1:2
        % Shuffle #1
        permutedC{mod} = NaN(DC{mod},DP{mod},DE{mod},DT{mod}); %% CHANGE
        for co = 1:DC{mod}
            permutedC{mod}(co, :, :, :) = squeeze(CondMatrix{mod}(co, randperm(DP{mod}),:,:)); %% CHANGE
        end
        
        if shrink_time
            [permutedC{mod},Q] = sub_time(permutedC{mod}, 4);
        else
            Q = DT{1};
        end
       %% 6). Use multivariate noise normalization if set
        if use_mvnn
            permutedC{mod} = MVNN(permutedC{mod}, DP{mod});
        end
        % make supertrials
        for co = 1:DC{mod}
            for step= 1:K %average by steps
                    pseudo_trialD{mod}(co,step,:,:)= squeeze(mean(permutedC{mod}(co,edges{mod}(step):edges{mod}(step)+H{mod}(step)-1,:,:),2)); %assign pseudo trial to pseudo_trial_D
            end
        end
    end
    DA = NaN(num_permutations,2,Q,Q);

    
    for condA=1:DC{mod}%M %loop through all conditions %% CONDITIONS
        for condB = condA+1:DC{mod}%M  %loop through all conditions >condA+1 (to get all pair-wise combinations
            %['condA:' num2str(condA) ' ' 'condB: ' num2str(condB)]
            for time_point1 =1:Q % all time points are independent
                for time_point2 = 1:Q
                
                    % from vis to aud
                    training_data=[squeeze(pseudo_trialD{1}(condA,:,:, time_point1)) ; squeeze(pseudo_trialD{1}(condB,:,:,time_point1))];
                    labels_train=[ones(1,K) 2*ones(1,K)];
                    model = libsvmtrain(labels_train', training_data,'-s 0 -t 0 -q');
                    testing_data=[squeeze(pseudo_trialD{2}(condA,:,:,time_point2)); squeeze(pseudo_trialD{2}(condB,:,:,time_point2))];
                    labels_test= [ones(1,K) 2*ones(1,K)];
                    [predicted_label, accuracy, decision_values] = libsvmpredict(labels_test', testing_data , model);
                    DA(perm, 1, time_point1, time_point2) = accuracy(1);
                    clear training_data labels_train model testing_data labels_test accuracy
                    
                    % from aud to vis
                    training_data=[squeeze(pseudo_trialD{2}(condA,:,:, time_point1)) ; squeeze(pseudo_trialD{2}(condB,:,:,time_point1))];
                    labels_train=[ones(1,K) 2*ones(1,K)];
                    model = libsvmtrain(labels_train', training_data,'-s 0 -t 0 -q');
                    testing_data=[squeeze(pseudo_trialD{1}(condA,:,:,time_point2)); squeeze(pseudo_trialD{1}(condB,:,:,time_point2))];
                    labels_test= [ones(1,K) 2*ones(1,K)];
                    [predicted_label, accuracy, decision_values] = libsvmpredict(labels_test', testing_data , model);
                    DA(perm, 2, time_point1, time_point2) = accuracy(1);
                    clear training_data labels_train model testing_data labels_test accuracy  
                end
           end % time point
        end % cond A
    end %cond B
toc
end % permutation
DA_end(KK, :,:,:)=squeeze(nanmean(DA,[1])); % dimensions: category, modality (train on b, train on a), time, time
clear DA
end
save([resdir, '/', 'cl_cat_' modality '_sub' num2str(subn, '%02d')], 'DA_end');

%figure(2); plot(1:200, DA_end); %average across conditionsA&B