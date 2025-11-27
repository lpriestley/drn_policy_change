%% Timecourse fMRI analysis %%
clearvars; close all; clc;

% Define repo root
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(mfilename('fullpath')), 'fmt'));


GLM = 'S5.3';
stat = 1;
textout = 0;
getpeak = 0;

timelock = 'offer-timelocked';
if strcmp(timelock, 'offer-timelocked')
    epochSuffix = '_epoched_offer';
elseif strcmp(timelock, 'action-timelocked')
    epochSuffix = '_epoched_action';
elseif strcmp(timelock, 'revision')
    epochSuffix = 'epoched';
end

timecourse_dir = fullfile(repo_root, 'data', 'fMRI', 'timecourse_raw', filesep);
tc_graph_dir   = fullfile(repo_root, 'data', 'fMRI', 'timecourse_outputs', filesep);

feat_name = 'GLM01';

write_data = 1;

%% Select subjects + ROIs
subjects = {'S103','S104','S105','S106','S107','S109','S110','S111','S112', 'S114', 'S115','S116','S117','S118','S119','S120','S121','S122','S123','S124','S125','S126', 'S128', 'S130', 'S131','S132','S133','S135'};

if ismember(GLM, {'4.2a', '4.2b', '4.2c', '4.2d', '4.2a_low', '4.2c_low', 'S5.3', 'S5.4'}) || ~isempty(regexp(GLM, '_win', 'once'))
    roi = {'DRN_custom'};
end

if ismember(GLM, {'4.1a', '4.1b','S7.1', '4.1c'})
    roi = {'DRN_custom', 'MBD', 'HB', 'BF', 'LC_Pauli'};
end

if ismember(GLM, {'S7.2'})
    roi = {'DRN_custom', 'BF'};
end

if ismember(GLM, {'4.5'})
    roi = {'DRN_custom', 'dACCsphere7', 'AIsphere7'};
end

if ismember(GLM, {'4.4'})
    roi = {'dACCsphere7', 'AIsphere7'};
end

if ismember(GLM, {'S6.1', 'S6.2'})
    roi = {'MBD'};
end

if ismember(GLM, {'4.3a'})
    seed_roi = {'dACCsphere7'}; roi = {'DRN_custom'};
end

if ismember(GLM, {'4.3b'})
    seed_roi = {'AIsphere7'}; roi = {'DRN_custom'};
end

if ismember(GLM, {'4.1a_supp'})
    roi = {'DRN_custom', 'fourth_ventricle', 'MRN'};
end

if ismember(GLM, {'4.2a_supp'})
    roi = {'DRN_custom', 'LC_Pauli'};
end

if ismember(GLM, {'4.2c_supp'})
    roi = {'DRN_custom', 'LC_Pauli'};
end

if ismember(GLM, {'4.2a_cortex', '4.2b_cortex', '4.2c_cortex', '4.2d_cortex'})
    roi = {'DRN_custom', 'dACCsphere7', 'AIsphere7'};
end

behavAll = csvread(fullfile(repo_root, 'data', 'behaviour', 'raw_lf.csv'), 1, 0);

%% Run GLM
for r = 1:numel(roi)

    z = 0;
    trial_data = [];

    subjects = unique(behavAll(:,1));

    for is=1:length(subjects)

        z = z + 1;
        s = subjects(is);

        % load data
        subStr = ['S', num2str(s)];

        % prepare behaviour
        dataBehav = behavAll(behavAll(:,1) == s, :);

        % remove trials in subjects that scanner crashed before the end
        if      s == 103 %remove the last trial
            dataBehav(end,:) = [];
        elseif  s == 104 %remove the last 7 trials
            dataBehav((end-6:end),:) = [];
        end

        % Specify generic regressors
        trial = dataBehav(:,2);
        offer = dataBehav(:,3);
        decision = dataBehav(:,4);
        trial_in_block = dataBehav(:,9);
        mu_val = dataBehav(:, 13);
        env_bin = dataBehav(:,28);
        prev_action = dataBehav(:, 46);
        prev_policy = dataBehav(:, 47);
        mu_val_pe = offer - mu_val;

        window_width = 5;

        % identify policy + action changes
        policy_change = NaN(1, length(trial))';
        for i=1:length(policy_change)
            if ~isnan(prev_policy(i))
                if prev_policy(i)==decision(i)
                    policy_change(i)=0;
                elseif prev_policy(i)~=decision(i)
                    policy_change(i)=1;
                end
            end
        end

        action_change = NaN(1, length(trial))';
        for i=1:length(action_change)
            if ~isnan(prev_action(i))
                if prev_action(i)==decision(i)
                    action_change(i)=0;
                elseif prev_action(i)~=decision(i)
                    action_change(i)=1;
                end
            end
        end

        pursue_change = NaN(1, length(trial))';
        for i=1:length(pursue_change)
            if ~isnan(prev_policy(i))
                if prev_policy(i)==0 & decision(i)==1
                    pursue_change(i)=1;
                else
                    pursue_change(i)=0;
                end
            end
        end

        reject_change = NaN(1, length(trial))';
        for i=1:length(reject_change)
            if ~isnan(prev_policy(i))
                if prev_policy(i)==1 & decision(i)==0
                    reject_change(i)=1;
                else
                    reject_change(i)=0;
                end
            end
        end

        congruent_change = NaN(1, length(trial))';
        for i=1:length(congruent_change)
            if ~isnan(prev_policy(i))
                if reject_change(i)==1 & env_bin(i)==1 % if they switch to reject in the rich environment
                    congruent_change(i)=1;
                elseif pursue_change(i)==1 & env_bin(i)==-1% if they switch to pursue in the poor environment
                    congruent_change(i)=1;
                else
                    congruent_change(i)=0;
                end
            end
        end

        incongruent_change = NaN(1, length(trial))';
        for i=1:length(incongruent_change)
            if ~isnan(prev_policy(i))
                if reject_change(i)==1 & env_bin(i)==-1 % if they switch to reject in the poor environment
                    incongruent_change(i)=1;
                elseif pursue_change(i)==1 & env_bin(i)==1% if they switch to pursue in the rich environment
                    incongruent_change(i)=1;
                else
                    incongruent_change(i)=0;
                end
            end
        end

        congruent_vs_incongruent = NaN(1, length(trial))';
        for i=1:length(congruent_vs_incongruent)
            if ~isnan(prev_policy(i))
                if congruent_change(i)==1
                    congruent_vs_incongruent(i)=1;
                elseif incongruent_change(i)==1
                    congruent_vs_incongruent(i)=-1;
                else
                    congruent_vs_incongruent(i)=0;
                end
            end
        end

        congruent_poor = NaN(1, length(trial))';
        for i=1:length(congruent_poor)
            if ~isnan(prev_policy(i))
                if congruent_change(i)==1 & offer(i)==10 & env_bin(i)==-1
                    congruent_poor(i)=1;
                else
                    congruent_poor(i)=0;
                end
            end
        end

        congruent_poor_low = NaN(1, length(trial))';
        for i=1:length(congruent_poor_low)
            if ~isnan(prev_policy(i))
                if congruent_change(i)==1 & offer(i)==5 & env_bin(i)==-1
                    congruent_poor_low(i)=1;
                else
                    congruent_poor_low(i)=0;
                end
            end
        end

        incongruent_poor = NaN(1, length(trial))';
        for i=1:length(incongruent_poor)
            if ~isnan(prev_policy(i))
                if incongruent_change(i)==1 & offer(i)==10 & env_bin(i)==-1
                    incongruent_poor(i)=1;
                else
                    incongruent_poor(i)=0;
                end
            end
        end

        congruent_rich = NaN(1, length(trial))';
        for i=1:length(congruent_rich)
            if ~isnan(prev_policy(i))
                if congruent_change(i)==1 & offer(i)==10 & env_bin(i)==1
                    congruent_rich(i)=1;
                else
                    congruent_rich(i)=0;
                end
            end
        end

        congruent_rich_low = NaN(1, length(trial))';
        for i=1:length(congruent_rich_low)
            if ~isnan(prev_policy(i))
                if congruent_change(i)==1 & offer(i)==5 & env_bin(i)==1
                    congruent_rich_low(i)=1;
                else
                    congruent_rich_low(i)=0;
                end
            end
        end

        incongruent_rich = NaN(1, length(trial))';
        for i=1:length(incongruent_rich)
            if ~isnan(prev_policy(i))
                if incongruent_change(i)==1 & offer(i)==10 & env_bin(i)==1
                    incongruent_rich(i)=1;
                else
                    incongruent_rich(i)=0;
                end
            end
        end

        congruent_rich_vs_congruent_poor = NaN(1, length(trial))';
        for i=1:length(congruent_rich_vs_congruent_poor)
            if ~isnan(prev_policy(i))
                if congruent_poor(i)==1
                    congruent_rich_vs_congruent_poor(i)=1;
                elseif congruent_rich(i)==1
                    congruent_rich_vs_congruent_poor(i)=-1;
                else
                    congruent_rich_vs_congruent_poor(i)=0;
                end
            end
        end

        % identify trials (option-agnostic) and encounters
        % (option-specific) after policy changes
        trial_after_policy_change = NaN(size(policy_change));
        trial_after_policy_change([false; ~isnan(policy_change(1:end-1))]) = 0;
        trial_after_policy_change([false; policy_change(1:end-1) == 1]) = 1;

        N = length(policy_change);
        encounter_after_policy_change = NaN(N,1);
        valid_mag = ~isnan(offer);
        encounter_after_policy_change(valid_mag) = 0;

        u = unique(offer(valid_mag));
        for k = 1:numel(u)
            x = u(k);
            idx = find(offer == x);
            pc_rows = idx(policy_change(idx) == 1);
            for t = 1:numel(pc_rows)
                p = pc_rows(t);
                nxt_pos = find(idx > p, 1, 'first');
                if ~isempty(nxt_pos)
                    encounter_after_policy_change(idx(nxt_pos)) = 1;
                end
            end
        end

        % Prepare regressors
        % main
        REG.trial = trial;
        REG.offer = offer;
        REG.decision = decision;
        REG.mu_val = mu_val;
        REG.env_bin = env_bin;
        REG.prev_action = prev_action;
        REG.prev_policy = prev_policy;
        REG.trial_in_block = trial_in_block;

        REG.policy_change = policy_change;
        REG.action_change = action_change;
        REG.trial_after_policy_change = trial_after_policy_change;
        REG.encounter_after_policy_change = encounter_after_policy_change;

        REG.congruent_poor = congruent_poor;
        REG.congruent_rich = congruent_rich;
        REG.incongruent_poor = incongruent_poor;
        REG.incongruent_rich = incongruent_rich;
        REG.congruent_poor_low = congruent_poor_low;
        REG.congruent_rich_low = congruent_rich_low;
        REG.congruent_rich_vs_congruent_poor = congruent_rich_vs_congruent_poor;

        REG.congruent_change = congruent_change;
        REG.incongruent_change = incongruent_change;
        REG.mu_val_pe = mu_val_pe;

        REG.constant = ones(length(REG.trial), 1);

        % load time-series data
        load(fullfile(timecourse_dir, subStr, feat_name, ['tc_' roi{r} epochSuffix]));


        %% Run GLMs

        % GLM4.1: Effect of policy-change
        if  strcmp(GLM, '4.1a')|strcmp(GLM, '4.1a_supp')

            % regressor design matrix
            dmat =  [REG.policy_change, REG.trial, REG.constant];

            dmat(isnan(REG.policy_change),:)=[];
            trial_data(isnan(REG.policy_change),:)=[];

            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));

            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % write effect size
            if  textout
                filename = [timecourse_dir,s,'/',session,'/',feat_name,'/',([roi{r},'_effectSize'])];
                save(filename,'betaOut');
            end
        end

        % GLM4.1b: Effect of action-change
        if  strcmp(GLM, '4.1b')

            % regressor design matrix
            dmat =  [REG.action_change, REG.trial, REG.constant];

            dmat(isnan(REG.action_change),:)=[];
            trial_data(isnan(REG.action_change),:)=[];

            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));

            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % write effect size
            if  textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename,'betaOut');
            end
        end

        % GLM S5: Effect of trial after policy-change
        if  strcmp(GLM, 'S5.3')

            % regressor design matrix
            dmat =  [REG.trial_after_policy_change, REG.trial, REG.constant];

            dmat(isnan(REG.trial_after_policy_change),:)=[];
            trial_data(isnan(REG.trial_after_policy_change),:)=[];

            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));

            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % write effect size
            if  textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename,'betaOut');
            end
        end

        % GLM SD: Effect of encounter after policy-change
        if  strcmp(GLM, 'S5.4')

            % regressor design matrix
            dmat =  [REG.encounter_after_policy_change, REG.trial, REG.constant];

            dmat(isnan(REG.encounter_after_policy_change),:)=[];
            trial_data(isnan(REG.encounter_after_policy_change),:)=[];

            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));

            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % write effect size
            if  textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename,'betaOut');
            end
        end

        % GLM4.2a: effect of congruent pursue change
        if strcmp(GLM, '4.2a')|strcmp(GLM, '4.2a_supp')|strcmp(GLM, '4.2a_cortex')

            % regressor design matrix (don't forget constant!)
            dmat =  [REG.congruent_poor, REG.trial, REG.constant];

            dmat(isnan(REG.congruent_poor),:)=[];
            trial_data(isnan(REG.congruent_poor),:)=[];

            % filter mid-opt trials only
            mag_index = offer;
            mag_index(isnan(REG.congruent_poor),:)=[];

            % filter according to environment
            dmat = dmat(mag_index==10,:);
            trial_data = trial_data(mag_index==10,:);

            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));

            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % write effect size
            if  textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename,'betaOut');
            end
        end

        % GLM4.2b: effect of incongruent pursue change
        if  strcmp(GLM, '4.2b')|strcmp(GLM, '4.2b_cortex')

            % regressor design matrix
            dmat =  [REG.incongruent_rich, REG.trial, REG.constant];

            dmat(isnan(REG.incongruent_rich),:)=[];
            trial_data(isnan(REG.incongruent_rich),:)=[];

            % filter mid-opt trials only
            mag_index = offer;
            mag_index(isnan(REG.incongruent_rich),:)=[];

            % filter according to environment
            dmat = dmat(mag_index==10,:);
            trial_data = trial_data(mag_index==10,:);

            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));

            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % write effect size
            if  textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename,'betaOut');
            end
        end

        % GLM4.2c: effect of congruent reject change
        if  strcmp(GLM, '4.2c')|strcmp(GLM, '4.2c_supp')|strcmp(GLM, '4.2c_cortex')

            % regressor design matrix (don't forget constant!)
            dmat =  [REG.congruent_rich, REG.trial, REG.constant];

            dmat(isnan(REG.congruent_rich),:)=[];
            trial_data(isnan(REG.congruent_rich),:)=[];

            % filter mid-opt trials only
            mag_index = offer;
            mag_index(isnan(REG.congruent_rich),:)=[];
            dmat = dmat(mag_index==10,:);
            trial_data = trial_data(mag_index==10,:);

            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));

            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % write effect size
            if  textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename,'betaOut');
            end
        end

        % GLM4.2d: effect of incongruent reject change
        if strcmp(GLM, '4.2d')|strcmp(GLM, '4.2d_cortex')

            % regressor design matrix (don't forget constant!)
            dmat =  [REG.incongruent_poor, REG.trial, REG.constant];

            dmat(isnan(REG.incongruent_poor),:)=[];
            trial_data(isnan(REG.incongruent_poor),:)=[];

            % filter mid-opt trials only
            mag_index = offer;
            mag_index(isnan(REG.incongruent_poor),:)=[];
            dmat = dmat(mag_index==10,:);
            trial_data = trial_data(mag_index==10,:);

            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));

            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % write effect size
            if  textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename,'betaOut');
            end
        end

        % GLM4.2a_low: effect of congruent pursue change in 5-pt options
        if  strcmp(GLM, '4.2a_low')

            % regressor design matrix (don't forget constant!)
            dmat =  [REG.congruent_poor_low, REG.trial, REG.constant];

            dmat(isnan(REG.congruent_poor_low),:)=[];
            trial_data(isnan(REG.congruent_poor_low),:)=[];

            % filter low-opt trials only
            mag_index = offer;
            mag_index(isnan(REG.congruent_poor_low),:)=[];
            dmat = dmat(mag_index==5,:);
            trial_data = trial_data(mag_index==5,:);

            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));

            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % write effect size
            if  textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename,'betaOut');
            end
        end

        % GLM4.2c: effect of congruent reject change in 5-pt options
        if strcmp(GLM, '4.2c_low')

            % regressor design matrix (don't forget constant!)
            dmat =  [REG.congruent_rich_low, REG.trial, REG.constant];

            dmat(isnan(REG.congruent_rich_low),:)=[];
            trial_data(isnan(REG.congruent_rich_low),:)=[];

            % filter low-opt trials only
            mag_index = offer;
            mag_index(isnan(REG.congruent_rich_low),:)=[];
            dmat = dmat(mag_index==5,:);
            trial_data = trial_data(mag_index==5,:);

            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));

            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % write effect size
            if  textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename,'betaOut');
            end
        end

        % GLM4.3: PPI as a function of policy-switch
        if  strcmp(GLM, '4.3a')||strcmp(GLM, '4.3b') % PPI as a function of policy-switch

            % load time-series data
            seed_TC = load([timecourse_dir,subStr,'/',feat_name,'/',(['tc_' seed_roi{1} epochSuffix])]);
            roi_TC  = load([timecourse_dir,subStr,'/',feat_name,'/',(['tc_' roi{r} epochSuffix])]);
            seed_TC.trial_data(isnan(REG.env_bin),:)=[];

            % run GLM
            for i = 1:size(seed_TC.trial_data,2)

                % load ROI time-series data
                REG.TC  = roi_TC.trial_data(:,i);
                dmat = [REG.env_bin, REG.trial, REG.constant];
                dmat = [REG.TC, dmat];

                % remove trials with no response
                dmat(isnan(REG.env_bin),:)=[];

                % normalise data
                dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

                % create PPI regressor
                REG.PPI = zscore (dmat(:,1) .* dmat(:,2));
                dmat = [REG.PPI, dmat];

                % contrast design matrix
                contrasts = diag(ones(size(dmat,2),1));

                % beta X time output
                betaOut(:,i) = ols(seed_TC.trial_data(:,i),dmat,contrasts);
                clear dmat contrasts REG.TC REG.PPI
            end
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;
        end

        if strcmp(GLM, '4.4') % effect of rich-vs-poor environment

            % regressor design matrix (don't forget constant!)
            dmat =  [REG.env_bin, REG.trial, REG.constant];

            dmat(isnan(REG.env_bin),:)=[];
            trial_data(isnan(REG.env_bin),:)=[];

            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));

            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % write effect size
            if  textout
                filename = [timecourse_dir,s,'/',session,'/',feat_name,'/',([roi{r},'_effectSize'])];
                save(filename,'betaOut');
            end
        end

        if strcmp(GLM, '4.5') % environment-specific coding of congruent changes

            % regressor design matrix 
            dmat =  [REG.congruent_rich_vs_congruent_poor, REG.trial, REG.constant];

            dmat(isnan(REG.congruent_rich_vs_congruent_poor),:)=[];
            trial_data(isnan(REG.congruent_rich_vs_congruent_poor),:)=[];

            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));

            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % write effect size
            if  textout
                filename = [timecourse_dir,s,'/',session,'/',feat_name,'/',([roi{r},'_effectSize'])];
                save(filename,'betaOut');
            end
        end

        % GLM S6a: effect of response
        if  strcmp(GLM, 'S6.1')

            % regressor design matrix
            dmat =  [REG.decision, REG.trial, REG.constant];

            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));

            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % write effect size
            if  textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename,'betaOut');
            end
        end

        % GLM S6b: effect of value difference [rw.mag(t) - ave.value(t)]
        if  strcmp(GLM, 'S6.2')

            % regressor design matrix
            dmat =  [REG.mu_val_pe, REG.decision, REG.trial, REG.constant];

            dmat(isnan(REG.mu_val_pe),:)=[];
            trial_data(isnan(REG.mu_val_pe),:)=[];

            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));

            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % write effect size
            if  textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename,'betaOut');
            end
        end

        % GLM S7a: effect of incongruent (i.e. exploratory) policy-siwtches
        if  strcmp (GLM, 'S7.1')

            % regressor design matrix (don't forget constant!)
            dmat =  [REG.incongruent_change, REG.congruent_change, REG.trial, REG.constant];

            dmat(isnan(REG.incongruent_change),:)=[];
            trial_data(isnan(REG.incongruent_change),:)=[];

            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));

            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % write effect size
            if  textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename,'betaOut');
            end
        end

        % GLM S7b: effect of incongruent vs environment-specific switches
        if  strcmp (GLM, 'S7.2')

            % regressor design matrix (don't forget constant!)
            dmat =  [REG.incongruent_change, REG.congruent_rich_vs_congruent_poor, REG.trial, REG.constant];

            dmat(isnan(REG.incongruent_change)|isnan(REG.congruent_rich_vs_congruent_poor),:)=[];
            trial_data(isnan(REG.incongruent_change)|isnan(REG.congruent_rich_vs_congruent_poor),:)=[];

            contrasts = diag(ones(size(dmat,2),1));

            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % write effect size
            if  textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename,'betaOut');
            end
        end

        % GLM S5.1 window 1: effect of congruent pursue switch in trials 1-6
        if strcmp(GLM, 'S5.1_win1')

            win_start = 1;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_poor, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_poor = REG.congruent_poor(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_poor), :) = [];
            trial_data(isnan(congruent_poor), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.1 window 2: effect of congruent pursue switch in trials 2-7
        if strcmp(GLM, 'S5.1_win2')

            win_start = 2;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_poor, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_poor = REG.congruent_poor(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_poor), :) = [];
            trial_data(isnan(congruent_poor), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.1 window 3: effect of congruent pursue switch in trials 3-8
        if strcmp(GLM, 'S5.1_win3')

            win_start = 3;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_poor, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_poor = REG.congruent_poor(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_poor), :) = [];
            trial_data(isnan(congruent_poor), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.1 window 4: effect of congruent pursue switch in trials 4-9
        if strcmp(GLM, 'S5.1_win4')

            win_start = 4;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_poor, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_poor = REG.congruent_poor(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_poor), :) = [];
            trial_data(isnan(congruent_poor), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.1 window 5: effect of congruent pursue switch in trials 5-10
        if strcmp(GLM, 'S5.1_win5')

            % --- Define window across all blocks: trials 1–10 within each block ---
            win_start = 5;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_poor, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_poor = REG.congruent_poor(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_poor), :) = [];
            trial_data(isnan(congruent_poor), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.1 window 6: effect of congruent pursue switch in trials 6-11
        if strcmp(GLM, 'S5.1_win6')

            win_start = 6;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_poor, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_poor = REG.congruent_poor(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_poor), :) = [];
            trial_data(isnan(congruent_poor), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.1 window 7: effect of congruent pursue switch in trials 7-12
        if strcmp(GLM, 'S5.1_win7')

            win_start = 7;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_poor, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_poor = REG.congruent_poor(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_poor), :) = [];
            trial_data(isnan(congruent_poor), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.1 window 8: effect of congruent pursue switch in trials 8-13
        if strcmp(GLM, 'S5.1_win8')

            win_start = 8;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_poor, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_poor = REG.congruent_poor(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_poor), :) = [];
            trial_data(isnan(congruent_poor), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.1 window 9: effect of congruent pursue switch in trials
        % 9-14
        if strcmp(GLM, 'S5.1_win9')

            win_start = 9;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_poor, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_poor = REG.congruent_poor(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_poor), :) = [];
            trial_data(isnan(congruent_poor), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.1 window 10: effect of congruent pursue switch in trials
        % 10-15
        if strcmp(GLM, 'S5.1_win10')

            win_start = 10;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_poor, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_poor = REG.congruent_poor(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_poor), :) = [];
            trial_data(isnan(congruent_poor), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.2 window 1: effect of congruent reject switch in trials 1-5
        if strcmp(GLM, 'S5.2_win1')

            win_start = 1;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_rich, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_rich = REG.congruent_rich(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_rich), :) = [];
            trial_data(isnan(congruent_rich), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.2 window 2: effect of congruent reject switch in trials
        % 2-6
        if strcmp(GLM, 'S5.2_win2')

            win_start = 2;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_rich, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_rich = REG.congruent_rich(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_rich), :) = [];
            trial_data(isnan(congruent_rich), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.2 window 3: effect of congruent reject switch in trials
        % 3-8
        if strcmp(GLM, 'S5.2_win3')

            win_start = 3;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_rich, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_rich = REG.congruent_rich(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_rich), :) = [];
            trial_data(isnan(congruent_rich), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.2 window 4: effect of congruent reject switch in trials
        % 4-9
        if strcmp(GLM, 'S5.2_win4')

            win_start = 4;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_rich, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_rich = REG.congruent_rich(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_rich), :) = [];
            trial_data(isnan(congruent_rich), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.2 window 5: effect of congruent reject switch in trials
        % 5-10
        if strcmp(GLM, 'S5.2_win5')

            win_start = 5;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_rich, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_rich = REG.congruent_rich(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_rich), :) = [];
            trial_data(isnan(congruent_rich), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.2 window 6: effect of congruent reject switch in trials
        % 6-11
        if strcmp(GLM, 'S5.2_win6')

            win_start = 6;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_rich, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_rich = REG.congruent_rich(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_rich), :) = [];
            trial_data(isnan(congruent_rich), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.2 window 7: effect of congruent reject switch in trials
        % 7-12
        if strcmp(GLM, 'S5.2_win7')

            win_start = 7;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_rich, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_rich = REG.congruent_rich(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_rich), :) = [];
            trial_data(isnan(congruent_rich), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.2 window 8: effect of congruent reject switch in trials
        % 8-13
        if strcmp(GLM, 'S5.2_win8')

            win_start = 8;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_rich, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_rich = REG.congruent_rich(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_rich), :) = [];
            trial_data(isnan(congruent_rich), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.2 window 9: effect of congruent reject switch in trials
        % 9-14
        if strcmp(GLM, 'S5.2_win9')

            win_start = 9;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_rich, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_rich = REG.congruent_rich(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_rich), :) = [];
            trial_data(isnan(congruent_rich), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

        % GLM S5.2 window 10: effect of congruent reject switch in trials
        % 10-15
        if strcmp(GLM, 'S5.2_win10')

            % --- Define window across all blocks: trials 1–10 within each block ---
            win_start = 10;
            win_end = win_start + window_width;
            in_window = (REG.trial_in_block >= win_start) & (REG.trial_in_block <= win_end);

            % --- Remove all other trials *after* constructing REG but before GLM ---
            dmat =  [REG.congruent_rich, REG.trial, REG.constant];
            trial_data = trial_data(in_window,:);
            dmat = dmat(in_window,:);
            congruent_rich = REG.congruent_rich(in_window);

            % --- Remove NaNs in policy_change ---
            dmat(isnan(congruent_rich), :) = [];
            trial_data(isnan(congruent_rich), :) = [];

            % --- Contrast matrix ---
            contrasts = diag(ones(size(dmat,2),1));

            % --- Normalise regressors (excluding constant) ---
            dmat(:,1:end-1) = zscore(dmat(:,1:end-1));

            % --- Run OLS regression ---
            betaOut = ols(trial_data, dmat, contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % --- Save effect size if needed ---
            if textout
                filename = fullfile(timecourse_dir, subStr, feat_name, [roi{r} '_effectSize']);
                save(filename, 'betaOut');
            end
        end

    end
end
%% write data out for graphs
contrast_to_plot = 1;

for r = 1:numel(roi)

    beta = squeeze(allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(contrast_to_plot,:,:));

    if write_data
        tc_data = {};
        tc_data(:, 1) = num2cell(mean(beta'));
        tc_data(:, 2) = num2cell(std(beta')/sqrt(size(beta',1)));
        tc_data(:, 3) = roi(r);
        tc_data(:, 4) = {['GLM_0', num2str(GLM)]};
        tc_data(:, 5) = {['contr_', num2str(contrast_to_plot)]};
        if strcmp(GLM, '4.3')
            filename = fullfile(tc_graph_dir, ['GLM_0' num2str(GLM) '_contr_' num2str(contrast_to_plot) '_' seed_roi{1} '_to_' roi{r} '_PPI_timecourse.csv']);
        else
            filename = fullfile(tc_graph_dir, ['GLM_0' num2str(GLM) '_contr_' num2str(contrast_to_plot) '_' roi{r} '_' timelock '_timecourse.csv']);

        end
        writetable(cell2table(tc_data), filename);
        clear tc_data
    end

end


%% stats

if stat

    clearvars window LOOT peak
    numsession = 1:numel(subjects);
    pre_win = 2;
    post_win = 8;
    upsample = 20;
    TR = 1.775;
    nsamples = round(((pre_win+post_win)./TR)*upsample);
    start = 40;
    finish = max(nsamples);
    reg = contrast_to_plot;

    % find peak in a specified window
    for r = 1:numel(roi)

        for s = 1:length(numsession)

            numsession(s) = [];
            window  = mean(squeeze(allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(reg,:,numsession)),2);

            [m,i]  = max(abs(window(start:finish)));
            LocPeak(s,r)=(i*10)/length(window);
            i      = i + (start-1);
            peak(s,r) = allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(reg,i,s);
            numsession = 1:numel(subjects);
            clear window fw

        end
        [h,p(r),ci,stats]= ttest(peak(:,r));
        t(r) = stats.tstat;
    end

    [bonf_p, bonf_h] = bonf_holm(p)

    if write_data
        tb = num2cell(peak);
        tb(:, size(tb,2)+1) = {['GLM_0', num2str(GLM)]};
        tb(:, size(tb,2)+1) = {['contr', num2str(contrast_to_plot)]};
        tb(:, size(tb,2)+1) = {timelock};
        tb = cell2table(tb, 'VariableNames', [roi, 'GLM', 'contrast', 'timelock']);
        writetable(tb, fullfile(tc_graph_dir,['GLM_0' num2str(GLM) '_contr_' num2str(contrast_to_plot) '_' timelock '_peaks.csv']));
        clear tc_data
    end
end
