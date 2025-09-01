% CLUSTER_1WAY - run hierarchical models on iEEG timeseries data with one
% within-channel fixed effect (2 levels), with the random effects of 
% subject and channel nested in subject. Perform cluster-based corrections 
% for multiple comparisons across timepoints, preserving the random effects
% structure.
%
% Download sample data file: https://drive.google.com/file/d/1i170RbUZuRQzY8gzEC99XBdo4uyCbsl4/view?usp=sharing
% Contains the following variables:
%   sid = subject ID
%   ch = channel label
%   hit_miss = hit (1) or miss (0) within-channel variable
%   data = timeseries data with one row per channel per level for the 
%       within-channel variable (186 rows x 139 timepoints)
%
% Copyright (c) 2025
% EL Johnson, PhD & AJO Dede, PhD

% set directories
datdir = pwd; % replace with your path to the downloaded data
savdir = datdir; % replace with your save path

% add path to subfunctions
addpath(fullfile(pwd, 'subfunctions'));

% load sample data file
data = load(fullfile(datdir, 'lme_data_1way'));

% initialize structure for model outputs
lme = [];
lme.F = nan(1, size(data.data, 2));
lme.p = lme.F;

% run model per timepoint, looping through timepoints
disp(' ');
for t = 1:size(data.data, 2)
    disp(['Running model on datapoint ' num2str(t) '/' num2str(size(data.data, 2)) '...']);

    % create table for model
    tmp_dat = table(data.sid, data.ch, data.hit_miss, data.data(:,t), ...
        'VariableNames', {'sid', 'ch', 'hit_miss', 'data'});
    
    % run model
    tmp_lme = fitlme(tmp_dat, 'data ~ hit_miss + (1|sid) + (1|sid:ch)');
    tmp_lme = anova(tmp_lme);

    % extract model outputs
    lme.F(t) = tmp_lme{2,2};
    lme.p(t) = tmp_lme{2,5};

    clear tmp*
end

% create null distributions for cluster testing
nperm = 1000; % number of permutations

hit_miss_main = lme.F'; % model stats

sid_ch = strcat(data.sid, '_', data.ch); % channel-nested-in-subject IDs
uid = unique(sid_ch); % unique IDs

% initialize null distribution
hit_miss_main_null = squeeze(zeros([size(hit_miss_main), nperm]));

% loop through permutations
disp(' ');
for p = 1:nperm
    disp(['Shuffling the data for permutation ' num2str(p) '/' num2str(nperm) '...']);

    hit_miss_shuff = zeros(length(data.hit_miss), 1);
    
    % loop through unique IDs
    for u = 1:length(uid)
        u_idx = find(ismember(sid_ch, uid{u}));
        
        % randomly shuffle within-channel variable labels
        if rand() > 0.5 % coin flip
            hit_miss_shuff(u_idx(1)) = 1; % set 1st row of pair to condition 1 (hit)
        else
            hit_miss_shuff(u_idx(2)) = 1; % set 2nd row of pair to condition 1 (hit)
        end
    end
    hit_miss_shuff = logical(hit_miss_shuff); % convert to boolean
    
    % run model per timepoint with shuffled labels, looping through timepoints
    for t = 1:size(data.data, 2)
        
        % create table for model
        tmp_dat = table(data.sid, data.ch, hit_miss_shuff, data.data(:,t), ...
            'VariableNames', {'sid', 'ch', 'hit_miss', 'data'});
        
        % run model
        tmp_lme = fitlme(tmp_dat, 'data ~ hit_miss + (1|sid) + (1|sid:ch)');
        tmp_lme = anova(tmp_lme);

        % extract model outputs
        hit_miss_main_null(t,p) = tmp_lme{2,2};
        
        clear tmp*
    end
end

clear data

% run cluster tests

% main effect of within-channel variable (hit/miss)
disp(' ');
[~,p,~] = cluster_test(hit_miss_main, hit_miss_main_null);
lme.p_clust = p(:)';

% save
save(fullfile(savdir, 'lme_clust_1way'), 'lme');
