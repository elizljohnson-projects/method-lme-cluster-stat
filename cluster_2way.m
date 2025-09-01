% CLUSTER_2WAY - run hierarchical models on iEEG timeseries data with one
% within-channel fixed effect (2 levels) and one between-channel fixed 
% effect (2 levels), with the random effects of subject and channel nested 
% in subject. Perform cluster-based corrections for multiple comparisons 
% across timepoints, preserving the random effects structure.
%
% Download sample data file: https://drive.google.com/file/d/1AYk6KAdVkZbDSJkTKwI4XJCmdEofD9rV/view?usp=sharing
% Contains the following variables:
%   sid = subject ID
%   ch = channel label
%   hit_miss = hit (1) or miss (0) within-channel variable
%   region = brain region (phrc or hc) between-channel variable
%   data = timeseries data with one row per channel per level for the 
%       within-channel variable (338 rows x 139 timepoints)
%
% Copyright (c) 2025
% EL Johnson, PhD & AJO Dede, PhD

% set directories
datdir = pwd; % replace with your path to the downloaded data
savdir = datdir; % replace with your save path

% add path to subfunctions
addpath(fullfile(pwd, 'subfunctions'));

% load sample data file
data = load(fullfile(datdir, 'lme_data_2way'));

% initialize structure for model outputs
lme = [];
lme.hit_miss = []; % hit/miss main effect
lme.hit_miss.F = nan(1, size(data.data, 2));
lme.hit_miss.p = lme.hit_miss.F;
lme.region = lme.hit_miss; % region main effect
lme.int = lme.hit_miss; % hit/miss*region interaction

% run model per timepoint, looping through timepoints
disp(' ');
for t = 1:size(data.data, 2)
    disp(['Running model on datapoint ' num2str(t) '/' num2str(size(data.data, 2)) '...']);

    % create table for model
    tmp_dat = table(data.sid, data.ch, data.hit_miss, data.region, data.data(:,t), ...
        'VariableNames', {'sid', 'ch', 'hit_miss', 'region', 'data'});
    
    % run model
    tmp_lme = fitlme(tmp_dat, 'data ~ hit_miss * region + (1|sid) + (1|sid:ch)');
    tmp_lme = anova(tmp_lme);

    % extract model outputs
    lme.hit_miss.F(t) = tmp_lme{2,2};
    lme.hit_miss.p(t) = tmp_lme{2,5};
    lme.region.F(t) = tmp_lme{3,2};
    lme.region.p(t) = tmp_lme{3,5};
    lme.int.F(t) = tmp_lme{4,2};
    lme.int.p(t) = tmp_lme{4,5};

    clear tmp*
end

% create null distributions for cluster testing
nperm = 1000; % number of permutations

hit_miss_main = lme.hit_miss.F'; % model stats
region_main = lme.region.F';
int = lme.int.F';

sid_ch = strcat(data.sid, '_', data.ch); % channel-nested-in-subject IDs
uid = unique(sid_ch); % unique IDs

regions = unique(data.region); % region labels

% initialize null distributions
hit_miss_main_null = squeeze(zeros([size(hit_miss_main), nperm]));
region_main_null = squeeze(zeros([size(region_main), nperm]));
int_null = squeeze(zeros([size(int), nperm]));

% loop through permutations
disp(' ');
for p = 1:nperm
    disp(['Shuffling the data for permutation ' num2str(p) '/' num2str(nperm) '...']);

    hit_miss_shuff = zeros(length(data.hit_miss), 1);
    region_shuff = cell(length(data.region), 1);
    n_region1 = sum(ismember(data.region, regions{1})) / length(unique(data.hit_miss));
    n_region2 = sum(ismember(data.region, regions{2})) / length(unique(data.hit_miss));
    
    % loop through unique IDs
    for u = 1:length(uid)
        u_idx = find(ismember(sid_ch, uid{u}));
        
        % randomly shuffle within-channel variable labels
        if rand() > 0.5 % coin flip
            hit_miss_shuff(u_idx(1)) = 1; % set 1st row of pair to condition 1 (hit)
        else
            hit_miss_shuff(u_idx(2)) = 1; % set 2nd row of pair to condition 1 (hit)
        end
        
        % randomly shuffle between-channel variable labels
        if rand() > (n_region1 / (n_region1 + n_region2))
            region_shuff{u_idx(1)} = regions{1}; % set both rows to region 1
            region_shuff{u_idx(2)} = regions{1};
            n_region2 = n_region2 - 1;
        else
            region_shuff{u_idx(1)} = regions{2}; % set both rows to region 2
            region_shuff{u_idx(2)} = regions{2};
            n_region1= n_region1 - 1;
        end
    end
    hit_miss_shuff = logical(hit_miss_shuff); % convert to boolean
    
    % run model per timepoint with shuffled labels, looping through timepoints
    for t = 1:size(data.data, 2)
        
        % create table for model
        tmp_dat = table(data.sid, data.ch, hit_miss_shuff, region_shuff, data.data(:,t), ...
            'VariableNames', {'sid', 'ch', 'hit_miss', 'region', 'data'});
        
        % run model
        tmp_lme = fitlme(tmp_dat, 'data ~ hit_miss * region + (1|sid) + (1|sid:ch)');
        tmp_lme = anova(tmp_lme);

        % extract model outputs
        hit_miss_main_null(t,p) = tmp_lme{2,2};
        region_main_null(t,p) = tmp_lme{3,2};
        int_null(t,p) = tmp_lme{4,2};
        
        clear tmp*
    end
end

clear data

% run cluster tests

% main effect of within-channel variable (hit/miss)
disp(' ');
[~,p,~] = cluster_test(hit_miss_main, hit_miss_main_null);
lme.hit_miss.p_clust = p(:)';

% main effect of between-channel variable (region)
disp(' ');
[~,p,~] = cluster_test(region_main, region_main_null);
lme.region.p_clust = p(:)';

% interaction (hit/miss*region)
disp(' ');
[~,p,~] = cluster_test(int, int_null);
lme.int.p_clust = p(:)';

% save
save(fullfile(savdir, 'lme_clust_2way'), 'lme');
