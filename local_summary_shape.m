% clear; clc; close all;
addpath('/Users/yejin/Documents/LocalSummaryShape/local_summary_shape/srvf');
addpath('/Users/yejin/Documents/LocalSummaryShape/local_summary_shape/functions');

folder_name = strcat('simul');
file_path   = strcat(folder_name,'/');
if ~exist(file_path, 'dir')
    mkdir(file_path); % create a folder to save results 
end

% import spatial shape data
data_path = '/Users/yejin/Documents/LocalSummaryShape/local_summary_shape/data/';
data_file_name = 'simul_data_50_range_1.mat';
data = load(strcat(data_path,data_file_name)).data;
data_matrix = data.data_matrix{1};
xy = data.xy{1};
x_min = data.spatial_domain(1);
x_max = data.spatial_domain(2);
y_min = data.spatial_domain(3);
y_max = data.spatial_domain(4);
range = (x_max - x_min)/2;
N_subregion = 4;
smoothing_param = 0.5;

% visualize spatial shape data
[fig,data_matrix_center_adj] = spatial_plot_shape(data_matrix, xy);
ylim([y_min, y_max]); 
xlim([x_min, x_max]);
xticks = get(gca, 'XTick');
yticks = get(gca, 'YTick');
set(gca, 'FontSize', 25);
set(gca, 'XTick', xticks);
set(gca, 'YTick', yticks);
tightfig;
saveas(fig, strcat(file_path,'spatial_shape_plot_',folder_name,'.png'));

% define local regions
subregions = split_square(x_min, x_max, y_min, y_max, sqrt(N_subregion));
ref        = compute_centers(subregions);
xs         = unique([subregions(:,1); subregions(:,2)]);
ys         = unique([subregions(:,3); subregions(:,4)]);
numRefs    = size(ref,1);
yline(ys,'--','Color','r','Linewidth',3);
xline(xs,'--','Color','r','Linewidth',3);
xy_local          = cell(numRefs,1);
data_matrix_local = cell(numRefs,1);
for i = 1:numRefs
    ind_local            = xy(:,1) >= subregions(i, 1) &...
                           xy(:,1) <= subregions(i, 2) &...
                           xy(:,2) >= subregions(i, 3) &...
                           xy(:,2) <= subregions(i, 4);
    xy_local{i}          = xy(ind_local,:);
    data_matrix_local{i} = data_matrix(:,:,ind_local);
end

ref_out        = cell(numRefs, 1);
X_s0_out       = cell(numRefs, 1);
q_s0_out       = cell(numRefs, 1);
lambda_out     = cell(numRefs, 1);
q_out          = cell(numRefs, 1);
delta_norm_out = cell(numRefs, 1);
psamp_out      = cell(numRefs, 1);

% get local summary
parfor i = 1:numRefs
    disp(i);
    s0         = ref(i, :);
    variog_file_path_name = fullfile(file_path, ...
       ['fitted_variog_at_', num2str(i)]);

    tic;
    [X_s0, q_s0, lambda_mat, q, delta_norm, ~] = ...
        Local_Summary0(s0, data_matrix_local{i}, ...
        xy_local{i}, variog_file_path_name, smoothing_param);
    toc;
    lambda_val     = lambda_mat(:,end);
    ref_out{i}     = s0;
    X_s0_out{i}    = X_s0;
    q_s0_out{i}    = q_s0;
    lambda_out{i}  = lambda_val;
    q_out{i}       = q;
    delta_norm_out{i} = delta_norm;
    
    psamp_out_pc1{i}   = PCdirection(q_s0, q, ...
        lambda_val,...
        1, strcat(folder_name,'/PC1_',num2str(i),'.png'));
    psamp_out_pc2{i}   = PCdirection(q_s0, q, ...
        lambda_val,...
        2, strcat(folder_name,'/PC2_',num2str(i),'.png'));
    psamp_out_pc3{i}   = PCdirection(q_s0, q, ...
        lambda_val,...
        3, strcat(folder_name,'/PC3_',num2str(i),'.png'));
end

% for i = 1:4
%     [~, K_out{i}] = PCdirection(out.q_s0{i},out.q{i}, ...
%         out.lambda{i},...
%         1, strcat('tmp.png'));
% end
% 
% for i = 1:4
%     [U,S,~]=svd(K_out{i});
%     s{i} = diag(S);
% end
% for i = 1:4
%     s{i} = s{i}./sum(s{i});
% end

out.ref         = ref_out;
out.X_s0        = X_s0_out;
out.q_s0        = q_s0_out;
out.lambda      = lambda_out;
out.q           = q_out;
out.delta_norm  = delta_norm_out;
out.psamp_pc1   = psamp_out_pc1;
out.psamp_pc2   = psamp_out_pc2;
out.psamp_pc3   = psamp_out_pc3;
out.xy          = xy;
save(fullfile(file_path, 'out.mat'), 'out');

psamp_mean = zeros(size(out.q{1},1),size(out.q{1},2),length(out.q));
midpoint = floor(size(out.psamp_pc1{1},3)/2) + 1;
for i = 1:size(psamp_mean,3)
psamp_mean(:,:,i) = out.psamp_pc1{i}(:,:,midpoint);
end
% obtain optimal rotation of psamp_mean wrt the mean
[~,~,~,~,~,Rbest,~] = FindElasticMean_extrinsic(psamp_mean);

[psamp_mean_aligned1, psamp_std_aligned1] = align_shapes(out.psamp_pc1, Rbest);
[psamp_mean_aligned2, psamp_std_aligned2] = align_shapes(out.psamp_pc2, Rbest);
[psamp_mean_aligned3, psamp_std_aligned3] = align_shapes(out.psamp_pc3, Rbest);

[fig1,min_dev1,max_dev1] = plot_local_summary_magnitude(ref, [x_min, x_max, y_min, y_max], xs, ys, ...
    psamp_mean_aligned1, psamp_std_aligned1, strcat(file_path,'local_smry1_pc1_',folder_name));
[fig2,min_dev2,max_dev2] = plot_local_summary_magnitude(ref, [x_min, x_max, y_min, y_max], xs, ys, ...
    psamp_mean_aligned2, psamp_std_aligned2, strcat(file_path,'local_smry1_pc2_',folder_name));
[fig3,min_dev3,max_dev3] = plot_local_summary_magnitude(ref, [x_min, x_max, y_min, y_max], xs, ys, ...
    psamp_mean_aligned3, psamp_std_aligned3, strcat(file_path,'local_smry1_pc3_',folder_name));

fig4 = plot_local_summary_direction(ref, [x_min, x_max, y_min, y_max], xs, ys, ...
    psamp_mean_aligned1, psamp_std_aligned1, strcat(file_path,'local_smry2_pc1_',folder_name));
fig5 = plot_local_summary_direction(ref, [x_min, x_max, y_min, y_max], xs, ys, ...
    psamp_mean_aligned2, psamp_std_aligned2, strcat(file_path,'local_smry2_pc2_',folder_name));
fig6 = plot_local_summary_direction(ref, [x_min, x_max, y_min, y_max], xs, ys, ...
    psamp_mean_aligned3, psamp_std_aligned3, strcat(file_path,'local_smry2_pc3_',folder_name));

min_min_dev = min([min_dev1,min_dev2,min_dev3]);
max_max_dev = max([max_dev1,max_dev2,max_dev3]);
clim([min_min_dev,max_max_dev]);