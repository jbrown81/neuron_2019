load('brown_neuron_2019_data.mat');
% dxs_152: diagnosis for scan, 0=bvftd, 1=svppa
% interscan_interval_152: elapsed time in days between baseline scan and followup scan

% nii file with all 246 Brainnetome cortical/subcortical ROIs + (unused) 27 SUIT cerbellar ROIs:
% brainnetome_SUIT_273_rois.nii

% excluded nodes (lower fMRI SNR + cerebellar)
out_242=[101 102 117 118 247:273];
% included nodes
in_242=setdiff(1:273,out_242);
% indices for 242 includes nodes from set of 273
nodes_242_273_lut=[in_242' (1:242)'];
nodes_273_242_lut=[];

% epicenters
epicenters=epicenters_152;
epicenters_corr=epicenters_corr_152;
epicenters_corr_scaled=epicenters_corr./repmat(epicenters_corr(:,1),[1 194]);
top_epicenters=epicenters(:,1);

% euclidean distance matrix: 
% centers_dist_mat
inv_centers_dist_mat=1./centers_dist_mat;
% network distance (shortest path) matrix:
% distance_mat_242_thr100
% derived from connect_metrics_242.m

% functional connectivity matrix:
% connect_mat_242_thr100

% in paper, criteria for multiple epicenters:
% atrophy / seed FC correlations that are within 95% of top atrophy / FC seed correlation
epis_top_n={};
for i=1:152
    cur_epi_top_n=epicenters(i,find(epicenters_corr_scaled(i,:)>.95));
    epis_top_n{i}=cur_epi_top_n;
end

multi_epi=false; % allow multiple epicenters for a subject
hazard_euclidean=true; % scale nodal hazard by Euclidean distance
n_neighbors_keep=5; % in paper, used top 5 neighbors for hazard
n_regions=242;

scan_num=1; % choose any value from 1-152
subj=subjs_152(scan_num); % corresponding subject ID
cur_epi_top_n=epis_top_n{scan_num};
cur_epi=epicenters(scan_num,1);
baseline=baseline_all(:,i); % baseline atrophy (W-scores)
followup=followup_all(:,i); % followup atrophy (W-scores)

for k=1:242
    clear neighbors neighbor_weights n_neighbors_keep_cur top_neighbors
    neighbors=find(connect_mat_242_thr100(k,:));
    neighbors_count(k)=length(neighbors);
    neighbor_weights=connect_mat_242_thr100(k,neighbors);
    [s1 i1]=sort(neighbor_weights,'descend');
    n_neighbors_keep_cur=min(length(neighbors),n_neighbors_keep);
    top_neighbors=neighbors(i1(1:n_neighbors_keep_cur));
    % nodal hazard calculation
    if hazard_euclidean
        hazard(k)=sum((connect_mat_242_thr100(k,top_neighbors).*baseline(top_neighbors)')./centers_dist_mat(k,top_neighbors));
    else
        hazard(k)=sum((connect_mat_242_thr100(k,top_neighbors).*baseline(top_neighbors)'));
    end
    
    % spatial hazard calculation
    spatial_neighbors=find(inv_centers_dist_mat(k,:));
    spatial_neighbor_weights=inv_centers_dist_mat(k,spatial_neighbors);
    [s2 i2]=sort(spatial_neighbor_weights,'descend');
    spatial_n_neighbors_keep_cur=min(length(spatial_neighbors),n_neighbors_keep);
    spatial_top_neighbors=spatial_neighbors(i2(2:spatial_n_neighbors_keep_cur+1));
    spatial_hazard(k)=sum((inv_centers_dist_mat(k,spatial_top_neighbors).*baseline(spatial_top_neighbors)'));

    % shortest path to epicenter calculation
    % based on pathlength distance matrix 
    if multi_epi
        top_n_dists=distance_mat_242_thr100(cur_epi_top_n,k);
        cur_epi_local=cur_epi_top_n(find(top_n_dists==min(top_n_dists)));
        if length(cur_epi_local)>1 % if two epicenters tie for their distance from the current node, choose one epicenter at random
            cur_r=ceil(rand(1)*length(cur_epi_local));
            cur_epi_local=cur_epi_local(cur_r);
        end
        spe(k)=distance_mat_242_thr100(cur_epi_local,k);
    else
        spe(k)=distance_mat_242_thr100(cur_epi,k);
    end
    
    % euclidean distance to epicenter
    if multi_epi
        eucdist(k)=centers_dist_mat(cur_epi_local,k);
    else
        eucdist(k)=centers_dist_mat(cur_epi,k);
    end
    
    % atrophy change from baseline to followup
    change(k)=followup(i)-baseline(i);
end