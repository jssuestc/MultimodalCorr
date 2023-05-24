clear;clc
n_rot=5000;
corr_type='spearman';
nroi_cortical=210; nroi_subcortical=36;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vertex_coords_sp, faces]= read_surf('\fsaverage\surf\lh.sphere.reg');
[vertices, label_annot, colortable] = read_annotation('\fsaverage\label\lh.BN_Atlas.annot');
ind=0;          % iteration counter
lh_centroid=[];    % initialisation of centroid array
for ic=1:colortable.numEntries      % loop over parcellated structures
    if isempty(strfind(colortable.struct_names{ic},'Unknown')) && isempty(strfind(colortable.struct_names{ic},'corpus')) % exclude "unknown" structures and corpus callosum from the parcellation
        ind=ind+1;                                                  % increment counter for every valid region
        label=colortable.table(ic,5);                               % ID of current parcel
        vertices_reg=find(label_annot==label);                      % vertices with given ID
        lh_centroid(ind,:)=mean(vertex_coords_sp(vertices_reg,:));
    end% average coordinates of all vertices within the current parcel to generate the centroid
end
lh_centroid=lh_centroid(all(~isnan(lh_centroid),2),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vertex_coords_sp, faces]= read_surf('\fsaverage\surf\rh.sphere.reg');
[vertices, label_annot, colortable] = read_annotation('\fsaverage\label\rh.BN_Atlas.annot');
ind=0;          % iteration counter
rh_centroid=[];    % initialisation of centroid array
for ic=1:colortable.numEntries      % loop over parcellated structures
    if isempty(strfind(colortable.struct_names{ic},'Unknown')) && isempty(strfind(colortable.struct_names{ic},'corpus')) % exclude "unknown" structures and corpus callosum from the parcellation
        ind=ind+1;                                                  % increment counter for every valid region
        label=colortable.table(ic,5);                               % ID of current parcel
        vertices_reg=find(label_annot==label);                      % vertices with given ID
        rh_centroid(ind,:)=mean(vertex_coords_sp(vertices_reg,:));
    end% average coordinates of all vertices within the current parcel to generate the centroid
end
rh_centroid=rh_centroid(all(~isnan(rh_centroid),2),:); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
perm_id_cortical = rotate_parcellation(lh_centroid, rh_centroid, n_rot); % rotate cortical parcellations
% permutated subcortical parcellations by randoming shuffling order
for w=1:n_rot
    perm_id_subcortical(:,w)=randperm(nroi_subcortical)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('\demo\input\map1.mat') % map1 for spatial correlation
load('\demo\input\map2.mat') % map2 for spatial correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[p_spin, r_dist,x_perm,y_perm]= perm_sphere_p(map1, map2, perm_id_cortical, perm_id_subcortical,corr_type);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[p_spin, r_dist,x_perm,y_perm] = spin_test(map1, map2,surface_name,parcellation_name,n_rot,correlation_type);






