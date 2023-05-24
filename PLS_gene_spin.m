clear;clc
% load predictor vairables, regional gene expresseion 
gene_regional_expression = importdata('\demo\input\gene_regional_expression_zscored.mat');
Y_data = importdata('\demo\input\Ydata.mat'); % load response variables 
Y=Y_data; % Response variable 
% z-score:
Y=zscore(Y);
outputDir = '\demo\output';
if ~isdir(outputDir)
    mkdir(outputDir);
end
% Note that here we use the left hemisphere only
X=gene_regional_expression(:,2:end); % Predictors
X(isnan(X)) = 0;
X=zscore(X);
n_rot=5000; % number of rotations
nroi_cortical=105; nroi_subcortical=18; % number of nodes
%
[vertex_coords_sp, faces]= read_surf('\fsaverage\surf\lh.sphere.reg');
[vertices, label_annot, colortable] = read_annotation('\fsaverage\label\lh.BN_Atlas.annot');
ind=0;          % iteration counter
lh_centroid_cortical=[];    % initialisation of centroid array
for ic=1:colortable.numEntries      % loop over parcellated structures
    if isempty(strfind(colortable.struct_names{ic},'Unknown')) && isempty(strfind(colortable.struct_names{ic},'corpus')) % exclude "unknown" structures and corpus callosum from the parcellation
        ind=ind+1;                                                  % increment counter for every valid region
        label=colortable.table(ic,5);                               % ID of current parcel
        vertices_reg=find(label_annot==label);                      % vertices with given ID
        lh_centroid_cortical(ind,:)=mean(vertex_coords_sp(vertices_reg,:));
    end% average coordinates of all vertices within the current parcel to generate the centroid
end
lh_centroid_cortical=lh_centroid_cortical(all(~isnan(lh_centroid_cortical),2),:);
perm_id_cortical = spin_map_lh(lh_centroid_cortical, n_rot); % rotate cortical parcellations of the left hemisphere
% % permutated subcortical parcellations of the left hemisphere
for w=1:n_rot
    perm_id_subcortical(:,w)=randperm(nroi_subcortical)';
end
%
%perform full PLS and plot variance in Y explained by top 15 components
%typically top 2 or 3 components will explain a large part of the variance
%(hopefully!)
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y);
save('MCM_fun_pls_result', 'XL', 'YL', 'XS', 'YS', 'BETA', 'PCTVAR', 'MSE', 'stats');
dim=15;
plot(1:dim,cumsum(100*PCTVAR(2,1:dim)),'-o','LineWidth',1.5,'Color',[140/255,0,0]);
set(gca,'Fontsize',14)
xlabel('Number of PLS components','FontSize',14);
ylabel('Percent Variance Explained in Y','FontSize',14);
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_spin=[];
for dim=1:10
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
    temp=cumsum(100*PCTVAR(2,1:dim));
    Rsquared = temp(dim);
    for j=1:n_rot    
        Yspin_cortical=Y(1:nroi_cortical,:); Yspin_cortical=Yspin_cortical(perm_id_cortical(:,j));
        Yperm_subcortical=Y(nroi_cortical+1:end,:); Yperm_subcortical=Yperm_subcortical(perm_id_subcortical(:,j));
        Yspin=[Yspin_cortical;Yperm_subcortical];
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Yspin,dim);
        temp=cumsum(100*PCTVAR(2,1:dim));
        Rsq(j) = temp(dim);
    end
    p_spin(dim)=length(find(Rsq>=Rsquared))/n_rot;
end
figure
plot(1:dim, p_spin,'ok','MarkerSize',8,'MarkerFaceColor','r');
xlabel('Number of PLS components','FontSize',14);
ylabel('p-value','FontSize',14);
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('\demo\input\ProbeID.mat')
genes = genes;
geneindex = 1:10027;

%%% 导入Excel 表格
FileProbes = '\demo\input\Probes.xlsx';
ProbeTable = readtable(FileProbes);
ProbeIDAll = ProbeTable.probe_id;
%%% �存储
GeneID = ProbeTable.gene_id;
GeneSymbol = ProbeTable.gene_symbol;
GeneName = ProbeTable.gene_name;
EntrezID = ProbeTable.entrez_id;
Chromosome = ProbeTable.chromosome;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y);
XL1 = XL(:, 1);
%store regions' IDs and weights in descending order of weight for both components:
[R1,p1]=corr([XS(:,1),XS(:,2)],Y);
%align PLS components with desired direction for interpretability
if R1(1,1)<0  %this is specific to the data shape we were using - will need ammending
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0 %this is specific to the data shape we were using - will need ammending
    stats.W(:,2)=-1*stats.W(:,2);
    XS(:,2)=-1*XS(:,2);
end
[PLS1w,x1] = sort(stats.W(:,1),'descend');
PLS1ids=genes(x1);
XL1 = XL1(x1);
geneindex1=geneindex(x1);
[PLS2w,x2] = sort(stats.W(:,2),'descend');
PLS2ids=genes(x2);
geneindex2=geneindex(x2);

%print out results
%     csvwrite('PLS1_ROIscores.csv',XS(:,1));
%     csvwrite('PLS2_ROIscores.csv',XS(:,2));

%define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
PLS2weights=[];
%start bootstrap
bootnum = 5000;
for i=1:bootnum
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled subjects
    Yr=Y(myresample,:); % define X for resampled subjects
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr); %perform PLS for resampled data    
    temp=stats.W(:,1);%extract PLS1 weights
    newW=temp(x1); %order the newly obtained weights the same way as initial PLS
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights=[PLS1weights,newW];%store (ordered) weights from this bootstrap run    
    %     temp=stats.W(:,2);%extract PLS2 weights
    %     newW=temp(x2); %order the newly obtained weights the same way as initial PLS
    %     if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
    %         newW=-1*newW;
    %     end
    %     PLS2weights=[PLS2weights,newW]; %store (ordered) weights from this bootstrap run
end

%get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');
% PLS2sw=std(PLS2weights');
%get bootstrap weights
temp1=PLS1w./PLS1sw';
% temp2=PLS2w./PLS2sw';
%order bootstrap weights (Z) and names of regions
[Z1 ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
XL1 = XL1(ind1);
geneindex1=geneindex1(ind1);
% [Z2 ind2]=sort(temp2,'descend');
% PLS2=PLS2ids(ind2);
% geneindex2=geneindex2(ind2);
a = 1- normcdf(abs(Z1));
[n_signif, index]  = fdr(a, 0.05);% 
Z1 = Z1(index);
PLS1 = PLS1(index);
XL1 = XL1(index);
geneindex1 = geneindex1(index);
for i = 1:length(PLS1)
    ind = find(ProbeIDAll == PLS1(i));
    Gene_info_all{i, 1} = GeneID(ind);
    Gene_info_all{i, 2} = char(GeneSymbol(ind));
    Gene_info_all{i, 3} = char(GeneName(ind));
    Gene_info_all{i, 4} = EntrezID(ind);
    if ~isnan(Chromosome(ind))
        Gene_info_all{i, 5} = Chromosome(ind);
    else
        Gene_info_all{i, 5} = 'X';
    end
    Gene_info_all{i, 6} = PLS1(i);
    Gene_info_all{i, 7} = geneindex1(i);
    Gene_info_all{i, 8} = Z1(i);
    Gene_info_all{i, 9} = XL1(i);
end
title = {'gene_id', 'gene_symbol', 'gene_name', 'entrez_id', 'chromosome', 'ProbeID', 'geneindex', 'PLS_W_Z', 'XL'};
Gene_info_all_1 = [title; Gene_info_all];
xlswrite([outputDir, filesep, 'Significant_gene_info.xls'], Gene_info_all_1);
clear Gene_info_all_1
% end
