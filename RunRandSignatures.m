function [rand_stats]=RunRandSignatures(exprData,geneTissue,signature,residuals,numRands)
% RunRandSignatures - Run dose resgression based on the signature
%   Inputs: exprData - the imputed expression data matrix
%         geneTissue - the names of the gene-tissue names (feature names)
%         signature - the trained signature
%         residuals - the residuals (dependant variable)
%   Outputs: rand_stats - the statistics of the regression repeats (R2)

rand_stats=zeros(numRands,1);
warning('off','stats:regress:RankDefDesignMat');
rng('shuffle');
features_str=cell(size(geneTissue,1),1);
for k=1:size(geneTissue,1)
    features_str{k}=[geneTissue{k,1},geneTissue{k,2}];
end

sigFeatures_str=cell(size(signature,1),1);
for k=1:size(signature,1)
    sigFeatures_str{k}=[signature{k,1},signature{k,2}];
end
relFeatures_str=cell(size(geneTissue,1),1);
for k=1:size(geneTissue,1)
    relFeatures_str{k}=[geneTissue{k,1},geneTissue{k,2}];
end
top_features_idx=ismember(features_str,sigFeatures_str);
rel_features_idx=ismember(features_str,relFeatures_str);

remaining_features_idx=find((rel_features_idx-top_features_idx)>0);

if (sum(top_features_idx)>0)
   
    for r=1:numRands
        p=randperm(length(remaining_features_idx));
        rand_features_idx=remaining_features_idx(p(1:size(signature,1)),:);

        test_data=exprData(:,rand_features_idx);
        nonconst_idx=std(test_data)>eps;
        test_data=zscore(test_data(:,nonconst_idx));
        [~,~,~,~,stats]=regress(residuals,[ones(size(residuals,1),1),test_data]);
        rand_stats(r,1)=stats(1);
    end
end

warning('on','stats:regress:RankDefDesignMat');

