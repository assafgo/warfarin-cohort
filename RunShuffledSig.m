function [rand_stats]=RunShuffledSig(exprData,geneTissue,signature,residuals,numRands)
% RunRandSignatures - Run dose resgression based on the signature
%   Inputs: exprData - the imputed expression data matrix
%         geneTissue - the names of the gene-tissue names (feature names)
%         signature - the trained signature
%         residuals - the residuals (dependant variable)
%   Outputs: rand_stats - the statistics of the regression repeats (R2)


rand_stats=zeros(numRands,1);
warning('off','stats:regress:RankDefDesignMat');
rng('shuffle');
% find the remaining not in the signature

features_str=cell(size(geneTissue,1),1);
for k=1:size(geneTissue,1)
    features_str{k}=[geneTissue{k,1},geneTissue{k,2}];
end

sigFeatures_str=cell(size(signature,1),1);
for k=1:size(signature,1)
    sigFeatures_str{k}=[signature{k,1},signature{k,2}];
end
top_features_idx=(ismember(features_str,sigFeatures_str));

if (sum(top_features_idx)==0)
    return
end

test_data=exprData(:,top_features_idx);
nonconst_idx=std(test_data)>eps;
test_data=zscore(test_data(:,nonconst_idx));
for r=1:numRands
    p=randperm(size(exprData,1));
    rand_residuals=residuals(p);

    [~,~,~,~,stats]=regress(rand_residuals,[ones(size(residuals,1),1),test_data]);
    rand_stats(r,1)=stats(1);
end

warning('on','stats:regress:RankDefDesignMat');
