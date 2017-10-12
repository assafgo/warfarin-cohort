function [ regression_stats] = RegressSig( exprData,geneTissue,signature,residuals )
% RegressSig - Run dose resgression based on the signature
%   Inputs: exprData - the imputed expression data matrix
%         geneTissue - the names of the gene-tissue names (feature names)
%         signature - the trained signature
%         residuals - the residuals (dependant variable)
%   Outputs: regression_stats - the statistics of the regression (R2)

features_str=cell(size(geneTissue,1),1);
for k=1:size(geneTissue,1)
    features_str{k}=[geneTissue{k,1},geneTissue{k,2}];
end

sigFeatures_str=cell(size(signature,1),1);
for k=1:size(signature,1)
    sigFeatures_str{k}=[signature{k,1},signature{k,2}];
end
top_features_idx=(ismember(features_str,sigFeatures_str));

selectedData=exprData(:,top_features_idx);
nonconst_idx=std(selectedData)>eps;
selectedData=zscore(selectedData(:,nonconst_idx));

[~,~,~,~,stats]=regress(residuals,[ones(size(residuals,1),1),selectedData]);
regression_stats=stats(1);

end

