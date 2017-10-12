function [signature]=SelectRobustSignature(signatureSets,selectedFeatures,cutoff)
% SelectRobustSignature selects a robust signature from multiple LASSO runs. 
% Input: signatureSets - the set of signatures returned from multiple runs
%        selectedFeatures - the selectedFeatures in any of the runs
%        cutoff - the minimal percentage of signatures a gene-tissue needs
%        to be in, in order to be considered robust
% Output: signature - a cell array with the robust signature

% create a set of gene-tissue strings

feature_appears=zeros(size(selectedFeatures,1),1);
features_str=cell(size(selectedFeatures,1),1);
for k=1:size(selectedFeatures,1)
    features_str{k}=[selectedFeatures{k,1},selectedFeatures{k,2}];
end

% count the number of appearances
for subset=1:size(signatureSets,1)
    feature_set=signatureSets{subset,1};
    if (~isempty(feature_set))
        feature_set_str=cell(size(feature_set,1),1);
        for k=1:size(feature_set,1)
            feature_set_str{k}=[feature_set{k,1},feature_set{k,2}];
        end
        idxFeature=ismember(features_str,feature_set_str);
        
        feature_appears(idxFeature)=feature_appears(idxFeature)+1;
    end
end
feature_appears=feature_appears/size(signatureSets,1);
selected_features=(feature_appears>=cutoff);
signature=selectedFeatures(selected_features,:);

end
