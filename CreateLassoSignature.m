function [allSigsGenes,featuresSets]=CreateLassoSignature(data,phen,featureNames,SUBSETS)
% CreateLassoSignature uses lasso with preconditioning and bootstrapping to create a
% signature. 
% Input: data - an MxN matrix, where M is the number of sampels and N is
% the number of features
%        phen - the dependant variable of size M 
%        featureNames - the names of the features (N-sized cell array)
%        SUBSETS - the optional number of subsets to run. default is 100
% Output: allSigsGenes - a cell array with the union of all the signaure
% genes. Is used to produce the signature
%         featuresSets - a cell array with SUBSETS number of signatures


if nargin<4
    SUBSETS=100;
end

s = RandStream('mrg32k3a','Seed','shuffle');
RandStream.setGlobalStream(s);

% remove features with zero variance
nz_var_features=(var(data)>eps);
data=data(:,nz_var_features);
featureNames=featureNames(nz_var_features,:);

featuresSets=cell(SUBSETS,1);
signatures_minMSE=[];
tic;
opt=struct('UseParallel','true','UseSubstreams','true','Streams',s);
for subset=1:SUBSETS
    if (mod(subset,1)==0)
        fprintf('%ld,',subset);
    end
    [B,lasso_fit]=lasso(data,phen,'NumLambda',100,'Alpha',1,'CV',5,'Options',opt);% 'Options',opt=struct('UseParallel','false','UseSubstreams','false');
    indexMinMSE=lasso_fit.IndexMinMSE;
    
    top_features_idx=(B(:,indexMinMSE)~=0);
    selected_features=featureNames(top_features_idx,:);
    if (~isempty(selected_features))
        featuresSets{subset,1} =  selected_features;
        signatures_minMSE=[signatures_minMSE;selected_features];
    end
end

toc;
signatures_minMSE=uniqueRowsCA(signatures_minMSE,'rows');
allSigsGenes=cell(1,1);
allSigsGenes{1}=signatures_minMSE;

end

function [signature,feature_appears]=SelectSignature(feature_sets,features,index,cutoff)
feature_appears=zeros(size(features,1),1);
features_str=cell(size(features,1),1);
for k=1:size(features,1)
    features_str{k}=[features{k,1},features{k,2}];
end
for subset=1:size(feature_sets,1)
    feature_set=feature_sets{subset,index};
    if (~isempty(feature_set))
        feature_set_str=cell(size(feature_set,1),1);
        for k=1:size(feature_set,1)
            feature_set_str{k}=[feature_set{k,1},feature_set{k,2}];
        end
        idxFeature=ismember(features_str,feature_set_str);
        
        feature_appears(idxFeature)=feature_appears(idxFeature)+1;
    end
end
feature_appears=feature_appears/size(feature_sets,1);
selected_features=(feature_appears>=cutoff);
signature=features(selected_features,:);
feature_appears=feature_appears(selected_features,:);
end




