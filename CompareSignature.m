function stats = CompareSignature( dirName,isAAorEUR, signature )
%CalcStats Calculate the statistics of using the signatures
%   computes for each of the signatures and each of the datasets
% Input: dirName - the name of the directory of the training files
%        isAAorEUR - train an AA signature (=0) or EUR (=1) 
%        signature - the signature learned on the training set
% Output: stats - the statistics of the signature, with R2, p-value
% compared to random signatures and p-value compared to shuffled signatures

if (isAAorEUR==0)
    fileName=[dirName,'/AA_validation.txt'];
elseif (isAAorEUR==1)
    fileName=[dirName,'/EUR_validation.txt'];
else
    error('isAAorEUR can take the values 0 or 1 only'); 
end

[exprData,residuals,geneTissue]=LoadData(fileName);

NUM_RAND=10000;

stats=zeros(1,3);
rand_stats=zeros(NUM_RAND,2);

% stats
stats(1,1) = RegressSig(exprData,geneTissue,signature,residuals);
% random signatures statistics
rand_stats(:,1) = RunRandSignatures(exprData,geneTissue,signature,residuals,NUM_RAND);
% shuffled signatures statistics
rand_stats(:,2) = RunShuffledSig(exprData,geneTissue,signature,residuals,NUM_RAND);

stats(1,2)=sum(stats(1,1)<rand_stats(:,1))/NUM_RAND;
stats(1,3)=sum(stats(1,1)<rand_stats(:,2))/NUM_RAND;
    
end

