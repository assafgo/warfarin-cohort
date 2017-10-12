function signature=runCreateSignature (dirName,isAAorEUR)
% runCreateSignature loads the relavant files and creates a robust signature via LASSO. 
% Input: dirName - the name of the directory of the training files
%        isAAorEUR - train an AA signature (=0) or EUR (=1) 
% Output: signature - a cell array with the robust signature

if (isAAorEUR==0)
    fileName=[dirName,'/AA_train.txt'];
elseif (isAAorEUR==1)
    fileName=[dirName,'/EUR_train.txt'];
else
    error('isAAorEUR can take the values 0 or 1 only'); 
end

NUMBER_OF_SUBSETS=100;
SUBSET_CUTOFF = 0.5;

[exprData,residuals,geneTissue]=LoadData(fileName);

[allGenesSelected,features_sets]=CreateLassoSignature(exprData,residuals,geneTissue,NUMBER_OF_SUBSETS);

signature=SelectRobustSignature(features_sets,allGenesSelected{1,1},SUBSET_CUTOFF);

end
