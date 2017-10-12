function [exprData,residuals,geneTissue]=LoadData(fileName)
% LoadData loads the relavant files. 
% Input: fileName - the name of the file to load
% Output: exprData - the imputed expression data matrix
%         residuals - the residuals (dependant variable)
%         geneTissue - the names of the gene-tissue names (feature names)

fid=fopen(fileName,'r');
tissues = textscan(fid,'%s',540,'Delimiter','\t');
geneNames = textscan(fid,'%s',540,'Delimiter','\t');
geneTissue=[tissues{1,1}(2:end),geneNames{1,1}(2:end)];

fclose(fid);

exprData=dlmread(fileName,'\t',2,0);
residuals=exprData(:,1);
exprData=exprData(:,2:end);

