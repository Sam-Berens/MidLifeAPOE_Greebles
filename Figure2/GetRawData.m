function [DataTable] = GetRawData()

warning('off','MATLAB:table:ModifiedAndSavedVarnames');
DataTable = readtable('..\Greebles_SummaryData.csv');
warning('on','MATLAB:table:ModifiedAndSavedVarnames');

DataTable.Genotype = categorical(DataTable.APOE_haplotype);

return