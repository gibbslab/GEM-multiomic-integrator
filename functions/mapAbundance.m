function abundanceTable = mapAbundance (base_model, abundance, proteinECNumbers)

abundanceTable = table();
abundance = addvars(abundance, ProteinECNumbers, 'After', 1);
abundance.protein = [];
ID = abundance(:, 2);
column_names = abundance.Properties.VariableNames;
parfor k=3:length(column_names)
    value = abundance(:, k);
    abundanceToMap =  [ID, value];
    abundanceToMap.Properties.VariableNames = {'id' 'value'};
    
    % Remove missing values
    abundanceToMap = rmmissing(abundanceToMap);
    abundanceToMap = abundanceToMap(~cellfun('isempty', abundanceToMap.id), :);
    
    [abundanceRxns, parsedPR] = mapAbundanceToReactions(base_model, abundanceToMap);
    abundanceRxns_table = table(abundanceRxns, 'VariableNames', column_names(k));
    abundanceTable = [abundanceTable, abundanceRxns_table];
end

end