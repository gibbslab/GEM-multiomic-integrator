function expressionTable = mapExpression(base_model, expresion)

Recon3D.genes = regexprep(Recon3D.genes,"\.[0-9]*","");

expresion_columns = expresion.Properties.VariableNames;
expressionRxns1 = cell(length(expresion_columns) - 2, 1);
expressionTable = table();
for k=2:length(expresion_columns)
    expressionToMap =  expresion(:, [1 k]);
    expressionToMap.Properties.VariableNames = {'gene' 'value'};
    [expressionRxns, parsedGPR] = mapExpressionToReactions(base_model, expressionToMap);
    expressionTable = addvars(expressionTable, expressionRxns);
end
expressionTable.Properties.VariableNames = expresion.Properties.VariableNames(2:end);