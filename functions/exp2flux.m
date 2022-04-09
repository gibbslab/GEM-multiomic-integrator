function model = exp2flux(model, expressionRxns)
    % Save boundries
    lb = model.lb;
    ub = model.ub;
    % Find exchanges
    modelExchanges = findExchangeReactions(model);


    % Modify boundries acording to expression
    model.lb = -1 * expressionRxns;
    model.ub = expressionRxns;

    model.lb(lb >= 0) = lb(lb >= 0); % keep irreversible reactions
    % Keep boundries in exchange reactions
    model.lb(modelExchanges) = lb(modelExchanges);
    model.ub(modelExchanges) = ub(modelExchanges);
end