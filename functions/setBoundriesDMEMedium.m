function model = setBoundriesDMEMedium(model)

modelExchanges = findExchangeReactions(model);

model.lb(ismember(model.rxns,model.rxns(modelExchanges)))=0;


% Open just the DMEM's exchange reactions
model = changeExBound(model, 'gly', -0.030, modelExchanges);
model = changeExBound(model, 'gln_L', -0.584, modelExchanges);
model = changeExBound(model, 'ile_L', -0.105, modelExchanges);
model = changeExBound(model, 'leu_L', -0.105, modelExchanges);
model = changeExBound(model, 'met_L', -0.030, modelExchanges);
model = changeExBound(model, 'phe_L', -0.066, modelExchanges);
model = changeExBound(model, 'ser_L', -0.042, modelExchanges);
model = changeExBound(model, 'thr_L', -0.095, modelExchanges);
% model = changeExBound(model, 'trp_L', -0.016, modelExchanges);
model = changeExBound(model, 'val_L', -0.094, modelExchanges);
% this is a derivate of Folic acid which is the original metabolite in the DMEM
model = changeExBound(model, 'dhf', -0.004, modelExchanges);
model = changeExBound(model, 'nad', -0.004, modelExchanges);
% model = changeExBound(model, 'ribfiv', -0.0004, modelExchanges);
model = changeExBound(model, 'inost', -0.0072, modelExchanges);
model = changeExBound(model, 'glc_D', -4.5, modelExchanges);

% Open the FBS' exchange realctions
model = changeExBound(model, 'inost', -69.9 * 10^-3, modelExchanges);
model = changeExBound(model, 'thymd', -1.5 * 10^-3, modelExchanges);
model = changeExBound(model, 'ribflv', -581.9 * 10^-6, modelExchanges);
model = changeExBound(model, 'pydx', -150.7 * 10^-6, modelExchanges);
model = changeExBound(model, 'cys_L', -100.2 * 10^-3, modelExchanges);
model = changeExBound(model, 'asn_L', -56.8 * 10^-3, modelExchanges);
model = changeExBound(model, 'trp_L', -44.2 * 10^-3, modelExchanges);
model = changeExBound(model, 'ser_L', -249.8 * 10^-3, modelExchanges);
model = changeExBound(model, 'ala_L', -50 * 10^-3, modelExchanges);
model = changeExBound(model, 'arg_L', -846.7 * 10^-3, modelExchanges);
model = changeExBound(model, 'val_L', -451.3 * 10^-3, modelExchanges);
model = changeExBound(model, 'ile_L', -415.2 * 10^-3, modelExchanges);
model = changeExBound(model, 'leu_L', -450.1 * 10^-3, modelExchanges);
model = changeExBound(model, 'lys_L', -624.1 * 10^-3, modelExchanges);
model = changeExBound(model, 'pro_L', -149.9 * 10^-3, modelExchanges);
model = changeExBound(model, 'thr_L', -448.8 * 10^-3, modelExchanges);
model = changeExBound(model, 'Lcystin', -130.2 * 10^-3, modelExchanges);
model = changeExBound(model, 'gln_L', -2.5, modelExchanges);
model = changeExBound(model, 'hxan', -15.4 * 10^-3, modelExchanges);
model = changeExBound(model, 'lipoate', -508.9 * 10^-3, modelExchanges);

%% FROM: https://www.researchgate.net/figure/Lipid-composition-of-FBS-and-FBSLess-mM-sera_fig13_49820996
% FBS' Triglycerides
model = changeExBound(model, 'tag_hs', -1.26, modelExchanges);
% FBS' phospholipids
% model = changeExBound(model, 'pe_hs', -0.72, modelExchanges);
model = changeExBound(model, 'ps_hs', -0.72, modelExchanges);
% FBS' cholesterol
model = changeExBound(model, 'chsterol', -1.1, modelExchanges); % HDL

%% Exchange gases
% O2, CO2, H2, NH4 and H2S b
model = changeExBound(model, 'o2', -10, modelExchanges);
model = changeExBound(model, 'no', -10, modelExchanges);

%% SOURCE: https://www.researchgate.net/publication/15827024_Purine_base_and_nucleoside_cytidine_and_uridine_concentrations_in_foetal_calf_and_other_sera
model = changeExBound(model, 'cytd', -1.3 * 10^-3, modelExchanges);        % Cytidine
model = changeExBound(model, 'uri', -5.1 * 10^-3, modelExchanges);        % Uridine
model = changeExBound(model, 'urate', -130 * 10^-3, modelExchanges);      % Urate
model = changeExBound(model, 'hxan', -74.7 * 10^-3, modelExchanges);      % Hypoxanthine


model = changeExBound(model, 'hmcr', -30.1 * 10^-3, modelExchanges);        % Homocitrulline

model = changeRxnBounds(model, 'EX_hmcr[e]'               , -0.0301, 'l');
model = changeRxnBounds(model, 'sink_chol[c]'             , -0.056811, 'l');
model = changeRxnBounds(model, 'sink_glygn2[c]'           , -0.52, 'l');
model = changeRxnBounds(model, 'sink_fad[c]'              , -0.5, 'l');


model = changeRxnBounds(model, 'EX_co2[e]', 0.530, 'u');
model = changeRxnBounds(model, 'EX_co2[e]', 0.515, 'l');


end


function model = changeExBound(model, metName, value, modelExchanges)
% Search for reactions containing the metabolite
metList = model.mets(strcmp(model.metReconMap, metName));
rxnsWithMet = findRxnsFromMets(model, metList);

% See which exchanges have the metabolite and chage de boundrie
rxnList = intersect(model.rxns(modelExchanges), rxnsWithMet);
model = changeRxnBounds(model, rxnList, value, 'l');

end