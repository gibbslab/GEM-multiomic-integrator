function result = mapAbundanceToReactions(model, abundance)
%MAPABUNDANCETOREACTIONS This function maps the protein abundances to 
% the reactions in the model.


%% 1. Get the list of ProteinECNumbers
% Maybe we can check if the proteins are already downloaded to not do that
% again
if isfile('ProteinECNumbers.mat')
    load('ProteinECNumbers.mat');
else
    ProteinECNumbers = getProteinCodes(abundance);
    save('ProteinECNumbers.mat', 'ProteinECNumbers');
end

%% 2. Take rxnECNumbers and format correctly 
% .rxnECNumbers ; = & or = | .- = all in the subgroup
Recon3D.rxnECNumbers = strrep(Recon3D.rxnECNumbers,";"," & ");
Recon3D.rxnECNumbers = strrep(Recon3D.rxnECNumbers,"or","|");
Recon3D.rxnECNumbers = strrep(Recon3D.rxnECNumbers,".-","all");

%% 3. Find matches between model and abundance
ECNumbers = regexp(model.rxnECNumbers, '\<(?!EC:|^\>)([0-9.\-]*)','match');
ECNumbersLogic = regexp(model.rxnECNumbers, '( or | ;)','match');

rxn_n = numel(ECNumbers);
result = cell (rxn_n, 1);
for i=1:rxn_n
    ReactionAbundance = 0;
    for j=1:numel(ECNumbers{i})
        % look for the EC number in the proteome
        matches = cellfun(@(x) strcmp(x, ECNumbers{i}{j}), ProteinECNumbers, 'UniformOutput', false);
        FirstIntex = cellfun(@(c) any(c(:)), matches);
        IDs = {abundance.ID{FirstIntex}};
        Abundances = str2double({abundance.Control_NHA1_veh_tech2{FirstIntex}});
        [AbundanceSumary, index] = max(Abundances);
        AbundaceIndex = FirstIntex(index);
        
%% 4 and = min,  or = max
        % Now we have the protein abundance
        % We need to check wether the logic is and/or to use min/max
        ReactionAbundance = max(ReactionAbundance, AbundanceSumary);
    end
    result{i} =  ReactionAbundance;
end

end



function ProteinECNumbers = getProteinCodes (abundance)
import matlab.net.*
import matlab.net.http.*

tic;
IDs = abundance.ID;
IDs_length = numel(abundance.ID);
ProteinECNumbers = cell(IDs_length, 1);

f = waitbar(0,'Please wait...');

parfor i=1:IDs_length
    UniProtURL = ['https://www.uniprot.org/uniprot/' IDs{i} '.txt'];
    txt = webread(UniProtURL); 
    
    % Get the Enzyme Commission numbers
    ECNumbers = regexp(txt, 'DE[ ]*EC=(?<ECNumbers>([0-9]+.([0-9]+.)*([0-9]|-)+))','names');
    % Get the Transporter Classification numbers 
    TCDBNumbers = regexp(txt, 'DR[ ]*TCDB; (?<TCDBNumbers>([0-9A-Z]+.([0-9A-Z]+.)*([0-9A-Z]|-)+))','names');
    

   ProteinECNumbers{i} = [
       reshape(struct2cell(ECNumbers), numel(ECNumbers), 1)
       reshape(struct2cell(TCDBNumbers), numel(TCDBNumbers), 1)
       ];
    
    waitbar(i/IDs_length, f,['Downloading (' num2str(i) '/' num2str(IDs_length) ')...']);
end
close(f)
toc;


end



