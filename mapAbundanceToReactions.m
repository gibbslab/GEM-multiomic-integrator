function ProteinECNumbers = mapAbundanceToReactions(model, abundance)
%MAPABUNDANCETOREACTIONS This function maps the protein abundances to 
% the reactions in the model.


%% 1. Get the list of ProteinECNumbers
ProteinECNumbers = getProteinCodes(abundance);




%% 3. Tomar rxnECNumbers formatear bien 
% .rxnECNumbers ; = & or = | .- = all in the subgroup
% añadir parentesis al rededor de cada uno x(number)


%% 4 and = min,  or = max


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



