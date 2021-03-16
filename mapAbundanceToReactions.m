function ProteinECNumbers = mapAbundanceToReactions(model, abundance)
%MAPABUNDANCETOREACTIONS This function maps the protein abundances to 
% the reactions in the model.


%% 1. Get the list of ProteinECNumbers
ProteinECNumbers = getProteinECNumbers(abundance);




%% 3. Tomar rxnECNumbers formatear bien 
% .rxnECNumbers ; = & or = | .- = all in the subgroup
% añadir parentesis al rededor de cada uno x(number)


%% 4 and = min,  or = max


end



function ProteinECNumbers = getProteinECNumbers (abundance)
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
    
    UniProtData = regexp(txt, 'DE[ ]*EC=(?<ECNumbers>([0-9]+.([0-9]+.)*([0-9]|-)+))','names');
    
    if not(isempty(UniProtData))
       ProteinECNumbers{i} = reshape(struct2cell(UniProtData), numel(UniProtData), 1);
    else
       ProteinECNumbers{i} = {''};
    end
    
    waitbar(i/IDs_length, f,'Downloading ...');
end
close(f)
toc;


end



