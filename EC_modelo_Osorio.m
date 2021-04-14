modelo=readCbModel('Astrocyte_Osorio2019.xml');
reglas=modelo.rules;
ECnumbers_osorio=readtable('ECnumbers_Osorio');

Var_rem=table2array(ECnumbers_osorio(:,'ECNumber'));


for i=1:2405
    Secuencia=strcat(['x(',num2str(i),')'],'');
    modelo.rules =  strrep(modelo.rules, Secuencia,Var_rem(i));
end


%remp = '( )|';
%modelo.rules =  strrep(modelo.rules, remp ,'');
%modelo.rules =  strrep(modelo.rules, '\)','');
