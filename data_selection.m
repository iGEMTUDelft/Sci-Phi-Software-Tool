close all 
clear 
clc
%%
tic
[~,~,Data_raw]=xlsread([pwd '/o576245-Refseq_species.xlsx']);
[~,~,Organisms]=xlsread([pwd '/organisms.xlsx']);

%%
Taxids = Data_raw(2:end,3);
Taxids= cell2mat(Taxids);
Organisms = cell2mat(Organisms);

codon_data_per_organism = [];
for i = 1:length(Organisms)
    pos = find(Taxids == Organisms(i));
    
    %get average if multiple
    if length(pos) > 1
        all = cell2mat(Data_raw(pos+2,13:end));
        average = sum(all,1)./length(pos);
        codon_data_per_organism = [codon_data_per_organism; average]
    else
        codon_data_per_organism = [codon_data_per_organism; cell2mat(Data_raw(pos+2,13:end))]
    end 
end

%%
%make top row for referencing 
top_row = [];
for i = 1:76
    if i >= 13
        top_row = [top_row aminolookup(nt2aa(Data_raw(1,i)))];
    end
end
top_row = cellstr(top_row);

second_row = cellstr(Data_raw(1,13:end));
second_row = ["Taxid" second_row];

combined = [" " top_row; second_row; Organisms codon_data_per_organism];
combined = cellstr(combined);

combined=cell2table(combined);
writetable(combined,'data_formatted.xlsx')
toc