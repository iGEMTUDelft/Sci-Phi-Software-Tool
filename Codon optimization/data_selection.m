function selecting_data = data_selection()
    load_data = false;
    [~,~,Organisms]=xlsread([pwd '/organisms.xlsx']);
    Organisms = cell2mat(Organisms); 
    fileID = fopen('previous_selected_orgs.txt');
    textfile = textscan(fileID,'%q');
    fclose(fileID);
    textfile = strsplit(string(textfile), ',');
    previous_orgs = zeros(length(textfile),1);
    for i = 1:length(textfile)
        previous_orgs(i) = str2num(textfile(i));
    end
    
    Organisms_sort = sort(Organisms);
    previous_orgs = sort(previous_orgs);
    
     if length(previous_orgs) == length(Organisms_sort)
         if ~isequal(previous_orgs,Organisms_sort)
             load_data = true;
         end
     else
         load_data = true;
     end 
    
    if load_data
        [~,~,Data_raw]=xlsread([pwd '/o576245-Refseq_species.xlsx']);

        Taxids = Data_raw(2:end,3);
        Taxids= cell2mat(Taxids);


        codon_data_per_organism = [];
        for i = 1:length(Organisms)
            pos = find(Taxids == Organisms(i));

            %get average if multiple
            if length(pos) > 1
                all = cell2mat(Data_raw(pos+2,13:end));
                average = sum(all,1)./length(pos);
                codon_data_per_organism = [codon_data_per_organism; average];
            else
                codon_data_per_organism = [codon_data_per_organism; cell2mat(Data_raw(pos+2,13:end))];
            end 
        end

        top_row = [];
        for i = 1:76
            if i >= 13
                top_row = [top_row aminolookup(nt2aa(Data_raw(1,i), 'AlternativeStartCodons', false))];
            end
        end
        top_row = cellstr(top_row);

        second_row = cellstr(Data_raw(1,13:end));
        second_row = ["Taxid" second_row];

        combined = [" " top_row; second_row; Organisms codon_data_per_organism];
        combined = cellstr(combined);

        combined=cell2table(combined);
        writetable(combined,'data_formatted.xlsx')

        GC_data_per_organism = [];
        for i = 1:length(Organisms)
            pos = find(Taxids == Organisms(i));

            %get average if multiple
            if length(pos) > 1
                all = cell2mat(Data_raw(pos+2,9));
                average = sum(all,1)./length(pos);
                GC_data_per_organism = [GC_data_per_organism; average];
            else
                GC_data_per_organism = [GC_data_per_organism; cell2mat(Data_raw(pos+2,9))];
            end 
        end
        GC_data_per_organism = array2table(GC_data_per_organism);
        writetable(GC_data_per_organism,'GC_content.xlsx')
     % write the organisms to remember
     dlmwrite('previous_selected_orgs.txt',Organisms','precision',10);    
    end
end