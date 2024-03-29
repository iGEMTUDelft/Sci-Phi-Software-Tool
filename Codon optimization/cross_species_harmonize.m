function cross_species_harmonized = cross_species_harmonize(input,RFC,TypeIIs, fileName)
    fileName = char(fileName);
    fileName = [fileName '.txt'];
    startpos = 1;

    input_CDS = upper(input(startpos:end));
    cross_species_harmonized = transform_seq(input, startpos);
    fid=fopen(fileName,'w');
    fprintf(fid, cross_species_harmonized);
    fclose(fid);
    if nt2aa(input_CDS,'AlternativeStartCodons', false) == nt2aa(cross_species_harmonized,'AlternativeStartCodons', false)
        disp("Optimized")
    end

    input_codon = [];
    output_codon = [];
    for i = 1:length(input_CDS)
        if mod(i,3) == 0
            input_codon = [input_codon upper(convertCharsToStrings([input_CDS(i-2) input_CDS(i-1) input_CDS(i)]))];
            output_codon = [output_codon convertCharsToStrings([cross_species_harmonized(i-2) cross_species_harmonized(i-1) cross_species_harmonized(i)])];
        end
    end
    input_CDS = upper(input_CDS);
    cross_species_harmonized
    disp(["Codons changed: " length(input_codon) - sum(input_codon == output_codon) " Percentage changed: " 100 - sum(input_codon == output_codon)/length(input_codon)*100]);
    disp(["Nucleotides changed: " length(input_CDS) - sum(input_CDS == cross_species_harmonized) " Percentage changed: " 100 - sum(input_CDS == cross_species_harmonized)/length(input_CDS)*100]);

    function char_seq_codon = transform_seq(input, startpos)

        %get CDS from full mRNA
        input_CDS = input(startpos:end);
        CDS_codons = [];
        for i = 1:length(input_CDS)
            if mod(i,3) == 0
                CDS_codons = ([CDS_codons convertCharsToStrings([input_CDS(i-2) input_CDS(i-1) input_CDS(i)])]);
            end
        end

        %process the data
        [Data_raw_codons, AA_ref, Codon_ref, AA_list, n_orgs, restriction_data] = process_data();

        %make a list of preferred codons
        [variances AA_ref_new, Codon_ref_new] = get_variance(Data_raw_codons, AA_ref, AA_list, n_orgs, Codon_ref);
        [ordered_codons, ordered_variances] = get_ordered_list(variances, AA_ref, AA_list, Codon_ref, Codon_ref_new);
        Codon_ref_new
        %make initial optimized sequence
        seq_codon = [];
        for i = 1:length(CDS_codons)
            codon = CDS_codons(i);
            %get the ordered codon list for this aminoacid
            syn_codon = ordered_codons(find(Codon_ref_new(1,:) == upper(codon)));
            codon = syn_codon(1); %grab the first one

            seq_codon = [seq_codon codon];
        end

        char_seq_codon = char(strjoin(seq_codon));
        char_seq_codon = char_seq_codon(~isspace(char_seq_codon))

        disp(char_seq_codon);

        %now we will iterate over restriction sites and grab different codons
        char_seq_codon = eliminate_restrictions(seq_codon, char_seq_codon, Codon_ref_new, ordered_codons, ordered_variances, restriction_data);
    end

    function [Data_raw_codons, AA_ref, Codon_ref, AA_list, n_orgs, restriction_data] = process_data()
        %import data
        [~,~,Data_raw_codons]=xlsread([pwd '/data_formatted.xlsx']);
        Data_raw_codons = Data_raw_codons(2:end,:);
        [~,~,restriction_data_raw]=xlsread([pwd '/restriction_enzyme_database.xlsx']);
        if RFC == 1 && TypeIIs == 0
            fileID = fopen('RFC.txt');
        elseif RFC == 0 && TypeIIs == 1
            fileID = fopen('TypeIIs.txt');
        else
            fileID = fopen('both.txt');
        end
        textfile = textscan(fileID,'%q');
        fclose(fileID);

        textfile = strsplit(string(textfile), ',');

        Restriction_Ref = string(restriction_data_raw(2:end,1));
        restriction_data = [];
        for i = 1:length(textfile)
            restriction_data = [restriction_data restriction_data_raw(find(Restriction_Ref == textfile(i)) + 1, 3)];
        end

        restriction_data = [textfile' restriction_data'];

        %make reference lists
        AA_list = ["phe" "leu" "ile" "met" "val" "tyr" "*" "his" "gln" "asn" "lys" "asp" "glu" "ser" "pro" "thr" "ala" "cys" "trp" "arg" "gly"];
        AA_raw = lower(string(Data_raw_codons(1,2:end)));
        AA_raw(find(AA_raw == "end")) = "*";
        Codon_raw = string(Data_raw_codons(2,2:end));

        AA_ref = string([1,length(AA_raw)]);
        Codon_ref = string([1,length(AA_raw)]);
        startpos = 1;
        for i = 1:length(AA_list)
            pos = find(AA_raw == AA_list(i));
            leng = length(pos);
            AA_ref(startpos:startpos+leng-1) = AA_list(i);
            Codon_ref(startpos:startpos+leng-1) = Codon_raw(pos);
            startpos = startpos + leng;
        end
        n_orgs = length(Data_raw_codons(3:end,1));
    end

    function [variances, AA_ref_new, Codon_ref_new] = get_variance(Data_raw_codons, AA_ref, AA_list, n_orgs, Codon_ref)
        %get the data for frequency
        freqs = zeros(n_orgs, length(AA_ref));

        for i = 1:n_orgs
            data = Data_raw_codons(2+i,2:end);
            freqs(i,:) = str2double(data);
        end

        %calculate percentage
        percentage = [];
        for i = 1:length(AA_list)
           aa = AA_list(i);
           pos = find(AA_ref == aa);
           aa_freq = freqs(:,pos);
           sum_rows_aa= sum(aa_freq, 2);

           freq_aa = aa_freq;
           N = length(sum_rows_aa);     
           for ii = 1:N
               freq_aa(ii,:)=freq_aa(ii,:)./sum_rows_aa(ii,1);
           end
           percentage = [percentage freq_aa];
        end

        %find the length of the new reference lists
        total_length = 0;
        for i = 1:length(AA_list)
           leng = length(find(AA_ref == AA_list(i)));
           leng = leng*leng;
           total_length = total_length + leng;
        end

        %define variances and new reference lists
        variances = zeros(1, total_length);
        AA_ref_new = strings([1,total_length]);
        Codon_ref_new = strings([2,total_length]);

        %make the new reference lists
        startpos_AA = 1;
        startpos_codons_r1 = 1;
        c = 1;
        for i = 1:length(AA_list)
            %new AA ref list
            leng = length(find(AA_ref == AA_list(i)));
            leng = leng*leng;
            AA_ref_new(startpos_AA:startpos_AA+leng-1) = AA_list(i);

            %new codon ref list
            codon = Codon_ref(c);
            leng_codons = leng/length(find(AA_ref == AA_list(i))); %divide up the length in the segments for each codon
            for ii = 1:leng
                Codon_ref_new(1,startpos_codons_r1:startpos_codons_r1+leng_codons-1) = codon;
                Codon_ref_new(2,startpos_codons_r1:startpos_codons_r1+length(find(AA_ref == AA_list(i)))-1) = Codon_ref(find(AA_ref == AA_list(i)));
                if mod(ii,leng_codons) == 0
                    c = c + 1;
                    if c <= length(Codon_ref)
                        codon = Codon_ref(c);
                    end
                    startpos_codons_r1 = startpos_codons_r1 + leng_codons;
                end
            end
            startpos_AA = startpos_AA+leng;  
        end

        %calculate the variances
        startpos_AA = 0;
        for i = 1:length(AA_list)
           pos = find(AA_ref == AA_list(i));
           freqs = percentage(:,pos);
           [m, n] = size(freqs);
           for j = 1:n
              for k = 1:n
                 column = [freqs(1,j); freqs(2:end,k)];
                 variances(1,startpos_AA+k) = var(column);        
              end
              startpos_AA = startpos_AA + n;
           end
        end    
    end

    function [ordered_codons, ordered_variances] = get_ordered_list(variances, ~, ~, Codon_ref, Codon_ref_new)
        ordered_codons = []; %string(size(Codon_ref_new));
        ordered_variances = [];  %string(size(variance));

        for i = 1:length(Codon_ref)
            %get the positions of this codon
            codon = Codon_ref(i);
            if codon == 'TTG'
                disp(find(Codon_ref_new(1,:) == codon))
            end
            codon_pos = find(Codon_ref_new(1,:) == codon);

            %get the list of synonymous codons and the variances
            syn_codons = Codon_ref_new(2,codon_pos);
            syn_variance = variances(codon_pos);

            %sort the variances
            sorted_list = sort(syn_variance);

            %use the ordered variance list to order the codons
            for ii = 1:length(syn_codons)
                pos = find(syn_variance  == sorted_list(ii));
                codon = syn_codons(pos);
                ordered_codons = [ordered_codons codon];
            end
            ordered_variances = [ordered_variances sorted_list];
        end
    end

    function output = eliminate_restrictions(seq_codon, char_seq_codon, Codon_ref_new, ordered_codons, ordered_variances, restrictions)
        %now we will iterate over restriction sites and grab different codons
        empty = false;
        restriction_pos_old = [];
        while ~empty
            %see if its truly empty of restriction sites, this is the best way my tired mind could
            %think of for now..
            c = 0;
            all_pos = [];
            for j = 1:length(restrictions(:,1))
                [~, restriction_pos] = rebasecuts(char_seq_codon, {char(restrictions(j, 1))});

                disp(char(restrictions(j, 1)))
                disp(restriction_pos)
                if ~isempty(restriction_pos)
                    c = c+1;
                    all_pos = [all_pos restriction_pos];                
                end
            end
            if ~isequal(restriction_pos_old,all_pos) 
                restriction_pos_old = all_pos;
            else
                if ~isempty(restriction_pos_old)
                    disp('This sequence could not be solved.');
                     empty = true;
                end              
            end
            
            if c == 0
                empty = true;
            end
            
            for i = 1:length(restrictions(:,1))
                %find the restriction site
                i
                [~, restriction_pos] = rebasecuts(char_seq_codon, {char(restrictions(i, 1))});
                cant_fix_restr_pos = false(1, length(restriction_pos));
                if ~isempty(restriction_pos)

                    %determine the length of the recognition site left and right of
                    %the recognition site the rebasecuts gives
                    recognition_seq = char(restrictions(i, 2));
                    length_right = length(recognition_seq) - find(recognition_seq == '''');
                    length_left = length(recognition_seq) - length_right - 1;

                    for ii = 1:length(restriction_pos)
                        if ~cant_fix_restr_pos(ii)
                            if mod(restriction_pos(ii), 3) == 0 
                                %end of codon
                                if mod(length_right, 3) == 0
                                    span_right = length_right/3;
                                else
                                    span_right = fix(length_right/3) + 1;
                                end 
                                if length_left > 2
                                    length_left = length_left - 2;
                                    if mod(length_left, 3) == 0
                                        span_left = length_left/3;
                                    else
                                        span_left = fix(length_left/3) + 1;
                                    end
                                else
                                    span_left = 0;
                                end
                                codon_start = restriction_pos(ii)/3 - span_left;
                            elseif mod(restriction_pos(ii), 3) == 2 
                                %middle of codon 
                                if length_left > 1
                                    length_left = length_left - 1;
                                    if mod(length_left, 3) == 0
                                        span_left = length_left/3;
                                    else
                                        span_left = fix(length_left/3) + 1;
                                    end 
                                else
                                    span_left = 0;
                                end
                                if length_right > 1
                                    length_right = length_right - 1;
                                    if mod(length_right, 3) == 0
                                        span_right = length_right/3;
                                    else
                                        span_right = fix(length_right/3) + 1;
                                    end
                                else
                                    span_right = 0;
                                end
                                codon_start = restriction_pos(ii)/3 + 1/3 - span_left;
                            elseif mod(restriction_pos(ii), 3) == 1 
                                %beginning of codon
                                if mod(length_left, 3) == 0
                                    span_left = length_left/3;
                                else
                                    span_left = fix(length_left/3) + 1;
                                end 

                                if length_right > 2
                                    length_right = length_right - 2;
                                    if mod(length_right, 3) == 0
                                        span_right = length_right/3;
                                    else
                                        span_right = fix(length_right/3) + 1;
                                    end
                                else
                                    span_right = 0;
                                end
                                codon_start = restriction_pos(ii)/3 - 1/3 - span_left;
                            end

                            %get all codons
                            total_span = span_right + span_left;
                            codon_start = uint8(codon_start);
                            total_span = uint8(total_span);

                            codons = seq_codon(codon_start:codon_start+total_span);

                            %define the next best choices
                            choices = 2.*ones(1,length(codons));
                            consider = true(1,length(codons));

                            choices(find(codons == 'ATG')) = 1;

                            fixed = false;
                            while ~fixed    
                                %get the codon we need to change and what to change it to
                                [alt_codon, pos, cant_fix] = get_alt_codon(codons, Codon_ref_new,  ordered_codons, ordered_variances, choices, consider);

                                seq_codon(codon_start+pos-1) = alt_codon;
                                char_seq_codon = char(strjoin(seq_codon));
                                char_seq_codon = char_seq_codon(~isspace(char_seq_codon));
                                disp(i)
                                [~, check] = rebasecuts(char_seq_codon, {char(restrictions(i, 1))});

                                if isempty(find(check == restriction_pos(ii)))
                                    %this site is now gone
                                    fixed = true;
                                else
                                    if cant_fix == true
                                        cant_fix_restr_pos = true;
                                        disp('Couldn''t fix it')
                                    else                                                         
                                        %not gone so change the choice of pos
                                        consider(pos) = false; %(pos) = choices(pos) + 1;
                                        if sum(consider) == 0
                                            %all of them have been considered so
                                            %now go to next choice
                                            if choices(pos) >= length(codons)
                                                fixed = true; %can't fix it so let go of loop
                                            else
                                                choices = choices + 1;
                                                choices(find(codons == 'ATG')) = 1;
                                                consider = true(1,length(codons));
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end  
            end
        end
        char_seq_codon = char(strjoin(seq_codon));
        char_seq_codon = char_seq_codon(~isspace(char_seq_codon));
        output = char_seq_codon;
    end

    function [alt_codon, pos, cant_fix] = get_alt_codon(codons, Codon_ref_new, ordered_codons, ordered_variances, choices, consider)
        %make list of delta_variance
        counter = 0;
        delta_variances = zeros(1, length(codons(1,:)));
        for j = 1:length(delta_variances)
            if consider(j)
                codon = codons(j);
                codon_pos = find(Codon_ref_new(1,:) == codon);

                codon_variances = ordered_variances(codon_pos);

                %only calculate difference if still possible otherwise just set it
                %at super high so it won't be considered for best option      
                if length(codon_variances) > choices(j) 
                    delta_variances(j) = codon_variances(choices(j) - 1) - codon_variances(choices(j));
                else
                    delta_variances(j) = 10^9;
                    choices(j) = length(codon_variances);
                end
            else
                %just make a super large value so it won't be considered
                delta_variances(j) = 10^10;
            end
        end

        %if all of them have been put up on way too high then its not fixable
        if counter == length(delta_variances)
            cant_fix = true;
        else
            cant_fix = false;
        end

        %get the one with the smallest difference
        pos = find(delta_variances == min(delta_variances));

        if length(pos) > 1
            %get the first one if there are more
            pos = pos(1);
        end
        codon_to_change = codons(pos);
        syn_codons_pos = find(Codon_ref_new(1,:) == codon_to_change);
        alt_codons = ordered_codons(syn_codons_pos);
        alt_codon = alt_codons(choices(pos));
    end
end