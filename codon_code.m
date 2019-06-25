clc, clear all
tic
input = 'atgtctaccaccatgggaattcaacctcctaaaaagaaacgtaaagttaatactattaatattgctaaaaatgacttctcagatattgaattagcagccattccatttaatacattagcagatcactatggtgaacgtttagcacgtgaacagttagcattagaacatgaatcatatgaaatgggtgaagcacgttttcgtaagatgttcgagcgtcagttaaaagcaggtgaagttgcagataatgcagcagccaaacctttaattactacattattacctaaaatgattgctcgtattaacgattggtttgaagaggttaaagcaaagcgtggtaaacgtcctacagcatttcagttcttacaagaaatcaaacctgaagcagttgcatatattactattaaaacaacattagcatgtttaacatcagcagataatacaacagttcaagcagttgcatcagcaattggtcgtgcaattgaagatgaagcacgttttggtcgtattcgtgatttagaagccaaacattttaaaaaaaatgttgaagaacagttaaacaaacgtgttggtcatgtttataaaaaagcatttatgcaggttgttgaagcagatatgttatcaaaaggtttattaggtggtgaagcatggtcatcatggcataaagaagattcaattcatgttggtgttcgttgtattgaaatgttaattgaatcaacaggtatggtttcattacatcgtcagaatgcaggtgttgttggtcaagattcagaaacaattgaattagcacctgaatatgcagaagcaattgcaacacgtgcaggtgcattagcaggtatttcaccaatgtttcagccttgtgttgttcctcctaaaccttggacaggtattacaggtggtggttattgggcaaatggtcgtcgtcctttagcattagttcgtacacattcaaaaaaagcattaatgcgttatgaagatgtttacatgcctgaagtttataaagccattaatattgcacagaatactgcatggaaaatcaacaagaaagttttagcagttgcaaatgttattactaaatggaaacattgtcctgttgaagatattcctgcaattgaacgtgaagaattaccaatgaaacctgaagatattgatatgaatcctgaagcattaacagcatggaaacgtgcagcagcagccgtttatcgtaaagataaagcacgtaaatcacgtcgtatttcattagagttcatgttagagcaagcaaacaagttcgcaaatcataaagccatttggtttccttataatatggattggcgtggtcgtgtttatgcagtttcaatgtttaatcctcaaggtaatgatatgaccaaaggtttattaaccttagctaaaggtaaacctattggtaaagaaggttattattggttaaaaatccatggtgcaaattgtgcaggtgttgataaagttccttttccagaacgtattaagttcattgaagaaaatcatgaaaatattatggcatgtgctaaatcaccattagaaaatacatggtgggcagaacaagattcacctttttgttttttagccttttgttttgaatatgcaggtgttcaacatcatggtttatcatataattgttcattaccattagcatttgatggttcatgttcaggtattcagcatttttcagcaatgttacgtgatgaagttggtggtcgtgcagttaacttattaccttcagaaacagttcaagatatctatggtattgttgctaaaaaagttaatgaaatcttacaggcagatgccattaatggtacagataatgaagttgttacagttacagatgaaaatacaggtgaaatatcagaaaaagttaaattaggtaccaaagcattagcaggtcaatggttagcatatggtgttacacgttcagttactaagcgttcagttatgacattagcatatggttcaaaagagttcggttttcgtcaacaggttttagaagataccattcaacctgcaattgattcaggtaaaggtttaatgtttactcaacctaatcaagcagcaggttatatggccaaattaatatgggaatcagtttcagttacagttgttgcagcagttgaagcaatgaattggttaaaatcagcagccaaattattagcagcagaagttaaagataaaaaaacaggtgaaatcttacgtaagcgttgtgcagttcattgggttacacctgatggttttcctgtttggcaagaatataaaaaacctattcaaacacgtttaaacttaatgtttttaggtcagttccgtttacaacctactattaatacaaacaaagattcagaaattgatgcacataaacaagaatcaggtattgcacctaacttcgttcattcacaagatggttcacatttacgtaaaacagttgtttgggcacatgaaaaatatggtattgaatcatttgcattaattcacgattcatttggtaccattccagcagatgcagcaaacttattcaaagcagttcgtgaaacaatggttgatacatatgaatcatgtgatgttttagcagacttctatgatcagttcgcagatcagttacatgaatcacagttagataaaatgcctgcattacctgccaaaggtaacttaaacttacgtgatattttagaatcagacttcgcatttgcctaa';

startpos = 1;

input_CDS = upper(input(startpos:end));
output_CDS = transform_seq(input, startpos);
%propInput = oligoprop(input_CDS)
%propOutput = oligoprop(output_CDS)

if nt2aa(input_CDS) == nt2aa(output_CDS)
    disp("Optimized")
    r_fold_input = rnafold(input_CDS);
    r_fold_output = rnafold(output_CDS);

    rnaplot(r_fold_input, 'seq', input_CDS, 'format', 'dot');
    title('Input')

    rnaplot(r_fold_output, 'seq', output_CDS, 'format', 'dot');
    title('Output')
    toc
end

input_codon = [];
output_codon = [];
for i = 1:length(input_CDS)
    if mod(i,3) == 0
        input_codon = [input_codon convertCharsToStrings([input_CDS(i-2) input_CDS(i-1) input_CDS(i)])];
        output_codon = [output_codon convertCharsToStrings([output_CDS(i-2) output_CDS(i-1) output_CDS(i)])];
    end
end

disp(["Codons changed: " sum(input_codon == output_codon) " Percentage: " sum(input_codon == output_codon)/length(input_codon)*100]);
disp(["Nucleotides changed: " sum(input_CDS == output_CDS) " Percentage: " sum(input_CDS == output_CDS)/length(input_CDS)*100]);

function char_seq_codon = transform_seq(input, startpos)
    
    %get CDS from full mRNA
    CDS = input(startpos:end);
    AA_CDS = nt2aa(CDS, 'AlternativeStartCodons', false);
    
    %process the data
    [Data_raw_codons, AA_ref, Codon_ref, AA_list, n_orgs, restriction_data] = process_data();
    
    %make a list of preferred codons
    variance = get_variance(Data_raw_codons, AA_ref, AA_list, n_orgs);
    [ordered_codons, ordered_variances] = get_ordered_list(variance, AA_ref, AA_list, Codon_ref);
    
    %restructure AA_ref (because ordering codons gives different structure)
    for i = 1:length(AA_ref)
        AA_ref(i) = nt2aa(char(ordered_codons(i)), 'AlternativeStartCodons', false);
    end
    
    %make initial optimized sequence
    seq_codon = [];
    for i = 1:length(AA_CDS)
        aa_1 = AA_CDS(i);
        
        %matlab sometimes uses END and sometimes uses * for stop codons,
        %force it to use * at all times
        if aa_1 == 'end'
            aa_1 = '*';
        end
        
        %get the ordered codon list for this aminoacid
        aa_codons = ordered_codons(find(AA_ref == aa_1));
        codon_1 = aa_codons(1); %grab the first one
        
        seq_codon = [seq_codon codon_1];
    end

    char_seq_codon = char(strjoin(seq_codon));
    char_seq_codon = char_seq_codon(~isspace(char_seq_codon));
    
    disp(char_seq_codon);
    
    %now we will iterate over restriction sites and grab different codons
    char_seq_codon = eliminate_restrictions(seq_codon, char_seq_codon, AA_ref, ordered_codons, ordered_variances, restriction_data);
    
end

function [Data_raw_codons, AA_ref, Codon_ref, AA_list, n_orgs, restriction_data] = process_data()
    %import data
    [~,~,Data_raw_codons]=xlsread([pwd '/data_proof.xlsx']);
    [~,~,restriction_data_raw]=xlsread([pwd '/restriction_enzyme_database.xlsx']);
    fileID = fopen('restriction_sites.txt');
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
    AA_ref = string(Data_raw_codons(1,13:end));
    AA_ref(find(AA_ref == "stop")) = "*";
    Codon_ref = string(Data_raw_codons(2,13:end));
    AA_list = ["phe" "leu" "ile" "met" "val" "tyr" "*" "his" "gln" "asn" "lys" "asp" "glu" "ser" "pro" "thr" "ala" "cys" "trp" "arg" "gly"];
    n_orgs = length(Data_raw_codons(3:end,1));
end

function variance = get_variance(Data_raw, AA_ref, AA_list, n_orgs)
        
    freqs = zeros(n_orgs, length(AA_ref));

    for i = 1:n_orgs
        data = Data_raw(2+i,13:end); 
        freqs(i,:) = cell2mat(data);
    end

    variance = [];
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
       var_aa = var(freq_aa);
       variance = [variance var_aa];
    end
    
end

function [ordered_codons, ordered_variances] = get_ordered_list(variance, AA_ref, AA_list, Codon_ref)
    ordered_codons = [];
    ordered_variances = [];

    for i = 1:length(AA_list)
        %get the positions of this aminoacid
        aa = AA_list(i);
        aa_pos = find(AA_ref == aa);
        
        %get the list of codons and the variances
        aa_codons = Codon_ref(aa_pos);
        aa_variance = variance(aa_pos);
        
        %sort the variances
        sorted_list = sort(aa_variance);
        
        %use the ordered variance list to order the codons
        for ii = 1:length(aa_codons)
            pos = find(aa_variance  == sorted_list(ii));
            codon = aa_codons(pos);
            ordered_codons = [ordered_codons codon];
        end
        ordered_variances = [ordered_variances sorted_list];
    end
end

function output = eliminate_restrictions(seq_codon, char_seq_codon, AA_ref, ordered_codons, ordered_variances, restrictions)
    %now we will iterate over restriction sites and grab different codons
    
    empty = false;
    
    while ~empty
        %see if its truly empty of restriction sites, this is the best way my tired mind could
        %think of for now..
        c = 0;
        
        for j = 1:length(restrictions(:,1))
            [~, restriction_pos] = rebasecuts(char_seq_codon, {char(restrictions(j, 1))});
            if ~isempty(restriction_pos)
                c = c+1;
            end
        end
        if c == 0
            empty = true;
        end
        for i = 1:length(restrictions(:,1))
            %find the restriction site
            [~, restriction_pos] = rebasecuts(char_seq_codon, {char(restrictions(i, 1))});
            
            if ~isempty(restriction_pos)

                %determine the length of the recognition site left and right of
                %the recognition site the rebasecuts gives
                recognition_seq = char(restrictions(i, 2));
                length_right = length(recognition_seq) - find(recognition_seq == '''');
                length_left = length(recognition_seq) - length_right - 1;

                for ii = 1:length(restriction_pos)
                     
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
                    codons = seq_codon(codon_start:codon_start+total_span);

                    %define the next best choices
                    choices = 2.*ones(1,length(codons));

                    fixed = false;
                    while ~fixed    
                        %get the codon we need to change and what to change it to
                        [alt_codon, pos, cant_fix] = get_alt_codon(codons, AA_ref,  ordered_codons, ordered_variances, choices);

                        seq_codon(codon_start+pos-1) = alt_codon;
                        char_seq_codon = char(strjoin(seq_codon));
                        char_seq_codon = char_seq_codon(~isspace(char_seq_codon));

                        [~, check] = rebasecuts(char_seq_codon, {char(restrictions(i, 1))});

                        if isempty(find(check == restriction_pos(ii)))
                            %this site is now gone
                            fixed = true;
                        else
                            if cant_fix == true
                                fixed = true;
                                disp('Couldn''t fix it')
                            else
                                %not gone so change the choice of pos
                                choices(pos) = choices(pos) + 1;
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

function [alt_codon, pos, cant_fix] = get_alt_codon(codons, AA_ref,  ordered_codons, ordered_variances, choices)
    %get list of the amino acids
    AAs = [];
    for i = 1:length(codons)
        AAs = [AAs nt2aa(char(codons(i)))];
    end

    %make list of delta_variance
    counter = 0;
    delta_variances = zeros(1, length(AAs));
    for i = 1:length(delta_variances)
        aa = AAs(i);
        aa_pos = find(AA_ref == aa);

        codon_variances = ordered_variances(aa_pos);
        
        %only calculate difference if still possible otherwise just set it
        %at super high so it won't be considered for best option
        if length(codon_variances) < choices(i)
            delta_variances(i) = codon_variances(choices(i) - 1) - codon_variances(choices(i));
        else
            delta_variances(i) = 10^9;
            counter = counter + 1;
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
    syn_codons_pos = find(AA_ref == nt2aa(char(codon_to_change))); 
    alt_codons = ordered_codons(syn_codons_pos);
    alt_codon = alt_codons(choices(pos));
end












