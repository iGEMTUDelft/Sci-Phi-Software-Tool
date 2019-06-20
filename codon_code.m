
input = 'gagacctatcagtgatagaggatccaaccgagcccaatataggactcagggtgccaaaaaaatgtctaccaccatgggaattcaacctcctaaaaagaaacgtaaagttaatactattaatattgctaaaaatgacttctcagatattgaattagcagccattccatttaatacattagcagatcactatggtgaacgtttagcacgtgaacagttagcattagaacatgaatcatatgaaatgggtgaagcacgttttcgtaagatgttcgagcgtcagttaaaagcaggtgaagttgcagataatgcagcagccaaacctttaattactacattattacctaaaatgattgctcgtattaacgattggtttgaagaggttaaagcaaagcgtggtaaacgtcctacagcatttcagttcttacaagaaatcaaacctgaagcagttgcatatattactattaaaacaacattagcatgtttaacatcagcagataatacaacagttcaagcagttgcatcagcaattggtcgtgcaattgaagatgaagcacgttttggtcgtattcgtgatttagaagccaaacattttaaaaaaaatgttgaagaacagttaaacaaacgtgttggtcatgtttataaaaaagcatttatgcaggttgttgaagcagatatgttatcaaaaggtttattaggtggtgaagcatggtcatcatggcataaagaagattcaattcatgttggtgttcgttgtattgaaatgttaattgaatcaacaggtatggtttcattacatcgtcagaatgcaggtgttgttggtcaagattcagaaacaattgaattagcacctgaatatgcagaagcaattgcaacacgtgcaggtgcattagcaggtatttcaccaatgtttcagccttgtgttgttcctcctaaaccttggacaggtattacaggtggtggttattgggcaaatggtcgtcgtcctttagcattagttcgtacacattcaaaaaaagcattaatgcgttatgaagatgtttacatgcctgaagtttataaagccattaatattgcacagaatactgcatggaaaatcaacaagaaagttttagcagttgcaaatgttattactaaatggaaacattgtcctgttgaagatattcctgcaattgaacgtgaagaattaccaatgaaacctgaagatattgatatgaatcctgaagcattaacagcatggaaacgtgcagcagcagccgtttatcgtaaagataaagcacgtaaatcacgtcgtatttcattagagttcatgttagagcaagcaaacaagttcgcaaatcataaagccatttggtttccttataatatggattggcgtggtcgtgtttatgcagtttcaatgtttaatcctcaaggtaatgatatgaccaaaggtttattaaccttagctaaaggtaaacctattggtaaagaaggttattattggttaaaaatccatggtgcaaattgtgcaggtgttgataaagttccttttccagaacgtattaagttcattgaagaaaatcatgaaaatattatggcatgtgctaaatcaccattagaaaatacatggtgggcagaacaagattcacctttttgttttttagccttttgttttgaatatgcaggtgttcaacatcatggtttatcatataattgttcattaccattagcatttgatggttcatgttcaggtattcagcatttttcagcaatgttacgtgatgaagttggtggtcgtgcagttaacttattaccttcagaaacagttcaagatatctatggtattgttgctaaaaaagttaatgaaatcttacaggcagatgccattaatggtacagataatgaagttgttacagttacagatgaaaatacaggtgaaatatcagaaaaagttaaattaggtaccaaagcattagcaggtcaatggttagcatatggtgttacacgttcagttactaagcgttcagttatgacattagcatatggttcaaaagagttcggttttcgtcaacaggttttagaagataccattcaacctgcaattgattcaggtaaaggtttaatgtttactcaacctaatcaagcagcaggttatatggccaaattaatatgggaatcagtttcagttacagttgttgcagcagttgaagcaatgaattggttaaaatcagcagccaaattattagcagcagaagttaaagataaaaaaacaggtgaaatcttacgtaagcgttgtgcagttcattgggttacacctgatggttttcctgtttggcaagaatataaaaaacctattcaaacacgtttaaacttaatgtttttaggtcagttccgtttacaacctactattaatacaaacaaagattcagaaattgatgcacataaacaagaatcaggtattgcacctaacttcgttcattcacaagatggttcacatttacgtaaaacagttgtttgggcacatgaaaaatatggtattgaatcatttgcattaattcacgattcatttggtaccattccagcagatgcagcaaacttattcaaagcagttcgtgaaacaatggttgatacatatgaatcatgtgatgttttagcagacttctatgatcagttcgcagatcagttacatgaatcacagttagataaaatgcctgcattacctgccaaaggtaacttaaacttacgtgatattttagaatcagacttcgcatttgcctaataaggcgcgccccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctc';

startpos = 62;

%propInput = oligoprop(input)
%propOutput = oligoprop(output)
input_CDS = upper(input(startpos:end));
output_CDS = transform_seq(input, startpos);

function char_seq_codon = transform_seq(input, startpos)
    
    %get CDS from full mRNA
    CDS = input(startpos:end);
    AA_CDS = nt2aa(CDS, 'AlternativeStartCodons', false);
    
    %process the data
    [Data_raw, AA_ref, Codon_ref, AA_list, n_orgs] = process_data();
    
    %make a list of preferred codons
    variance = get_variance(Data_raw, AA_ref, AA_list, n_orgs);
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
    
    %now we will iterate over restriction sites and grab different codons
    char_seq_codon = eliminate_restrictions(seq_codon, char_seq_codon, AA_ref, ordered_codons)
end

function [Data_raw, AA_ref, Codon_ref, AA_list, n_orgs] = process_data()
    %import data
    [~,~,Data_raw]=xlsread([pwd '/data_proof.xlsx']);
    
    %make reference lists
    AA_ref = string(Data_raw(1,13:end));
    AA_ref(find(AA_ref == "stop")) = "*";
    Codon_ref = string(Data_raw(2,13:end));
    AA_list = ["phe" "leu" "ile" "met" "val" "tyr" "*" "his" "gln" "asn" "lys" "asp" "glu" "ser" "pro" "thr" "ala" "cys" "trp" "arg" "gly"];
    n_orgs = length(Data_raw(3:end,1));
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

function output = eliminate_restrictions(seq_codon, char_seq_codon, AA_ref, ordered_codons)
    %now we will iterate over restriction sites and grab different codons
    empty = false;
    while ~empty
        enzyme_list = ["CfoI"]; %TODO import enzyme list
        for i = 1:length(enzyme_list)
            [~, sites] = rebasecuts(char_seq_codon, {char(enzyme_list(i))});
            for ii = 1:length(sites)
                choice = 2;        
                fixed = false;
                while ~fixed
                    if mod(sites(ii), 3) == 0 %end of codon so change either this one or the next one
                        codon = char(seq_codon(sites(i)/3));
                        aa = nt2aa(codon, 'AlternativeStartCodons', false);
                        alt_codons = ordered_codons(find(AA_ref == aa));  

                        %change the codon for the next option
                        seq_codon(sites(ii)/3) = alt_codons(choice);

                        char_seq_codon = char(strjoin(seq_codon));
                        char_seq_codon = char_seq_codon(~isspace(char_seq_codon));

                        [~, sites_new] = rebasecuts(char_seq_codon, {char(enzyme_list(i))});
                        check = find(sites_new == sites(ii));
                        if isempty(check)
                            fixed = true;
                        else
                            choice = choice + 1;
                            if choice > length(alt_codons)
                                fixed = true;
                                disp('Couldnt find better codon')
                            end
                        end

                    elseif mod(sites(ii), 3) == 2 %middle of codon so only change this one
                        codon = char(seq_codon(sites(ii)/3 + 1/3));
                        aa = nt2aa(codon, 'AlternativeStartCodons', false);
                        alt_codons = ordered_codons(find(AA_ref == aa));  

                        %change the codon for the next option
                        seq_codon(sites(ii)/3 + 1/3) = alt_codons(choice);

                        char_seq_codon = char(strjoin(seq_codon));
                        char_seq_codon = char_seq_codon(~isspace(char_seq_codon));

                        [~, sites_new] = rebasecuts(char_seq_codon, {char(enzyme_list(i))});
                        check = find(sites_new == sites(ii));
                        if isempty(check)
                            fixed = true;
                        else
                            choice = choice + 1;
                            if choice > length(alt_codons)
                                fixed = true;
                                disp('Couldnt find better codon')
                            end                   
                        end
                    elseif mod(sites(ii), 3) == 1 %beginning of codon so either change this one or the previous one
                        codon = char(seq_codon(sites(ii)/3 - 1/3));
                        aa = nt2aa(codon, 'AlternativeStartCodons', false);
                        alt_codons = ordered_codons(find(AA_ref == aa));  

                        %change the codon for the next option
                        seq_codon(sites(ii)/3 - 1/3) = alt_codons(choice);

                        char_seq_codon = char(strjoin(seq_codon));
                        char_seq_codon = char_seq_codon(~isspace(char_seq_codon));
                        
                        %check presence of the restriction site
                        [~, sites_new] = rebasecuts(char_seq_codon, {char(enzyme_list(i))});
                        check = find(sites_new == sites(ii));
                        if isempty(check)
                            fixed = true;
                        else
                            choice = choice + 1;
                            if choice > length(alt_codons)
                                fixed = true;
                                disp('Couldnt find better codon')
                            end
                        end
                    end
                end
            end
        end
        
        %final check to see if nothing is left
        [~, sites] = rebasecuts(char_seq_codon, {char(enzyme_list(i))});
        disp(sites)
        if isempty(sites)
            empty = true;
        end
    end
    output = char_seq_codon;
end