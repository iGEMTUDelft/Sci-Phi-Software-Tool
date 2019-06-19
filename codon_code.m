clc, clear all
input = 'gagacctatcagtgatagaggatccaaccgagcccaatataggactcagggtgccaaaaaaatgtctaccaccatgggaattcaacctcctaaaaagaaacgtaaagttaatactattaatattgctaaaaatgacttctcagatattgaattagcagccattccatttaatacattagcagatcactatggtgaacgtttagcacgtgaacagttagcattagaacatgaatcatatgaaatgggtgaagcacgttttcgtaagatgttcgagcgtcagttaaaagcaggtgaagttgcagataatgcagcagccaaacctttaattactacattattacctaaaatgattgctcgtattaacgattggtttgaagaggttaaagcaaagcgtggtaaacgtcctacagcatttcagttcttacaagaaatcaaacctgaagcagttgcatatattactattaaaacaacattagcatgtttaacatcagcagataatacaacagttcaagcagttgcatcagcaattggtcgtgcaattgaagatgaagcacgttttggtcgtattcgtgatttagaagccaaacattttaaaaaaaatgttgaagaacagttaaacaaacgtgttggtcatgtttataaaaaagcatttatgcaggttgttgaagcagatatgttatcaaaaggtttattaggtggtgaagcatggtcatcatggcataaagaagattcaattcatgttggtgttcgttgtattgaaatgttaattgaatcaacaggtatggtttcattacatcgtcagaatgcaggtgttgttggtcaagattcagaaacaattgaattagcacctgaatatgcagaagcaattgcaacacgtgcaggtgcattagcaggtatttcaccaatgtttcagccttgtgttgttcctcctaaaccttggacaggtattacaggtggtggttattgggcaaatggtcgtcgtcctttagcattagttcgtacacattcaaaaaaagcattaatgcgttatgaagatgtttacatgcctgaagtttataaagccattaatattgcacagaatactgcatggaaaatcaacaagaaagttttagcagttgcaaatgttattactaaatggaaacattgtcctgttgaagatattcctgcaattgaacgtgaagaattaccaatgaaacctgaagatattgatatgaatcctgaagcattaacagcatggaaacgtgcagcagcagccgtttatcgtaaagataaagcacgtaaatcacgtcgtatttcattagagttcatgttagagcaagcaaacaagttcgcaaatcataaagccatttggtttccttataatatggattggcgtggtcgtgtttatgcagtttcaatgtttaatcctcaaggtaatgatatgaccaaaggtttattaaccttagctaaaggtaaacctattggtaaagaaggttattattggttaaaaatccatggtgcaaattgtgcaggtgttgataaagttccttttccagaacgtattaagttcattgaagaaaatcatgaaaatattatggcatgtgctaaatcaccattagaaaatacatggtgggcagaacaagattcacctttttgttttttagccttttgttttgaatatgcaggtgttcaacatcatggtttatcatataattgttcattaccattagcatttgatggttcatgttcaggtattcagcatttttcagcaatgttacgtgatgaagttggtggtcgtgcagttaacttattaccttcagaaacagttcaagatatctatggtattgttgctaaaaaagttaatgaaatcttacaggcagatgccattaatggtacagataatgaagttgttacagttacagatgaaaatacaggtgaaatatcagaaaaagttaaattaggtaccaaagcattagcaggtcaatggttagcatatggtgttacacgttcagttactaagcgttcagttatgacattagcatatggttcaaaagagttcggttttcgtcaacaggttttagaagataccattcaacctgcaattgattcaggtaaaggtttaatgtttactcaacctaatcaagcagcaggttatatggccaaattaatatgggaatcagtttcagttacagttgttgcagcagttgaagcaatgaattggttaaaatcagcagccaaattattagcagcagaagttaaagataaaaaaacaggtgaaatcttacgtaagcgttgtgcagttcattgggttacacctgatggttttcctgtttggcaagaatataaaaaacctattcaaacacgtttaaacttaatgtttttaggtcagttccgtttacaacctactattaatacaaacaaagattcagaaattgatgcacataaacaagaatcaggtattgcacctaacttcgttcattcacaagatggttcacatttacgtaaaacagttgtttgggcacatgaaaaatatggtattgaatcatttgcattaattcacgattcatttggtaccattccagcagatgcagcaaacttattcaaagcagttcgtgaaacaatggttgatacatatgaatcatgtgatgttttagcagacttctatgatcagttcgcagatcagttacatgaatcacagttagataaaatgcctgcattacctgccaaaggtaacttaaacttacgtgatattttagaatcagacttcgcatttgcctaataaggcgcgccccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctc';

startpos = 62;

%propInput = oligoprop(input)
%propOutput = oligoprop(output)
input_CDS = upper(input(startpos:end))
output_CDS = transform_seq(input, startpos)
    

function new_seq = transform_seq(input, startpos)
    
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
      
    seq_codon = [];
    for i = 1:length(AA_CDS)
        aa = AA_CDS(i);
        
        %matlab sometimes uses END and sometimes uses * for stop codons,
        %force it to use * at all times
        if aa == 'end'
            aa = '*';
        end
        
        %get the ordered codon list for this aminoacid
        aa_codons = ordered_codons(find(AA_ref == aa));
        codon = aa_codons(1); %grab the first one
        
        seq_codon = [seq_codon codon];
    end
    new_seq = char(strjoin(seq_codon));
    new_seq = new_seq(~isspace(new_seq));
    
    %verify that the new sequence has the same AA sequence
    if nt2aa(new_seq, 'AlternativeStartCodons', false) == AA_CDS
        disp("Optimized");
    end
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
