clc, clear all
disp(transform_seq("MATLA", "aa"));

function new_seq = transform_seq(seq, type)

    if type == "nt"
        seq = nt2aa(seq);
    end
    
    if ~ischar(seq)
        seq = char(seq);
    end

    [~,~,Data_raw]=xlsread('/Users/Dennis/Documents/MATLAB/Sci-Phi-Software-Tool/data_proof.xlsx');

    AA_ref = string(Data_raw(1,13:end));
    AA_ref(find(AA_ref == "stop")) = "*";
    Triplet_ref = string(Data_raw(2,13:end));
    AA_list = ["phe" "leu" "ile" "met" "val" "tyr" "*" "his" "gln" "asn" "lys" "asp" "glu" "ser" "pro" "thr" "ala" "cys" "trp" "arg" "gly"];
    n_orgs = length(Data_raw(3:end,1));

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

    seq_codon = [];
    for i = 1:length(seq)
        aa = lower(aminolookup(seq(i)));
        min_variance = min(variance(find(AA_ref == aa)));
        min_variance_pos = find(variance == min_variance);

        %it might occur that the min value also occurs somewhere else in the
        %variance array, we should take the one from the right AA
        if length(min_variance_pos) ~= 1
            for i = 1:length(min_variance_pos)
               if AA_ref(min_variance_pos(i)) == aa
                  codon = Triplet_ref(min_variance_pos(i));
               end
            end
        else 
            codon = Triplet_ref(variance == min_variance);
        end
        seq_codon = [seq_codon codon];
    end
    new_seq = strjoin(seq_codon);
end

