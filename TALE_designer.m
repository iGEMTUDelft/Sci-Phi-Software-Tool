
function TALE_nt = TALE_designer(target_recog)
    %optimized spacers
    spacers_aa = ["GGKQALETVQRLLPVLCQDHGLTPEQVVAIAS" "GGKQALETVQRLLPVLCQAHGLTPDQVVAIAS" "GGKQALETVQRLLPVLCQAHGLTPAQVVAIAS" "GGKQALETVQRLLPVLCQDHGLTPDQVVAIAS"]; 

    prefix_nt = 'ATGGTAGATTTAAGAACTTTAGGATATTCACAGCAGCAACAGGAAAAGATCAAGCCCAAAGTTAGGTCGACAGTCGCGCAGCATCACGAAGCGCTGGTTGGTCATGGGTTTACACATGCCCACATCGTAGCCTTATCGCAGCACCCTGCAGCCCTTGGCACGGTCGCCGTCAAGTACCAGGACATGATTGCGGCGTTGCCGGAAGCCACACATGAGGCGATCGTCGGTGTGGGGAAACAGTGGAGCGGAGCCCGAGCGCTTGAGGCCCTGTTGACGGTCGCGGGAGAGCTGAGAGGGCCTCCCCTTCAGCTGGACACGGGCCAGTTGCTGAAGATCGCGAAGCGGGGAGGAGTCACGGCGGTCGAGGCGGTGCACGCGTGGCGCAATGCGCTCACGGGAGCACCCCTCAATCTGACCCCGGATCAGGTGGTTGCGATCGCCAGT';
    suffix_nt = 'AACAACGGGGGCAGACCCGCACTGGAGTCAATCGTGGCCCAGCTTTCGAGGCCGGACCCCGCGCTGGCCGCACTCACTAATGATCATCTTGTAGCGCTGGCCTGCCTCGGCGGACGACCCGCCTTGGATGCGGTGAAGAAGGGGCTCCCGCACGCGCCTGCATTGATTAAGCGGACCAACAGAAGGATTCCCGAGAGGACATCACATCGAGTGGCAGATCACGCGCAAGTGGTCCGCGTGCTCGGATTCTTCCAGTGTCACTCCCACCCCGCACAAGCGTTCGATGACGCCATGACTCAATTTGGTATGTCGAGACACGGACTGCTGCAGCTCTTTCGTAGAGTCGGTGTCACAGAACTCGAGGCCCGCTCGGGCACACTGCCTCCCGCCTCCCAGCGGTGGGACAGGATTCTCCAAGCGAGCGGTATGAAACGCGCGAAGCCTTCACCTACGTCAACTCAGACACCTGACCAGGCGAGCCTTCATGCGTTCGCAGACTCGCTGGAGAGGGATTTGGACGCGCCCTCGCCCATGCATGAAGGGGACCAAACTCGCGCGTCATAA';

    prefix_aa = 'MVDLRTLGYSQQQQEKIKPKVRSTVAQHHEALVGHGFTHAHIVALSQHPAALGTVAVKYQDMIAALPEATHEAIVGVGKQWSGARALEALLTVAGELRGPPLQLDTGQLLKIAKRGGVTAVEAVHAWRNALTGAPLNLTPDQVVAIAS';
    suffix_aa = 'NNGGRPALESIVAQLSRPDPALAALTNDHLVALACLGGRPALDAVKKGLPHAPALIKRTNRRIPERTSHRVADHAQVVRVLGFFQCHSHPAQAFDDAMTQFGMSRHGLLQLFRRVGVTELEARSGTLPPASQRWDRILQASGMKRAKPSPTSTQTPDQASLHAFADSLERDLDAPSPMHEGDQTRAS*';

    TALE_rules = ["NI" "NG" "NN" "HD"; "A" "T" "G" "C"];

    %first we use two times spacer 1 then spacer 2 etc
    RVD_AA = [];
    c = 1;
    for ii = 1:length(target_recog)
        spacer_AA = spacers_aa(c);
        if ii == length(target_recog)
            spacer_AA = spacers_aa(1);
        elseif ii == 1
            spacer_AA = spacers_aa(1);
        end
        AA = TALE_rules(1, find(TALE_rules(2,:) == target_recog(ii)));

        RVD_AA = [RVD_AA AA spacer_AA];

        if c>=4
            c = 1;
        else
            c = c+1;
        end
    end
    RVD_AA = char(strjoin(RVD_AA));
    RVD_AA = RVD_AA(~isspace(RVD_AA));

    TALE_nt = [prefix_nt make_RVD_nt(RVD_AA) suffix_nt];
    TALE_nt = char(strjoin(TALE_nt));
    TALE_nt = TALE_nt(~isspace(TALE_nt));

    TALE_AA = [prefix_aa string(RVD_AA) suffix_aa];
    TALE_AA = char(strjoin(TALE_AA));
    TALE_AA = TALE_AA(~isspace(TALE_AA));

    function RVD_nt = make_RVD_nt(RVD_AA)

        %reference component consists of TALEsp1 
        reference_RVD_nt = 'AATGGCGGAGGCAAACAAGCACTCGAAACTGTACAGCGCCTCCTGCCGGTACTGTGCCAAGATCATGGCTTGACGCCTGAGCAAGTCGTAGCTATTGCATCAAACATCGGTGGTAAACAGGCGCTGGAAACCGTACAACGATTACTCCCTGTCTTATGCCAGGCACACGGTCTGACCCCCGACCAAGTTGTAGCCATTGCGTCTAATGGCGGCGGGAAACAGGCCCTGGAAACGGTCCAACGTCTGTTACCCGTTCTGTGTCAGGCTCACGGTCTGACCCCTGCCCAGGTAGTTGCAATTGCCAGCAACATCGGCGGGAAACAAGCGCTGGAGACTGTGCAGCGTCTGCTCCCTGTGTTATGCCAAGATCATGGGCTCACTCCGGATCAGGTGGTGGCCATCGCTTCCAATATTGGCGGTAAACAGGCGCTGGAGACAGTGCAACGACTTTTACCTGTTCTCTGCCAGGATCATGGTCTAACTCCCGAGCAGGTCGTCGCCATCGCCTCTCATGACGGCGGGAAACAAGCGTTGGAAACTGTCCAGCGACTCCTGCCGGTTTTGTGCCAGGCCCACGGGCTTACTCCTGACCAGGTAGTTGCGATCGCGTCAAATGGGGGTGGCAAACAAGCCCTCGAAACCGTGCAACGCCTGCTGCCCGTCTTGTGCCAAGCTCATGGGCTGACTCCGGCGCAAGTAGTCGCGATTGCGAGCCACGATGGCGGTAAGCAGGCACTGGAAACGGTTCAGCGCCTGCTCCCGGTTCTATGCCAGGATCACGGCCTGACCCCGGACCAGGTCGTCGCGATCGCGTCAAATATCGGTGGCAAACAAGCTTTGGAGACAGTACAGCGCCTGTTACCAGTGCTTTGCCAGGACCATGGTCTGACCCCTGAGCAAGTAGTGGCGATCGCTTCTAATATTGGGGGCAAACAAGCGCTGGAAACAGTACAGCGTCTGTTACCGGTCCTATGCCAGGCACATGGCCTGACCCCTGATCAGGTGGTAGCCATTGCCAGTCATGATGGCGGTAAACAGGCGCTTGAGACTGTCCAACGTCTGCTGCCGGTCCTCTGTCAGGCTCATGGCCTGACGCCAGCTCAAGTCGTGGCTATCGCTTCGCATGATGGCGGAAAACAGGCACTGGAGACTGTGCAGCGACTGTTGCCAGTTCTGTGTCAGGATCACGGTTTAACTCCGGACCAGGTGGTCGCTATTGCGTCGAATGGGGGCGGTAAACAAGCGCTGGAAACTGTGCAACGTTTGCTCCCAGTTCTGTGCCAGGACCATGGGCTGACTCCGGAACAGGTAGTGGCCATTGCTTCTAATATTGGTGGGAAACAGGCGCTGGAAACCGTGCAGCGCCTGCTTCCAGTGCTTTGCCAGGCCCATGGCCTGACGCCAGATCAGGTGGTTGCTATAGCCAGCAATGGCGGCGGTAAACAGGCCCTCGAAACCGTCCAGCGCCTGCTCCCTGTGCTGTGCCAGGCCCATGGGCTTACCCCAGCGCAAGTAGTGGCGATTGCGTCTAATATTGGTGGTAAACAGGCGTTGGAGACTGTACAACGCCTGCTGCCAGTTTTATGCCAAGATCATGGTCTGACCCCTGAGCAGGTAGTGGCTATTGCATCC';
        reference_RVD_AA = nt2aa(reference_RVD_nt,'AlternativeStartCodons', false);

        reference_RVD_codons = [];

        for i = 1:length(reference_RVD_nt)
            if mod(i,3) == 0
                reference_RVD_codon = string([reference_RVD_nt(i-2) reference_RVD_nt(i-1) reference_RVD_nt(i)]);
                reference_RVD_codons = [reference_RVD_codons reference_RVD_codon];
            end
        end

        %process the codon usage data
        [Data_raw_codons, AA_ref, Codon_ref, AA_list] = process_data();
        [differences, AA_ref_new, Codon_ref_new] = get_difference(Data_raw_codons, AA_ref, AA_list, Codon_ref);

        %we will now look at the reference to make a decision on the right
        %codon
        RVD_codons = string([1, length(RVD_AA)]);
        for j = 1:length(RVD_AA)
            %if they are the same AA, just use the same codon
            if reference_RVD_AA(j) == RVD_AA(j);
                RVD_codons(j) = reference_RVD_codons(j);
            else
                %pos, so we know where to look for alternative codons
                relative_pos = find(AA_ref == lower(aminolookup(RVD_AA(j))));
                reference_pos = find(Codon_ref_new(1,:) == reference_RVD_codons(j));

                pos = relative_pos+reference_pos(1);

                diff = differences(pos);

                alt_pos = pos(find(diff == min(diff)));
                RVD_codons(j) = Codon_ref_new(2,alt_pos);

            end
        end
        RVD_nt = char(strjoin(RVD_codons));
        RVD_nt = string(RVD_nt(~isspace(RVD_nt)));
    end

    function [Data_raw_codons, AA_ref, Codon_ref, AA_list] = process_data()
        %import data
        [~,~,Data_raw_codons]=xlsread([pwd '/data_proof.xlsx']);

        %make reference lists
        AA_list = ["phe" "leu" "ile" "met" "val" "tyr" "*" "his" "gln" "asn" "lys" "asp" "glu" "ser" "pro" "thr" "ala" "cys" "trp" "arg" "gly"];
        AA_raw = string(Data_raw_codons(1,13:end));
        AA_raw(find(AA_raw == "stop")) = "*";
        Codon_raw = string(Data_raw_codons(2,13:end));

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
    end

    function [differences, AA_ref_new, Codon_ref_new] = get_difference(Data_raw_codons, AA_ref, AA_list, Codon_ref)

        %get the data for frequency
        data = Data_raw_codons(3,13:end); 
        freqs = cell2mat(data);

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

        %we will calculate the difference for each codon with each codon, 
        %we need a new reference list for that
        total_length = 64.^2;

        %define difference and new reference lists
        differences = zeros(1, total_length);
        AA_ref_new = strings([1,total_length]);
        Codon_ref_new = strings([2,total_length]);

        %make the new reference lists
        startpos_AA = 1;
        startpos_codons = 1;
        c = 1;
        for i = 1:length(AA_list)
            %new AA ref list
            leng = length(find(AA_ref == AA_list(i)));
            leng = leng*64;
            AA_ref_new(startpos_AA:startpos_AA+leng-1) = AA_list(i);

            %new codon ref list
            codon = Codon_ref(c);
            leng_codons = 64; 
            for ii = 1:leng
                Codon_ref_new(1,startpos_codons:startpos_codons+leng_codons-1) = codon;
                Codon_ref_new(2,startpos_codons:startpos_codons+leng_codons-1) = Codon_ref;
                if mod(ii,leng_codons) == 0
                    c = c + 1;
                    if c <= length(Codon_ref)
                        codon = Codon_ref(c);
                    end
                    startpos_codons = startpos_codons + leng_codons;
                end
            end
            startpos_AA = startpos_AA+leng;  
        end
        

        %calculate the differences
        startpos = 0;
        for i = 1:length(Codon_ref)
           freq_codon = percentage(i);
           for ii = 1:length(percentage)
               differences(ii+startpos) = abs(percentage(ii) - freq_codon);
           end
           startpos = startpos+64;
        end    
    end
end
