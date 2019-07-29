function recognition_seq = TALE_checker(TALE_CDS)

    TALE_rules = ["NI" "NG" "NN" "HD"; "A" "T" "G" "C"];
    TALE_AA_char = nt2aa(TALE_CDS,'AlternativeStartCodons', false);
    TALE_AA_str  = string(TALE_AA_char);

    %%TALEsp1 and TALEsp2
    spacers = ["GGKQALETVQRLLPVLCQDHGLTPEQVVAIAS" "GGKQALETVQRLLPVLCQAHGLTPDQVVAIAS"]; 

    %%TALE1
    %spacers = ["GGKQALETVQRLLPVLCQAHGLTPEQVVAIAS" "GGKQALETVQRLLPVLCQAHGLTPEQVVAIAS"]; 

    spacer_1_pos = strfind(TALE_AA_str, spacers(1));
    spacer_2_pos = strfind(TALE_AA_str, spacers(2));
    
    if spacer_1_pos(1) < spacer_2_pos(1)
        startpos = spacer_1_pos(1);
    else
        startpos = spacer_2_pos(1);
    end

    RVD = [];
    spacer_left = true;
    while spacer_left
        RVD = [RVD string([TALE_AA_char(startpos-2) TALE_AA_char(startpos-1)])];
        startpos = startpos + 34;
        if isempty(strfind(string(TALE_AA_char(startpos:end)), spacers(1))) && isempty(strfind(string(TALE_AA_char(startpos:end)), spacers(2)))
            spacer_left = false;
        end
    end

    recognition_seq = [];
    for i = 1:length(RVD)
        nt = TALE_rules(2, find(TALE_rules(1,:) == RVD(i)));
        recognition_seq = [recognition_seq nt];
    end
    recognition_seq = char(strjoin(recognition_seq));
    recognition_seq = recognition_seq(~isspace(recognition_seq));
end