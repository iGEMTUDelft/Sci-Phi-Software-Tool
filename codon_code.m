close all
clear all
clc
%%

[~,~,Data__raw]=xlsread('/Users/dennis/Documents/MATLAB/iGEM/Sci-Phi-Software-Tool/data_proof.xlsx');
final_variance= zeros(1,64);
%% phe
A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
phe={'phe'};
reduced_Data = ismember(A,phe);
a = [vector_nul  reduced_Data];
Data_phe = Data__raw(:,a > 0);
Alleen_waardes_phe = Data_phe(3:end,2:end);
Alleen_waardes_phe = cell2mat(Alleen_waardes_phe);
sum_rows_phe= sum(Alleen_waardes_phe, 2);


freq_phe = Alleen_waardes_phe;
N = length(sum_rows_phe);     
for i= 1: N
    freq_phe(i,:)=freq_phe(i,:)./sum_rows_phe(i,1);
end

phe_var = var(freq_phe);

for i=1:numel(phe_var)
    final_variance(i)=phe_var(i);
end

%% leu
A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
leu={'leu'};
reduced_Data = ismember(A,leu);
a = [vector_nul  reduced_Data];
Data_leu = Data__raw(:,a > 0);
Alleen_waardes_leu = Data_leu(3:end,2:end);
Alleen_waardes_leu = cell2mat(Alleen_waardes_leu);
sum_rows_leu= sum(Alleen_waardes_leu, 2);


freq_leu = Alleen_waardes_leu;
N = length(sum_rows_leu);     
for i= 1: N
    freq_leu(i,:)=freq_leu(i,:)./sum_rows_leu(i,1);
end

leu_var = var(freq_leu);

for i=1:numel(leu_var)
    final_variance(numel(phe_var)+i)=leu_var(i);
end

%% ile
A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
ile={'ile'};
reduced_Data = ismember(A,ile);
a = [vector_nul  reduced_Data];
Data_ile = Data__raw(:,a > 0);
Alleen_waardes_ile = Data_ile(3:end,2:end);
Alleen_waardes_ile = cell2mat(Alleen_waardes_ile);
sum_rows_ile= sum(Alleen_waardes_ile, 2);


freq_ile = Alleen_waardes_ile;
N = length(sum_rows_ile);     
for i= 1: N
    freq_ile(i,:)=freq_ile(i,:)./sum_rows_ile(i,1);
end

ile_var = var(freq_ile);

for i=1:numel(ile_var)
    final_variance(numel(phe_var)+numel(leu_var)+i)=ile_var(i);
end

%% met
A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
met={'met'};
reduced_Data = ismember(A,met);
a = [vector_nul  reduced_Data];
Data_met = Data__raw(:,a > 0);
Alleen_waardes_met = Data_met(3:end,2:end);
Alleen_waardes_met = cell2mat(Alleen_waardes_met);
sum_rows_met= sum(Alleen_waardes_met, 2);


freq_met = sum_rows_met;
N = length(sum_rows_met);     
for i= 1: N
    freq_met(i,:)=freq_met(i,:)./sum_rows_met(i,1);
end

met_var = var(freq_met);

for i=1:numel(met_var)
    final_variance(numel(phe_var)+numel(leu_var)+numel(ile_var)+i)=met_var(i);
end
 
%% val
A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
val={'val'};
reduced_Data = ismember(A,val);
a = [vector_nul reduced_Data];
Data_val = Data__raw(:,a > 0);
Alleen_waardes_val = Data_val(3:end,2:end);
Alleen_waardes_val = cell2mat(Alleen_waardes_val);
sum_rows_val= sum(Alleen_waardes_val, 2);



freq_val = Alleen_waardes_val;
N = length(sum_rows_val);     
for i= 1: N
    freq_val(i,:)=freq_val(i,:)./sum_rows_val(i,1);
end

val_var = var(freq_val);

for i=1:numel(val_var)
    final_variance(numel(phe_var)+numel(leu_var)+numel(ile_var)+numel(met_var)+i)=val_var(i);
end
%% tyr
A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
tyr={'tyr'};
reduced_Data = ismember(A,tyr);
a = [vector_nul reduced_Data];
Data_tyr = Data__raw(:,a > 0);
Alleen_waardes_tyr = Data_tyr(3:end,2:end);
Alleen_waardes_tyr = cell2mat(Alleen_waardes_tyr);
sum_rows_tyr= sum(Alleen_waardes_tyr, 2);


freq_tyr = Alleen_waardes_tyr;
N = length(sum_rows_tyr);     
for i= 1: N
    freq_tyr(i,:)=freq_tyr(i,:)./sum_rows_tyr(i,1);
end

tyr_var = var(freq_tyr);

for i=1:numel(tyr_var)
    final_variance(numel(phe_var)+numel(leu_var)+numel(ile_var)+numel(met_var)+numel(val_var)+i)= tyr_var(i);
end
%% stop

A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
stop={'stop'};
reduced_Data = ismember(A,stop);
a = [vector_nul reduced_Data];
Data_stop = Data__raw(:,a > 0);
Alleen_waardes_stop = Data_stop(3:end,2:end);
Alleen_waardes_stop = cell2mat(Alleen_waardes_stop);
sum_rows_stop= sum(Alleen_waardes_stop, 2);


freq_stop = Alleen_waardes_stop;
N = length(sum_rows_stop);     
for i= 1: N
    freq_stop(i,:)=freq_stop(i,:)./sum_rows_stop(i,1);
end

stop_var = var(freq_stop);

for i=1:numel(stop_var)
    final_variance(numel(phe_var)+numel(leu_var)+numel(ile_var)+numel(met_var)+numel(val_var)+numel(tyr_var)+i)= stop_var(i);
end
%% his

A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
his={'his'};
reduced_Data = ismember(A,his);
a = [vector_nul reduced_Data];
Data_his = Data__raw(:,a > 0);
Alleen_waardes_his = Data_his(3:end,2:end);
Alleen_waardes_his = cell2mat(Alleen_waardes_his);
sum_rows_his= sum(Alleen_waardes_his, 2);


freq_his = Alleen_waardes_his;
N = length(sum_rows_his);     
for i= 1: N
    freq_his(i,:)=freq_his(i,:)./sum_rows_his(i,1);
end

his_var = var(freq_his);

for i=1:numel(his_var)
    final_variance(numel(phe_var)+numel(leu_var)+numel(ile_var)+numel(met_var)+numel(val_var)+numel(tyr_var)+numel(stop_var)+i)= his_var(i);
end

%% gln


A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
gln={'gln'};
reduced_Data = ismember(A,gln);
a = [vector_nul reduced_Data];
Data_gln = Data__raw(:,a > 0);
Alleen_waardes_gln = Data_gln(3:end,2:end);
Alleen_waardes_gln = cell2mat(Alleen_waardes_gln);
sum_rows_gln= sum(Alleen_waardes_gln, 2);


freq_gln = Alleen_waardes_gln;
N = length(sum_rows_gln);     
for i= 1: N
    freq_gln(i,:)=freq_gln(i,:)./sum_rows_gln(i,1);
end

gln_var = var(freq_gln);

for i=1:numel(gln_var)
    final_variance(numel(phe_var)+numel(leu_var)+numel(ile_var)+numel(met_var)+numel(val_var)+numel(tyr_var)+numel(stop_var)+numel(his_var)+i)= gln_var(i);
end

%% asn
A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
asn={'asn'};
reduced_Data = ismember(A,asn);
a = [vector_nul reduced_Data];
Data_asn = Data__raw(:,a > 0);
Alleen_waardes_asn = Data_asn(3:end,2:end);
Alleen_waardes_asn = cell2mat(Alleen_waardes_asn);
sum_rows_asn= sum(Alleen_waardes_asn, 2);



freq_asn = Alleen_waardes_asn;
N = length(sum_rows_asn);     
for i= 1: N
    freq_asn(i,:)=freq_asn(i,:)./sum_rows_asn(i,1);
end

asn_var = var(freq_asn);
 
 for i=1:numel(asn_var)
     final_variance(numel(phe_var)+numel(leu_var)+numel(ile_var)+numel(met_var)+numel(val_var)+numel(tyr_var)+numel(stop_var)+numel(his_var)+numel(gln_var)+i)= asn_var(i);
end


%% lys

A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
lys={'lys'};
reduced_Data = ismember(A,lys);
a = [vector_nul reduced_Data];
Data_lys = Data__raw(:,a > 0);
Alleen_waardes_lys = Data_lys(3:end,2:end);
Alleen_waardes_lys = cell2mat(Alleen_waardes_lys);
sum_rows_lys= sum(Alleen_waardes_lys, 2);

freq_lys = Alleen_waardes_lys;
N = length(sum_rows_lys);     
for i= 1: N
    freq_lys(i,:)=freq_lys(i,:)./sum_rows_lys(i,1);
end

lys_var = var(freq_lys);
 
 for i=1:numel(lys_var)
     final_variance(numel(phe_var)+numel(leu_var)+numel(ile_var)+numel(met_var)+numel(val_var)+numel(tyr_var)+numel(stop_var)+numel(his_var)+numel(gln_var)+numel(asn_var)+i)= lys_var(i);
 end
%% asp

A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
asp={'asp'};
reduced_Data = ismember(A,asp);
a = [vector_nul reduced_Data];
Data_asp = Data__raw(:,a > 0);
Alleen_waardes_asp = Data_asp(3:end,2:end);
Alleen_waardes_asp = cell2mat(Alleen_waardes_asp);
sum_rows_asp= sum(Alleen_waardes_asp, 2);



freq_asp = Alleen_waardes_asp;
N = length(sum_rows_asp);     
for i= 1: N
    freq_asp(i,:)=freq_asp(i,:)./sum_rows_asp(i,1);
end

asp_var = var(freq_asp);
 
 for i=1:numel(asp_var)
     final_variance(numel(phe_var)+numel(leu_var)+numel(ile_var)+numel(met_var)+numel(val_var)+numel(tyr_var)+numel(stop_var)+numel(his_var)+numel(gln_var)+numel(asn_var)+numel(lys_var)+i)= asp_var(i);
 end
%% glu
A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
glu={'glu'};
reduced_Data = ismember(A,glu);
a = [vector_nul reduced_Data];
Data_glu = Data__raw(:,a > 0);
Alleen_waardes_glu = Data_glu(3:end,2:end);
Alleen_waardes_glu = cell2mat(Alleen_waardes_glu);
sum_rows_glu= sum(Alleen_waardes_glu, 2);



freq_glu = Alleen_waardes_glu;
N = length(sum_rows_glu);     
for i= 1: N
    freq_glu(i,:)=freq_glu(i,:)./sum_rows_glu(i,1);
end

glu_var = var(freq_glu);
 
 for i=1:numel(glu_var)
     final_variance(numel(phe_var)+numel(leu_var)+numel(ile_var)+numel(met_var)+numel(val_var)+numel(tyr_var)+numel(stop_var)+numel(his_var)+numel(gln_var)+numel(asn_var)+numel(lys_var)+numel(asp_var)+i)= glu_var(i);
 end
%% ser
A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
ser={'ser'};
reduced_Data = ismember(A,ser);
a = [vector_nul reduced_Data];
Data_ser = Data__raw(:,a > 0);
Alleen_waardes_ser = Data_ser(3:end,2:end);
Alleen_waardes_ser = cell2mat(Alleen_waardes_ser);
sum_rows_ser= sum(Alleen_waardes_ser, 2);


freq_ser = Alleen_waardes_ser;
N = length(sum_rows_ser);     
for i= 1: N
    freq_ser(i,:)=freq_ser(i,:)./sum_rows_ser(i,1);
end

ser_var = var(freq_ser);
 
 for i=1:numel(ser_var)
     final_variance(numel(phe_var)+numel(leu_var)+numel(ile_var)+numel(met_var)+numel(val_var)+numel(tyr_var)+numel(stop_var)+numel(his_var)+numel(gln_var)+numel(asn_var)+numel(lys_var)+numel(asp_var)+numel(glu_var)+i)= ser_var(i);
 end
%% pro

A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
pro={'pro'};
reduced_Data = ismember(A,pro);
a = [vector_nul reduced_Data];
Data_pro = Data__raw(:,a > 0);
Alleen_waardes_pro = Data_pro(3:end,2:end);
Alleen_waardes_pro = cell2mat(Alleen_waardes_pro);
sum_rows_pro= sum(Alleen_waardes_pro, 2);


freq_pro = Alleen_waardes_pro;
N = length(sum_rows_pro);     
for i= 1: N
    freq_pro(i,:)=freq_pro(i,:)./sum_rows_pro(i,1);
end

pro_var = var(freq_pro);
 
 for i=1:numel(pro_var)
     final_variance(numel(phe_var)+numel(leu_var)+numel(ile_var)+numel(met_var)+numel(val_var)+numel(tyr_var)+numel(stop_var)+numel(his_var)+numel(gln_var)+numel(asn_var)+numel(lys_var)+numel(asp_var)+numel(glu_var)+numel(ser_var)+i)= pro_var(i);
 end

%%
A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
thr={'thr'};
reduced_Data = ismember(A,thr);
a = [vector_nul reduced_Data];
Data_thr = Data__raw(:,a > 0);
Alleen_waardes_thr = Data_thr(3:end,2:end);
Alleen_waardes_thr = cell2mat(Alleen_waardes_thr);
sum_rows_thr= sum(Alleen_waardes_thr, 2);
 
 
 
freq_thr = Alleen_waardes_thr;
N = length(sum_rows_thr);     
for i= 1: N
    freq_thr(i,:)=freq_thr(i,:)./sum_rows_thr(i,1);
end

thr_var = var(freq_thr);
 
 for i=1:numel(thr_var)
     final_variance(numel(phe_var)+numel(leu_var)+numel(ile_var)+numel(met_var)+numel(val_var)+numel(tyr_var)+numel(stop_var)+numel(his_var)+numel(gln_var)+numel(asn_var)+numel(lys_var)+numel(asp_var)+numel(glu_var)+numel(ser_var)+numel(pro_var)+i)= thr_var(i);
 end
 
 
%% ala
A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
ala={'ala'};
reduced_Data = ismember(A,ala);
a = [vector_nul  reduced_Data];
Data_ala = Data__raw(:,a > 0);
Alleen_waardes_ala = Data_ala(3:end,2:end);
Alleen_waardes_ala = cell2mat(Alleen_waardes_ala);
sum_rows_ala= sum(Alleen_waardes_ala, 2);


freq_ala = Alleen_waardes_ala;
N = length(sum_rows_ala);     
for i= 1: N
    freq_ala(i,:)=freq_ala(i,:)./sum_rows_ala(i,1);
end


ala_var = var(freq_ala);

 
 for i=1:numel(ala_var)
     final_variance(numel(phe_var)+numel(leu_var)+numel(ile_var)+numel(met_var)+numel(val_var)+numel(tyr_var)+numel(stop_var)+numel(his_var)+numel(gln_var)+numel(asn_var)+numel(lys_var)+numel(asp_var)+numel(glu_var)+numel(ser_var)+numel(pro_var)+numel(thr_var)+i)= ala_var(i);
 end
%% cys
A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
cys={'cys'};
reduced_Data = ismember(A,cys);
a = [vector_nul  reduced_Data];
Data_cys = Data__raw(:,a > 0);
Alleen_waardes_cys = Data_cys(3:end,2:end);
Alleen_waardes_cys = cell2mat(Alleen_waardes_cys);
sum_rows_cys= sum(Alleen_waardes_cys, 2);


freq_cys = Alleen_waardes_cys;
N = length(sum_rows_cys);     
for i= 1: N
    freq_cys(i,:)=freq_cys(i,:)./sum_rows_cys(i,1);
end


cys_var = var(freq_cys);

 
 for i=1:numel(cys_var)
     final_variance(numel(phe_var)+numel(leu_var)+numel(ile_var)+numel(met_var)+numel(val_var)+numel(tyr_var)+numel(stop_var)+numel(his_var)+numel(gln_var)+numel(asn_var)+numel(lys_var)+numel(asp_var)+numel(glu_var)+numel(ser_var)+numel(pro_var)+numel(thr_var)+numel(ala_var)+i)= cys_var(i);
 end
 
 %% trp
A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
trp={'trp'};
reduced_Data = ismember(A,trp);
a = [vector_nul  reduced_Data];
Data_trp = Data__raw(:,a > 0);
Alleen_waardes_trp = Data_trp(3:end,2:end);
Alleen_waardes_trp = cell2mat(Alleen_waardes_trp);
sum_rows_trp= sum(Alleen_waardes_trp, 2);


freq_trp = Alleen_waardes_trp;
N = length(sum_rows_trp);     
for i= 1: N
    freq_trp(i,:)=freq_trp(i,:)./sum_rows_trp(i,1);
end


trp_var = var(freq_trp);

 
 for i=1:numel(trp_var)
     final_variance(numel(phe_var)+numel(leu_var)+numel(ile_var)+numel(met_var)+numel(val_var)+numel(tyr_var)+numel(stop_var)+numel(his_var)+numel(gln_var)+numel(asn_var)+numel(lys_var)+numel(asp_var)+numel(glu_var)+numel(ser_var)+numel(pro_var)+numel(thr_var)+numel(ala_var)+numel(cys_var)+i)= trp_var(i);
 end
 
 %% arg
A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
arg={'arg'};
reduced_Data = ismember(A,arg);
a = [vector_nul  reduced_Data];
Data_arg = Data__raw(:,a > 0);
Alleen_waardes_arg = Data_arg(3:end,2:end);
Alleen_waardes_arg= cell2mat(Alleen_waardes_arg);
sum_rows_arg= sum(Alleen_waardes_arg, 2);


freq_arg = Alleen_waardes_arg;
N = length(sum_rows_arg);
for i= 1: N
    freq_arg(i,:)=freq_arg(i,:)./sum_rows_arg(i,1);
end


arg_var = var(freq_arg);

 
 for i=1:numel(arg_var)
     final_variance(numel(phe_var)+numel(leu_var)+numel(ile_var)+numel(met_var)+numel(val_var)+numel(tyr_var)+numel(stop_var)+numel(his_var)+numel(gln_var)+numel(asn_var)+numel(lys_var)+numel(asp_var)+numel(glu_var)+numel(ser_var)+numel(pro_var)+numel(thr_var)+numel(ala_var)+numel(cys_var)+numel(trp_var)+i)= arg_var(i);
 end
 
 
 
 
 %% gly
A = Data__raw(1,13:end);
vector_nul = zeros(1,12);
vector_nul(1,4)= 1;
gly={'gly'};
reduced_Data = ismember(A,gly);
a = [vector_nul  reduced_Data];
Data_gly = Data__raw(:,a > 0);
Alleen_waardes_gly = Data_gly(3:end,2:end);
Alleen_waardes_gly= cell2mat(Alleen_waardes_gly);
sum_rows_gly= sum(Alleen_waardes_gly, 2);


freq_gly = Alleen_waardes_gly;
N = length(sum_rows_gly);     
for i= 1: N
    freq_gly(i,:)=freq_gly(i,:)./sum_rows_gly(i,1);
end


gly_var = var(freq_gly);

 
 for i=1:numel(gly_var)
     final_variance(numel(phe_var)+numel(leu_var)+numel(ile_var)+numel(met_var)+numel(val_var)+numel(tyr_var)+numel(stop_var)+numel(his_var)+numel(gln_var)+numel(asn_var)+numel(lys_var)+numel(asp_var)+numel(glu_var)+numel(ser_var)+numel(pro_var)+numel(thr_var)+numel(ala_var)+numel(cys_var)+numel(trp_var)+numel(arg_var)+i)= gly_var(i);
 end

%%
final_variance;
final_AA = ["Phe" "Phe" "Leu" "Leu" "Leu" "Leu" "Leu" "Leu" "Ile" "Ile" "Ile" "Met" "Val" "Val" "Val" "Val" "Tyr" "Tyr" "*" "*" "*" "His" "His" "Gln" "Gln" "Asn" "Asn" "Lys" "Lys" "Asp" "Asp" "Glu" "Glu" "Ser" "Ser" "Ser" "Ser" "Ser" "Ser" "Pro" "Pro" "Pro" "Pro" "Thr" "Thr" "Thr" "Thr" "Ala" "Ala" "Ala" "Ala" "Cys" "Cys" "Trp" "Arg" "Arg" "Arg" "Arg" "Arg" "Arg" "Gly" "Gly" "Gly" "Gly"];
final_triplet = ["TTT" "TTC" "TTA" "TTG" "CTT" "CTC" "CTA" "CTG" "ATT" "ATC" "ATA" "ATG" "GTT" "GTC" "GTA" "GTG" "TCT" "TAC" "TCA" "TCG" "TGA" "CCT" "CCC" "CCA" "CCG" "ACT" "ACC" "ACA" "ACG" "GCT" "GCC" "GCA" "GGG" "TAT" "TAC" "TAA" "TAG" "AGT" "AGC" "CAT" "CAC" "CAA" "CAG" "AAT" "AAC" "AAA" "AAG" "GAT" "GAC" "GAA" "GAG" "TGT" "TGC" "TGG" "CGT" "CGC" "CGA" "CGG" "AGA" "AGG" "GGT" "GGC" "GGA" "GGG"];

seq_AA = 'MATLA';
seq_codon = [];
for i = 1:length(seq_AA)
    AA = aminolookup(seq_AA(i));
    min_variance = min(final_variance(find(final_AA == AA)));
    codon = final_triplet(final_variance == min_variance);
    seq_codon = [seq_codon codon];
end
seq_codon = strjoin(seq_codon)

