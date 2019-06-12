% Hi guys. 
% In order to work with this script you need to do a few steps
% before it works. First, you need to load in the data by opening this 
% Matlab code in the same map as the data base. Secondly, you need to 
% choose with which GC content concentration you want to work. In case you 
% want to work with all the bacteria in this databes, then you disable all
% the GC-content subparts. 

close all
clear all
clc
 
[~,~,Data__raw]=xlsread('/Users/osmanesen/Documents/MATLAB/optimum_codon_databasa.xlsx');
%[~,~,Data_raw]=xlsread('optimum_codon_databasa.xlsx');%

% %%
Data_Names = Data__raw(2:end,1);
[~,~,Data__search]=xlsread('/Users/osmanesen/Documents/MATLAB/organism.xlsx');
Data__search = Data__search'; 
Data__search = Data__search';
reduced_Data = ismember(Data_Names,Data__search);
Data_reduced = Data__raw(reduced_Data > 0,:); 
Data__raw = Data_reduced;
%%
GC_content = cell2mat(Data__raw(2:end,2));
%% GC content < 30 %             Run this part in order to obtain the histogram for this group of bacteria and disable the other GC content parts  
% GC_Low_30 = Data__raw(GC_content < 30,:);
% Data__raw = GC_Low_30;

%% GC content > 60 %             Run this part in order to obtain the histogram for this group of bacteria and disable the other GC content parts
 %GC_High_60 = Data__raw(GC_content > 60, :);
 %Data__raw = GC_High_60;

%% GC content  > 50  & < 60      Run this part in order to obtain the histogram for this group of bacteria and disable the other GC content parts
% GC_Low_60 = Data__raw(GC_content < 60, :);
%GC_content_60 = cell2mat(GC_Low_60(2:end,2));
%GC_Betw_50_60 = GC_Low_60(GC_content_60 > 50, :);

%Data__raw = GC_Betw_50_60;

%% GC content  > 40  & < 50      Run this part in order to obtain the histogram for this group of bacteria and disable the other GC content parts
%GC_Low_50 = Data__raw(GC_content < 50, :);
%GC_content_50 = cell2mat(GC_Low_50(2:end,2));
%GC_Betw_40_50 = GC_Low_50(GC_content_50 > 40, :);

%Data__raw = GC_Betw_40_50;

%% GC content  > 30  & < 40      Run this part in order to obtain the histogram for this group of bacteria and disable the other GC content parts
%GC_Low_40 = Data__raw(GC_content < 40, :);
%GC_content_40 = cell2mat(GC_Low_40(2:end,2));
%GC_Betw_30_40 = GC_Low_40(GC_content_40 > 30, :);

%Data__raw = GC_Betw_30_40;
%% This part will create the data for each amino acid 

data_Ala = Data__raw(2:end,3);
data_Arg = Data__raw(2:end,4);
data_Asn = Data__raw(2:end,5);
data_Asp = Data__raw(2:end,6);
data_Cys = Data__raw(2:end,7);
data_Gln = Data__raw(2:end,8);
data_Glu = Data__raw(2:end,9);
data_Gly = Data__raw(2:end,10);
data_His = Data__raw(2:end,11);
data_Ile = Data__raw(2:end,12);
data_Lue = Data__raw(2:end,13);
data_Lys = Data__raw(2:end,14);
data_Phe = Data__raw(2:end,15);
data_Pro = Data__raw(2:end,16);
data_Ser = Data__raw(2:end,17);
data_Thr = Data__raw(2:end,18);
data_Tyr = Data__raw(2:end,19);
data_Val = Data__raw(2:end,20);

%% Plotting best Ala  
 Codon_Ala = unique(data_Ala);
 N = numel(Codon_Ala);
 
 for i = 1 : N
 B = count(data_Ala,Codon_Ala(i,1));
 Count_Ala(:,i) = sum(B);
 end
 
 Count_Ala = Count_Ala';
 
overview_Ala = table( Codon_Ala, Count_Ala)

figure ;
bar(Count_Ala)
set(gca,'XTickLabel',{'GCA','GCC','GCG','GCT','None'})
title('\fontsize{20} Ala')
ylabel('Frequency of most favorable codon')

%% Plotting best Arg
 Codon_Arg = unique(data_Arg);
 N = numel(Codon_Arg);
 
 for i = 1: N
 B = count(data_Arg,Codon_Arg(i,1));
 Count_Arg(:,i) = sum(B);
 end
 
 Count_Arg = Count_Arg';
 
overview_Arg = table( Codon_Arg, Count_Arg)

figure ;
bar(Count_Arg)
set(gca,'XTickLabel',{'AGA','AGG','CGC','CGG','CGT','None'})
title('\fontsize{20}Arg')
ylabel('Frequency of most favorable codon')

%% Plotting best Asn 
 Codon_Asn = unique(data_Asn);
 N = numel(Codon_Asn);
 
 for i = 1: N
 B = count(data_Asn,Codon_Asn(i,1));
 Count_Asn(:,i) = sum(B);
 end
 
 Count_Asn = Count_Asn';
 
overview_Asn = table( Codon_Asn, Count_Asn)

figure ;
bar(Count_Asn)
set(gca,'XTickLabel',{'ACC','AAT','None'})
title('\fontsize{20}Arg')
ylabel('Frequency of most favorable codon')

%% Plotting best Asp
 Codon_Asp = unique(data_Asp);
 N = numel(Codon_Asp);
 
 for i = 1: N
 B = count(data_Asp,Codon_Asp(i,1));
 Count_Asp(:,i) = sum(B);
 end
 
 Count_Asp = Count_Asp';
 
overview_Asp = table( Codon_Asp, Count_Asp)

figure ;
bar(Count_Asp)
set(gca,'XTickLabel',{'GAC','GAT','None'})
title('\fontsize{20}Asp')
ylabel('Frequency of most favorable codon')

%% Plotting best Cys
 Codon_Cys = unique(data_Cys);
 N = numel(Codon_Cys);
 
 for i = 1: N
 B = count(data_Cys,Codon_Cys(i,1));
 Count_Cys(:,i) = sum(B);
 end
 
 Count_Cys = Count_Cys';
 
overview_Cys = table( Codon_Cys, Count_Cys)


figure ;
bar(Count_Cys)
set(gca,'XTickLabel',{'TGC','TGT','None'})
title('\fontsize{20}Cys')
ylabel('Frequency of most favorable codon')

%% Plotting best Gln
 Codon_Gln = unique(data_Gln);
 N = numel(Codon_Gln);
 
 for i = 1: N
 B = count(data_Gln,Codon_Gln(i,1));
 Count_Gln(:,i) = sum(B);
 end
 
 Count_Gln = Count_Gln';
 
overview_Gln = table( Codon_Gln, Count_Gln)


figure ;
bar(Count_Gln)
set(gca,'XTickLabel',{'CAA','CAG','None'})
title('\fontsize{20}Gln')
ylabel('Frequency of most favorable codon')
%% Plotting best Glu
 Codon_Glu = unique(data_Glu);
 N = numel(Codon_Glu);
 
 for i = 1: N
 B = count(data_Glu,Codon_Glu(i,1));
 Count_Glu(:,i) = sum(B);
 end
 
 Count_Glu = Count_Glu';
 
overview_Glu = table( Codon_Glu, Count_Glu)

figure ;
bar(Count_Glu)
set(gca,'XTickLabel',{'GAA','GAG','None'})
title('\fontsize{20}Glu')
ylabel('Frequency of most favorable codon')

%% Plotting best Gly
 Codon_Gly = unique(data_Gly);
 N = numel(Codon_Gly);
 
 for i = 1: N
 B = count(data_Gly,Codon_Gly(i,1));
 Count_Gly(:,i) = sum(B);
 end
 
 Count_Gly = Count_Gly';
 
overview_Gly = table( Codon_Gly, Count_Gly)

figure ;
bar(Count_Gly)
set(gca,'XTickLabel',{'GGA','GGC','GGG','GGT','None'})
title('\fontsize{20}Gly')
ylabel('Frequency of most favorable codon')

%% Plotting best His
 Codon_His = unique(data_His);
 N = numel(Codon_His);
 
 for i = 1: N
 B = count(data_His,Codon_His(i,1));
 Count_His(:,i) = sum(B);
 end
 
 Count_His = Count_His';
 
overview_His = table( Codon_His, Count_His)

figure ;
bar(Count_His)
set(gca,'XTickLabel',{'CAC','CAT','None'})
title('\fontsize{20}His')
ylabel('Frequency of most favorable codon')

%% Plotting best Lue
 Codon_Lue = unique(data_Lue);
 N = numel(Codon_Lue);
 
 for i = 1: N
 B = count(data_Lue,Codon_Lue(i,1));
 Count_Lue(:,i) = sum(B);
 end
 
 Count_Lue = Count_Lue';
 
overview_Lue = table( Codon_Lue, Count_Lue)

figure ;
bar(Count_Lue)
set(gca,'XTickLabel',{'CTA','CTC','CTG','CTT','TTA','TTG','None'})
title('\fontsize{20}Lue')
ylabel('Frequency of most favorable codon')
%% Plotting best Ile
 Codon_Ile = unique(data_Ile);
 N = numel(Codon_Ile);
 
 for i = 1: N
 B = count(data_Ile,Codon_Ile(i,1));
 Count_Ile(:,i) = sum(B);
 end
 
 Count_Ile = Count_Ile';
 
overview_Ile = table( Codon_Ile, Count_Ile)

figure ;
bar(Count_Ile)
set(gca,'XTickLabel',{'ATA','ATC','ATT','None'})
title('\fontsize{20}Ile')
ylabel('Frequency of most favorable codon')
%% Plotting best Lys
 Codon_Lys = unique(data_Lys);
 N = numel(Codon_Lys);
 
 for i = 1: N
 B = count(data_Lys,Codon_Lys(i,1));
 Count_Lys(:,i) = sum(B);
 end
 
 Count_Lys = Count_Lys';
 
overview_Lys = table( Codon_Lys, Count_Lys)

figure ;
bar(Count_Lys)
set(gca,'XTickLabel',{'AAA','AAG','None'})
title('\fontsize{20}Lys')
ylabel('Frequency of most favorable codon')

%% Plotting best Phe
 Codon_Phe = unique(data_Phe);
 N = numel(Codon_Phe);
 
 for i = 1: N
 B = count(data_Phe,Codon_Phe(i,1));
 Count_Phe(:,i) = sum(B);
 end
 
 Count_Phe = Count_Phe';
 
overview_Phe = table( Codon_Phe, Count_Phe)

figure ;
bar(Count_Phe)
set(gca,'XTickLabel',{'TTC','TTT','GGG','None'})
title('\fontsize{20}Phe')
ylabel('Frequency of most favorable codon')

%% Plotting best Pro
 Codon_Pro = unique(data_Pro);
 N = numel(Codon_Pro);
 
 for i = 1: N
 B = count(data_Pro,Codon_Pro(i,1));
 Count_Pro(:,i) = sum(B);
 end
 
 Count_Pro = Count_Pro';
 
overview_Pro = table( Codon_Pro, Count_Pro)

figure ;
bar(Count_Pro)
set(gca,'XTickLabel',{'CCA','CCC','CCG','CCT','None'})
title('\fontsize{20}Pro')
ylabel('Frequency of most favorable codon')
%% Plotting best Ser
 Codon_Ser = unique(data_Ser);
 N = numel(Codon_Ser);
 
 for i = 1: N
 B = count(data_Ser,Codon_Ser(i,1));
 Count_Ser(:,i) = sum(B);
 end
 
 Count_Ser = Count_Ser';
 
overview_Ser = table( Codon_Ser, Count_Ser)

figure ;
bar(Count_Ser)
set(gca,'XTickLabel',{'AGC','AGT','TCA','TCC','TCG','TCT','None'})
title('\fontsize{20}Ser')
ylabel('Frequency of most favorable codon')

%% Plotting best Thr
 Codon_Thr = unique(data_Thr);
 N = numel(Codon_Thr);
 
 for i = 1: N
 B = count(data_Thr,Codon_Thr(i,1));
 Count_Thr(:,i) = sum(B);
 end
 
 Count_Thr = Count_Thr';
 
overview_Thr = table( Codon_Thr, Count_Thr)

figure ;
bar(Count_Thr)
set(gca,'XTickLabel',{'ACA','ACC','ACG','ACT','None'})
title('\fontsize{20}Thr')
ylabel('Frequency of most favorable codon')

%% Plotting best Tyr
 Codon_Tyr = unique(data_Tyr);
 N = numel(Codon_Tyr);
 
 for i = 1: N
 B = count(data_Tyr,Codon_Tyr(i,1));
 Count_Tyr(:,i) = sum(B);
 end
 
 Count_Tyr = Count_Tyr';
 
overview_Tyr = table( Codon_Tyr, Count_Tyr)

figure ;
bar(Count_Tyr)
set(gca,'XTickLabel',{'TAC','TAT','None'})
title('\fontsize{20}Tyr')
ylabel('Frequency of most favorable codon')

%% Plotting best Val
 Codon_Val = unique(data_Val);
 N = numel(Codon_Val);
 
 for i = 1: N
 B = count(data_Val,Codon_Val(i,1));
 Count_Val(:,i) = sum(B);
 end
 
 Count_Val = Count_Val';
 
overview_Val = table( Codon_Val, Count_Val)

figure ;
bar(Count_Val)
set(gca,'XTickLabel',{'GTA','GTC','GTT','None'})
title('\fontsize{20}Val')
ylabel('Frequency of most favorable codon')


%% Overview of all the histograms
% 
% figure;
% subplot(6,3,1)
% bar(Count_Ala)
% title('\fontsize{20} Ala')
% ylabel('Freq codon')
% 
% subplot(6,3,2)
% bar(Count_Arg)
% title('\fontsize{20}Arg')
% ylabel('Freq codon')
% 
% subplot(6,3,3)
% bar(Count_Asn)
% title('\fontsize{20}Arg')
% ylabel('Freq codon')
% 
% subplot(6,3,4)
% bar(Count_Asp)
% title('\fontsize{20}Asp')
% ylabel('Freq codon')
% 
% subplot(6,3,5)
% bar(Count_Cys)
% title('\fontsize{20}Cys')
% ylabel('Freq codon')
% 
% subplot(6,3,6)
% bar(Count_Gln)
% title('\fontsize{20}Gln')
% ylabel('Freq codon')
% 
% subplot(6,3,7)
% bar(Count_Glu)
% title('\fontsize{20}Glu')
% ylabel('Freq codon')
% 
% subplot(6,3,8)
% bar(Count_Gly)
% title('\fontsize{20}Gly')
% ylabel('Freq codon')
% 
% subplot(6,3,9)
% bar(Count_His)
% title('\fontsize{20}His')
% ylabel('Freq codon')
% 
% subplot(6,3,10)
% bar(Count_Lue)
% title('\fontsize{20}Lue')
% ylabel('Freq codon')
% 
% subplot(6,3,11)
% bar(Count_Ile)
% title('\fontsize{20}Ile')
% ylabel('Freq codon')
% 
% subplot(6,3,12)
% bar(Count_Lys)
% title('\fontsize{20}Lys')
% ylabel('Freq codon')
% 
% subplot(6,3,13)
% bar(Count_Phe)
% title('\fontsize{20}Phe')
% ylabel('Freq codon')
% 
% subplot(6,3,14)
% bar(Count_Pro)
% title('\fontsize{20}Pro')
% ylabel('Freq codon')
% 
% subplot(6,3,15)
% bar(Count_Ser)
% title('\fontsize{20}Ser')
% ylabel('Freq codon')
% 
% subplot(6,3,16)
% bar(Count_Thr)
% title('\fontsize{20}Thr')
% ylabel('Freq codon')
% 
% subplot(6,3,17)
% bar(Count_Tyr)
% title('\fontsize{20}Tyr')
% ylabel('Freq codon')
% 
% subplot(6,3,18)
% bar(Count_Val)
% title('\fontsize{20}Val')
% ylabel('Freq codon')
 %% Overview of all the histograms

figure;
subplot(6,3,1)
bar(Count_Ala)
set(gca,'XTickLabel',{'GCC','GCG','GCT','None'})
title('\fontsize{20} Ala')
ylabel('Freq codon')

subplot(6,3,2)
bar(Count_Arg)
set(gca,'XTickLabel',{'AGA','AGG','CGC','CGG','CGT','None'})
title('\fontsize{20}Arg')
ylabel('Freq codon')

subplot(6,3,3)
bar(Count_Asn)
set(gca,'XTickLabel',{'ACC','AAT','None'})
title('\fontsize{20}Arg')
ylabel('Freq codon')

subplot(6,3,4)
bar(Count_Asp)
set(gca,'XTickLabel',{'GAC','GAT','None'})
title('\fontsize{20}Asp')
ylabel('Freq codon')

subplot(6,3,5)
bar(Count_Cys)
set(gca,'XTickLabel',{'TGC','TGT','None'})
title('\fontsize{20}Cys')
ylabel('Freq codon')

subplot(6,3,6)
bar(Count_Gln)
set(gca,'XTickLabel',{'CAA','CAG','None'})
title('\fontsize{20}Gln')
ylabel('Freq codon')

subplot(6,3,7)
bar(Count_Glu)
set(gca,'XTickLabel',{'GAA','GAG','None'})
title('\fontsize{20}Glu')
ylabel('Freq codon')

subplot(6,3,8)
bar(Count_Gly)
set(gca,'XTickLabel',{'GGA','GGC','GGG','GGT','None'})
title('\fontsize{20}Gly')
ylabel('Freq codon')

subplot(6,3,9)
bar(Count_His)
set(gca,'XTickLabel',{'CAC','CAT','None'})
title('\fontsize{20}His')
ylabel('Freq codon')

subplot(6,3,10)
bar(Count_Lue)
set(gca,'XTickLabel',{'CTA','CTC','CTG','CTT','TTA','TTG','None'})
title('\fontsize{20}Lue')
ylabel('Freq codon')

subplot(6,3,11)
bar(Count_Ile)
set(gca,'XTickLabel',{'ATA','ATC','ATT','None'})
title('\fontsize{20}Ile')
ylabel('Freq codon')

subplot(6,3,12)
bar(Count_Lys)
set(gca,'XTickLabel',{'AAA','AAG','None'})
title('\fontsize{20}Lys')
ylabel('Freq codon')

subplot(6,3,13)
bar(Count_Phe)
set(gca,'XTickLabel',{'TTC','TTT','GGG','None'})
title('\fontsize{20}Phe')
ylabel('Freq codon')

subplot(6,3,14)
bar(Count_Pro)
set(gca,'XTickLabel',{'CCA','CCC','CCG','CCT','None'})
title('\fontsize{20}Pro')
ylabel('Freq codon')

subplot(6,3,15)
bar(Count_Ser)
set(gca,'XTickLabel',{'AGC','TCA','TCC','TCG','TCT','None'})
title('\fontsize{20}Ser')
ylabel('Freq codon')

subplot(6,3,16)
bar(Count_Thr)
set(gca,'XTickLabel',{'ACA','ACC','ACG','ACT','None'})
title('\fontsize{20}Thr')
ylabel('Freq codon')

subplot(6,3,17)
bar(Count_Tyr)
set(gca,'XTickLabel',{'TAC','TAT','None'})
title('\fontsize{20}Tyr')
ylabel('Freq codon')

subplot(6,3,18)
bar(Count_Val)
set(gca,'XTickLabel',{'GTA','GTC','GTG','GTT','None'})
title('\fontsize{20}Val')
ylabel('Freq codon')