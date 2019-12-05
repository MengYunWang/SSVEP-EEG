%% Transform the time domain into frequency domain

% Created by M.-Y. Wang
% 05-12-2019

%% FFT for each condition
clear all
clc

%--------------------------------------------------------------------Condition 1 Neutral
data1_name = dir ('F:\face-random\Preprocessing\Conditions\Condition1_Neutral\*.set');
EEG = pop_loadset('filename',data1_name(1).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition1_Neutral');
%define parameters
time2ft = dsearchn (EEG.times',[300;2100]);
data_n = size (EEG.CSD(:,time2ft(1):time2ft(2),:),2);
data_hz = linspace (0,EEG.srate/2,(EEG.srate/2)/0.1+1); %(EEG.srate/2)/0.2 = data_n/2, 0.2 is the resolution we want
data_N = (EEG.srate/2)/0.1*2;
%initialization 
Neutral_amp = zeros (EEG.nbchan,length(data_hz),length(data1_name)); %channels*data_points*sub_numbs
Neutral_pow = zeros (EEG.nbchan,length(data_hz),length(data1_name));
Happy_amp = zeros (EEG.nbchan,length(data_hz),length(data1_name));
Happy_pow = zeros (EEG.nbchan,length(data_hz),length(data1_name));
N2H_amp = zeros (EEG.nbchan,length(data_hz),length(data1_name));
N2H_pow = zeros (EEG.nbchan,length(data_hz),length(data1_name));
H2N_amp = zeros (EEG.nbchan,length(data_hz),length(data1_name));
H2N_pow = zeros (EEG.nbchan,length(data_hz),length(data1_name));
Static_amp = zeros (EEG.nbchan,length(data_hz),length(data1_name));
Static_pow = zeros (EEG.nbchan,length(data_hz),length(data1_name));
Dynamic_amp = zeros (EEG.nbchan,length(data_hz),length(data1_name));
Dynamic_pow = zeros (EEG.nbchan,length(data_hz),length(data1_name));

for ii = 1:length(data1_name);
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition1_Neutral\');
    temp_data = squeeze (nanmean(EEG.CSD(:,time2ft(1):time2ft(2),:),3));
    data1_ft = fft (temp_data,data_N,2);
    data1_apt1 = abs (data1_ft)./data_N;
    data1_apt2 = data1_apt1 (:,1:data_N/2+1);
    data1_apt2(:,2:end-1) = 2.*data1_apt2(:,2:end-1);
    data1_pow = data1_apt2.^2;
    Neutral_amp(:,:,ii) = data1_apt2;
    Neutral_pow(:,:,ii) = data1_pow;
end
    Neutral_amp(:,:,[10,15,17,19,21]) = nan;
    Neutral_pow(:,:,[10,15,17,19,21]) = nan;

%--------------------------------------------------------------------Condition 2 Happy
data2_name = dir ('F:\face-random\Preprocessing\Conditions\Condition2_Happy\*.set');
for ii = 1:length(data2_name);
    EEG = pop_loadset('filename',data2_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition2_Happy\');
    temp_data = squeeze (nanmean(EEG.CSD(:,time2ft(1):time2ft(2),:),3));
    data2_ft = fft (temp_data,data_N,2);
    data2_apt1 = abs (data2_ft)./data_N;
    data2_apt2 = data2_apt1 (:,1:data_N/2+1);
    data2_apt2(:,2:end-1) = 2.*data2_apt2(:,2:end-1);
    data2_pow = data2_apt2.^2;
    Happy_amp(:,:,ii) = data2_apt2;
    Happy_pow(:,:,ii) = data2_pow;
end
    Happy_amp(:,:,[10,15,17,19,21]) = nan;
    Happy_pow(:,:,[10,15,17,19,21]) = nan;

%--------------------------------------------------------------------Condition 3 N2H
data3_name = dir ('F:\face-random\Preprocessing\Conditions\Condition3_N2H\*.set');
for ii = 1:length(data3_name);
    EEG = pop_loadset('filename',data3_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition3_N2H\');
    temp_data = squeeze (nanmean(EEG.CSD(:,time2ft(1):time2ft(2),:),3));
    data3_ft = fft (temp_data,data_N,2);
    data3_apt1 = abs (data3_ft)./data_N;
    data3_apt2 = data3_apt1 (:,1:data_N/2+1);
    data3_apt2(:,2:end-1) = 2.*data3_apt2(:,2:end-1);
    data3_pow = data3_apt2.^2;
    N2H_amp(:,:,ii) = data3_apt2;
    N2H_pow(:,:,ii) = data3_pow;
end
    N2H_amp(:,:,[10,15,17,19,21]) = nan;
    N2H_pow(:,:,[10,15,17,19,21]) = nan;

%--------------------------------------------------------------------Condition 4 H2N
data4_name = dir ('F:\face-random\Preprocessing\Conditions\Condition4_H2N\*.set');
for ii = 1:length(data4_name);
    EEG = pop_loadset('filename',data4_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition4_H2N\');
    temp_data = squeeze (nanmean(EEG.CSD(:,time2ft(1):time2ft(2),:),3));
    data4_ft = fft (temp_data,data_N,2);
    data4_apt1 = abs (data4_ft)./data_N;
    data4_apt2 = data4_apt1 (:,1:data_N/2+1);
    data4_apt2(:,2:end-1) = 2.*data4_apt2(:,2:end-1);
    data4_pow = data4_apt2.^2;
    H2N_amp(:,:,ii) = data4_apt2;
    H2N_pow(:,:,ii) = data4_pow;
end
    H2N_amp(:,:,[10,15,17,19,21]) = nan;
    H2N_pow(:,:,[10,15,17,19,21]) = nan;

%--------------------------------------------------------------------Static
for ii = 1:length(data2_name);
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition1_Neutral\');
    temp_data1 = squeeze (nanmean(EEG.CSD(:,time2ft(1):time2ft(2),:),3));
    EEG = pop_loadset('filename',data2_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition2_Happy\');
    temp_data2 = squeeze (nanmean(EEG.CSD(:,time2ft(1):time2ft(2),:),3));
    temp_data = (temp_data1 + temp_data2)./2;
    data2_ft = fft (temp_data,data_N,2);
    data2_apt1 = abs (data2_ft)./data_N;
    data2_apt2 = data2_apt1 (:,1:data_N/2+1);
    data2_apt2(:,2:end-1) = 2.*data2_apt2(:,2:end-1);
    data2_pow = data2_apt2.^2;
    Static_amp(:,:,ii) = data2_apt2;
    Static_pow(:,:,ii) = data2_pow;
end
    Static_amp(:,:,[10,15,17,19,21]) = nan;
    Static_pow(:,:,[10,15,17,19,21]) = nan;    
    
%--------------------------------------------------------------------Dynamic
for ii = 1:length(data2_name);
    EEG = pop_loadset('filename',data3_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition3_N2H\');
    temp_data1 = squeeze (nanmean(EEG.CSD(:,time2ft(1):time2ft(2),:),3));
    EEG = pop_loadset('filename',data4_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition4_H2N\');
    temp_data2 = squeeze (nanmean(EEG.CSD(:,time2ft(1):time2ft(2),:),3));
    temp_data = (temp_data1 + temp_data2)./2;
    data2_ft = fft (temp_data,data_N,2);
    data2_apt1 = abs (data2_ft)./data_N;
    data2_apt2 = data2_apt1 (:,1:data_N/2+1);
    data2_apt2(:,2:end-1) = 2.*data2_apt2(:,2:end-1);
    data2_pow = data2_apt2.^2;
    Dynamic_amp(:,:,ii) = data2_apt2;
    Dynamic_pow(:,:,ii) = data2_pow;
end
    Dynamic_amp(:,:,[10,15,17,19,21]) = nan;
    Dynamic_pow(:,:,[10,15,17,19,21]) = nan;     
    
    
    
save fft_amp_pow Neutral_amp Neutral_pow Happy_amp  Happy_pow N2H_amp N2H_pow...
    H2N_amp H2N_pow Static_amp Static_pow Dynamic_amp Dynamic_pow data_hz EEG

%% Compute the SNR and Z-score
    Neutral_SNR = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));  Neutral_Z = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));
    Happy_SNR = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));    Happy_Z = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));
    N2H_SNR = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));      N2H_Z = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));
    H2N_SNR = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));      H2N_Z = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));
    Static_SNR = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));   Static_Z = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));
    Dynamic_SNR = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));  Dynamic_Z = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));
    for mm = 1:size (Neutral_amp,3);
        for ii = 11:length (data_hz)-11;
            neutral_temp = zeros (64,18);
            neutral_temp(:,1:9) = squeeze (Neutral_amp(:,ii-10:ii-2,mm)); neutral_temp (:,10:18) = squeeze (Neutral_amp(:,ii+2:ii+10,mm));
            [~,indx1] = max (neutral_temp,[],2); for jj = 1:64; neutral_temp (jj,indx1(jj)) = NaN; end
            [~,indx2] = max (neutral_temp,[],2); for jj = 1:64; neutral_temp (jj,indx2(jj)) = NaN; end
                        
            happy_temp = zeros (64,18);
            happy_temp(:,1:9) = squeeze (Happy_amp(:,ii-10:ii-2,mm)); happy_temp (:,10:18) = squeeze (Happy_amp(:,ii+2:ii+10,mm));
            [~,indx1] = max (happy_temp,[],2); for jj = 1:64; happy_temp (jj,indx1(jj)) = NaN; end
            [~,indx2] = max (happy_temp,[],2); for jj = 1:64; happy_temp (jj,indx2(jj)) = NaN; end
            
            base_temp1 = (neutral_temp + happy_temp)./2;
            Neutral_SNR (:,ii,mm) = squeeze(Neutral_amp(:,ii,mm)) ./ nanmean(base_temp1,2);
            Neutral_Z (:,ii,mm)= (squeeze(Neutral_amp(:,ii,mm)) - nanmean(base_temp1,2)) ./ nanstd(base_temp1,[],2);
            Happy_SNR (:,ii,mm) = squeeze (Happy_amp(:,ii,mm)) ./ nanmean(base_temp1,2);
            Happy_Z (:,ii,mm) = (squeeze (Happy_amp(:,ii,mm)) - nanmean(base_temp1,2)) ./ nanstd(base_temp1,[],2);
%----------------------------------------------------------------            
            N2H_temp = zeros (64,18);
            N2H_temp(:,1:9) = squeeze (N2H_amp(:,ii-10:ii-2,mm)); N2H_temp (:,10:18) = squeeze (N2H_amp(:,ii+2:ii+10,mm));
            [~,indx1] = max (N2H_temp,[],2);   for jj = 1:64; N2H_temp (jj,indx1(jj)) = NaN; end
            [~,indx2] = max (N2H_temp,[],2);   for jj = 1:64; N2H_temp (jj,indx2(jj)) = NaN; end;
            
            H2N_temp = zeros (64,18);
            H2N_temp(:,1:9) = squeeze (H2N_amp(:,ii-10:ii-2,mm)); H2N_temp (:,10:18) = squeeze (H2N_amp(:,ii+2:ii+10,mm));
            [~,indx1] = max (H2N_temp,[],2); for jj = 1:64; H2N_temp (jj,indx1(jj)) = NaN; end
            [~,indx2] = max (H2N_temp,[],2); for jj = 1:64; H2N_temp (jj,indx2(jj)) = NaN; end

            base_temp2 = (N2H_temp + H2N_temp)./2;
            N2H_SNR (:,ii,mm) = squeeze (N2H_amp(:,ii,mm)) ./ nanmean(base_temp2,2);
            N2H_Z (:,ii,mm) = (squeeze (N2H_amp(:,ii,mm)) - nanmean(base_temp2,2)) ./ nanstd(base_temp2,[],2);
            H2N_SNR (:,ii,mm) = squeeze (H2N_amp(:,ii,mm)) ./ nanmean(base_temp2,2);
            H2N_Z (:,ii,mm) = (squeeze (H2N_amp(:,ii,mm)) - nanmean(base_temp2,2)) ./ nanstd(base_temp2,[],2);
%---------------------------------------------------------------            
            Static_temp = zeros (64,18);
            Static_temp(:,1:9) = squeeze (Static_amp(:,ii-10:ii-2,mm)); Static_temp (:,10:18) = squeeze (Static_amp(:,ii+2:ii+10,mm));
            [~,indx1] = max (Static_temp,[],2); for jj = 1:64; Static_temp (jj,indx1(jj)) = NaN; end
            [~,indx2] = max (Static_temp,[],2); for jj = 1:64; Static_temp (jj,indx2(jj)) = NaN; end
            
            Dynamic_temp = zeros (64,18);
            Dynamic_temp(:,1:9) = squeeze (Dynamic_amp(:,ii-10:ii-2,mm)); Dynamic_temp (:,10:18) = squeeze (Dynamic_amp(:,ii+2:ii+10,mm));
            [~,indx1] = max (Dynamic_temp,[],2); for jj = 1:64; Dynamic_temp (jj,indx1(jj)) = NaN; end
            [~,indx2] = max (Dynamic_temp,[],2); for jj = 1:64; Dynamic_temp (jj,indx2(jj)) = NaN; end
            
            base_temp3 = (Static_temp + Dynamic_temp)./2;
            Static_SNR (:,ii,mm) = squeeze (Static_amp(:,ii,mm)) ./ nanmean(base_temp3,2);
            Static_Z (:,ii,mm) = (squeeze (Static_amp(:,ii,mm)) - nanmean(base_temp3,2)) ./ nanstd(base_temp3,[],2);
            Dynamic_SNR (:,ii,mm) = squeeze (Dynamic_amp(:,ii,mm)) ./ nanmean(base_temp3,2);
            Dynamic_Z (:,ii,mm) = (squeeze (Dynamic_amp(:,ii,mm)) - nanmean(base_temp3,2)) ./ nanstd(base_temp3,[],2);
        end
    end
save fft_SNRZ_new  Neutral_SNR  Neutral_Z  Happy_SNR Happy_Z  N2H_SNR N2H_Z...
     H2N_SNR H2N_Z Static_SNR Static_Z Dynamic_SNR Dynamic_Z data_hz EEG

