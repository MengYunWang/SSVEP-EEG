%% Transform the data into time-frequency domain
% define frequency parameter> define other wavelet parameters> initialize
% output TF data>  main function

% Created by M.-Y. Wang
% 12-10-2017

%%
clear all
clc
% -------------------------------------Initialize the parameter-------------

data1_name = dir ('F:\face-random\Preprocessing\Conditions\Condition1_Neutral\*.set');
data2_name = dir ('F:\face-random\Preprocessing\Conditions\Condition2_Happy\*.set');
data3_name = dir ('F:\face-random\Preprocessing\Conditions\Condition3_N2H\*.set');
data4_name = dir ('F:\face-random\Preprocessing\Conditions\Condition4_H2N\*.set');
EEG = pop_loadset('filename',data1_name(1).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition1_Neutral');

% frequency parameters
frex =  [10,20]; % define the lowest freq

% other wavelet parameters
range_cycles = [4 6]; % define the cycles; can use the fixed number 3 or 6
s = range_cycles./(2*pi*frex);
wavtime = -2:1/EEG.srate:2;
half_wave = (length(wavtime)-1)/2;
nWave = length(wavtime);

% --------------------------------------------------- Condition1-Neutral ------------------
cd F:\face-random\Preprocessing\Conditions\Condition1_Neutral
% initialize output time-frequency data

% Neutral_tfamp = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name)); % freq * chan * time * subs
Neutral_power = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name));
for ii = 1:length(data1_name);
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition1_Neutral\');
    ERP_data = squeeze(nanmean(EEG.CSD,3));
    % FFT parameters
    nConv = nWave + EEG.pnts - 1;
    nConv_pow2 = 2^nextpow2(nConv);
    alldata = ERP_data;   
    dataX   = fft(alldata,nConv_pow2,2);
    
    for fi=1:length(frex);
        wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
        waveletX = fft(wavelet,nConv_pow2);
        waveletX = waveletX./max(waveletX);
        as = ifft(bsxfun(@times, dataX, waveletX),nConv_pow2,2);
        as = as (:,1:nConv);
        as = as(:,half_wave+1:end-half_wave);
%         Neutral_tfamp (fi,:,:,ii) = abs(as);
        Neutral_power (fi,:,:,ii) = abs(as).^2;
    end
end
%         Neutral_tfamp (:,:,:,[10,15,17,19,21]) = nan;
        Neutral_power (:,:,:,[10,15,17,19,21]) = nan;
% --------------------------------------------------- Condition2-Happy ------------------
cd F:\face-random\Preprocessing\Conditions\Condition2_Happy

% Happy_tfamp = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name));
Happy_power = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name));
for ii = 1:length(data2_name);
    EEG = pop_loadset('filename',data2_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition2_Happy\');
    ERP_data = squeeze(nanmean(EEG.CSD,3));
    % FFT parameters
    nConv = nWave + EEG.pnts - 1;
    nConv_pow2 = 2^nextpow2(nConv);    
    alldata = ERP_data;   
    dataX   = fft(alldata,nConv_pow2,2);
    
    for fi=1:length(frex);
        wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
        waveletX = fft(wavelet,nConv_pow2);
        waveletX = waveletX./max(waveletX);
        as = ifft(bsxfun(@times, dataX, waveletX),nConv_pow2,2);
        as = as (:,1:nConv);
        as = as(:,half_wave+1:end-half_wave);
%         Happy_tfamp (fi,:,:,ii) = abs(as);
        Happy_power (fi,:,:,ii) = abs(as).^2;
    end
end
%         Happy_tfamp (:,:,:,[10,15,17,19,21]) = nan;
        Happy_power (:,:,:,[10,15,17,19,21]) = nan;
% --------------------------------------------------- Condition3-N2H ------------------
cd F:\face-random\Preprocessing\Conditions\Condition3_N2H

% N2H_tfamp = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name));
N2H_power = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name));
for ii = 1:length(data3_name);
    EEG = pop_loadset('filename',data3_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition3_N2H\');
    ERP_data = squeeze(nanmean(EEG.CSD,3)); 
    % FFT parameters
    nConv = nWave + EEG.pnts - 1;
    nConv_pow2 = 2^nextpow2(nConv); 
    alldata = ERP_data;   
    dataX   = fft(alldata,nConv_pow2,2);
    
    for fi=1:length(frex);
        wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
        waveletX = fft(wavelet,nConv_pow2);
        waveletX = waveletX./max(waveletX);
        as = ifft(bsxfun(@times, dataX, waveletX),nConv_pow2,2);
        as = as (:,1:nConv);
        as = as(:,half_wave+1:end-half_wave);
%         N2H_tfamp (fi,:,:,ii) = abs(as);
        N2H_power (fi,:,:,ii) = abs(as).^2;
    end
end
%         N2H_tfamp (:,:,:,[10,15,17,19,21]) = nan;
        N2H_power (:,:,:,[10,15,17,19,21]) = nan;
% --------------------------------------------------- Condition4-H2N------------------
cd F:\face-random\Preprocessing\Conditions\Condition4_H2N

% H2N_tfamp = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name));
H2N_power = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name));
for ii = 1:length(data4_name);
    EEG = pop_loadset('filename',data4_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition4_H2N\');
    ERP_data = squeeze(nanmean(EEG.CSD,3)); 
    % FFT parameters
    nConv = nWave + EEG.pnts - 1;
    nConv_pow2 = 2^nextpow2(nConv); 
    alldata = ERP_data;   
    dataX   = fft(alldata,nConv_pow2,2);
    
    for fi=1:length(frex);
        wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
        waveletX = fft(wavelet,nConv_pow2);
        waveletX = waveletX./max(waveletX);
        as = ifft(bsxfun(@times, dataX, waveletX),nConv_pow2,2);
        as = as (:,1:nConv);
        as = as(:,half_wave+1:end-half_wave);
%         H2N_tfamp (fi,:,:,ii) = abs(as);
        H2N_power (fi,:,:,ii) = abs(as).^2;
    end
end
%         H2N_tfamp (:,:,:,[10,15,17,19,21]) = nan;
        H2N_power (:,:,:,[10,15,17,19,21]) = nan;

        
%--------------------------------------------------------------Static
% Static_tfamp = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name)); % freq * chan * time * subs
Static_power = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name));
for ii = 1:length(data1_name);
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition1_Neutral\');
    ERP_temp1 = squeeze(nanmean(EEG.CSD,3));
    EEG = pop_loadset('filename',data2_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition2_Happy\');
    ERP_temp2 = squeeze(nanmean(EEG.CSD,3));
    
    ERP_data = (ERP_temp1 + ERP_temp2)./2;
    % FFT parameters
    nConv = nWave + EEG.pnts - 1;
    nConv_pow2 = 2^nextpow2(nConv);
    alldata = ERP_data;   
    dataX   = fft(alldata,nConv_pow2,2);
    
    for fi=1:length(frex);
        wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
        waveletX = fft(wavelet,nConv_pow2);
        waveletX = waveletX./max(waveletX);
        as = ifft(bsxfun(@times, dataX, waveletX),nConv_pow2,2);
        as = as (:,1:nConv);
        as = as(:,half_wave+1:end-half_wave);
%         Static_tfamp (fi,:,:,ii) = abs(as);
        Static_power (fi,:,:,ii) = abs(as).^2;
    end
end
%         Static_tfamp (:,:,:,[10,15,17,19,21]) = nan;
        Static_power (:,:,:,[10,15,17,19,21]) = nan;

%--------------------------------------------------------------Dynamic
% Dynamic_tfamp = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name)); % freq * chan * time * subs
Dynamic_power = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name));
for ii = 1:length(data1_name);
    EEG = pop_loadset('filename',data3_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition3_N2H\');
    ERP_temp1 = squeeze(nanmean(EEG.CSD,3));
    EEG = pop_loadset('filename',data4_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition4_H2N\');
    ERP_temp2 = squeeze(nanmean(EEG.CSD,3));
    
    ERP_data = (ERP_temp1 + ERP_temp2)./2;
    % FFT parameters
    nConv = nWave + EEG.pnts - 1;
    nConv_pow2 = 2^nextpow2(nConv);
    alldata = ERP_data;   
    dataX   = fft(alldata,nConv_pow2,2);
    
    for fi=1:length(frex);
        wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
        waveletX = fft(wavelet,nConv_pow2);
        waveletX = waveletX./max(waveletX);
        as = ifft(bsxfun(@times, dataX, waveletX),nConv_pow2,2);
        as = as (:,1:nConv);
        as = as(:,half_wave+1:end-half_wave);
%         Dynamic_tfamp (fi,:,:,ii) = abs(as);
        Dynamic_power (fi,:,:,ii) = abs(as).^2;
    end
end
%         Dynamic_tfamp (:,:,:,[10,15,17,19,21]) = nan;
        Dynamic_power (:,:,:,[10,15,17,19,21]) = nan;        
        
        time2save = EEG.times;
        
% ---------------------------------------- Compute the baseline corrected TF power
baseline_time = [-500, -300];
baseindex = dsearchn (time2save',baseline_time');

baseline_power_sta =( nanmean(Neutral_power(:,:,baseindex(1):baseindex(2),:),3) + nanmean(Happy_power(:,:,baseindex(1):baseindex(2),:),3))./2;
baseline_power_dyn = ( nanmean(N2H_power(:,:,baseindex(1):baseindex(2),:),3) + nanmean(H2N_power(:,:,baseindex(1):baseindex(2),:),3))./2;
baseline_power = (nanmean(Static_power(:,:,baseindex(1):baseindex(2),:),3) + nanmean(Dynamic_power(:,:,baseindex(1):baseindex(2),:),3))./2; 
%         baseline_power_neutral = nanmean (Neutral_power(:,:,baseindex(1):baseindex(2),:),3);
        Neutral_dB = 10.*log10((Neutral_power)./repmat(baseline_power_sta,1,1,length(time2save),1));
        
%         baseline_power_happy = nanmean(Happy_power(:,:,baseindex(1):baseindex(2),:),3);
        Happy_dB = 10.*log10((Happy_power)./repmat(baseline_power_sta,1,1,length(time2save),1));
        
%         baseline_power_n2h = nanmean(N2H_power(:,:,baseindex(1):baseindex(2),:),3);
        N2H_dB = 10.*log10((N2H_power)./repmat(baseline_power_dyn,1,1,length(time2save),1));
        
%         baseline_power_h2n = nanmean(H2N_power(:,:,baseindex(1):baseindex(2),:),3);
        H2N_dB = 10.*log10((H2N_power)./repmat(baseline_power_dyn,1,1,length(time2save),1));
        
%         baseline_power_static = nanmean(Static_power(:,:,baseindex(1):baseindex(2),:),3);
        Static_dB = 10.*log10((Static_power)./repmat(baseline_power,1,1,length(time2save),1));
        
%         baseline_power_dynamic = nanmean(Dynamic_power(:,:,baseindex(1):baseindex(2),:),3);
        Dynamic_dB = 10.*log10((Dynamic_power)./repmat(baseline_power,1,1,length(time2save),1));

save  ssVEP_TF_10+20 frex time2save Neutral_dB Happy_dB N2H_dB  H2N_dB Static_dB Dynamic_dB -v7.3

 