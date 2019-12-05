%% Granger Causality
% 2018-01-23
% Meng-yun Wang

%% ---------------------------------------------------------STATIC FACE
clear all
clc
%% --------------------------------------------------Condition 1--Neutral
data1_name = dir ('F:\face-random\Preprocessing\Conditions\Condition1_Neutral\*.set');
%---------------------------------------innitialize parameter
EEG = pop_loadset('filename',data1_name(1).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition1_Neutral\');
EEG = pop_resample( EEG, 200);
timewin = 2100; % Granger prediction parameters in ms
time_points = round(timewin/(1000/EEG.srate));% convert parameters to indices

times2start = 0; % 
times2startidx = dsearchn(EEG.times' ,times2start);% convert requested times to indices

% initialize
bic1 = nan(length(data1_name),50); % Bayes info criteria (hard-coded to order=50)
data_station1 = nan (length(data1_name),1);

for ii = 1 :length(data1_name);
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition1_Neutral\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    data2go = nan(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
    % %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),nanmean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    tempdata = eegdata(:,times2startidx:times2startidx+time_points,:);
    
    % detrend and zscore all data
    for triali=1:size(tempdata,3)
        tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
        tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
    end
    
    % check covariance stationarity
    unit_root = cca_check_cov_stat_mtrial(tempdata,EEG.trials,time_points+1,[]);% channels * trials
    data_station1 (ii) = sum (sum(unit_root));
    
    % reshape tempdata for armorf
    tempdata = reshape(tempdata,2,(time_points+1)*EEG.trials);
    
    % %----------------------------------------Get the order
    % test bic1 for optimal model order at each time point
    for bici=1:size(bic1,2)
        % run model
        [~,E] = armorf(tempdata,EEG.trials,time_points+1,bici);
        % compute Bayes Information Criteria
        bic1(ii,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
    end
end

%%
order = 100;% ms
order_points   = round(order/(1000/EEG.srate));
[x2y1,y2x1] = deal(nan(length(data1_name),1)); % the function deal assigns inputs to all outputs

for ii = 1 :length(data1_name);
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition1_Neutral\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);    
    data2go = nan(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital > x
    data2go(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal > y
 % %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),nanmean(data2go(:,:,:),3)); % remove ERP from selected electrodes to improve stationarity
   
        % data from all trials in this time window
        tempdata = eegdata(:,times2startidx:times2startidx + time_points,:);
        
        % detrend and zscore all data
        for triali=1:size(tempdata,3)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,(time_points+1)*EEG.trials);
  % %---------------------------------------main procedure-------------------------
        
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,(time_points+1),order_points);
        [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,(time_points+1),order_points);
        [Axy,E] = armorf(tempdata     ,EEG.trials,(time_points+1),order_points);
        
        % time-domain causal  estimate
        y2x1(ii,1)=log(Ex/E(1,1));
        x2y1(ii,1)=log(Ey/E(2,2));
end

%% ---------------------------------------------------------Condition 2

data2_name = dir ('F:\face-random\Preprocessing\Conditions\Condition2_Happy\*.set');
% initialize
bic2 = nan(length(data2_name),50); % Bayes info criteria (hard-coded to order=50)
data_station2 = nan (length(data2_name),1);
[x2y2,y2x2] = deal(nan(length(data2_name),1)); % the function deal assigns inputs to all outputs

for ii = 1:length(data2_name);
    EEG = pop_loadset('filename',data2_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition2_Happy\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);    
    data2go = nan(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (mean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (mean(EEG.CSDr([1,33,34],:,:))); %frontal
    %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),nanmean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
        
    % data from all trials in this time window
        tempdata = eegdata(:,times2startidx:times2startidx + time_points,:);
        
        % detrend and zscore all data
        for triali=1:size(tempdata,3)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        
        % check covariance stationarity
        unit_root = cca_check_cov_stat_mtrial(tempdata,EEG.trials,time_points+1,[]);
        data_station2 (ii,1) = sum (sum(unit_root));
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,(time_points+1)*EEG.trials);
        %----------------------------------------Get the order
        % test bic2 for optimal model order at each time point
        for bici=1:size(bic2,2)
            % run model
            [~,E] = armorf(tempdata,EEG.trials,(time_points+1),bici);
            % compute Bayes Information Criteria
            bic2(ii,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
        end
        %---------------------------------------main procedure-------------------------
        
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,(time_points+1),order_points);
        [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,(time_points+1),order_points);
        [Axy,E] = armorf(tempdata     ,EEG.trials,(time_points+1),order_points);
        
        % time-domain causal  estimate
        y2x2(ii,1)=log(Ex/E(1,1)); % Front to Occi
        x2y2(ii,1)=log(Ey/E(2,2)); % Occi to Front
end

 
%% ------------------------------------------------------------------Condition 3
data3_name = dir ('F:\face-random\Preprocessing\Conditions\Condition3_N2H\*.set');
% initialize
bic3 = nan(length(data3_name),50); % Bayes info criteria (hard-coded to order=50)
data_station3 = nan (length(data3_name),1);
[x2y3,y2x3] = deal(nan(length(data3_name),1)); % the function deal assigns inputs to all outputs

for ii = 1:length(data3_name);
    EEG = pop_loadset('filename',data3_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition3_N2H\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    data2go = nan(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
%---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),nanmean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
        % data from all trials in this time window
        tempdata = eegdata(:,times2startidx:times2startidx+time_points,:);
        
        % detrend and zscore all data
        for triali=1:size(tempdata,3)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
             % check covariance stationarity
        unit_root = cca_check_cov_stat_mtrial(tempdata,EEG.trials,time_points+1,[]);
        data_station3 (ii,1) = sum (sum(unit_root));
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,(time_points+1)*EEG.trials);
%----------------------------------------Get the order
        % test bic2 for optimal model order at each time point
        for bici=1:size(bic3,2)
            % run model
            [~,E] = armorf(tempdata,EEG.trials,(time_points+1),bici);
            % compute Bayes Information Criteria
            bic3(ii,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
        end
%---------------------------------------main procedure-------------------------
        
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,(time_points+1),order_points);
        [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,(time_points+1),order_points);
        [Axy,E] = armorf(tempdata     ,EEG.trials,(time_points+1),order_points);
        
        % time-domain causal  estimate
        y2x3(ii,1)=log(Ex/E(1,1));
        x2y3(ii,1)=log(Ey/E(2,2));
end
%% --------------------------------------------------------------------Condition 4
data4_name = dir ('F:\face-random\Preprocessing\Conditions\Condition4_H2N\*.set');
% initialize
bic4 = nan(length(data4_name),50); % Bayes info criteria (hard-coded to order=15)
data_station4 = nan (length(data4_name),1);
[x2y4,y2x4] = deal(nan(length(data4_name),1)); % the function deal assigns inputs to all outputs

for ii = 1:length(data4_name);
    EEG = pop_loadset('filename',data4_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition4_H2N\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    data2go = zeros(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
 %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),nanmean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
        % data from all trials in this time window
        tempdata = eegdata(:,times2startidx:times2startidx + time_points,:);
        
        % detrend and zscore all data
        for triali=1:size(tempdata,3)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        % check covariance stationarity
        unit_root = cca_check_cov_stat_mtrial(tempdata,EEG.trials,(time_points+1),[]);
        data_station4 (ii,1) = sum (sum(unit_root));
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,(time_points+1)*EEG.trials);
%----------------------------------------Get the order
        % test bic2 for optimal model order at each time point
        for bici=1:size(bic4,2)
            % run model
            [~,E] = armorf(tempdata,EEG.trials,(time_points+1),bici);
            % compute Bayes Information Criteria
            bic4(ii,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
        end
%----------------------------------------main procedure
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,(time_points+1),order_points);
        [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,(time_points+1),order_points);
        [Axy,E] = armorf(tempdata     ,EEG.trials,(time_points+1),order_points);

        % time-domain causal  estimate
        y2x4(ii,1)=log(Ex/E(1,1));
        x2y4(ii,1)=log(Ey/E(2,2));
end

%% ---------------------------------------------------Static Condition
data1_name = dir ('F:\face-random\Preprocessing\Conditions\Condition1_Neutral\*.set');
data2_name = dir ('F:\face-random\Preprocessing\Conditions\Condition2_Happy\*.set');
EEG = pop_loadset('filename',data1_name(1).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition1_Neutral\');

bic_s = nan (length(data1_name),50); % Bayes info criteria (hard-coded to order=50)
data_station_s = nan (length(data1_name),1);
[x2y_s,y2x_s] = deal(nan(length(data1_name),1)); % the function deal assigns inputs to all outputs

for ii = 1:length(data1_name);
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition1_Neutral\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    data2go1 = nan(2,EEG.pnts,EEG.trials);
    data2go1(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go1(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
    EEG = pop_loadset('filename',data2_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition2_Happy\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    data2go2 = nan(2,EEG.pnts,EEG.trials);
    data2go2(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go2(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
    if size (data2go1,3) > size (data2go2,3)
       data2go1(:,:,size(data2go2,3)+1:size(data2go1,3)) = [];
    elseif size (data2go1,3) < size (data2go2,3)
       data2go2(:,:,size(data2go1,3)+1:size(data2go2,3)) = [];
    end
    
    data2go = (data2go1 + data2go2)./2;
    %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),nanmean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
        % data from all trials in this time window
        tempdata = eegdata(:,times2startidx:times2startidx + time_points,:);
        
        % detrend and zscore all data
        for triali=1:size(tempdata,3)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        % check covariance stationarity
        unit_root = cca_check_cov_stat_mtrial(tempdata,size(data2go,3),(time_points+1),[]);
        data_station_s (ii,1) = sum (sum(unit_root));
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,(time_points+1)*size(data2go,3));
%----------------------------------------Get the order
        % test bic2 for optimal model order at each time point
        for bici=1:size(bic_s,2)
            % run model
            [~,E] = armorf(tempdata,size(data2go,3),(time_points+1),bici);
            % compute Bayes Information Criteria
            bic_s(ii,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
        end
%----------------------------------------main procedure
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),size(data2go,3),(time_points+1),order_points);
        [Ay,Ey] = armorf(tempdata(2,:),size(data2go,3),(time_points+1),order_points);
        [Axy,E] = armorf(tempdata     ,size(data2go,3),(time_points+1),order_points);

        % time-domain causal  estimate
        y2x_s(ii,1)=log(Ex/E(1,1));
        x2y_s(ii,1)=log(Ey/E(2,2));
end

%% --------------------------------------------Dynamic conditions

data1_name = dir ('F:\face-random\Preprocessing\Conditions\Condition1_Neutral\*.set');
data2_name = dir ('F:\face-random\Preprocessing\Conditions\Condition2_Happy\*.set');
data3_name = dir ('F:\face-random\Preprocessing\Conditions\Condition3_N2H\*.set');
data4_name = dir ('F:\face-random\Preprocessing\Conditions\Condition4_H2N\*.set');
% initialize
bic_d = nan(length(data3_name),50); % Bayes info criteria (hard-coded to order=50)
data_station_d = nan (length(data3_name),1);
[x2y_d,y2x_d] = deal(nan(length(data3_name),1)); % the function deal assigns inputs to all outputs

for ii = 1:length(data3_name);

    EEG = pop_loadset('filename',data3_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition3_N2H\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    data2go3 = nan(2,EEG.pnts,EEG.trials);
    data2go3(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go3(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
    EEG = pop_loadset('filename',data4_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition4_H2N\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    data2go4 = zeros(2,EEG.pnts,EEG.trials);    
    data2go4(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go4(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
    if size (data2go3,3) > size (data2go4,3)
       data2go3(:,:,size(data2go4,3)+1:size(data2go3,3)) = [];
    elseif size (data2go3,3) < size (data2go4,3)
       data2go4(:,:,size(data2go3,3)+1:size(data2go4,3)) = [];
    end
    data2go = (data2go3 + data2go4)./2;
    
    %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),nanmean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
        % data from all trials in this time window
        tempdata = eegdata(:,times2startidx:times2startidx + time_points,:);
        
        % detrend and zscore all data
        for triali=1:size(tempdata,3)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        % check covariance stationarity
        unit_root = cca_check_cov_stat_mtrial(tempdata,size(data2go,3),(time_points+1),[]);
        data_station_d (ii,1) = sum (sum(unit_root));
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,(time_points+1)*size(data2go,3));
%----------------------------------------Get the order
        % test bic2 for optimal model order at each time point
        for bici=1:size(bic_d,2)
            % run model
            [~,E] = armorf(tempdata,size(data2go,3),(time_points+1),bici);
            % compute Bayes Information Criteria
            bic_d(ii,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
        end
%----------------------------------------main procedure
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),size(data2go,3),(time_points+1),order_points);
        [Ay,Ey] = armorf(tempdata(2,:),size(data2go,3),(time_points+1),order_points);
        [Axy,E] = armorf(tempdata     ,size(data2go,3),(time_points+1),order_points);

        % time-domain causal  estimate
        y2x_d(ii,1)=log(Ex/E(1,1));
        x2y_d(ii,1)=log(Ey/E(2,2));
end
%% --------------------------------------------Resting state

% initialize
timewin_b = 700; % Granger prediction parameters in ms
time_points_b = round(timewin_b/(1000/EEG.srate));% convert parameters to indices

times2start_b = -700; % 
times2startidx_b = dsearchn(EEG.times' ,times2start_b);
[x2y_b,y2x_b] = deal(nan(length(data3_name),1)); % the function deal assigns inputs to all outputs

for ii = 1:length(data3_name);
   EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition1_Neutral\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    data2go1 = nan(2,EEG.pnts,EEG.trials);
    data2go1(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go1(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
    EEG = pop_loadset('filename',data2_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition2_Happy\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    data2go2 = nan(2,EEG.pnts,EEG.trials);
    data2go2(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go2(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
    EEG = pop_loadset('filename',data3_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition3_N2H\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    data2go3 = nan(2,EEG.pnts,EEG.trials);
    data2go3(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go3(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
    EEG = pop_loadset('filename',data4_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition4_H2N\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    data2go4 = zeros(2,EEG.pnts,EEG.trials);    
    data2go4(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go4(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
    [num, ~] = min ([size(data2go1,3), size(data2go2,3), size(data2go3,3),size(data2go4,3)]);
    
    data2go1_temp = data2go1(:,:,1:num);
    data2go2_temp = data2go2(:,:,1:num);
    data2go3_temp = data2go3(:,:,1:num);
    data2go4_temp = data2go4(:,:,1:num);
       
    data2go = (data2go1_temp + data2go2_temp + data2go3_temp + data2go4_temp)./4;
   
    %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),nanmean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
        % data from all trials in this time window
        tempdata = eegdata(:,times2startidx_b:times2startidx_b + time_points_b,:);
        
        % detrend and zscore all data
        for triali=1:size(tempdata,3)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
%----------------------------------------main procedure
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),size(data2go,3),(time_points_b+1),order_points);
        [Ay,Ey] = armorf(tempdata(2,:),size(data2go,3),(time_points_b+1),order_points);
        [Axy,E] = armorf(tempdata     ,size(data2go,3),(time_points_b+1),order_points);

        % time-domain causal  estimate
        y2x_b(ii,1)=log(Ex/E(1,1));
        x2y_b(ii,1)=log(Ey/E(2,2));
end
%%


bic1([10 15 17 19 21],:) = nan; bic2([10 15 17 19 21],:) = nan;
bic3([10 15 17 19 21],:) = nan; bic4([10 15 17 19 21],:) = nan;
bic_s([10 15 17 19 21],:) = nan; bic_d([10 15 17 19 21],:) = nan;

data_station1 ([10 15 17 19 21]) = nan; data_station2 ([10 15 17 19 21]) = nan; 
data_station3 ([10 15 17 19 21]) = nan; data_station4 ([10 15 17 19 21]) = nan; 
data_station_s ([10 15 17 19 21]) = nan; data_station_d ([10 15 17 19 21]) = nan; 


x2y1 ([10 15 17 19 21]) = nan; x2y2 ([10 15 17 19 21]) = nan; 
x2y3 ([10 15 17 19 21]) = nan; x2y4 ([10 15 17 19 21]) = nan; 
x2y_s ([10 15 17 19 21]) = nan; x2y_d ([10 15 17 19 21]) = nan; 
x2y_b ([10 15 17 19 21]) = nan;

x2y1_bc = x2y1 - x2y_b; x2y2_bc = x2y2 - x2y_b;
x2y3_bc = x2y3 - x2y_b; x2y4_bc = x2y4 - x2y_b;
x2y_s_bc = x2y_s - x2y_b; x2y_d_bc = x2y_d - x2y_b;

y2x1 ([10 15 17 19 21]) = nan; y2x2 ([10 15 17 19 21]) = nan;
y2x3 ([10 15 17 19 21]) = nan; y2x4 ([10 15 17 19 21]) = nan;
y2x_s ([10 15 17 19 21]) = nan; y2x_d ([10 15 17 19 21]) = nan;
y2x_b ([10 15 17 19 21]) = nan;

y2x1_bc = y2x1 - y2x_b; y2x2_bc = y2x2 - y2x_b;
y2x3_bc = y2x3 - y2x_b; y2x4_bc = y2x4 - y2x_b;
y2x_s_bc = y2x_s - y2x_b; y2x_d_bc = y2x_d - y2x_b;

save TGC_whole  bic1 bic2 bic3 bic4 bic_s bic_d ...
    data_station1 data_station2 data_station3 data_station4 data_station_s data_station_d...
    x2y1 y2x1 x2y2 y2x2 x2y3 y2x3 x2y4 y2x4 x2y_s y2x_s x2y_d y2x_d x2y_b y2x_b...
    x2y1_bc y2x1_bc x2y2_bc y2x2_bc x2y3_bc y2x3_bc x2y4_bc y2x4_bc x2y_s_bc y2x_s_bc x2y_d_bc y2x_d_bc 


