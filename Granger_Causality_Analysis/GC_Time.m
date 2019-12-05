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
timewin = 600; % Granger prediction parametersin ms
timewin_points = round(timewin/(1000/EEG.srate));% convert parameters to indices

times2save = -400:50:2200; % temporal down-sample results (but not data!)in ms
times2saveidx = dsearchn(EEG.times' ,times2save');% convert requested times to indices

% % initialize
% bic1 = nan(length(data1_name),length(times2save),50); % Bayes info criteria (hard-coded to order=40)
% data_station1 = nan (length(data1_name),length(times2save));
% 
% for ii = 1 :length(data1_name);
%     EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition1_Neutral\');
%     EEG = pop_resample( EEG, 200);
%     EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
%     data2go = zeros(2,EEG.pnts,EEG.trials);
%     data2go(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
%     data2go(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
%  % %---------------------------------------Preprocessing
%     eegdata = bsxfun(@minus,data2go(:,:,:),nanmean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
%     for timei=1:length(times2save);
%         % data from all trials in this time window
%         tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
%         
%         % detrend and zscore all data
%         for triali=1:size(tempdata,3)
%             tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
%             tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
%         end
%         
%         % check covariance stationarity
%         unit_root = cca_check_cov_stat_mtrial(tempdata,EEG.trials,timewin_points,[]);% channels * trials
%         data_station1 (ii,timei) = sum (sum(unit_root));
%         % reshape tempdata for armorf
%         tempdata = reshape(tempdata,2,(timewin_points)*EEG.trials);
%   % %----------------------------------------Get the order
%         % test bic1 for optimal model order at each time point
%         for bici=1:size(bic1,3)
%             % run model
%             [~,E] = armorf(tempdata,EEG.trials,timewin_points,bici);
%             % compute Bayes Information Criteria
%             bic1(ii,timei,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
%         end
%     end
% end
%
order = 100;% ms
order_points   = round(order/(1000/EEG.srate));
[x2y1,y2x1] = deal(nan(length(data1_name),length(times2save))); % the function deal assigns inputs to all outputs

for ii = 1 :length(data1_name);
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition1_Neutral\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);    
    data2go = nan(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
 % %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),nanmean(data2go(:,:,:),3)); % remove ERP from selected electrodes to improve stationarity
    for timei=1:length(times2save);
        % data from all trials in this time window
        tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        
        % detrend and zscore all data
        for triali=1:size(tempdata,3)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
  % %---------------------------------------main procedure-------------------------
        
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,timewin_points,order_points);
        [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,timewin_points,order_points);
        [Axy,E] = armorf(tempdata     ,EEG.trials,timewin_points,order_points);
        
        % time-domain causal  estimate
        y2x1(ii,timei)=log(Ex/E(1,1));
        x2y1(ii,timei)=log(Ey/E(2,2));
    end
end

%% ---------------------------------------------------------Condition 2

data2_name = dir ('F:\face-random\Preprocessing\Conditions\Condition2_Happy\*.set');
% initializeg
% bic2 = nan(length(data2_name),length(times2save),50); % Bayes info criteria (hard-coded to order=50)
% data_station2 = nan (length(data2_name),length(times2save));
[x2y2,y2x2] = deal(nan(length(data2_name),length(times2save))); % the function deal assigns inputs to all outputs

for ii = 1:length(data2_name);
    EEG = pop_loadset('filename',data2_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition2_Happy\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]); 
    data2go = zeros(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
    %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),nanmean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    for timei=1:length(times2save)
        % data from all trials in this time window
        tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        
        % detrend and zscore all data
        for triali=1:size(tempdata,3)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        
%         % check covariance stationarity
%         unit_root = cca_check_cov_stat_mtrial(tempdata,EEG.trials,timewin_points,[]);
%         data_station2 (ii,timei) = sum (sum(unit_root));
%         % reshape tempdata for armorf
%         tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
%         %----------------------------------------Get the order
%         % test bic2 for optimal model order at each time point
%         for bici=1:size(bic2,3)
%             % run model
%             [~,E] = armorf(tempdata,EEG.trials,timewin_points,bici);
%             % compute Bayes Information Criteria
%             bic2(ii,timei,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
%         end
        %---------------------------------------main procedure-------------------------
        
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,timewin_points,order_points);
        [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,timewin_points,order_points);
        [Axy,E] = armorf(tempdata     ,EEG.trials,timewin_points,order_points);
        
        % time-domain causal  estimate
        y2x2(ii,timei)=log(Ex/E(1,1)); % Front to Occi
        x2y2(ii,timei)=log(Ey/E(2,2)); % Occi to Front
    end
end

%% ------------------------------------------------------------------Condition 3
data3_name = dir ('F:\face-random\Preprocessing\Conditions\Condition3_N2H\*.set');
% initialize
% bic3 = nan(length(data3_name),length(times2save),50); % Bayes info criteria (hard-coded to order=50)
% data_station3 = nan (length(data3_name),length(times2save));
[x2y3,y2x3] = deal(nan(length(data3_name),length(times2save))); % the function deal assigns inputs to all outputs

for ii = 1:length(data3_name);
    EEG = pop_loadset('filename',data3_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition3_N2H\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]); 
    data2go = nan(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
%---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),nanmean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    for timei=1:length(times2save)
        % data from all trials in this time window
        tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        
        % detrend and zscore all data
        for triali=1:size(tempdata,3)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        
%         % check covariance stationarity
%         unit_root = cca_check_cov_stat_mtrial(tempdata,EEG.trials,timewin_points,[]);
%         data_station3 (ii,timei) = sum (sum(unit_root));
%         % reshape tempdata for armorf
%         tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
% %----------------------------------------Get the order
%         % test bic2 for optimal model order at each time point
%         for bici=1:size(bic3,3)
%             % run model
%             [~,E] = armorf(tempdata,EEG.trials,timewin_points,bici);
%             % compute Bayes Information Criteria
%             bic3(ii,timei,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
%         end
%---------------------------------------main procedure-------------------------
        
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,timewin_points,order_points);
        [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,timewin_points,order_points);
        [Axy,E] = armorf(tempdata     ,EEG.trials,timewin_points,order_points);
        
        % time-domain causal  estimate
        y2x3(ii,timei)=log(Ex/E(1,1));
        x2y3(ii,timei)=log(Ey/E(2,2));        
    end
end 
%% --------------------------------------------------------------------Condition 4
data4_name = dir ('F:\face-random\Preprocessing\Conditions\Condition4_H2N\*.set');
% initialize
% bic4 = nan(length(data4_name),length(times2save),50); % Bayes info criteria (hard-coded to order=50)
% data_station4 = nan (length(data4_name),length(times2save));
[x2y4,y2x4] = deal(nan(length(data4_name),length(times2save))); % the function deal assigns inputs to all outputs

for ii = 1:length(data4_name);
    EEG = pop_loadset('filename',data4_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition4_H2N\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]); 
    data2go = nan(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
 %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),nanmean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    for timei=1:length(times2save)
        % data from all trials in this time window
        tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        
        % detrend and zscore all data
        for triali=1:size(tempdata,3)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        
%         % check covariance stationarity
%         unit_root = cca_check_cov_stat_mtrial(tempdata,EEG.trials,timewin_points,[]);
%         data_station4 (ii,timei) = sum (sum(unit_root));
%         % reshape tempdata for armorf
%         tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
% %----------------------------------------Get the order
%         % test bic2 for optimal model order at each time point
%         for bici=1:size(bic4,3)
%             % run model
%             [~,E] = armorf(tempdata,EEG.trials,timewin_points,bici);
%             % compute Bayes Information Criteria
%             bic4(ii,timei,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
%         end
%----------------------------------------main procedure
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,timewin_points,order_points);
        [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,timewin_points,order_points);
        [Axy,E] = armorf(tempdata     ,EEG.trials,timewin_points,order_points);

        % time-domain causal  estimate
        y2x4(ii,timei)=log(Ex/E(1,1));
        x2y4(ii,timei)=log(Ey/E(2,2));
    end
end
%% --------------------------------------------------Static

% initialize
%     bic_s = nan(length(data1_name),length(times2save),50); % Bayes info criteria (hard-coded to order=50)
%     data_station_s = nan (length(data1_name),length(times2save));
    [x2y_s,y2x_s] = deal(nan(length(data1_name),length(times2save))); % the function deal assigns inputs to all outputs

for ii = 1 :length(data1_name);
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition1_Neutral\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]); 
    data2go1 = nan(2,EEG.pnts,EEG.trials);
    data2go1(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go1(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
    EEG = pop_loadset('filename',data2_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition2_Happy\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]); 
    data2go2 = zeros(2,EEG.pnts,EEG.trials);
    data2go2(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go2(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
    if size (data2go1,3) > size (data2go2,3)
       data2go1(:,:,size(data2go2,3)+1:size(data2go1,3)) = [];
    elseif size (data2go1,3) < size (data2go2,3)
       data2go2(:,:,size(data2go1,3)+1:size(data2go2,3)) = [];
    end
    
    data2go = (data2go1 + data2go2)./2;
    
 % %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),nanmean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    for timei=1:length(times2save);
        % data from all trials in this time window
        tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        
        % detrend and zscore all data
        for triali=1:size(tempdata,3)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        
%         % check covariance stationarity
%         unit_root = cca_check_cov_stat_mtrial(tempdata,size (data2go,3),timewin_points,[]);% channels * trials
%         data_station_s (ii,timei) = sum (sum(unit_root));
%         % reshape tempdata for armorf
%         tempdata = reshape(tempdata,2,timewin_points*size (data2go,3));
%   % %----------------------------------------Get the order
%         % test bic1 for optimal model order at each time point
%         for bici=1:size(bic_s,3)
%             % run model
%             [~,E] = armorf(tempdata,size(data2go,3),timewin_points,bici);
%             % compute Bayes Information Criteria
%             bic_s(ii,timei,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
%         end
     % %---------------------------------------main procedure-------------------------
        
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),size(data2go,3),timewin_points,order_points);
        [Ay,Ey] = armorf(tempdata(2,:),size(data2go,3),timewin_points,order_points);
        [Axy,E] = armorf(tempdata     ,size(data2go,3),timewin_points,order_points);
        
        % time-domain causal  estimate
        y2x_s(ii,timei)=log(Ex/E(1,1));
        x2y_s(ii,timei)=log(Ey/E(2,2));    
    end
end
%% ----------------------------------------------------Dynamic

% initialize
% bic_d = nan(length(data1_name),length(times2save),50); % Bayes info criteria (hard-coded to order=50)
% data_station_d = nan (length(data1_name),length(times2save));
[x2y_d,y2x_d] = deal(nan(length(data1_name),length(times2save))); % the function deal assigns inputs to all outputs

for ii = 1 :length(data1_name);
    EEG = pop_loadset('filename',data3_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition3_N2H\');
    EEG = pop_resample( EEG, 200);
    EEG.CSDr = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]); 
    data2go1 = nan(2,EEG.pnts,EEG.trials);
    data2go1(1,:,:) = squeeze (nanmean(EEG.CSDr([27:30,64],:,:))); %occipital
    data2go1(2,:,:) = squeeze (nanmean(EEG.CSDr([1,33,34],:,:))); %frontal
    EEG = pop_loadset('filename',data4_name(ii).name,'filepath','F:\face-random\Preprocessing\Conditions\Condition4_H2N\');
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
    
 % %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),nanmean(data2go(:,:,:),3)); % remove ERP from selected electrodes to improve stationarity
    for timei=1:length(times2save);
        % data from all trials in this time window
        tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        
        % detrend and zscore all data
        for triali=1:size(tempdata,3)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
%                % check covariance stationarity
%         unit_root = cca_check_cov_stat_mtrial(tempdata,size (data2go,3),timewin_points,[]);% channels * trials
%         data_station_d (ii,timei) = sum (sum(unit_root));
%         
%         % reshape tempdata for armorf
%         tempdata = reshape(tempdata,2,timewin_points*size(data2go,3));
%           % %----------------------------------------Get the order
%         % test bic1 for optimal model order at each time point
%         for bici=1:size(bic_d,3)
%             % run model
%             [~,E] = armorf(tempdata,size(data2go,3),timewin_points,bici);
%             % compute Bayes Information Criteria
%             bic_d(ii,timei,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
%         end
  % %---------------------------------------main procedure-------------------------
        
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),size(data2go,3),timewin_points,order_points);
        [Ay,Ey] = armorf(tempdata(2,:),size(data2go,3),timewin_points,order_points);
        [Axy,E] = armorf(tempdata     ,size(data2go,3),timewin_points,order_points);
        
        % time-domain causal  estimate
        y2x_d(ii,timei)=log(Ex/E(1,1));
        x2y_d(ii,timei)=log(Ey/E(2,2));
    end
end

x2y1 ([10 15 17 19 21],:) = nan; x2y2 ([10 15 17 19 21],:) = nan; 
x2y3 ([10 15 17 19 21],:) = nan; x2y4 ([10 15 17 19 21],:) = nan; 
x2y_s ([10 15 17 19 21],:) = nan; x2y_d ([10 15 17 19 21],:) = nan; 

y2x1 ([10 15 17 19 21],:) = nan; y2x2 ([10 15 17 19 21],:) = nan;
y2x3 ([10 15 17 19 21],:) = nan; y2x4 ([10 15 17 19 21],:) = nan;
y2x_s ([10 15 17 19 21],:) = nan; y2x_d ([10 15 17 19 21],:) = nan;

%% Baseline Correction
    baseline_wind = [-200 0];
    [~,baseline_idx1] = find (times2save == -200);
    [~,baseline_idx2] = find (times2save == 0);
% 
    base1 = (x2y1(:,baseline_idx1:baseline_idx2) + x2y2 (:,baseline_idx1:baseline_idx2))./2;
    granger_bs1 = nanmean (base1,2);
    
    base2 = (x2y3(:,baseline_idx1:baseline_idx2) + x2y4 (:,baseline_idx1:baseline_idx2))./2;
    granger_bs2 = nanmean (base2,2);
    
    base3 = (x2y_s(:,baseline_idx1:baseline_idx2) + x2y_d (:,baseline_idx1:baseline_idx2))./2;
    granger_bs3 = nanmean (base3,2);
    
    base4 = (y2x1(:,baseline_idx1:baseline_idx2) + y2x2 (:,baseline_idx1:baseline_idx2))./2;
    granger_bs4 = nanmean (base4,2);
    
    base5 = (y2x3(:,baseline_idx1:baseline_idx2) + y2x4 (:,baseline_idx1:baseline_idx2))./2;
    granger_bs5 = nanmean (base5,2);
    
    base6 = (y2x_s(:,baseline_idx1:baseline_idx2) + y2x_d (:,baseline_idx1:baseline_idx2))./2;
    granger_bs6 = nanmean (base6,2);
    
    x2y1_bs = bsxfun(@minus, x2y1, granger_bs1); x2y2_bs = bsxfun(@minus, x2y2, granger_bs1);
    x2y3_bs = bsxfun(@minus, x2y3, granger_bs2); x2y4_bs = bsxfun(@minus, x2y4, granger_bs2);
    x2y_s_bs = bsxfun(@minus, x2y_s, granger_bs3); x2y_d_bs = bsxfun(@minus, x2y_d, granger_bs3);
    
    y2x1_bs = bsxfun(@minus, y2x1, granger_bs4); y2x2_bs = bsxfun(@minus, y2x2, granger_bs4);
    y2x3_bs = bsxfun(@minus, y2x3, granger_bs5); y2x4_bs = bsxfun(@minus, y2x4, granger_bs5); 
    y2x_s_bs = bsxfun(@minus, y2x_s, granger_bs6); y2x_d_bs = bsxfun(@minus, y2x_d, granger_bs6); 
%%
    save TGC_100+600+50  times2save x2y1 y2x1 x2y2 y2x2 x2y3 y2x3 x2y4 y2x4 x2y_s y2x_s x2y_d y2x_d...
    x2y1_bs y2x1_bs x2y2_bs y2x2_bs x2y3_bs y2x3_bs x2y4_bs y2x4_bs x2y_s_bs y2x_s_bs x2y_d_bs y2x_d_bs
% bic1 bic2 bic3 bic4 bic_s bic_d ...
%         data_station1 data_station2 data_station3 data_station4 data_station_s data_station_d...