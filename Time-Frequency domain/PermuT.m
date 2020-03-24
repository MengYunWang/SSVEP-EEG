%¡¡permutation t tests
% Actual t values > swap two conditions for 1k times and compute the t
% values> get the largest and smallest t value for each time> get a t value
% distribution > decide whether the actual t value is significant or not
%Created by M.-Y. Wang
%30-11-2019

%% permutation t tests

clear all
clc

load ssVEP_TF_10+20.mat
%------------------------------------------------------Neutral vs Happy
    neutral_data = squeeze (nanmean (Neutral_dB (1,[27:30,64],:,:),2));
    happy_data = squeeze (nanmean (Happy_dB (1,[27:30,64],:,:),2));
    neutral_data (:, [10 15 17 19 21]) = [];
    happy_data (:, [10 15 17 19 21]) = [];
    data_differ = happy_data - neutral_data;
    sub_num = size (data_differ, 2);
    t_real1 = nanmean(data_differ,2) ./ (nanstd (data_differ,[],2)./sqrt(sub_num));
    
    npermutation = 1000;
%     t_perm = nan(length(time2save),,npermutation);
    t_max = nan(1000,1); t_min = nan(1000,1);
    for permi = 1:npermutation;
        all_data = cat(2,neutral_data,happy_data);
        rand_idx = randperm (size (all_data,2));
        neutral_perm = all_data(:,rand_idx(1:sub_num));
        happy_perm = all_data(:,rand_idx(sub_num+1:end));
        perm_differ = happy_perm - neutral_perm;
        t_perm = nanmean(perm_differ,2) ./ (nanstd (perm_differ,[],2) ./ sqrt(sub_num));
        t_max(permi) = max(t_perm(701:2801)); % t values of the 0-2100ms
        t_min(permi) = min(t_perm(701:2801));
    end
%     lower_t1 = prctile(t_perm, 2.5,2);
%     upper_t1 = prctile(t_perm,97.5,2);
    lower_t1 = prctile(t_min, 2.5);
    upper_t1 = prctile(t_max,97.5);
%------------------------------------------------------N2H vs H2N
    n2h_data = squeeze (nanmean (N2H_dB (1,[27:30,64],:,:),2));
    h2n_data = squeeze (nanmean (H2N_dB (1,[27:30,64],:,:),2));
    n2h_data (:, [10 15 17 19 21]) = [];
    h2n_data (:, [10 15 17 19 21]) = [];
    data_differ = n2h_data - h2n_data;
    sub_num = size (data_differ, 2);
    t_real2 = nanmean(data_differ,2) ./ (nanstd (data_differ,[],2)./sqrt(sub_num));

%     t_perm = nan(length(time2save),npermutation);
    t_max = nan(1000,1); t_min = nan(1000,1);
    for permi = 1:npermutation;
        all_data = cat(2,n2h_data,h2n_data);
        rand_idx = randperm (size (all_data,2));
        n2h_perm = all_data(:,rand_idx(1:sub_num));
        h2n_perm = all_data(:,rand_idx(sub_num+1:end));
        perm_differ = n2h_perm - h2n_perm;
        t_perm = nanmean(perm_differ,2) ./ (nanstd (perm_differ,[],2) ./ sqrt(sub_num));
        t_max(permi) = max(t_perm(701:2801));
        t_min(permi) = min(t_perm(701:2801));
    end
%     lower_t2 = prctile(t_perm, 2.5,2);
%     upper_t2 = prctile(t_perm,97.5,2);
    lower_t2 = prctile(t_min, 2.5);
    upper_t2 = prctile(t_max,97.5);

%------------------------------------------------------Dynamic vs Static
    static_data = squeeze (nanmean (Static_dB (1,[27:30,64],:,:),2));
    dynamic_data = squeeze (nanmean (Dynamic_dB (1,[27:30,64],:,:),2));
    static_data (:, [10 15 17 19 21]) = [];
    dynamic_data (:, [10 15 17 19 21]) = [];
    data_differ = dynamic_data - static_data;
    sub_num = size (data_differ, 2);
    t_real3 = nanmean(data_differ,2) ./ (nanstd (data_differ,[],2)./sqrt(sub_num));
    
%     t_perm = nan(length(time2save),npermutation);
    t_max = nan(1000,1); t_min = nan(1000,1);
    for permi = 1:npermutation;
        all_data = cat(2,static_data,dynamic_data);
        rand_idx = randperm (size (all_data,2));
        static_perm = all_data(:,rand_idx(1:sub_num));
        dynamic_perm = all_data(:,rand_idx(sub_num+1:end));
        perm_differ = dynamic_perm - static_perm;
%         t_perm (:,permi)= nanmean(perm_differ,2) ./ (nanstd (perm_differ,[],2) ./ sqrt(sub_num));
        t_perm = nanmean(perm_differ,2) ./ (nanstd (perm_differ,[],2) ./ sqrt(sub_num));
        t_max(permi) = max(t_perm(701:end-400));
        t_min(permi) = min(t_perm(701:end-400));
    end
%     lower_t3 = prctile(t_perm, 2.5,2);    
%     upper_t3 = prctile(t_perm,97.5,2);
    lower_t3 = prctile(t_min, 2.5);
    upper_t3 = prctile(t_max,97.5);                             
    
    save permu_t t_real1 t_real2 t_real3 lower_t1 lower_t2 lower_t3...
        upper_t1 upper_t2 upper_t3
%%
    figure, clf
    set (gcf,'color','w')
    plot (time2save, t_real1, '-k','linewidth',2.5)
    hold on
    plot ([time2save(1), time2save(end)],[lower_t1,lower_t1],'-b','linewidth',2)
    hold on
    plot ([time2save(1), time2save(end)],[upper_t1,upper_t1],'-r','linewidth',2)
    hold on
    line ([-100,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on
    line ([0,0],[-5,5],'color','k','linewidth',1,'linestyle','--');
    set (gca,'xlim',[-100 2500],'ylim',[-5,5],'ytick',-6:2:6,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('t-values','FontSize',28,'fontweight','bold','fontname','arial black')
    
    figure, clf
    set (gcf,'color','w')
    plot (time2save, t_real2, '-k','linewidth',2.5)
    hold on
    plot ([time2save(1), time2save(end)],[lower_t2, lower_t2],'-b','linewidth',2)
    hold on
    plot ([time2save(1), time2save(end)],[upper_t2, upper_t2],'-r','linewidth',2)
    hold on
    line ([-100,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on
    line ([0,0],[-5,5],'color','k','linewidth',1,'linestyle','--');
    set (gca,'xlim',[-100 2500],'ylim',[-5,5],'ytick',-6:2:6,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('t-values','FontSize',28,'fontweight','bold','fontname','arial black')

    figure, clf
    set (gcf,'color','w')
    plot (time2save, t_real3, '-k','linewidth',2.5)
    hold on
    plot ([time2save(1), time2save(end)],[lower_t3, lower_t3],'-b','linewidth',2)
    hold on
    plot ([time2save(1), time2save(end)],[upper_t3, upper_t3],'-r','linewidth',2)
    hold on
    line ([-100,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on
    line ([0,0],[-5,5],'color','k','linewidth',1,'linestyle','--');
    set (gca,'xlim',[-100 2500],'ylim',[-5,5],'ytick',-6:2:6,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('t-values','FontSize',28,'fontweight','bold','fontname','arial black')


