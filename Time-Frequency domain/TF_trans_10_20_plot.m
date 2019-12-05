%% Transform the data into time-frequency domain
% define frequency parameter> define other wavelet parameters> initialize
% output TF data>  main function

% Created by M.-Y. Wang
% 12-10-2017

%% --------------------------Plot an averaged ssVEP amplitude of Occipital electrodes
% % 27:30 O1 IZ OZ POZ 64 O2
% figure, clf
% set (gcf,'color','w')
%     plot (time2save, squeeze (nanmean(nanmean(Neutral_tfamp(1,[27:30,64],:,:),4),2)),'-k','linewidth',2.5);
%     hold on
%     plot (time2save, squeeze (nanmean(nanmean(Happy_tfamp(1,[27:30,64],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',2.5);
%     hold on
%     plot (time2save, squeeze (nanmean(nanmean(N2H_tfamp(1,[27:30,64],:,:),4),2)),'-b','linewidth',2.5);
%     hold on
%     plot (time2save, squeeze (nanmean(nanmean(H2N_tfamp(1,[27:30,64],:,:),4),2)),'-r','linewidth',2.5);
%     hold on
%     line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
%     hold on
%     line ([1700,1700],[0,2],'color','k','linewidth',1,'linestyle','--');
%     set (gca,'xlim',[-500 2500],'ylim',[0,11],'ytick',0:2:10,'linewidth',3);
%     set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
%     ylabel (['\mu','V'],'FontSize',28,'fontweight','bold','fontname','arial black')
% 
%     figure, clf
% set (gcf,'color','w')
%     plot (time2save, squeeze (nanmean(nanmean(Neutral_tfamp(2,[27:30,64],:,:),4),2)),'-k','linewidth',2.5);
%     hold on
%     plot (time2save, squeeze (nanmean(nanmean(Happy_tfamp(2,[27:30,64],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',2.5);
%     hold on
%     plot (time2save, squeeze (nanmean(nanmean(N2H_tfamp(2,[27:30,64],:,:),4),2)),'-b','linewidth',2.5);
%     hold on
%     plot (time2save, squeeze (nanmean(nanmean(H2N_tfamp(2,[27:30,64],:,:),4),2)),'-r','linewidth',2.5);
%     hold on
%     line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
%     hold on
% 
% %     set (gca,'xlim',[-50 2500],'ylim',[0.15,1],'ytick',0:0.2:1.4,'linewidth',3);
%     set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
%     ylabel (['\mu','V'],'FontSize',28,'fontweight','bold','fontname','arial black')
% %     legend Neutral  Happy N2H H2N

%% --------------------------Plot an averaged ssVEP (dB) of Occipital electrodes
% 27:30 O1 IZ OZ POZ 64 O2
    figure, clf
    set (gcf,'color','w')
    plot (time2save, squeeze (nanmean(nanmean(Neutral_dB(2,[27:30,64],:,:),4),2)),'--','color',[.6,.6,.6],'linewidth',2.5);
    hold on
    plot (time2save, squeeze (nanmean(nanmean(Happy_dB(2,[27:30,64],:,:),4),2)),'--k','linewidth',2.5);
    hold on
    line ([-100,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on
    line ([0,0],[-10,20],'color','k','linewidth',1,'linestyle','--');
    
    set (gca,'xlim',[-100 2500],'ylim',[-4,14],'ytick',-5:5:15,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')    

    figure, clf
    set (gcf,'color','w')
    plot (time2save, squeeze (nanmean(nanmean(N2H_dB(1,[27:30,64],:,:),4),2)),'-b','linewidth',2.5);
    hold on
    plot (time2save, squeeze (nanmean(nanmean(H2N_dB(1,[27:30,64],:,:),4),2)),'-r','linewidth',2.5);
    hold on
    line ([-100,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on
    line ([0,0],[-10,20],'color','k','linewidth',1,'linestyle','--');
    
    set (gca,'xlim',[-500 2500],'ylim',[-6,16],'ytick',-5:5:15,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')

    figure, clf
    set (gcf,'color','w')
    plot (time2save, squeeze (nanmean(nanmean(Static_dB(1,[27:30,64],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',2.5);
    hold on
    plot (time2save, squeeze (nanmean(nanmean(Dynamic_dB(1,[27:30,64],:,:),4),2)),'-k','linewidth',2.5);
    hold on
    line ([-100,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on
    line ([0,0],[-10,20],'color','k','linewidth',1,'linestyle','--');
    
    set (gca,'xlim',[-100 2500],'ylim',[-6,16],'ytick',-5:5:15,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')
    

% %     legend Neutral  Happy N2H H2N
%% --------------------------Plot an averaged ssVEP of Occipital electrodes across all conditions
% % 27:30 O1 IZ OZ POZ 64 O2
% figure, clf
%  tf1 = zeros (4,4000); tf2 = zeros (4,4000);
%  tf1(1,:) = squeeze (nanmean(nanmean(Neutral_tfamp(1,[27:30,64],:,:),4),2)); tf2(1,:) = squeeze (nanmean(nanmean(Neutral_tfamp(2,[27:30,64],:,:),4),2));
%  tf1(2,:) = squeeze (nanmean(nanmean(Happy_tfamp(1,[27:30,64],:,:),4),2)); tf2(2,:) = squeeze (nanmean(nanmean(Happy_tfamp(2,[27:30,64],:,:),4),2));
%  tf1(3,:)= squeeze (nanmean(nanmean(N2H_tfamp(1,[27:30,64],:,:),4),2));  tf2(3,:)= squeeze (nanmean(nanmean(N2H_tfamp(2,[27:30,64],:,:),4),2));
%  tf1(4,:) = squeeze (nanmean(nanmean(H2N_tfamp(1,[27:30,64],:,:),4),2)); tf2(4,:) = squeeze (nanmean(nanmean(H2N_tfamp(2,[27:30,64],:,:),4),2));
% set (gcf,'color','w')
%     plot (time2save, nanmean (tf1),'-k','linewidth',3);
%     hold on;
%     plot (time2save, nanmean (tf2),':k','linewidth',3)
% 
% %     set (gca,'xlim',[-200 2500],'ylim',[0.15,1.3],'ytick',0:0.2:1.4,'linewidth',3);
%     set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
%     ylabel (['\mu','V'],'FontSize',28,'fontweight','bold','fontname','arial black')
%  legend 10Hz 20Hz
 %% --------------------------Plot an averaged ssVEP (dB) of Occipital electrodes across all conditions
% 27:30 O1 IZ OZ POZ 64 O2
    figure, clf
    tf1 = zeros (4,800); tf2 = zeros (4,800);
    tf1(1,:) = squeeze (nanmean(nanmean(Neutral_dB(1,[27:30,64],:,:),4),2)); tf2(1,:) = squeeze (nanmean(nanmean(Neutral_dB(2,[27:30,64],:,:),4),2));
    tf1(2,:) = squeeze (nanmean(nanmean(Happy_dB(1,[27:30,64],:,:),4),2)); tf2(2,:) = squeeze (nanmean(nanmean(Happy_dB(2,[27:30,64],:,:),4),2));
    tf1(3,:)= squeeze (nanmean(nanmean(N2H_dB(1,[27:30,64],:,:),4),2));  tf2(3,:)= squeeze (nanmean(nanmean(N2H_dB(2,[27:30,64],:,:),4),2));
    tf1(4,:) = squeeze (nanmean(nanmean(H2N_dB(1,[27:30,64],:,:),4),2)); tf2(4,:) = squeeze (nanmean(nanmean(H2N_dB(2,[27:30,64],:,:),4),2));
    set (gcf,'color','w')
    plot (time2save, nanmean (tf1),'-k','linewidth',3);
    hold on;
    plot (time2save, nanmean (tf2),':k','linewidth',3)

    set (gca,'xlim',[-700 2500],'ylim',[-7,10],'ytick',-5:5:10,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')
%     legend 10Hz 20Hz
    
%% --------------------------Plot an averaged ssVEP amplitude of Frontal electrodes
% % 1 FP1 33 FPZ 34 FP2
% figure, clf
% set (gcf,'color','w')
%     plot (time2save, squeeze (nanmean(nanmean(Neutral_tfamp(1,[1,33,34],:,:),4),2)),'-k','linewidth',2.5);
%     hold on
%     plot (time2save, squeeze (nanmean(nanmean(Happy_tfamp(1,[1,33,34],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',2.5);
%     hold on
%     plot (time2save, squeeze (nanmean(nanmean(N2H_tfamp(1,[1,33,34],:,:),4),2)),'-b','linewidth',2.5);
%     hold on
%     plot (time2save, squeeze (nanmean(nanmean(H2N_tfamp(1,[1,33,34],:,:),4),2)),'-r','linewidth',2.5);
%     hold on
%     line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
% %     hold on
% %     line ([1700,1700],[0,2],'color','k','linewidth',1,'linestyle','--');
% %     set (gca,'xlim',[-50 2500],'ylim',[0.15,1.3],'ytick',0:0.2:1.4,'linewidth',3);
%     set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
%     ylabel (['\mu','V'],'FontSize',28,'fontweight','bold','fontname','arial black')
% 
%     figure, clf
% set (gcf,'color','w')
%     plot (time2save, squeeze (nanmean(nanmean(Neutral_tfamp(2,[1,33,34],:,:),4),2)),'-k','linewidth',2.5);
%     hold on
%     plot (time2save, squeeze (nanmean(nanmean(Happy_tfamp(2,[1,33,34],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',2.5);
%     hold on
%     plot (time2save, squeeze (nanmean(nanmean(N2H_tfamp(2,[1,33,34],:,:),4),2)),'-b','linewidth',2.5);
%     hold on
%     plot (time2save, squeeze (nanmean(nanmean(H2N_tfamp(2,[1,33,34],:,:),4),2)),'-r','linewidth',2.5);
%     hold on
%     line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
%     hold on
% 
% %     set (gca,'xlim',[-50 2500],'ylim',[0.15,1],'ytick',0:0.2:1.4,'linewidth',3);
%     set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
%     ylabel (['\mu','V'],'FontSize',28,'fontweight','bold','fontname','arial black')
% %     legend Neutral  Happy N2H H2N

%% --------------------------Plot an averaged ssVEP (dB) of Frontal electrodes
% 1 FP1 33 FPZ 34 FP2
figure, clf
set (gcf,'color','w')
    plot (time2save, squeeze (nanmean(nanmean(Neutral_dB(1,[1,33,34],:,:),4),2)),'-k','linewidth',2.5);
    hold on
    plot (time2save, squeeze (nanmean(nanmean(Happy_dB(1,[1,33,34],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',2.5);
    hold on
    plot (time2save, squeeze (nanmean(nanmean(N2H_dB(1,[1,33,34],:,:),4),2)),'-b','linewidth',2.5);
    hold on
    plot (time2save, squeeze (nanmean(nanmean(H2N_dB(1,[1,33,34],:,:),4),2)),'-r','linewidth',2.5);
    hold on
    line ([-100,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
%     hold on
%     line ([1700,1700],[0,2],'color','k','linewidth',1,'linestyle','--');
%     set (gca,'xlim',[0 2500],'ylim',[-4,16],'ytick',-5:5:15,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')

    figure, clf
set (gcf,'color','w')
    plot (time2save, squeeze (nanmean(nanmean(Static_dB(1,[1,33,34],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',2.5);
    hold on
    plot (time2save, squeeze (nanmean(nanmean(Dynamic_dB(1,[1,33,34],:,:),4),2)),'-k','linewidth',2.5);
    hold on
    line ([-100,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on

%     set (gca,'xlim',[0 2500],'ylim',[-4,16],'ytick',-5:5:15,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')
%     legend Neutral  Happy N2H H2N
%% --------------------------Plot an averaged ssVEP of Frontal electrodes across all conditions
% % 1 FP1 33 FPZ 34 FP2
% figure, clf
%  tf1 = zeros (4,4000); tf2 = zeros (4,4000);
%  tf1(1,:) = squeeze (nanmean(nanmean(Neutral_tfamp(1,[1,33,34],:,:),4),2)); tf2(1,:) = squeeze (nanmean(nanmean(Neutral_tfamp(2,[1,33,34],:,:),4),2));
%  tf1(2,:) = squeeze (nanmean(nanmean(Happy_tfamp(1,[1,33,34],:,:),4),2)); tf2(2,:) = squeeze (nanmean(nanmean(Happy_tfamp(2,[1,33,34],:,:),4),2));
%  tf1(3,:)= squeeze (nanmean(nanmean(N2H_tfamp(1,[1,33,34],:,:),4),2));  tf2(3,:)= squeeze (nanmean(nanmean(N2H_tfamp(2,[1,33,34],:,:),4),2));
%  tf1(4,:) = squeeze (nanmean(nanmean(H2N_tfamp(1,[1,33,34],:,:),4),2)); tf2(4,:) = squeeze (nanmean(nanmean(H2N_tfamp(2,[1,33,34],:,:),4),2));
% set (gcf,'color','w')
%     plot (time2save, nanmean (tf1),'-k','linewidth',3);
%     hold on;
%     plot (time2save, nanmean (tf2),':k','linewidth',3)
% 
% %     set (gca,'xlim',[-200 2500],'ylim',[0.15,1.3],'ytick',0:0.2:1.4,'linewidth',3);
%     set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
%     ylabel (['\mu','V'],'FontSize',28,'fontweight','bold','fontname','arial black')
%  legend 10Hz 20Hz
 %% --------------------------Plot an averaged ssVEP (dB) of Frontal electrodes across all conditions
% 1 FP1 33 FPZ 34 FP2
    figure, clf
    tf1 = zeros (4,4000); tf2 = zeros (4,4000);
    tf1(1,:) = squeeze (nanmean(nanmean(Neutral_dB(1,[1,33,34],:,:),4),2)); tf2(1,:) = squeeze (nanmean(nanmean(Neutral_dB(2,[1,33,34],:,:),4),2));
    tf1(2,:) = squeeze (nanmean(nanmean(Happy_dB(1,[1,33,34],:,:),4),2)); tf2(2,:) = squeeze (nanmean(nanmean(Happy_dB(2,[1,33,34],:,:),4),2));
    tf1(3,:)= squeeze (nanmean(nanmean(N2H_dB(1,[1,33,34],:,:),4),2));  tf2(3,:)= squeeze (nanmean(nanmean(N2H_dB(2,[1,33,34],:,:),4),2));
    tf1(4,:) = squeeze (nanmean(nanmean(H2N_dB(1,[1,33,34],:,:),4),2)); tf2(4,:) = squeeze (nanmean(nanmean(H2N_dB(2,[1,33,34],:,:),4),2));
    set (gcf,'color','w')
    plot (time2save, nanmean (tf1),'-k','linewidth',3);
    hold on;
    plot (time2save, nanmean (tf2),':k','linewidth',3)

%     set (gca,'xlim',[-200 2500],'ylim',[0.15,1.3],'ytick',0:0.2:1.4,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')
    legend 10Hz 20Hz    