%% Transform the data into time-frequency domain
% define frequency parameter> define other wavelet parameters> initialize
% output TF data>  main function

% Created by M.-Y. Wang
% 05-12-2019

%% --------------------------Plot an averaged ssVEP (dB) of Occipital electrodes
% 27:30 O1 IZ OZ POZ 64 O2
    figure, clf
    set (gcf,'color','w')
    plot (time2save, squeeze (nanmean(nanmean(Neutral_dB([27:30,64],:,:),4))),'--','color',[.6,.6,.6],'linewidth',2.5);
    hold on
    plot (time2save, squeeze (nanmean(nanmean(Happy_dB([27:30,64],:,:),4))),'--k','linewidth',2.5);
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
    plot (time2save, squeeze (nanmean(nanmean(N2H_dB([27:30,64],:,:),4))),'-b','linewidth',2.5);
    hold on
    plot (time2save, squeeze (nanmean(nanmean(H2N_dB([27:30,64],:,:),4))),'-r','linewidth',2.5);
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
    plot (time2save, squeeze (nanmean(nanmean(Static_dB([27:30,64],:,:),4))),'-','color',[.6,.6,.6],'linewidth',2.5);
    hold on
    plot (time2save, squeeze (nanmean(nanmean(Dynamic_dB([27:30,64],:,:),4))),'-k','linewidth',2.5);
    hold on
    line ([-100,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on
    line ([0,0],[-10,20],'color','k','linewidth',1,'linestyle','--');
    
    set (gca,'xlim',[-100 2500],'ylim',[-6,16],'ytick',-5:5:15,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')
      
