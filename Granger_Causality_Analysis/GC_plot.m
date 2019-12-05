%% Plot GC-time domain
clear all
clc
load TGC_100+600+50

%--------------------------------From occipital to frontal
figure,clf
set (gcf,'color','w')
value = spcrv([times2save;squeeze(nanmean(x2y_d))]);
plot(value(1,:),value(2,:),'k','linewidth',3)
hold on;
value = spcrv([times2save;squeeze(nanmean(x2y_s))]);
plot(value(1,:),value(2,:),'color',[0.6 0.6 0.6],'linewidth',3)
% plot (times2save,squeeze(mean(x2y1)))
% plot (times2save,squeeze(mean(y2x1)))
set (gca,'xlim',[-200 2100],'xtick',-500:500:2000,'ylim',[0.055 0.085],'ytick',0:0.01:0.1,'linewidth',3);
set (gca,'FontSize',20,'fontweight','bold','fontname','arial black')


figure,clf
set (gcf,'color','w')
value = spcrv([times2save;squeeze(nanmean(x2y2))]);
plot(value(1,:),value(2,:),'--k','linewidth',3)
hold on
value = spcrv([times2save;squeeze(nanmean(x2y1))]); % smooth the line----function 'spcrv'
plot(value(1,:),value(2,:),'--','color',[0.6 0.6 0.6],'linewidth',3)
hold on;
set (gca,'xlim',[-200 2100],'xtick',-500:500:2000,'ylim',[0.055 0.085],'ytick',0:0.01:0.1,'linewidth',3);
set (gca,'FontSize',20,'fontweight','bold','fontname','arial black')

figure,clf
set (gcf,'color','w')
value = spcrv([times2save;squeeze(nanmean(x2y3))]);
plot(value(1,:),value(2,:),'b','linewidth',3)
hold on;
value = spcrv([times2save;squeeze(nanmean(x2y4))]);
plot(value(1,:),value(2,:),'r','linewidth',3)
% plot (times2save,squeeze(mean(x2y1)))
% plot (times2save,squeeze(mean(y2x1)))
set (gca,'xlim',[-200 2100],'xtick',-500:500:2000,'ylim',[0.035 0.055],'ytick',0:0.01:0.06,'linewidth',3);
set (gca,'FontSize',20,'fontweight','bold','fontname','arial black')


%-----------------------------------Frontal to occipital

figure,clf
set (gcf,'color','w')
value = spcrv([times2save;squeeze(nanmean(y2x_d))]);
plot(value(1,:),value(2,:),'k','linewidth',3)
hold on;
value = spcrv([times2save;squeeze(nanmean(y2x_s))]);
plot(value(1,:),value(2,:),'color',[0.6 0.6 0.6],'linewidth',3)
% plot (times2save,squeeze(mean(x2y1)))
% plot (times2save,squeeze(mean(y2x1)))
set (gca,'xlim',[-200 2100],'xtick',-500:500:2000,'ylim',[0.005 0.025],'ytick',0:0.01:0.03,'linewidth',3);
set (gca,'FontSize',20,'fontweight','bold','fontname','arial black')

figure,clf
set (gcf,'color','w')
value = spcrv([times2save;squeeze(nanmean(y2x2))]);
plot(value(1,:),value(2,:),'--k','linewidth',3)
hold on
value = spcrv([times2save;squeeze(nanmean(y2x1))]);
plot(value(1,:),value(2,:),'--','color',[0.6 0.6 0.6],'linewidth',3)
set (gca,'xlim',[-200 2100],'xtick',-500:500:2000,'ylim',[0.005 0.025],'ytick',0:0.01:0.03,'linewidth',3);
set (gca,'FontSize',20,'fontweight','bold','fontname','arial black')


figure,clf
set (gcf,'color','w')
value = spcrv([times2save;squeeze(nanmean(y2x3))]);
plot(value(1,:),value(2,:),'b','linewidth',3)
hold on
 value = spcrv([times2save;squeeze(nanmean(y2x4))]);
plot(value(1,:),value(2,:),'r','linewidth',3)
set (gca,'xlim',[-200 2100],'xtick',-500:500:2000,'ylim',[0.005 0.025],'ytick',0:0.01:0.03,'linewidth',3);
set (gca,'FontSize',20,'fontweight','bold','fontname','arial black')

