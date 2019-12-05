%% plot the fft results

% Created by M.-Y. Wang
% 25-10-2017

%% Plot amplitude for each condition
clear all
clc

load fft_SNRZ
figure, clf, 
set (gcf,'color','w')
for chani = 1:16;
    subplot (4, 4,chani)
    plot (data_hz, squeeze (nanmean(Neutral_SNR(chani,:,:),3)),'--k','linewidth',2.5);
    hold on
    plot (data_hz, squeeze (nanmean(Happy_SNR(chani,:,:),3)),'-.','color',[.5,.5,.5],'linewidth',2.5);
    hold on
    plot (data_hz, squeeze (nanmean(N2H_SNR(chani,:,:),3)),'-b','linewidth',2.5');
    hold on
    plot (data_hz, squeeze (nanmean(H2N_SNR(chani,:,:),3)),'-r','linewidth',2.5);
        set (gca,'xlim',[5 35],'xtick',5:5:35,'ylim',[1.1 1.3])
    title ([num2str(chani), EEG.chanlocs(chani).labels ])
end
legend Neutral Happy N2H H2N

figure, clf
set (gcf,'color','w')
for chani = 17:32;
    subplot (4, 4,chani-16)
    plot (data_hz, squeeze (nanmean(Neutral_SNR(chani,:,:),3)),'--k','linewidth',2.5);
    hold on
    plot (data_hz, squeeze (nanmean(Happy_SNR(chani,:,:),3)),'-.','color',[.5,.5,.5],'linewidth',2.5);
    hold on
    plot (data_hz, squeeze (nanmean(N2H_SNR(chani,:,:),3)),'-b','linewidth',2.5');
    hold on
    plot (data_hz, squeeze (nanmean(H2N_SNR(chani,:,:),3)),'-r','linewidth',2.5);
    set (gca,'xlim',[5 35],'xtick',5:5:35,'ylim',[1.8 3])
    title ([num2str(chani), EEG.chanlocs(chani).labels ])
end
legend Neutral Happy N2H H2N

figure, clf
set (gcf,'color','w')
for chani = 33:48;
    subplot (4, 4,chani-32)
    plot (data_hz, squeeze (nanmean(Neutral_SNR(chani,:,:),3)),'--k','linewidth',2.5);
    hold on
    plot (data_hz, squeeze (nanmean(Happy_SNR(chani,:,:),3)),'-.','color',[.5,.5,.5],'linewidth',2.5);
    hold on
    plot (data_hz, squeeze (nanmean(N2H_SNR(chani,:,:),3)),'-b','linewidth',2.5');
    hold on
    plot (data_hz, squeeze (nanmean(H2N_SNR(chani,:,:),3)),'-r','linewidth',2.5);
    set (gca,'xlim',[5 35],'xtick',5:5:35,'ylim',[1.1 1.3])
    title ([num2str(chani), EEG.chanlocs(chani).labels ])
end
legend Neutral Happy N2H H2N

figure, clf
set (gcf,'color','w')
for chani = 49:64;
    subplot (4, 4,chani-48)
    plot (data_hz, squeeze (nanmean(Neutral_SNR(chani,:,:),3)),'--k','linewidth',2.5);
    hold on
    plot (data_hz, squeeze (nanmean(Happy_SNR(chani,:,:),3)),'-.','color',[.5,.5,.5],'linewidth',2.5);
    hold on
    plot (data_hz, squeeze (nanmean(N2H_SNR(chani,:,:),3)),'-b','linewidth',2.5');
    hold on
    plot (data_hz, squeeze (nanmean(H2N_SNR(chani,:,:),3)),'-r','linewidth',2.5);
    set (gca,'xlim',[5 35],'xtick',5:5:35,'ylim',[1.5 3])
    title ([num2str(chani), EEG.chanlocs(chani).labels ])
end
legend Neutral Happy N2H H2N
%% plot brian topology across all conditions
temptdata = ( Neutral_SNR + Happy_SNR  + N2H_SNR + H2N_SNR )./4;
data2plot = squeeze (nanmean (temptdata,3));

temptdata_s = ( Neutral_SNR + Happy_SNR )./2;
data2plot_s = squeeze (nanmean (temptdata_s,3));
temptdata_d = ( N2H_SNR + H2N_SNR )./2;
data2plot_d = squeeze (nanmean (temptdata_d,3));

figure , clf
set (gcf,'color','w')
%clo = othercolor ('Paired7');
%topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
% topoplot (data2plot (:,19),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
% topoplot (data2plot_d (:,19),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
topoplot (data2plot_s (:,101),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
set(gca,'clim',[0 3.5])
colorbar([0.85 0.2 0.03 0.5],'ytick',0:1:4,'YTickLabel',{'0','1','2','3','4'},'fontsize',18,'fontweight','bold','fontname','arial black')

% figure, clf
% set (gcf,'color','w')
% % topoplot (data2plot (:,37),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
% % topoplot (data2plot_d (:,37),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
% topoplot (data2plot_s (:,201),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
% set(gca,'clim',[0 2.5])
% colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')
%% plot 10 and 20 brian topology under different conditions
%---------------------------------------10 hz under different conditions
figure , clf
set (gcf,'color','w')
%clo = othercolor ('Paired7');
%topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
topoplot (squeeze (nanmean (Neutral_SNR (:,101,:),3)),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
set(gca,'clim',[0 3.5])
% colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')

figure , clf
set (gcf,'color','w')
%clo = othercolor ('Paired7');
%topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
topoplot (squeeze (nanmean (Happy_SNR (:,101,:),3)),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
set(gca,'clim',[0 3.5])
% colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')

figure , clf
set (gcf,'color','w')
%clo = othercolor ('Paired7');
%topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
topoplot (squeeze (nanmean (N2H_SNR (:,101,:),3)),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
set(gca,'clim',[0 3.5])
% colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')

figure , clf
set (gcf,'color','w')
%clo = othercolor ('Paired7');
%topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
topoplot (squeeze (nanmean (H2N_SNR (:,101,:),3)),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
set(gca,'clim',[0 3.5])
% colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')

figure , clf
set (gcf,'color','w')
%clo = othercolor ('Paired7');
%topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
topoplot (squeeze (nanmean (Static_SNR (:,101,:),3)),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
set(gca,'clim',[0 3.5])
% colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')

figure , clf
set (gcf,'color','w')
%clo = othercolor ('Paired7');
%topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
topoplot (squeeze (nanmean (Dynamic_SNR (:,101,:),3)),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
set(gca,'clim',[0 3.5])
% colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')

% %---------------------------------------20 hz under different conditions
% figure , clf
% set (gcf,'color','w')
% %clo = othercolor ('Paired7');
% %topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
% topoplot (squeeze (nanmean (Neutral_SNR (:,101,:),3)),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
% set(gca,'clim',[0 2.5])
% % colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')
% 
% figure , clf
% set (gcf,'color','w')
% %clo = othercolor ('Paired7');
% %topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
% topoplot (squeeze (nanmean (Happy_SNR (:,101,:),3)),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
% set(gca,'clim',[0 2.5])
% % colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')
% 
% figure , clf
% set (gcf,'color','w')
% %clo = othercolor ('Paired7');
% %topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
% topoplot (squeeze (nanmean (N2H_SNR (:,101,:),3)),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
% set(gca,'clim',[0 2.5])
% % colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')
% 
% figure , clf
% set (gcf,'color','w')
% %clo = othercolor ('Paired7');
% %topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
% topoplot (squeeze (nanmean (H2N_SNR (:,101,:),3)),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
% set(gca,'clim',[0 2.5])
% % colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')

%% Plot SNR of Occipital across all conditions and subjects
figure, clf
 SNR_all = zeros (4,5001);
 SNR_all(1,:) = squeeze (nanmean(nanmean(Neutral_SNR([26:30,63,64],:,:),3),1));
 SNR_all(2,:) = squeeze (nanmean(nanmean(Happy_SNR([26:30,63,64],:,:),3),1));
 SNR_all(3,:)= squeeze (nanmean(nanmean(N2H_SNR([26:30,63,64],:,:),3),1));
 SNR_all(4,:) = squeeze (nanmean(nanmean(H2N_SNR([26:30,63,64],:,:),3),1));
 
 set (gcf,'color','w')
 plot (data_hz,nanmean(SNR_all),'-k','linewidth',3);
 set (gca,'xlim',[5 35],'xtick',0:10:40,'ylim',[0.15 4.85],'ytick',0:1:4,'linewidth',3)
 set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
 title ('Occipital','FontSize',28,'fontweight','bold','fontname','arial black')
%  xlabel ('frequency (Hz)','FontSize',28,'fontweight','bold','fontname','arial black')
 ylabel ('SNR','FontSize',28,'fontweight','bold','fontname','arial black')

%% Plot SNR of Occipital under different conditions
figure, clf
set (gcf,'color','w')
plot (data_hz, squeeze (nanmean(nanmean(Neutral_SNR([26:30,63,64],:,:),3),1)),'-k','linewidth',2.5);
hold on
plot (data_hz, squeeze (nanmean(nanmean(Happy_SNR([26:30,63,64],:,:),3),1)),'-','color',[.6,.6,.6],'linewidth',2.5);
hold on
plot (data_hz, squeeze (nanmean(nanmean(N2H_SNR([26:30,63,64],:,:),3),1)),'-b','linewidth',2.5);
hold on
plot (data_hz, squeeze (nanmean(nanmean(H2N_SNR([26:30,63,64],:,:),3),1)),'-r','linewidth',2.5);
 set (gca,'xlim',[9.5 10.5],'xtick',0:10:40,'ylim',[1.6 2],'ytick',0:1:3,'linewidth',2.5,'box','off')
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title ('Occipital','FontSize',16,'fontweight','bold','fontname','arial black')
% xlabel ('frequency (Hz)','FontSize',16,'fontweight','bold','fontname','arial black')
% ylabel (['Amplitude(','\mu','V)'],'FontSize',16,'fontweight','bold','fontname','arial black')
% legend Neutral Happy N2H H2N
 %% Plot SNR of Frontal across all conditions and subjects
figure, clf
 SNR_all = zeros (4,5001);
 SNR_all(1,:) = squeeze (nanmean(nanmean(Neutral_SNR([1,33,34],:,:),3),1));
 SNR_all(2,:) = squeeze (nanmean(nanmean(Happy_SNR([1,33,34],:,:),3),1));
 SNR_all(3,:)= squeeze (nanmean(nanmean(N2H_SNR([1,33,34],:,:),3),1));
 SNR_all(4,:) = squeeze (nanmean(nanmean(H2N_SNR([1,33,34],:,:),3),1));
 
 set (gcf,'color','w')
 plot (data_hz,nanmean(SNR_all),'-k','linewidth',3);
 set (gca,'xlim',[5 35],'xtick',0:10:30,'ylim',[0.15 4.85],'ytick',0:1:4,'linewidth',3)
 set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
 title ('Frontal','FontSize',28,'fontweight','bold','fontname','arial black')
%  xlabel ('frequency (Hz)','FontSize',28,'fontweight','bold','fontname','arial black')
 ylabel ('SNR','FontSize',28,'fontweight','bold','fontname','arial black')
 %% Plot SNR of frontal under different conditions
figure, clf
set (gcf,'color','w')
plot (data_hz, squeeze (nanmean(nanmean(Neutral_SNR([1,33,34],:,:),3),1)),'-k','linewidth',2.5);
hold on
plot (data_hz, squeeze (nanmean(nanmean(Happy_SNR([1,33,34],:,:),3),1)),'-','color',[.6,.6,.6],'linewidth',2.5);
hold on
plot (data_hz, squeeze (nanmean(nanmean(N2H_SNR([1,33,34],:,:),3),1)),'-b','linewidth',2.5);
hold on
plot (data_hz, squeeze (nanmean(nanmean(H2N_SNR([1,33,34],:,:),3),1)),'-r','linewidth',2.5);
 set (gca,'xlim',[5 25],'xtick',4:2:26,'ylim',[0.75 1.5],'ytick',0:.25:1.5,'linewidth',2.5)
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
% title ('Frontal','FontSize',16,'fontweight','bold','fontname','arial black')
% xlabel ('frequency (Hz)','FontSize',16,'fontweight','bold','fontname','arial black')
% ylabel (['Amplitude(','\mu','V)'],'FontSize',16,'fontweight','bold','fontname','arial black')
% legend Neutral Happy N2H H2N