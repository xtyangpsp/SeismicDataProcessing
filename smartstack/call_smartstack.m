close all;
load('datafortest.mat');
dt=abs(timeflag(2)-timeflag(1));
din=egfraw(:,1:100);
tmin=-300;
tmax=300;
itmin=round((tmin - timeflag(1))/dt);
itmax=round((tmax - timeflag(1))/dt);


[nsamp,ndata]=size(din);
scale=0.95;
reftraceid=77; % the id of one trace with clear signal

figure('position',[400 400 1050 650]); 
subplot(2,2,[1,3]); hold on;
for i=1:ndata
    plot(timeflag(itmin:itmax),scale*din(itmin:itmax,i)/max(abs(din(itmin:itmax,i))) + i,'k');
%     dtemp=hilbert(din(:,i));
%     plot(timeflag(itmin:itmax),scale*abs(dtemp(itmin:itmax))/max(abs(dtemp(itmin:itmax))) + i,'r');
%     drawnow;
end

h0=plot(timeflag(itmin:itmax),scale*din(itmin:itmax,reftraceid)/max(abs(din(itmin:itmax,reftraceid))) + reftraceid,'g','linewidth',2);

dmean=nanmean(din,2);
h1=plot(timeflag(itmin:itmax),scale*dmean(itmin:itmax)/max(abs(dmean(itmin:itmax))) + ndata + 2,'b','linewidth',1);

dmedian=nanmedian(din,2);
h2=plot(timeflag(itmin:itmax),scale*dmedian(itmin:itmax)/max(abs(dmedian(itmin:itmax))) + ndata + 4,'r','linewidth',1);

[drobust,stackweight]=smartstack(din,'stacktype','robust');
h3=plot(timeflag(itmin:itmax),scale*drobust(itmin:itmax)/max(abs(drobust(itmin:itmax))) + ndata + 6,'m','linewidth',1);
xlim([tmin tmax])
hold off;
ylim([0 ndata+8])
legend([h0,h1,h2,h3],'reference','Mean','Median','Robust: entire','location','southeast');
xlabel('cross-correlation time (seconds)')
box on;
axis on;
set(gca,'TickDir','out','fontsize',14)

subplot(2,2,2); hold on;
h0=plot(timeflag(itmin:itmax),scale*din(itmin:itmax,reftraceid)/max(abs(din(itmin:itmax,reftraceid))) + 2,'g','linewidth',2);

h1=plot(timeflag(itmin:itmax),scale*dmean(itmin:itmax)/max(abs(dmean(itmin:itmax))) +4,'b','linewidth',1);

h2=plot(timeflag(itmin:itmax),scale*dmedian(itmin:itmax)/max(abs(dmedian(itmin:itmax))) + 6,'r','linewidth',1);

h3=plot(timeflag(itmin:itmax),scale*drobust(itmin:itmax)/max(abs(drobust(itmin:itmax))) + 8,'m','linewidth',1);
xlim([tmin tmax])
hold off;
ylim([0 10])
% legend([h0,h1,h2,h3],'reference','Mean','Median','Robust: entire','location','northwest');
xlabel('cross-correlation time (seconds)')
box on;
axis on;
set(gca,'TickDir','out','fontsize',14)


subplot(2,2,4); hold on;
plot(timeflag(itmin:itmax),scale*drobust(itmin:itmax)/max(abs(drobust(itmin:itmax))) + 2,'k','linewidth',1);
drawnow;

robusttmin = 0;
robusttmax = 200;
robustitmin=round((robusttmin - timeflag(1))/dt);
robustitmax=round((robusttmax - timeflag(1))/dt);

[drobust,~]=smartstack(din,'stackwindow',[robustitmin,robustitmax],'stacktype','robust');
plot(timeflag(itmin:itmax),scale*drobust(itmin:itmax)/max(abs(drobust(itmin:itmax))) + 4,'b','linewidth',1);
drawnow;
robusttmin = 50;
robusttmax = 150;
robustitmin=round((robusttmin - timeflag(1))/dt);
robustitmax=round((robusttmax - timeflag(1))/dt);

[drobust,~]=smartstack(din,'stackwindow',[robustitmin,robustitmax],'stacktype','robust');
plot(timeflag(itmin:itmax),scale*drobust(itmin:itmax)/max(abs(drobust(itmin:itmax))) + 6,'r','linewidth',1);
drawnow;
xlim([tmin tmax])
hold off;
ylim([0 8])
legend('Rubust: entire','Robust: 0-200','Robust: 50-150','location','northwest');
xlabel('cross-correlation time (seconds)')
box on;
axis on;
set(gca,'TickDir','out','fontsize',14)
saveas(gca,'smartstack_example.png','png')

%% plot spectrum of the traces
myylim=[0.01 2.5];
figure('position',[400 400 950 600])
subplot(2,2,1)
dref=din(:,reftraceid)/max(abs(din(:,reftraceid)));
spectrogram(dref,300,150,600,5,'yaxis');ax=gca;ax.YScale = 'log';
ylim(myylim)
title('Reference trace');
set(gca,'FontSize',14)

subplot(2,2,2)
dmeantemp=dmean/max(abs(dmean));
spectrogram(dmean,300,150,300,5,'yaxis');ax=gca;ax.YScale = 'log';
title('Mean stack');
ylim(myylim)
set(gca,'FontSize',14)

subplot(2,2,3)
dmediantemp=dmedian/max(abs(dmedian));
spectrogram(dmediantemp,300,150,300,5,'yaxis');ax=gca;ax.YScale = 'log';
title('Median');
ylim(myylim)
set(gca,'FontSize',14)

subplot(2,2,4)
drobusttemp=drobust/max(abs(drobust));
spectrogram(drobusttemp,300,150,600,5,'yaxis');ax=gca;ax.YScale = 'log';
title('Robust stack');
ylim(myylim)
set(gca,'FontSize',14)

saveas(gca,'smartstack_example_spectrograms.png','png')