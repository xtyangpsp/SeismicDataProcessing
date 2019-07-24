close all;
load('datafortest.mat');
dt=abs(timeflag(2)-timeflag(1));
din=egfraw(:,1:100);
tmin=-300;
tmax=300;
itmin=round((tmin - timeflag(1))/dt);
itmax=round((tmax - timeflag(1))/dt);


[nsamp,ndata]=size(din);
scale=0.75;

figure('position',[400 400 1050 650]); 
subplot(2,2,[1,3]); hold on;
for i=1:ndata
    plot(timeflag(itmin:itmax),scale*din(itmin:itmax,i)/max(abs(din(itmin:itmax,i))) + i,'k');
%     dtemp=hilbert(din(:,i));
%     plot(timeflag(itmin:itmax),scale*abs(dtemp(itmin:itmax))/max(abs(dtemp(itmin:itmax))) + i,'r');
%     drawnow;
end
dmean=nanmean(din,2);
h1=plot(timeflag(itmin:itmax),scale*dmean(itmin:itmax)/max(abs(dmean(itmin:itmax))) + ndata + 2,'b','linewidth',1);

dmedian=nanmedian(din,2);
h2=plot(timeflag(itmin:itmax),scale*dmedian(itmin:itmax)/max(abs(dmedian(itmin:itmax))) + ndata + 4,'r','linewidth',1);

[drobust,stackweight]=smartstack(din,'stacktype','robust');
h3=plot(timeflag(itmin:itmax),scale*drobust(itmin:itmax)/max(abs(drobust(itmin:itmax))) + ndata + 6,'m','linewidth',1);
xlim([tmin tmax])
hold off;
ylim([0 ndata+8])
legend([h1,h2,h3],'mean','median','robust: entire','location','northwest');
xlabel('cross-correlation time (seconds)')
box on;
axis on;
set(gca,'TickDir','out','fontsize',14)

subplot(2,2,2); hold on;
h1=plot(timeflag(itmin:itmax),scale*dmean(itmin:itmax)/max(abs(dmean(itmin:itmax))) + 2,'b','linewidth',1);

h2=plot(timeflag(itmin:itmax),scale*dmedian(itmin:itmax)/max(abs(dmedian(itmin:itmax))) + 4,'r','linewidth',1);

h3=plot(timeflag(itmin:itmax),scale*drobust(itmin:itmax)/max(abs(drobust(itmin:itmax))) + 6,'m','linewidth',1);
xlim([tmin tmax])
hold off;
ylim([0 8])
legend([h1,h2,h3],'mean','median','robust: entire','location','northwest');
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
legend('rubust: entire','robust: 0-200','robust: 50-150','location','northwest');
xlabel('cross-correlation time (seconds)')
box on;
axis on;
set(gca,'TickDir','out','fontsize',14)