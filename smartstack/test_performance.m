%
close all;
load('datafortest.mat');
dt=abs(timeflag(2)-timeflag(1));
din=egfraw(:,1:end);
tmin=-300;
tmax=300;
itmin=round((tmin - timeflag(1))/dt);
itmax=round((tmax - timeflag(1))/dt);

[nsamp,ndata]=size(din);
scale=0.95;
reftraceid=77;

% form good pool
pool_good=egfraw(:,[54:71,76:81,122,126:138,146:153]);
pool_bad=egfraw(:,[1:47,83:115,124,140,143,155]);

figure('position',[400 400 1050 650]); 
subplot(2,2,[1 3])
hold on;
for i=1:ndata
    plot(timeflag(itmin:itmax),scale*din(itmin:itmax,i)/max(abs(din(itmin:itmax,i))) + i,'k');
%     dtemp=hilbert(din(:,i));
%     plot(timeflag(itmin:itmax),scale*abs(dtemp(itmin:itmax))/max(abs(dtemp(itmin:itmax))) + i,'r');
%     drawnow;
end
xlim([tmin tmax])
hold off;
ylim([0 ndata+1])
title('all');
xlabel('cross-correlation time (seconds)')
box on;
axis on;
set(gca,'TickDir','out','fontsize',14)

subplot(2,2,2)
hold on;
for i=1:size(pool_good,2)
    plot(timeflag(itmin:itmax),scale*pool_good(itmin:itmax,i)/max(abs(pool_good(itmin:itmax,i))) + i,'k');
%     dtemp=hilbert(din(:,i));
%     plot(timeflag(itmin:itmax),scale*abs(dtemp(itmin:itmax))/max(abs(dtemp(itmin:itmax))) + i,'r');
%     drawnow;
end
xlim([tmin tmax])
hold off;
ylim([0 size(pool_good,2)+1])
title('good');
xlabel('cross-correlation time (seconds)')
box on;
axis on;
set(gca,'TickDir','out','fontsize',14)

subplot(2,2,4)
hold on;
for i=1:size(pool_bad,2)
    plot(timeflag(itmin:itmax),scale*pool_bad(itmin:itmax,i)/max(abs(pool_bad(itmin:itmax,i))) + i,'k');
%     dtemp=hilbert(din(:,i));
%     plot(timeflag(itmin:itmax),scale*abs(dtemp(itmin:itmax))/max(abs(dtemp(itmin:itmax))) + i,'r');
%     drawnow;
end
xlim([tmin tmax])
hold off;
ylim([0 size(pool_bad,2)+1])
title('bad');
xlabel('cross-correlation time (seconds)')
box on;
axis on;
set(gca,'TickDir','out','fontsize',14)

saveas(gca,'smartstack_example_pooled.png','png')
%% 
testpoolsize=size(pool_good,2);
figure('position',[400 400 1050 850]); 
ccmean=nan(size(pool_good,2),1);
ccmedian=nan(size(pool_good,2),1);
ccrobust=nan(size(pool_good,2),1);
stackweightall=nan(testpoolsize,size(pool_good,2));
refstack=smartstack(pool_good);
% loop to generate test pool
for i=1:size(pool_good,2)
    clear idx_good idx_bad
    idx_good=randperm(size(pool_good,2),i);
    clear testpool drobust stackweight
    if i< size(pool_good,2)
        idx_bad=randperm(size(pool_bad,2),testpoolsize - i);
        testpool=[pool_good(:,idx_good),pool_bad(:,idx_bad)];
    else
        testpool=pool_good(:,idx_good);
    end
    
    subplot(3,3,[1 4]);
    hold on;
    dstack=smartstack(testpool,'stacktype','mean');
    ctemp=corrcoef(refstack,dstack);
    ccmean(i)=ctemp(1,2);
    if i==1
        title('mean');
        xlim([tmin tmax])
        ylim([0 testpoolsize+1])
        xlabel('cross-correlation time (seconds)')
        ylabel('number of good trace')
        box on;
        axis on;
        set(gca,'TickDir','out','fontsize',14)
    end
    plot(timeflag(itmin:itmax),scale*dstack(itmin:itmax)/max(abs(dstack(itmin:itmax))) + i,'k','linewidth',1);
    if i==size(pool_good,2);hold off;end
    
    subplot(3,3,[2 5]);
    hold on;
    dstack=smartstack(testpool,'stacktype','median');
    ctemp=corrcoef(refstack,dstack);
    ccmedian(i)=ctemp(1,2);
    if i==1
        title('median')
        xlim([tmin tmax])
        ylim([0 testpoolsize+1])
        xlabel('cross-correlation time (seconds)')
        ylabel('number of good trace')
        box on;
        axis on;
        set(gca,'TickDir','out','fontsize',14)
    end
    plot(timeflag(itmin:itmax),scale*dstack(itmin:itmax)/max(abs(dstack(itmin:itmax))) + i,'k','linewidth',1);
    if i==size(pool_good,2);hold off;end
    
    subplot(3,3,[3 6]);
    hold on;
    [dstack, stackweight]=smartstack(testpool,'stacktype','robust');
    stackweightall(:,i)=sort(stackweight);
    ctemp=corrcoef(refstack,dstack);
    ccrobust(i)=ctemp(1,2);
    if i==1
        title('robust stack')
        xlim([tmin tmax])
        ylim([0 testpoolsize+1])
        xlabel('cross-correlation time (seconds)')
        ylabel('number of good trace')
        box on;
        axis on;
        set(gca,'TickDir','out','fontsize',14)
    end
    plot(timeflag(itmin:itmax),scale*dstack(itmin:itmax)/max(abs(dstack(itmin:itmax))) + i,'k','linewidth',1);
    if rem(i,2)==0
        text(tmax-0.2*abs(tmax),i,num2str(max(abs(dstack(itmin:itmax)))),'color','r','fontsize',10);
    else
        text(tmin+0.05*abs(tmin),i,num2str(max(abs(dstack(itmin:itmax)))),'color','r','fontsize',10);
    end
    if i==size(pool_good,2);hold off;end
    
end

%plot stack weight
subplot(3,3,7)
hold on;
title('sorted stack weight');
plot(stackweightall(:,[1 5 15 30 46]));
hold off;
xlim([0 testpoolsize+1])
ylabel('stack weight')
xlabel('sorted order')
legend('1/46','5/46','15/46','30/46','46/46');
box on;
axis on;
set(gca,'TickDir','out','fontsize',14)

%plot convergence rate
subplot(3,3,[8 9])
hold on;
title('convergence rate');
plot(ccmean,'ko-');
plot(ccmedian,'b^-');
plot(ccrobust,'rs-');
hold off;
xlim([0 testpoolsize+1])
ylabel('corrcoef with refstack')
xlabel('number of good trace')
box on;
axis on;
legend('mean','median','robust');
set(gca,'TickDir','out','fontsize',14)
saveas(gca,'smartstack_example_pooltest.png','png')