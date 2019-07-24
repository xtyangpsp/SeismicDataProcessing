function seisplot(t,d,base,fillpolarity,colorlist)
%Modified from: http://image.sciencenet.cn/olddata/kexue.com.cn/upload/blog/file/2010/4/2010422214637890138.pdf
%
%t: time
%d: data.
%base: fill baseline value
%
%fillpolarity: +, -, or +-, or 'N'. This determines which side to fill. If both 
%   this argument and base argument are not given or given as 'N', then the plot 
%   will not fill the wiggles. If base value is specified, then, by default,
%   the program fills positive side.
%
%colorlist: string cell list colors for each side. E.g.,{'r','b'} for
%   fillpolarity of +-. Note that the size of colorlist mush be at least the same as the
%   number of sides to fill, 1 or 2. Default is {'r','b'} for +-.

if(nargin<2)
    error('Not enough arguments.');
elseif nargin==2
    base=nan; %no basevalue
    fillpolarity='N';
elseif nargin==3
    if ~isnan(base)
        fillpolarity='+';
        colorlist={'r'};
    else
        fillpolarity='N';
    end
elseif nargin==4
    if strcmp(fillpolarity,'+')
        colorlist={'r'};
    elseif strcmp(fillpolarity,'-')
        colorlist={'b'};
    elseif strcmp(fillpolarity,'+-')
        colorlist={'r','b'};
    end 
end

n=length(d); 
z1=zeros(1,n);
% baseline=nan(n,1);
% d=d/max(abs(d)); %
% base=0; %
if ~isnan(base) && ~strcmp(fillpolarity,'N')
    z1(1:n)=d+base;
    % baseline(z1==base)=base;
    w1=nan(n,1);
    w2=nan(n,1);
    for i=1:n
        w1(i)=base; 
        w2(i)=base;
        if z1(i) > base
            w1(i)=z1(i); 
        elseif z1(i) < base
            w2(i)=z1(i); 
        end
    end
    w1(1)=base; w1(n)=base; %
    hold on;
    if strcmp(fillpolarity,'+')
        area(t,w1,base,'facecolor',colorlist{1},'edgecolor',colorlist{1},'linewidth',.5);
        plot(t,w2,'-k','linewidth',.25);
    elseif strcmp(fillpolarity,'-')
        area(t,w2,base,'facecolor',colorlist{1},'edgecolor',colorlist{1},'linewidth',.5);
        plot(t,w1,'-k','linewidth',.25);
    elseif strcmp(fillpolarity,'+-') || strcmp(fillpolarity,'-+')
        area(t,w1,base,'facecolor',colorlist{1},'edgecolor',colorlist{1},'linewidth',.5);
        area(t,w2,base,'facecolor',colorlist{2},'edgecolor',colorlist{2},'linewidth',.5);
    else
%         hold off;
        error('Wrong fillpolarity string: should be +, -, or +-.');
    end
    plot(t,zeros(n,1)+base,'-k','linewidth',.5);
    plot(t,z1,'-k','linewidth',.25);
%     hold off;
%     grid on;
else
    hold on;
    plot(t,d,'-k','linewidth',1);
%     hold off;
end
end