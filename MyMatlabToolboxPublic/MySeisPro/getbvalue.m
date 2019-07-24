function [Mc b sigma]=getbvalue(mag_data,dM,flag_plot,flag_Mc,method)
% USAGE: [Mc b sigma]=getbvalue(mag_data,dM,flag_plot, flag_Mc,'method')
%
% Estimate b-value for frequency-magnitude data using maximum likelihood
% method. 
%
% b-value is calculated following formula 3-9 in Marzocchi and Sandri (ANNALS OF GEOPHYSICS, 2003).
% It is originally proposed by Tinti and Mulargia (1987).
%
% Xiaotao Yang  Mar 20, 2013  @ Indiana University Bloomington
%
% IN:   mag_data        -> magnitude vector data.
%       dM              -> step length for magnitude.
%       flag_plot       -> 0, mute plot; 1, plot data and fit.
%       flag_Mc         -> 0<= Mc <=10.0, fixed/given Mc; -99, find the best-fit Mc.
%       method          -> choose method used for estimating b-value
%           'AU1965':       no bin,   b=1/(log(10)*(mean_mag - Mc))  => Aki, 1965;
%                            Utsu, 1965
%           'UTSU1966':     with bin,   b=1/(log(10)*(mean_mag - (Mc-dM/2)))
%                            =>Utsu,1966
%           'BENDER1983':   with bin,   b=              => Bender, 1983   %
%           NOT READY YET!
%           'TM1987':       with bin,    b=(1+(dM/(mean_mag -
%                           Mc)))/(log(10)*dM)  => Tinti and Mulargia, 1987
% OUT:  Mc              -> Magnitude of completeness, use the best-fit
%                       result.
%       b               -> b-value for the FMD.
%       sigma           -> uncertainties of b-value.


%data=load('dceri_mag.d');
%data=mag_cata_nut_oiink;
%Mc=1.3;
%dM=0.1; %dMc=0.01;
data=mag_data;

switch method
        case {'AU1965'}
            disp( 'Method: Aki & Utsu, 1965');
        case {'UTSU1966'}
            disp('Method: Utsu, 1966');
        case {'TM1987'}
            disp('Method: Tinti & Mulargia, 1987');
        otherwise
            error('***Unknown method, please check input!');
end
    
M=min(data):dM:max(data);
% flag_Mc
if (flag_Mc == -99)

%Mc0=0.1:.01:5;
Mc0=M; %min(data):dM:max(data);

% figure;
% plot(M,log10(cumulative_frequency(data,M)),'^'); hold on
for i=1:length(Mc0)
    
    for j=1:length(M)
        if(M(j) >= Mc0(i))
            start=j;
            break;
        end
    end
    k=1;
    for j=1:length(data)
        if (data(j)>=Mc0(i))
            dcut(k)=data(j);
            k=k+1;
        end
    end
    sm=mean(dcut);          % sample mean of cut data.
    Mcut=M(start:length(M));
%     Mc0(i)
    switch method
        case {'AU1965'}
            b_hat0=1/(log(10)*(sm - Mc0(i)));
            sigma0=b_hat0/sqrt(length(dcut));
        case {'UTSU1966'}
            b_hat0=1/(log(10)*(sm - (Mc0(i)-dM/2)));
            sigma0=b_hat0/sqrt(length(dcut));
        case {'TM1987'}
            p=1+(dM/(sm-Mc0(i)));
            %p0(i)=p;
            b_hat0=log(p)/(log(10)*dM);
            sigma0=(1-p)/(log(10)*dM*(sqrt(length(dcut)*p)));
        otherwise
            error('***Unknown method, please check input!');
    end
    b_hat(i)=b_hat0;
    sigma2(i)=abs(sigma0);
    %
    if (b_hat(i) <0); i= i-1;break;end
    dobs=log10(cumulative_frequency(dcut,Mcut));
    
    max_dobs(i)=dobs(1);
    dpred0=-b_hat(i)*M;
    dpred(i,:)=dpred0+(max_dobs(i)-dpred0(start));
    %dpred=10.^dpred;
    
    %R(i)=abs(100-100*(sum(abs(10.^dobs-10.^dpred))/length(data)));
    R(i)=(dobs-dpred(i,start:length(M)))*(dobs-dpred(i,start:length(M)))';
    clear dcut;
    
%     plot(M,dpred(i,:),'r');
end
% hold off;

xx=Mc0(1:i);
% x_fine=linspace(min(xx),max(xx),1000);
% R_fine=spline(xx,R,x_fine);
%plot(x_fine,R_fine,'.'); 


yy=diff([R R(i)]);
 
temp2=find(yy>=0);
Mc=(xx(temp2(1)-1)+xx(temp2(1)))/2;
%Mc=xx(temp(1));
b=(b_hat(temp2(1)-1)+b_hat(temp2(1)))/2;
% plot(xx,yy,'r'); 

%yy=diff([y y(i)]);
% [V I]=max(yy);

I=temp2(1);

sigma=sigma2(I);

dpred2=-b*M;
dpred_ok=dpred2+(max_dobs(I)-dpred2(I));

if (flag_plot==1)
    figure;
    subplot(2,1,1);
    plot(xx,yy,'o'); hold on;
    plot(xx,yy,'r'); 
    plot(xx,zeros(length(yy),1));hold off;

    subplot(2,1,2);
    plot(M,log10(cumulative_frequency(data,M)),'^'); hold on
    plot(M,dpred_ok,'r'); 
    plot([Mc Mc], [min(dpred_ok) max(dpred_ok)],'k');hold off;
    title(['Mc=',num2str(Mc),'; b=',num2str(b)]);
    
elseif (flag_plot==0)
else
    error('***Wrong plot flag: 0, mute; 1, display the plot. SEE USAGE: help getbvalue');
end

elseif (flag_Mc>=min(data) && flag_Mc<=max(data))  %for fixed Mc given
    Mc=flag_Mc;
    
    for j=1:length(M)
        if(M(j) >= Mc)
            start=j;
            break;
        end
    end
    k=1;
    for j=1:length(data)
        if (data(j)>=Mc)
            dcut0(k)=data(j);
            k=k+1;
        end
    end
    dcut=zeros(k-1);
    dcut=dcut0(1:k-1);
    
    Mcut=M(start:length(M));
    sm=mean(dcut);
    switch method
        case {'AU1965'}
            b_hat=1/(log(10)*(sm - Mc));
            sigma0=b_hat/sqrt(length(dcut));
        case {'UTSU1966'}
            b_hat=1/(log(10)*(sm - (Mc-dM/2)));
            sigma0=b_hat/sqrt(length(dcut));
        case {'TM1987'}
            p=1+(dM/(sm-Mc));
            b_hat=log(p)/(log(10)*dM);
            sigma0=(1-p)/(log(10)*dM*(sqrt(length(dcut)*p)));
        otherwise
            error('***Unknown method, please check input!');
    end
    %
    
    if (b_hat <0); error('***Negative b-value. Check Mc and the data set!');end
    dobs=log10(cumulative_frequency(data,M));
    
    Mc_dobs=interp1(M,dobs,Mc);
    dpred0=-b_hat*M;
    Mc_dpred=interp1(M,dpred0,Mc);
    
    Mc_diff=Mc_dobs-Mc_dpred;
    dpred=dpred0+Mc_diff;
    
    %R(i)=abs(100-100*(sum(abs(10.^dobs-10.^dpred))/length(data)));
    %R=(dobs-dpred(start:length(M)))*(dobs-dpred(start:length(M)))';
    
b=b_hat;
sigma=abs(sigma0);

if (flag_plot==1)
    figure;
    plot(M,log10(cumulative_frequency(data,M)),'^'); hold on
    plot(M,dpred,'r'); 
    plot([Mc Mc],[min(dpred) max(dpred)],'k');hold off;
    title(['Mc=',num2str(Mc),'; b=',num2str(b)]);

elseif (flag_plot==0)
else
    error('***Wrong plot flag: 0, mute; 1, display the plot. SEE USAGE: help getbvalue');
end

else
    error('***Give flag_Mc must be -99 or a value inside the data magnitude range.');
end

return;
end

% Other formulas for b-value calculation. b_hat_2_3 means from formula 2.3 in paper Marzocchi and Sandri (ANNALS OF GEOPHYSICS, 2003).
% b_hat_2_3=1/(log(10)*(mean(data) - Mc))
% 
% b_hat_3_1=1/(log(10)*(mean(data) - (Mc-dM/2)))
% 
% p=1+(dM/(mean(data)-Mc));
% 
% b_hat_3_9=log(p)/(log(10)*dM)

% clear variables.