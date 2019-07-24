function [origin Mout Mout_datetime]=locatebysa(stat_arrivals,ttmodel,T0,SAitn,MCMCitn,bounds_lon,bounds_lat,bounds_depth)
%Locating earthquake in spherical earth using simulated annealing optimization algorithm.
%Calls:
%   [origin Mout Mout_datetime]=locatebysa(stat_arrivals,model,T0,TN,SAitn,MCMCitn,bounds_lat,bounds_lat,bounds_depth)
%Inputs:
%   stat_arrivals:
%               A structure contains a series of station parameters: lon
%               lat elevation date time error. Here we use the first
%               P arrivals.
%   model:
%               Astructure of the earth velocity model. You can create this
%               model using mkreadnd.m from TTBOX (see subroutine descriptions).
%   T0:
%               Initial/starting temperature.
%   SAitn, MCMCitn:
%               Number of iterations for SA loop and for MCMC inner loop.
%   bounds_lon,lat,depth:
%               Bounds for those three parameters.
%
%Outputs:
%   origin:
%               A structure of the earthquake origin: origdate origtime  lon lat
%               depth.
%
%Subroutines used:
%   mkttime.m:
%               Calculate the seismic wave travel times in spherical earth 
%               by giving earth model.From TTBOX25102012 (http://www.dr-knapmeyer.de/downloads/)
%               by Martin Knapmeyer. [tt,p]=mkttime(phase,delta,h,model).
%   mktp.m:
%               Calculate the traveltime in a layer. Used here for the
%               station elevation corrections.
%Xiaotao Yang   @ Indiana University Bloomington
%Versions:
%   Nov 22, 2012 : created.
%   Dec 9, 2012 : use Antelope ttcalc instead of mkttime.
%--------------------------------------------------------------------------

%%%%%%======== Get parameters from the arguments.
% If call is correct"

% If arguments are enough and correct?

%%%%%%======== Set values for constants.
phase='P';
ttmethod='tttaup';
%ttmodel='iasp91';
R_earth=6371;
secondsperday=3600*24;
km2deg=180/(R_earth*pi);
Nmax=200;  % max steps used to find optimized origin time.
dttmax=10^-2;   % max traveltime difference between two optimized points used to find the origin time.
delta_origintime=0.1;  % in days. deviation of origin time used to give the bounds for finding the best origin time.
% This value is used to generate the initial two points for the
% optimization process.
alpha=0.3;   % used in golden section optimization process.
Tdecay=0.7;  %T=T*Tdecay.

L=length(stat_arrivals.lon);
%L_step_lon=(max(bounds_lon) - min(bounds_lon))/2;         % step length used in Metropolis-Hasting random walk.
%L_step_lat=(max(bounds_lat) - min(bounds_lat))/2;
%L_step_depth=(max(bounds_depth) - min(bounds_depth))/2;   % for crustal and upper mantle earthquakes.

%%%%%%======== Initialize variables.
initial_origdate=stat_arrivals.date{1};
initial_origtime=stat_arrivals.time{1};
% initial_depth=rand*(max(bounds_depth)-min(bounds_depth))+min(bounds_depth);
% initial_lon=rand*(max(bounds_lon)-min(bounds_lon))+min(bounds_lon);
% initial_lat=rand*(max(bounds_lat)-min(bounds_lat))+min(bounds_lat);

initial_depth=random('unif',min(bounds_depth),max(bounds_depth));
initial_lon=random('unif',min(bounds_lon),max(bounds_lon));
initial_lat=random('unif',min(bounds_lat),max(bounds_lat));

origin=struct('origdate',...
              'origtime',...
              'lon',[],...
              'lat',[],...
              'depth',[]);

origin.origdate=initial_origdate;
origin.origtime=initial_origtime;
origin.lon=initial_lon;
origin.lat=initial_lat;
origin.depth=initial_depth;


%% Simulated Annealing algorithm.
% (1) cooling schedule. Billings 1994.
%for i=1:SAitn   %this one comes from Kaj's class handout.
%    T(i)=T0-(T0-TN)(SAitn+1)/SAitn+(T0-TN)(SAitn+1)/(SAitn*(i+1));
%end
%T=T0*log10(2)./log10((1:SAitn)+1);
T=T0*Tdecay.^(0:(SAitn - 1));
factor=log10(2)./log10((1:SAitn)+1);

% (2) begin to cool down.
cov=diag(stat_arrivals.error.^2);
invcov=inv(cov);   % inverse of covariance matrix with 1/sigma^2.
Mout=zeros(SAitn,3);
M0=[initial_lon initial_lat initial_depth];

tmp=sprintf('%s,%s',initial_origdate, initial_origtime);
Mout_datetime=zeros(1,SAitn);
M0_datetime=datenum(tmp); % convert date and time to a number.

arrivals_obs=zeros(1,L);
for i=1:L
    tmp=sprintf('%s,%s',stat_arrivals.date{i},stat_arrivals.time{i});
    arrivals_obs(i)=datenum(tmp);
end

for i=1:SAitn
    
    %%%%%% calculate first arrival travel times for all stations.
    fprintf('==>> Cooling step: %d, T= %f\n',i, T(i));
    fprintf('    Initial model: %f  %f  %f  %f\n',M0,M0_datetime);
    
    tmin0=zeros(1,L);
    for j = 1:L
        %fprintf('Processing station: %d\n',j);
        delta=km2deg*geodist(stat_arrivals.lon(j),stat_arrivals.lat(j),M0(1),M0(2));
        %[tt,p]=mkttime(phase,delta,M0(3),model);
        %[tmin0(j) I]=min(tt);
        tmin0(j)=subtocallttcalcsh(ttmethod, ttmodel, phase, delta, M0(3));
        %pmin0(j)=p(I);
    end
    
    % convert travel time to arrivals in UTC time.
    arrivals_pred0=tmin0/secondsperday + M0_datetime;
    
    rr0=secondsperday*(arrivals_obs - arrivals_pred0);

    loglike0=-0.5*rr0*invcov*rr0'/T(i);
    %loglike0=-0.5*rr0*rr0'/T(i);
    %intialize radnom number generator
    %rand('state', sum(100*clock)); %reset random number generator
    
    %define 'previous values'
    Mprev=M0;
    Mprev_datetime=M0_datetime;
    loglikeprev=loglike0;
    %fprintf('loglikeprev out of loop: %f\n',loglikeprev);
    
    NumSamples=MCMCitn;
    
    M_all=zeros(length(M0),NumSamples);
    M_all_datetime=zeros(length(M0_datetime),NumSamples);
    all_loglike=zeros(1,NumSamples);
    
    %fprintf('--> Metropolis-Hasting sampling ...\n');
    loop=1;
    
    num_good=0; num_bad=0;num_bonus=0;
    while loop<=NumSamples
        fprintf('    -->Loop: %d - %d\n',i,loop);
        %take a random step in model space    
%         M = Mprev + rand(size(Mprev))*2*diag([L_step_lon L_step_lat L_step_depth]) ...
%             - [L_step_lon L_step_lat L_step_depth];
        M(1) = Mprev(1) + random('unif',-5*factor(i),5*factor(i));  %for lat and lon: -1 to 1.
        M(2) = Mprev(2) + random('unif',-5*factor(i),5*factor(i));  %for lat and lon: -1 to 1.
        
        M(3)=abs(Mprev(3) + random('unif',-1,1));  % for depth -1 to 1
        
        if M(1)>bounds_lon(2) || M(1)<bounds_lon(1) ||...
            M(2)>bounds_lat(2)|| M(2)<bounds_lat(1) ||...
            M(3)>bounds_depth(2)|| M(3)<bounds_depth(1)
            fprintf('      Model is out of bound. Dropped!\n');
            
        else
        
            tmin=zeros(1,L);
            for j = 1:L
                delta=km2deg*geodist(stat_arrivals.lon(j),stat_arrivals.lat(j),M(1),M(2));
            %[tt,p]=mkttime(phase,delta,M(3),model);
            %[tmin(j) I]=min(tt);
                tmin(j)=subtocallttcalcsh(ttmethod, ttmodel, phase, delta, M(3));
            %pmin(j)=p(I);
            end
        
                                                    %% ++++++++++++=========== find the best origin time for the
        %%%%%%%%current location: lon, lat, depth
        %fprintf('     Finding origin time ...\n');
            x=[Mprev_datetime - delta_origintime; Mprev_datetime + delta_origintime];
            df=abs(dtt_fitness(arrivals_obs,tmin,x(1)) - dtt_fitness(arrivals_obs,tmin,x(2)));
            kk=1;
            while df > dttmax
                if kk <= Nmax
                    x1=(x(2) - x(1))*alpha + x(1);
                    x2=(x(2) - x(1))*(1 - alpha) + x(1);
           
                    if dtt_fitness(arrivals_obs,tmin,x1) < dtt_fitness(arrivals_obs,tmin,x(1)), x(1)=x1;end
                    if dtt_fitness(arrivals_obs,tmin,x2) < dtt_fitness(arrivals_obs,tmin,x(2)), x(2)=x2;end
           
                    df=abs(dtt_fitness(arrivals_obs,tmin,x(1)) - dtt_fitness(arrivals_obs,tmin,x(2)));

                    kk = kk + 1;
                else
                    break;
                end
           
            end
        %df
            M_datetime=mean(x);
            %end of finding best fit origin time.
            
            fprintf('      Testing model: %f  %f  %f %f\n',M, M_datetime);
        %fprintf('    --> Best fit origin time: %f\n',M_datetime);
%         x_temp=datevec(x_hat);
%         M_datetime(1)=sprintf('%d/%d/%d',x_temp(2),x_temp(3),x_temp(1));
%         M_datetime(2)=sprintf('%d:%d:%f',x_temp(4),x_temp(5),x_temp(6));
        
        % convert travel time to arrivals in UTC time.
            arrivals_pred=tmin/secondsperday + M_datetime;
        
            rr=secondsperday*(arrivals_obs - arrivals_pred);

            loglike=-0.5*rr*invcov*rr'/T(i);
            
            %fprintf('loglikeprev: %f; loglike: %f\n',loglikeprev,loglike);
            
            %loglike=-0.5*rr*rr'/T(i);
        %use Metropolis rule to decide whether or not to accept model
            rat=exp(loglike-loglikeprev);
        

            if loglike > loglikeprev
                accept=1;
                disp('Good one!');
                num_good=num_good + 1;
            else
                r=rand;
%                 rat
                if r<rat
                    accept=1;
                    disp('Bonus!!!');
                    num_bonus=num_bonus + 1;
                else
                    disp('Bad one. Dropped!');
                    accept=0;
                    num_bad=num_bad+1;
                end
            end

    
            if accept==1;
                loglikeprev = loglike;
                Mprev = M;
                Mprev_datetime = M_datetime;
            else
                loglike = loglikeprev;
                M = Mprev;
                M_datetime = Mprev_datetime;
            end

            M_all(:,loop)=M;
            M_all_datetime(:,loop)=M_datetime;
        
            all_loglike(loop)=loglike;
            loop=loop+1;
        end
        
    end %loop for MCMC
    
    fprintf('Number of good: %d\nNumber of bad: %d\nNumber of bonus: %d\n\n',num_good,num_bad,num_bonus);
    
%     [max_loglike Ind]=max(all_loglike);
%     M0=M_all(:,Ind);% update the initial value.
%     M0_datetime=M_all_datetime(:,Ind);
    
    %if mod(i,2)==0
    subplot(3,4,i);
    plot(M_all(1,:),M_all(2,:),'.');
    %xlim(bounds_lon);
    %ylim(bounds_lat);
    title(['itn=',int2str(i),'; T=',num2str(T(i))]);
%     end
    
    M0=M;
    M0_datetime=M_datetime;
    
    Mout(i,:)=M0;
    Mout_datetime(i)=M0_datetime;
end

% get ready for returning values.
fprintf('Getting ready to return values ...\n');
tmp=datevec(M0_datetime);
origin.origdate=sprintf('%d/%d/%d',tmp(2),tmp(3),tmp(1));
origin.origtime=sprintf('%d:%d:%5f',tmp(4),tmp(5),tmp(6));
origin.lon=M0(1);
origin.lat=M0(2);
origin.depth=M0(3);

%remove temporary files.
fprintf('\nRemoving temporary files ...\n');
!rm par.tmp error.tmp tmp.out tt.out

%clear temp variables.
%clear df x1 x2 x tmp i j kk loop Mprev M0 M M0_datetime M_datetime M_all M_all_datetime loglike loglike0 loglikeprev;

return;
end  % end of function locatebysa.