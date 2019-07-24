function [trdata,trhdr]=readrf(ftrace, mdversion, nheaderlines)
%This subroutine reads the receiver function trace in text 
%file, which is output of RFeditor3+ or newer with speficied metadata version.
%See the first line of output file for metadata version specifier.
%
%USAGE: [trdata,trhdr]=readrf(ftrace)
%		[trdata,trhdr]=readrf(ftrace, mdversion)
%       [trdata,trhdr]=readrf(ftrace, mdversion, nheaderlines)
%Modified on: 12/31/2015
%
%trace header lines
% station            : CCM
% start_time         :  9/09/2011   2:43:45.970
% evid               :        117
% samples            :       6001
% dt                 :   0.025000
% t0                 :   0.000000
% stack_weight       :  0.0111
% RT_xcorcoe         :  0.4915
% RFQualityIndex     :  0.3438
% DeconSuccessIndex  :  0.4995
% niteration         :         81
% nspike             :         57
% epsilon            :   29.01390
% peakamp            :    0.24970
% averamp            :    0.00692
% rawsnr             :    0.60000

%to make the program compatible with earlier usage, default mdversion is 1.
if(nargin==1) 
    mdversion=1;
    nheaderlines=0;
elseif(nargin==2)
    nheaderlines=1;
end
  
fid=fopen(ftrace);
fprintf('Reading trace from: %s \n',ftrace);

if(mdversion==1)
	trhdr_default=struct('sta','-',...
             'UTCstart','-',...
             'evid',-9999,...
             'nsamp',-9999,...
             'dt',-9999.9,...
             't0',-9999.9,...
             'stackweight',-9999.9,...
             'rt_xcorcoe',-9999.9,...
             'rfqi',-9999.9,...
             'dsi',-9999.9,...
             'niteration',-9999,...
             'nspike',-9999,...
             'epsilon',-9999.9,...
             'peakamp',-9999.9,...
             'averamp',-9999.9,...
             'rawsnr',-9999.9);    %trace header
	
	trhdr=trhdr_default;
    for k=1:nheaderlines
        temp=fgetl(fid);
    end
    temp=fgetl(fid);
    
    k=strfind(temp,':');
    trhdr.sta=sscanf(temp((k+1):length(temp)),'%s');
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.UTCstart=sscanf(temp((k+1):length(temp)),'%s');
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.evid=sscanf(temp((k+1):length(temp)),'%d');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.nsamp=sscanf(temp((k+1):length(temp)),'%d');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.dt=sscanf(temp((k+1):length(temp)),'%g');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.t0=sscanf(temp((k+1):length(temp)),'%g');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.stackweight=sscanf(temp((k+1):length(temp)),'%g');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.rt_xcorcoe=sscanf(temp((k+1):length(temp)),'%g');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.rfqi=sscanf(temp((k+1):length(temp)),'%g');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.dsi=sscanf(temp((k+1):length(temp)),'%g');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.niteration=sscanf(temp((k+1):length(temp)),'%d');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.nspike=sscanf(temp((k+1):length(temp)),'%d');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.epsilon=sscanf(temp((k+1):length(temp)),'%g');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.peakamp=sscanf(temp((k+1):length(temp)),'%g');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.averamp=sscanf(temp((k+1):length(temp)),'%g');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.rawsnr=sscanf(temp((k+1):length(temp)),'%g');

elseif(mdversion==2)
	trhdr_default=struct('sta','-',...
             'UTCstart','-',...
             'evid',-9999,...
             'magnitude',-9999.9,...
             'magtype','-',...
             'nsamp',-9999,...
             'dt',-9999.9,...
             't0',-9999.9,...
             'stackweight',-9999.9,...
             'rt_xcorcoe',-9999.9,...
             'rfqi',-9999.9,...
             'dsi',-9999.9,...
             'niteration',-9999,...
             'nspike',-9999,...
             'epsilon',-9999.9,...
             'peakamp',-9999.9,...
             'averamp',-9999.9,...
             'rawsnr',-9999.9);    %trace header  

	trhdr=trhdr_default;
    for k=1:nheaderlines
        temp=fgetl(fid);
    end
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.sta=sscanf(temp((k+1):length(temp)),'%s');
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.UTCstart=sscanf(temp((k+1):length(temp)),'%s');
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.evid=sscanf(temp((k+1):length(temp)),'%d');
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.magnitude=sscanf(temp((k+1):length(temp)),'%g');
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.magtype=sscanf(temp((k+1):length(temp)),'%s');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.nsamp=sscanf(temp((k+1):length(temp)),'%d');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.dt=sscanf(temp((k+1):length(temp)),'%g');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.t0=sscanf(temp((k+1):length(temp)),'%g');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.stackweight=sscanf(temp((k+1):length(temp)),'%g');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.rt_xcorcoe=sscanf(temp((k+1):length(temp)),'%g');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.rfqi=sscanf(temp((k+1):length(temp)),'%g');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.dsi=sscanf(temp((k+1):length(temp)),'%g');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.niteration=sscanf(temp((k+1):length(temp)),'%d');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.nspike=sscanf(temp((k+1):length(temp)),'%d');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.epsilon=sscanf(temp((k+1):length(temp)),'%g');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.peakamp=sscanf(temp((k+1):length(temp)),'%g');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.averamp=sscanf(temp((k+1):length(temp)),'%g');
    
    temp=fgetl(fid);
    k=strfind(temp,':');
    trhdr.rawsnr=sscanf(temp((k+1):length(temp)),'%g');

else
	error('ERROR: wrong metadata version specifier, only use 1 or 2. Default is 1.');   
end

trdata=fscanf(fid,'%g',[1,inf]);

fclose(fid);

return;
end