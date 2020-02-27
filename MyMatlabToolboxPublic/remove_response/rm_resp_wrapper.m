%This is a wrapper function to remove the instrument response from multiple SAC files.
%This function is designed to process a large amount of data sets/sac
%files.
%It calls individual function to process each sac file.
%
function rm_resp_wrapper(datadir,respdir,outdir,freqmin,npoles,plotfig)
% datadir='testdata';
% respdir=strcat(datadir,'/IRISDMC');
% npoles=2;
% freqmin=0.01;
figflag=1;
if nargin < 5
    error('Not enough input arguments!');
elseif nargin==5
    figflag=plotfig;
end
if ~exist(outdir,'dir')
    system(['mkdir ' outdir]);
end
allfiles1=dir([datadir, '/*.sac']);
allfiles2=dir([datadir, '/*.SAC']);
allfiles_temp=[allfiles1;allfiles2];
allfiles=cell(size(allfiles_temp));
for j=1:size(allfiles_temp,1)
    allfiles{j}=allfiles_temp(j).name;
end
allfiles_unique=unique(allfiles);
nfiles=size(allfiles_unique,1);

for i=1:nfiles
    % read EGF file
    
    filename=allfiles_unique{i}; 
    
    intrace=readsac([datadir,'/',filename]);
    
    %demean and detrend
    dtemp=intrace.DATA1;
    intrace.DATA1=detrend(dtemp-nanmean(dtemp));
    outtrace=rm_resp_sac(intrace,freqmin,npoles,respdir,'SACPZ','.');
    
    outtrace.FILENAME = strcat(outdir,'/',intrace.FILENAME);
    
    if ~writesac(outtrace)
        error(['Failed to save: ' outtrace.FILENAME '. Make sure the directory exists!']);
    end
    
    if figflag
        dt=intrace.DELTA;
        timeflag=nan(intrace.NPTS,1);
        tt=intrace.B:dt:intrace.E;
        if length(tt)< intrace.NPTS
            timeflag(1:end-1)=tt;
            timeflag(end)=intrace.E+dt;
        else
            timeflag=tt;
        end
        figure('Position',[400 400 700 350]);
        hold on;
        plot(timeflag,intrace.DATA1/max(abs(intrace.DATA1)),'b-','color',[.5 .5 .5],'linewidth',2);
        plot(timeflag,outtrace.DATA1/max(abs(outtrace.DATA1)),'r-');
        legend('raw','rm_resp');
        xlabel('time (s)');
        hold off;
        axis on;
        grid on;
        box on;
        title(strcat(intrace.KNETWK,'.',intrace.KSTNM,'.',intrace.KCMPNM));
        set(gcf,'PaperPositionMode','auto');   
        eval(['print -dpng -r300 ' filename(1:end-4) '_rmresp.png']);
        close
    end
end

disp(['Removed responses from ',num2str(nfiles),' SAC files']);
end

