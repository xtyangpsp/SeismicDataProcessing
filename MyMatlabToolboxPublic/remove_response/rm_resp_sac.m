function outtrace = rm_resp_sac(intrace,lo_corner,npoles,pole_zero_dir,pole_zero_namebase,pole_zero_separator)
% function to remove instrument response of irisFetch data structure or SACPZ files
% written by Ge Jin, 2014/02/27, adapted by H Janiszewski
% Modified by Xiaotao Yang to process generic SAC file, without specifying
% event information.
%

pole_zero_file=strcat(pole_zero_namebase,pole_zero_separator,intrace.KNETWK,...
    pole_zero_separator,intrace.KSTNM,'*',pole_zero_separator,intrace.KCMPNM);
pole_zero_file_names=dir(fullfile(pole_zero_dir,pole_zero_file));

pzexist=0;
for nfile=1:length(pole_zero_file_names)
    pole_zero_file_name=pole_zero_file_names(nfile).name;
    
    pzfn_good=[pole_zero_dir '/' pole_zero_file_name];
    pzexist=1;
    break
end

isfigure = 0;

data = intrace.DATA1;
data = detrend(data);
data = flat_hanning_win(1:length(data),data,1,length(data),50);

N = intrace.NPTS;
delta = intrace.DELTA;
T = N*delta;

if mod(N,2)
     faxis = [0:(N-1)/2,-(N-1)/2:-1]*(1/T);
else
     faxis = [0:N/2,-N/2+1:-1]*(1/T);
end

%%%%%%%%%%
% My edits haj
%%%%%%%%%%

if pzexist == 1
    [zz,pp,constant] = read_sac_pole_zero(pzfn_good);

    zeros=zz;
    poles=pp;
    gain=constant;
end

w = faxis.*2*pi;
resp = ones(size(w));
if isempty(zeros) == 1
    zeros = 0;
end
for ip = 1:length(poles)
	resp = resp./(1i*w - poles(ip));
end
for ip = 1:length(zeros)
	resp = resp.*(1i*w - zeros(ip));
end
resp = resp*gain;
if isfigure
	figure;
	clf
	set(gcf,'position',[360   514   900   400]);
	hold on
	subplot(1,2,1)
	set(gca,'fontsize',18)
	semilogy(faxis,abs(resp),'rx');
	subplot(1,2,2)
	set(gca,'fontsize',18)
	plot(faxis,angle(resp),'rx');
end

lo_w=2*pi*lo_corner;
hpfiltfrq=( ((w./lo_w).^(2*npoles))./(1+(w./lo_w).^(2*npoles)) );
norm_trans=hpfiltfrq./resp;    % this is normalization transfer function
norm_trans(isnan(norm_trans)) = 0;

fftdata = fft(data);
fftdata = fftdata(:).*norm_trans(:);
data_cor = real(ifft(fftdata));

outtrace = intrace;
outtrace.DATA1 = data_cor;

disp(['Station: ',intrace.KSTNM,'.',intrace.KCMPNM,' deconv to ',num2str(intrace.IDEP)]);

return

