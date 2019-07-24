function extract_kernels(root_to_event,kernel_pair,fband,fbid)
% extract_kernels
% Modified by Xiaotao Yang from draw_kenels (by Wei Zhang).
% This function extracts and saves kernels for multiple periods to *.mat
% files for later use (e.g., plotting).
% fbid: indices of frequency bands for the current kernel pair.

% clear all

%MFILEROOT='/net/fs01/data/tibet/code/mfiles';
%MFILEROOT='/net/fs01/data/yang/easthemi/mfiles';
MFILEROOT='../../mfiles';
path([MFILEROOT '/fun-spool'],path);

%OUTPUT_ROOT=['./'];

% ----------------------- parameter -----------------------
% flag_surf=1;
% flag_pcolor=1;
% flag_overlap = 1;
% flag_print = 0;
% flag_jetwr = 1;
flag_km=0;

% minlon=283.8; maxlon=293.2; minlat=38.2;maxlat=47.8;
% topo
%load /net/fs01/data/yang/easthemi/misc/topogrid.mat
% load ../../misc/topogrid.mat
% 
% topo=qt;
% 
% % coastline
% !cp /home/yang/Proj/Shared/coastline/coastline_1to5m.dat coast.dat
% load coast.dat
% fidin=fopen('coast_NewEngland.dat','r');
% coast=textscan(fidin,'%f %f\n','CommentStyle','>');
% coastlon=coast{1};coastlat=coast{2};
% fclose(fidin);
% %indx=find(coastlon < 137 | coastlon > 146);
% indx=find(coastlon < -27.6 | coastlon > 155);
% coastlon(indx)=[];coastlat(indx)=[];clear indx;
% indx=find(coastlat < -52 | coastlat > 52);
% coastlon(indx)=[];coastlat(indx)=[];clear indx
% 
% clen=length(coastlat); cR(1:clen,1)=6381e3;
% [cy,cx,cz]=sph2cart(coastlon*pi/180,coastlat*pi/180,cR);

% topo lat long defined as the simulation grid
%RUN_ROOT='/net/fs01/data/yang/easthemi/ite_01/sim.station/skel/fx';
%RUN_ROOT='../sim.station/skel/fx';
RUN_ROOT=strcat(root_to_event,'/fz/');
fnm_conf   =[RUN_ROOT '/' 'SeisFD3D.conf'];
dir_coord  =[RUN_ROOT '/' 'input'];
% dir_metric =[RUN_ROOT '/' 'input'];
% dir_media  =[RUN_ROOT '/' 'input'];
% dir_source =[RUN_ROOT '/' 'input.src'];
% dir_station=[RUN_ROOT '/' 'input'];
% dir_out    =[RUN_ROOT '/' 'output'];

id = 2; % snapshot id (see SeisFD3d.conf) 
subs=[1,1,1];subc=[-1,-1,1];subt=[1,1,1];
[snapinfosurf]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);

[CLATSURF,LONSURF,RSURF]=gather_coord(snapinfosurf,'coorddir',dir_coord);
nxsurf=size(CLATSURF,1);nysurf=size(CLATSURF,2);nzsurf=size(CLATSURF,3);

LATSURF=pi/2-CLATSURF;
clatsurf=CLATSURF*180/pi;
lonsurf=LONSURF*180/pi;
latsurf=90-clatsurf;

%Zin=topo;
Zlat=LATSURF;Zlon=LONSURF;
% define area to be plotted
%flatmin=-45*pi/180;
flatmin=Zlat(end,1);flatmax=Zlat(1,1);
%flonmin=-30*pi/180;
flonmin=Zlon(1,1);flonmax=Zlon(1,end);
nindlat=find(Zlat<flatmin);
nindlon=find(Zlon<flonmin);
% Zin(nindlat)=NaN; Zin(nindlon)=NaN;
% phi_min=flatmin; phi_max=flatmax;
% theta_min=flonmin; theta_max=flonmax;

%% for kernels
%
%OUTPUT_ROOT='./';
%OUTPUT_ROOT='../sim.station/CN.GAC/fz/';
OUTPUT_ROOT=RUN_ROOT; %'../sim.station/TA.G61A/fz/';
fnm_conf=[OUTPUT_ROOT 'SeisFD3D.conf'];
%pnm_metric='/net/fs01/data/yang/easthemi/ite_01/sim.input/';
pnm_metric='../sim.input/';

% kernel_pair='TA.G61A/LD.TRNY';

% fbid=[2 3 4 5 6]; %frequency band IDs.
%fbid=[2];
% fband=[0.0067 0.01333;0.01 0.02;0.01333 0.0286;0.02 0.04;0.0286 0.0667; 0.05 .1; 0.0667 0.1333; .1 .2];
pband=flip(1./fband,2);
% [nfb, nc]=size(fband);

id=1;
for ifb=1:length(fbid)

    pnm_out=strcat(kernel_pair,'/BHZ/',num2str(fbid(ifb)),'/T1T2.P2/');

    %period = '50-100s'; %f2
    %period = '35-75s'; %f3
    %period = '25-50s'; %f4
    %period = '15-35s'; %f5
    %period = '10-20s'; %f6
    % period = '7.5-15s'; %f7
    % for vertical-vertical cross correlation
    % wavename = 'Rayleigh wave';
    subs=[ 1 1 1]; subc=[ -1 -1 -1]; subt=[ 1   1  1 ];

    var_list=[]; scl_caxis_list=[];
    var_list{end+1}='phase_Vp'; scl_caxis_list{end+1}=[-1.5e-14,1.5e-14];
    var_list{end+1}='phase_Vs'; scl_caxis_list{end+1}=[-1.5e-14,1.5e-14];

    % -------------------- load data --------------------------

    [snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);

    %-- get coord data
    [X,Y,Z]=gather_coord(snapinfo,'coorddir',pnm_metric);
    %X=90-X/pi*180; Y=Y/pi*180; Z=6371e3-Z; Z=Z/1e3;
    R=Z; CLAT=X; LON=Y;
    nx=size(X,1);ny=size(X,2);nz=size(X,3);

    clat=CLAT*180/pi;
    Klon=LON*180/pi;
    Klat=90-clat;

    [y,x,z]=sph2cart(LON,pi/2-CLAT,R);
    Kdepth=6371-Z/1000;
    str_unit='m';
    if flag_km
       x=x/1e3;y=y/1e3;z=z/1e3; str_unit='km';
    end

    % ----------------------- plot kernel -----------------------------------

    nvar=length(var_list);

    for n=1:nvar

        [V,~]=gather_dist(snapinfo,id,var_list{n},'outdir',pnm_out);

        V=double(V);
        kernelphase=var_list{n};
        thispband=pband(fbid(ifb),:);
        fstr=strcat(num2str(round(pband(fbid(ifb),1))),'-',num2str(round(pband(fbid(ifb),2))),'s');
        %save data to file.
        save(['kernel_',strrep(kernel_pair,'/','_'),fstr,'_',kernelphase,'.mat'],...
            'kernel_pair','kernelphase','Klon','Klat','Kdepth','V','thispband');
    end % nvar
end

end