% $Date: 2009-01-14 16:34:32 -0500 (Wed, 14 Jan 2009) $
% $Revision: 42 $
% $LastChangedBy: zhangw $

clear all

set_mfiles_path
%[fnm_conf,dir_coord,dir_metric,dir_media,dir_source, ...
%  dir_station,dir_out]=get_simul_path('root','./fz2');
%dir_coord='../prem.nowater/input'
%dir_media='../prem.nowater/input'

[fnm_conf,dir_coord,dir_metric,dir_media,dir_source, ...
  dir_station,dir_out]=get_simul_path('root','./fz.dr1.ricker');

dir_coord='./model.dr1'
dir_media='./model.dr1'
dir_station='./model.dr1.station'

% ----------------------- parameter -----------------------

%N=4; Fc=1/10;

cmp=[];cmp0=[];cmp1=[];dir_sgt=[];
cmp{end+1}= 'Uz'; cmp0{end+1}='1'; cmp1{end+1}='( 1)'; dir_sgt{end+1}='fz.dr1.ricker/output';
ncmp=numel(cmp);

%M=[ 0 0 0; 0 0 0; 0 1 0 ]; %000
%M=[ 1 0 0; 0 1 0; 0 0 1 ]; %000
%M=src_mech(0,90,0); %001
%M=src_mech(0,90,90); %002
%M=src_mech(0,90,180); %003
%M=src_mech(0,90,270); %004
%M=src_mech(0,90,45); %005
%M=src_mech(0,90,135); %006
%M=src_mech(0,90,225); %007
%M=src_mech(0,90,315); %008
%M=src_mech(0,45,45); %009
%M=src_mech(45,45,45); %010
%M=src_mech(45,45,90); %011
% 8km
%M=[ 1 0 0; 0 1 0; 0 0 1 ]; %012
%M=src_mech(45,45,90); %013
%M=src_mech(0,90,0); %014
%M=src_mech(0,90,45); %015
%M=src_mech(0,45,45); %016
%M=src_mech(45,45,45); %017
M=src_mech(45,45,90); %018

id=5;
% 25*(15+10)+1-51-50
%subs=[1,576-25*5,222,1];subc=[1,1,1,10000];subt=[1,1,1,1];  % 3
subs=[1,576-25*10,222,1];subc=[1,1,1,10000];subt=[1,1,1,1];  % 3

%subs=[1,576,10,1];subc=[1,1,1,10000];subt=[1,1,1,1];  % 2
%subs=[1,551,10,1];subc=[1,1,1,10000];subt=[1,1,1,1];  % 3
%subs=[1,501,10,1];subc=[1,1,1,10000];subt=[1,1,1,1];  % 5
%subs=[1,376,10,1];subc=[1,1,1,10000];subt=[1,1,1,1];  % 10
%subs=[1,251,10,1],subc=[1,1,1,10000];subt=[1,1,1,1];  % 15

spec_line='r';

% -------------------- load data --------------------------
[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);
[CLAT,LON,R]=gather_coord(snapinfo,'coorddir',dir_coord);
nx=size(CLAT,1);ny=size(CLAT,2);nz=size(CLAT,3);

clat=CLAT*180/pi;
lon=LON*180/pi;

for n=1:ncmp
    [U{n},T{n}]=retrieve_seismo_sgt(snapinfo,id,M,'outdir',dir_sgt{n},'mediadir',dir_media);
    U{n}=permute(U{n},[4 1 2 3]);
    nt=length(T{n});
    stept=T{n}(2)-T{n}(1);
end

if 0
% ------------------ convolv stf ------------------------
%S0=stf_gauss(T{1},5.575,18);D0=stf_gaussderiv(T{1},5.575,18);
S1=stf_ricker(T{1},0.05,30);D1=stf_rickerderiv(T{1},0.05,30);
for n=1:ncmp
for k=1:nz
for j=1:ny
for i=1:nx
    U1=conv(D1,U{n}(:,i,j,k));
    U{n}(2:nt,i,j,k)=U1(1:nt-1)*stept;
end
end
end
end

end

for n=1:ncmp
for k=1:nz
for j=1:ny
for i=1:nx
    V{n}(:,i,j,k)=gradient(U{n}(:,i,j,k),stept);
    A{n}(:,i,j,k)=gradient(V{n}(:,i,j,k),stept);
end
end
end
end

% -- filter --
if exist('Fc','var')
for n=1:ncmp
for k=1:nz
for j=1:ny
for i=1:nx
   Feff=1/stept/2; Wn=Fc/Feff; [b,a]=butter(N,Wn);
   U{n}(:,i,j,k)=filter(b,a,U{n}(:,i,j,k));
   V{n}(:,i,j,k)=filter(b,a,V{n}(:,i,j,k));
   A{n}(:,i,j,k)=filter(b,a,A{n}(:,i,j,k));
end
end
end
end
end

for n=1:ncmp
    Vi0(n)=max(max(max(max((V{n})))));
    Ui0(n)=max(max(max(max((U{n})))));
    Ai0(n)=max(max(max(max((A{n})))));
end
V0=max(Vi0);
U0=max(Ui0);
A0=max(Ai0);

% ------------------ plot figure ------------------------
for m=1:ncmp
    %figure(m);
    subplot(3,1,m)
    n=0;
for k=1:nz
for j=1:ny
for i=1:nx
    n=n+1;
    W=eval([cmp{m}(1) '{' num2str(n) '}(:,i,j,k)'])/eval(cmp0{m})*eval(cmp1{m});
    plot(T{n},W,spec_line);
    hold on;
end
end
end
    title(cmp{m});
end

if 0
subplot(3,1,1)
title('North')
xlim([0,500])
subplot(3,1,2)
title('East')
xlim([0,500])
subplot(3,1,3)
title('Radius')
xlim([0,500])
end

