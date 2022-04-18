%save velocity model to netcdf format.
close all;
ite_nm = 'ite_0.05deg_01';
velocitytag='S'; % 'P' for P velocities, 'PR' for Poisson's ratios.
idx = 0; % 0, absolute velocity; 1, velocity perturbation
savefigtag=1;
ncmodelfile='YangETAL_NACratonFWANT_Vs2022.nc';
outparfile=[ite_nm, '_dataforIRISEMC.cdl'];

modelid='NACraton_FWANT_Vs2022';
maskrayfile=['RayCoverOutline_',ite_nm,'_15-30s_cutoff10.mat'];
%
maparea.lon=[-97.5,-79.2]; % temporary 
maparea.lat=[34.5, 47.5];
depthrange=[0,150];
dlon=0.05;
dlat=0.05;

%read previous model
fnm_conf=['./SeisFD3D.conf_' ite_nm];
dir_coord=['./input_' ite_nm];
dir_media=['./updated_input_' ite_nm];
disp(['Read model... ' dir_media]);

largenumber=10e5;

id = 0; subs=[1,1,1];subc=[-1,-1,-1];subt=[1,1,1];
[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);
[XSIM,YSIM,ZSIM]=gather_coord(snapinfo,'coorddir',dir_coord);
% convert from radian to degrees
XSIM=90-XSIM*180/pi; %latitude
YSIM=YSIM*180/pi;

%define the area of plot (exclude pmls)
npml=12; %number of pml layers
minlat=XSIM(end,npml+1,end);maxlat=XSIM(1,1,end-npml);
minlon=YSIM(1,1,end-npml);maxlon=YSIM(npml+1,end,end);

mrh=gather_media(snapinfo,'rho','mediadir',dir_media);
mmu=gather_media(snapinfo,'mu','mediadir',dir_media);
mla=gather_media(snapinfo,'lambda','mediadir',dir_media);
mvp=((mla+2*mmu)./mrh).^0.5;
mvs=(mmu./mrh).^0.5;

x=YSIM(1,:,1)-360;
y=XSIM(:,1,1);
z=6371-squeeze(ZSIM(1,1,:))/1000;

idx_lon=find(x>=maparea.lon(1) & x<=maparea.lon(2));
idx_lat=find(y>=maparea.lat(1) & y<=maparea.lat(2));
idx_depth=find(z>=depthrange(1) & z<=depthrange(2));
outp=mvp/1000;
outs=mvs/1000;

%%
% outrh=mrh/1000;
vunit='kilometer sec^-1';
% rhunit='gram cm^-3';
fidout=fopen(outparfile,'w');
%
fprintf(fidout,'netcdf %s {\n',modelid);
fprintf(fidout,'dimensions: \n');
fprintf(fidout,'depth = %d;\n',length(idx_depth));
fprintf(fidout,'latitude = %d;\n',length(idx_lat));
fprintf(fidout,'longitude = %d;\n',length(idx_lon));

fprintf(fidout,'\nvariables: \n');
fprintf(fidout,'float depth(depth) ;\n');
fprintf(fidout,'depth:long_name = "depth below earth surface" ;\n');
fprintf(fidout,'depth:units = "kilometer ";\n');
fprintf(fidout,'depth:positive = "down ";\n');

fprintf(fidout,'float latitude(latitude) ;\n');
fprintf(fidout,'latitude:long_name = "Latitude; positive north" ;\n');
fprintf(fidout,'latitude:units = "degrees_north" ;\n');
fprintf(fidout,'latitude:standard_name = "latitude" ;\n');

fprintf(fidout,'float longitude(longitude) ;\n');
fprintf(fidout,'longitude:long_name = "Longitude; positive east" ;\n');
fprintf(fidout,'longitude:units = "degrees_east" ;\n');
fprintf(fidout,'longitude:standard_name = "longitude" ;\n');

fprintf(fidout,'float vs(depth, latitude, longitude) ;\n');
fprintf(fidout,'vs:long_name = "Shear-wave Velocity" ;\n');
fprintf(fidout,'vs:display_name = "Vs (km/s)" ;\n');
fprintf(fidout,'vs:units = "km.s-1" ;\n');
fprintf(fidout,'vs:_FillValue = 99999.f ;\n');

fprintf(fidout,'float vp(depth, latitude, longitude) ;\n');
fprintf(fidout,'vp:long_name = "P-wave Velocity" ;\n');
fprintf(fidout,'vp:display_name = "Vs (km/s)" ;\n');
fprintf(fidout,'vp:units = "km.s-1" ;\n');
fprintf(fidout,'vp:_FillValue = 99999.f ;\n');
% fprintf(fidout,'vs:valid_range = "-1000.f, 1000.f";\n');

%global attributes.
fprintf(fidout,'\n\n:title = "Shear wave velocity model for the North American craton" ;\n');
fprintf(fidout,':id = "%s" ;\n',modelid);
fprintf(fidout,':summary = "Shear wave velocity model for the top 150 km beneath North');
fprintf(fidout,' American craton from full-wave ambient noise tomogrpahy.";\n');

fprintf(fidout,':keywords = "shear wave velocity, ambient noise tomography, North American craton,');
fprintf(fidout,' s wave, velocity, Xiaotao Yang";\n');

fprintf(fidout,':Metadata_Conventions = "Unidata Dataset Discovery v1.0";\n');
fprintf(fidout,':creator_name = "Xiaotao Yang" ;\n');
fprintf(fidout,':creator_url = "https://sites.google.com/site/xiaotaoyanggeo" ;\n');
fprintf(fidout,':creator_email = "xtyang@purdue.edu OR stcyang@gmail.com" ;\n');
fprintf(fidout,':institution = "Purdue University" ;\n');

fprintf(fidout,':acknowledgment = "Model was provided by Xiaotao Yang at Purdue University. ');
fprintf(fidout,' This research was supported by the startup funding for Xiaotao Yang at Purdue');
fprintf(fidout,' University. All the seismic data used in this study were downloaded from IRIS ');
fprintf(fidout,' Data Management Center. The computational resources used in this study were ');
fprintf(fidout,'provided by the Purdue research computing." ;\n');

fprintf(fidout,':reference = "Yang, X. et al. (in preparation)" ;\n');
fprintf(fidout,':references = "http://www.iris.edu/dms/products/emc-references\n" ;\n');
fprintf(fidout,':history = "produced in 2022" ;\n');
fprintf(fidout,':comment = "model converted to netCDF by the author with template from IRIS DMC" ;\n');

%
fprintf(fidout,':geospatial_lat_min = %g;\n',min(y(idx_lat)));
fprintf(fidout,':geospatial_lat_max = %g;\n',max(y(idx_lat)));
fprintf(fidout,':geospatial_lat_units = "degrees_north";\n');
fprintf(fidout,':geospatial_lat_resolution = %g;\n',dlat);
fprintf(fidout,':geospatial_lon_min = %g;\n',min(x(idx_lon)));
fprintf(fidout,':geospatial_lon_max = %g;\n',max(x(idx_lon)));
fprintf(fidout,':geospatial_lon_units = "degrees_east";\n');
fprintf(fidout,':geospatial_lon_resolution = %g;\n',dlon);
fprintf(fidout,':geospatial_vertical_min = %g;\n',min(z(idx_depth)));
fprintf(fidout,':geospatial_vertical_max = %g;\n',max(z(idx_depth)));
fprintf(fidout,':geospatial_vertical_units = "kilometer";\n');
fprintf(fidout,':geospatial_vertical_positive = "down";\n');

fprintf(fidout,'data:\n\n');
fprintf(fidout,'depth = ');
for i=length(idx_depth):-1:1
    fprintf(fidout,'%g',z(idx_depth(i)));
    if i==1
        fprintf(fidout,' ;\n\n');
    else
        fprintf(fidout,',');
    end
end
fprintf(fidout,'latitude = ');
for i=length(idx_lat):-1:1
    fprintf(fidout,'%g',y(idx_lat(i)));
    if i==1
        fprintf(fidout,' ;\n\n');
    else
        fprintf(fidout,',');
    end
end
fprintf(fidout,'longitude = ');
for i=1:length(idx_lon)
    fprintf(fidout,'%g',x(idx_lon(i)));
    if i==length(idx_lon)
        fprintf(fidout,' ;\n\n');
    else
        fprintf(fidout,',');
    end
end
clear amask;
load(maskrayfile);
amask=nan(length(squeeze(XSIM(:,1,1))),length(squeeze(YSIM(1,:,1))));
for i=1:size(amask,1)
        clear id00;
        id00=inpolygon(YSIM(1,:,1)-360,XSIM(i,1,1)*ones(size(amask,2),1),raycover.data(:,1),...
                raycover.data(:,2));
        amask(i,id00)=1;
end
            
disp('writing Vs ...');
fprintf(fidout,'vs = \n');
for i=length(idx_depth):-1:1
    for j=length(idx_lat):-1:1
        for k=1:length(idx_lon)    
            maskval=amask(idx_lat(j),idx_lon(k));
            if i==1 && j==1 && k==length(idx_lon)
                if isnan(maskval)
                    fprintf(fidout,'NaN ;\n\n');
                else
                    fprintf(fidout,'%g ;\n\n',outs(idx_lat(j),idx_lon(k),idx_depth(i)));
                end
            else
                if isnan(maskval)
                    fprintf(fidout,'NaNf, ');
                else
                    fprintf(fidout,'%g, ',outs(idx_lat(j),idx_lon(k),idx_depth(i)));
                end
            end
        end
        fprintf(fidout,'\n');
    end
end
disp('writing Vp ...');
fprintf(fidout,'vp = \n');
for i=length(idx_depth):-1:1
    for j=length(idx_lat):-1:1
        for k=1:length(idx_lon)
            maskval=amask(idx_lat(j),idx_lon(k));
            if i==1 && j==1 && k==length(idx_lon)
                if isnan(maskval)
                    fprintf(fidout,'NaN ;\n\n');
                else
                    fprintf(fidout,'%g ;\n\n',outp(idx_lat(j),idx_lon(k),idx_depth(i)));
                end
            else
                if isnan(maskval)
                    fprintf(fidout,'NaNf, ');
                else
                    fprintf(fidout,'%g, ',outp(idx_lat(j),idx_lon(k),idx_depth(i)));
                end
            end
        end
        fprintf(fidout,'\n');
    end
end
fprintf(fidout,'}\n');
fclose(fidout);

%% convert to netcdf file.
disp(['Converting to netcdf file: ',ncmodelfile])
unix(['/usr/local/bin/ncgen3 ','-o ',ncmodelfile,' ',outparfile]);

%% check netCDF file.
disp('Loading netcdf file to check ...')
rawdepth=nc_varget(ncmodelfile,'depth');
rawlon=nc_varget(ncmodelfile,'longitude');
rawlat=nc_varget(ncmodelfile,'latitude');
vtag='vs';

lonind=find(rawlon>=maparea.lon(1) & rawlon<= maparea.lon(2));
latind=find(rawlat>=maparea.lat(1) & rawlat<= maparea.lat(2));
depthind=find(rawdepth>=depthrange(1) & rawdepth<= depthrange(2));

subsetlon=rawlon(lonind);
subsetlat=rawlat(latind);
subsetdepth=rawdepth(depthind);

subsetvs=nc_varget(ncmodelfile,vtag,[depthind(1)-1 latind(1)-1 lonind(1)-1],...
    [length(depthind) length(latind) length(lonind)]);
subsetvs(subsetvs> largenumber)=nan;
%%
state=[];
load ../00masterdatafiles/us_states.mat;
mapzlist=[13,22,31,40,51,63,74,97,109];
figure('Position',[400 400 1200 700]);
for i=1:length(mapzlist)
    z=mapzlist(i);
    [minz,mapz]=min(abs(subsetdepth - z));
%     mapz=mapzlist(i);
    subplot(3,3,i)
    pcolor(subsetlon, subsetlat,squeeze(subsetvs(mapz,:,:))); hold on;
    for sb=1:length(state)
            plot(state(sb).polygon(:,1), state(sb).polygon(:,2),'color',[.5 .5 .5],'LineWidth',1);
    end
    shading flat;
    colormap('jetwr');
    %caxis([4.3 4.7]);
    title(['Vs at ' num2str(round(subsetdepth(mapz))) ' km'],'FontSize',14);
    axis([maparea.lon(1) maparea.lon(2) maparea.lat(1) maparea.lat(2)]);
    daspect([1 cosd(mean(maparea.lat)) 1]);
    colorbar;
    set(gca,'FontSize',12)
    hold off;
end
