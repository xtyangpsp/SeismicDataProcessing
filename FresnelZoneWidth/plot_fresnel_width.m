function plot_fresnel_width(T,max_depth,savetofile,vp,vs,mod1dfile)
%plot_fresnel_width.
%usage: plot_fresnel_width(period, max_depth)
%Reference: Pavlis (2011). Three-dimensional, wavefield imaging of broadband 
%           seismic array data, Computers & Geosciences (37), 1054-1066

if(nargin==1) 
    max_depth=100;
    vp=8.1;
    vs=4.5;
    mod1dfile='iasp91.mod1d';
    savetofile='no';
elseif(nargin==2)
    vp=8.1;
    vs=4.5;
    mod1dfile='iasp91.mod1d';
    savetofile='no';
elseif(nargin==3)
    vp=8.1;
    vs=4.5;
    mod1dfile='iasp91.mod1d';
elseif(nargin==4)
    vs=4.5;
    mod1dfile='iasp91.mod1d';
elseif(nargin==5)
    mod1dfile='iasp91.mod1d';
end

if(isempty(num2str(vp))); vp=8.1;end
if(isempty(num2str(vs))); vp=4.5;end
if(isempty(mod1dfile)); mod1dfile='iasp91.mod1d';end
if(isempty(savetofile)); savetofile='no';end

dz=1.0;

system(['/opt/antelope/5.5/contrib/bin/getT2D ',num2str(mod1dfile),...
    ' 60 50 -rti ',num2str(dz),' -md ',num2str(max_depth),' >temp.d']);

td=load('temp.d');
%!rm temp.d
td_z2tau=(td(:,2)./td(:,5))';
td_z2tau=td_z2tau(2:length(td_z2tau));
%size(td_z2tau)
lat0=38;lon0=-88.5;
theta=0:0.1:2*pi;
R0=6371.0;

z=dz:dz:max_depth;
%size(z)
D=zeros(length(T),length(z));
for i=1:length(T)
    %D(i,:)=2*vs*sqrt((z./vs + T(i)/2).^2 - (z./vs).^2);
    D(i,:)=2*vs*sqrt((vp*T(i)*z)./(td_z2tau*(vp-vs)) + T(i)^2/4);
    
    if(strcmp(savetofile,'yes'))
        fid=fopen(['Fwidth_circle_',num2str(T(i)),'.dat'],'w');
        for j=1:length(z)
            x=lon0+rad2deg(cos(theta)*D(i,j)/(2*(R0-z(j))));
            y=lat0+rad2deg(sin(theta)*D(i,j)/(2*(R0-z(j))));
            for k=1:length(theta)
                fprintf(fid,'%g  %g  %g\n',x(k),y(k),z(j));
            end
            fprintf(fid,'>\n');

        end
        fclose(fid);
    end
    %plot(D,z);
    %hold on;
end    
%hold off;
plot(D,z);
xlabel('Fresnel zone width (km)');
ylabel('Depth (km)');

set(gca,'YDir','reverse');
grid on;


for j=1:length(T)
    LS{j}=num2str(T(j));
end

legend(LS);

end