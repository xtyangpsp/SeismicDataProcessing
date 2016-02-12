function [vp,vs,rou]=myreadmod1d(fmod)
%Read velocity and density information from 1d model with format of mod1d.
%output model for each component (vp, vs, rou) includes 3 columns of depth,
%velocity, gradient.
%USAGE: 
%Xiaotao Yang @ Indiana University
%2016.1.27

fid=fopen(fmod,'r');

%get number of lines
n=0;
tline = fgetl(fid);
while ischar(tline)
    n=n+1;
    tline = fgetl(fid);
end

vmodelp=zeros(n,3);
vmodels=zeros(n,3);
vmodelrou=zeros(n,3);

%go back to file start by rewinding.
frewind(fid);
ip=1;is=1;ir=1;
tline = fgetl(fid);
while ischar(tline)
   vpline = strfind(tline, 'Pvelocity');
   if(~isempty(vpline))
       vmodelp(ip,1:3)=sscanf(tline((vpline+9):length(tline)),'%g%g%g');
       tline = fgetl(fid);
       ip=ip+1;
   else
      vsline = strfind(tline, 'Svelocity');
      if(~isempty(vsline))
          vmodels(is,1:3)=sscanf(tline((vsline+9):length(tline)),'%g%g%g');
          tline = fgetl(fid);
          is=is+1;
      else
          vrouline = strfind(tline, 'Density');
          if(~isempty(vrouline))
              vmodelrou(ir,1:3)=sscanf(tline((vrouline+7):length(tline)),'%g%g%g');
              tline = fgetl(fid);
              ir=ir+1;
          end
      end
   end
end

vp=vmodelp(1:(ip-1),:);
vs=vmodels(1:(is-1),:);
rou=vmodelrou(1:(ir-1),:);
%vmodelrou(:,1)

fclose(fid);

end