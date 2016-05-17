function [vmodel,header]=get_layered_vmodel(vmodel_all,header_all,station)
%vmodel: model array (nlayers x 4), columns 1-4 for vp, vs, density, and thickness.

%set initial value;
[nsta,nl,s3]=size(vmodel_all);
vmodel=zeros(1,nl,4);
header=struct('station',station,'nlayers',nl);

for i=1:nsta
    if(strcmp(header_all(i).station, station))
        vmodel=vmodel_all(i,:,:);
        break;
    end
    
end

end