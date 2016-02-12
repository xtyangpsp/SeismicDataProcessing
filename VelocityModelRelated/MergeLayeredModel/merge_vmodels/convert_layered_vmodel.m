function [vmodel,nlayers]=convert_layered_vmodel(model,z)
%model: contails N rows and 3 columns of vp, vs, density
%z: corresponding depth
%USAGE: [vmodel,nlayers]=convert_layered_vmodel(model,z)
%Returns:
%   vmodel: format is vp, vs, density, thickness(km)
%   nlayers: number of layers, may include 0.0km thick bottom layer,
%   depending on the given model.
%
%Xiaotao Yang @IU
%2016.02.09
%
nd=length(z); %number of datapoints.
nlayers=1;
for i=1:(nd-1)
    
    if(z(i+1) - z(i)>0.001)
        vmodel(1,nlayers,1:4)=[model(i,:), z(i+1) - z(i)];
        nlayers=nlayers+1;
        if(i==(nd-1))
            nlayers=nlayers-1;
        end
    elseif(i==(nd-1))
        vmodel(1,nlayers,1:4)=[model(i,:), 0.0];
    end 

end

end