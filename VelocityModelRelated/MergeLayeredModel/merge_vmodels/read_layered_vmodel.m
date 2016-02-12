function [model,header]=read_layered_vmodel(modelfile)
%Reads layered model with the format for StaVariableLayerSynthetic used by
%migsimulation.
%USAGE: [model,header]=read_layered_vmodel(modelfile)
%Output:
%1. model: model array (nlayers, 4), columns 1-4 for vp, vs, density, and thickness.
%2. header: header struct with member station and nlayers. 
%*NOTE:
%In output model: all stations have the same number of layers. But sations
%with less than the maximum number of layers will be filled with NaN. Do
%remember to check for NaN when using the model data.
%
%Format example:
% LG19 8
%  5.5 2.7 2.65 0.803723
%  3.3 1.9 2.4 0.012796
%  5 2.7 2.65 0.450224
%  5 2.2 2.65 0.056131
%  6.5 3.5 2.75 1.9259
%  6 3.3 2.6 0.160257
%  5.5 3.2 2.65 0.133474
%  6.2 3.5 2.8 0
%columns: vp, vs, density, thickness(km)
%
%Xiaotao Yang @IU, 2016.02.09
%

%define default header structure and values.
header_default=struct('station','-','nlayers',0);

fprintf('readlayeredmodel: Reading on file: %s\n',modelfile);

fid=fopen(modelfile,'r');

tline=fgetl(fid);
i=1;
while ischar(tline)
    header(i)=header_default;
    header(i).station=sscanf(tline,'%s',[1 1]);
    header(i).nlayers=sscanf(tline,'%*s %d'); %use * to skip station string.

    %tmodel(i)=zeros(header.nlayers,4);
    %header(i).station
    tmodel=fscanf(fid,'%g %g %g %g\n',[4 header(i).nlayers]);
    %model(i,:,:)=tmodel';
    
    nlayer(i)=header(i).nlayers;
    
    tline=fgetl(fid);
    i=i+1;

end

model=nan(i-1,max(nlayer),4);

frewind(fid);

tline=fgetl(fid);
i=1;
while ischar(tline)

    tmodel=fscanf(fid,'%g %g %g %g\n',[4 header(i).nlayers]);
    model(i,1:header(i).nlayers,:)=tmodel';

    tline=fgetl(fid);
    i=i+1;

end

fprintf('readlayeredmodel: read %d stations.\n Maximum number of layers: %d\n',i-1,max(nlayer));

fclose(fid);

end