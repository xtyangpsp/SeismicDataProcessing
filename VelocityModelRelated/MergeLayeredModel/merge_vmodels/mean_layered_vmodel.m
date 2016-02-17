function meanmodel=mean_layered_vmodel(modelfile)
%Derive mean model of the input layered model. The modelfile should have
%the following format for each point:
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
%Xiaotao Yang @ IU
%2016.02.17
%

[model0,header0]=read_layered_vmodel(modelfile);

%error if the models have different number of layers.
nm=length(header0);

nl=header0(1).nlayers;
meanmodel=nan(1,nl,4);
for i=2:nm
    if(header0(i).nlayers ~= nl)
        error(['ERROR: inconsistent number of layers for model labeled as: ',header0(i).station]);
    end
    nl=header0(i).nlayers;
end

for i=1:nl
    meanmodel(1,i,:)=[mean(model0(:,i,1)),mean(model0(:,i,2)),mean(model0(:,i,3)),mean(model0(:,i,4))];
end

end
