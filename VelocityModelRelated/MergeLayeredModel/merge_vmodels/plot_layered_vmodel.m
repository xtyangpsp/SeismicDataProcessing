function [z,vmodel]=plot_layered_vmodel(model,plotornot)
%plotornot: default is yes, use tags other than 'yes','yesplot','plotyes'
%to suppress plot and only get the output. This option is designed for mass
%processing multiple stations to get z and vmodel data.
%
%model: is a layered model with format:vp  vs   density   thickness
%For example:
%  5.5 2.7 2.65 0.803723
%  3.3 1.9 2.4 0.012796
%  5 2.7 2.65 0.450224
%  5 2.2 2.65 0.056131
%  6.5 3.5 2.75 1.9259
%  6 3.3 2.6 0.160257
%  5.5 3.2 2.65 0.133474
%  6.2 3.5 2.8 0
%It is assumed that the first layer starts from the surface with depth =
%0.0km.
%
%Xiaotao Yang @IU
%2016.02.09

if(nargin==1)
    plotornot='yes';
end

[c1,c2,c3]=size(model);
%Generate depth vector
z=zeros(1,2*c2);
vmodel=(zeros(2*c2,3)); %for vp, vs, density;
z(2)=z(1)+model(1,1,4);
vmodel(1,:)=model(1,1,1:3);
vmodel(2,:)=vmodel(1,:);
% j=2;
for i=2:c2
    z(2*i-1)=z(2*i-2);
    z(2*i)=z(2*i-1)+model(1,i,4);
    
    vmodel(2*i-1,:)=model(1,i,1:3);
    vmodel(2*i,:)=model(1,i,1:3);
    
end

if(strcmp(plotornot,'yes') || strcmp(plotornot,'plotyes') || strcmp(plotornot,'yesplot'))
    plot(vmodel,z); axis ij;
    % minvp=min(vmodel(:,1));
    % minvs=min(vmodel(:,2));
    % minrho=min(vmodel(:,3));
    % maxvp=max(vmodel(:,1));
    % maxvs=max(vmodel(:,2));
    % maxrho=max(vmodel(:,3));

    xlim([min([min(vmodel(:,1)), min(vmodel(:,2)), min(vmodel(:,3))])-0.5, ...
        (max([max(vmodel(:,1)), max(vmodel(:,2)), max(vmodel(:,3))])+0.5)]);
    legend('vp','vs','density','Location','south');
    xlabel('Velocity (km/s)');ylabel('Depth (km)');
% else
%     disp('!!! Caution: Plot is suppressed.');
end


end