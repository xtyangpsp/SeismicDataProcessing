function merge_layered_vmodels(modelfile1,modelfile2,outfilename,plotornot,extractstation,yrange1,yrange2,yrange3)
%USAGE:
%merge_layered_vmodels(modelfile1,modelfile2,outfilename): plot all
%stations, not good for data with too many stations. Use '-' for outfile name to turn off save to file.
%
%merge_layered_vmodels(modelfile1,modelfile2,outfilename,plotornot): set
%flag plotornot to turn plot on or off. Use '-' for outfile name to turn off save to file.
%
%merge_layered_vmodels(modelfile1,modelfile2,outfilename,plotornot,extractstation):
%only merge specified station by extractstation. Use '-' for outfile name to turn off save to file.
%
%merge_layered_vmodels(modelfile1,modelfile2,outfilename,plotornot,extractstation,yrange)
%specify y range for plot, ignore when ploting option is off.
%
%1. Modelfile1 is dominant. only use modelfile2 when modelfile1 does not have
%coverage. for satations included in both model files, merge them.
%2. plotornot: default is yes, use tags other than 'yes','yesplot','plotyes'
%to suppress plot and only get the output. This option is designed for mass
%processing multiple stations.
%Output format example:
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
singlestation=0;
if(nargin==3)
    plotornot='yes';
    singlestation=0;
    yrangeauto1=1;
    yrangeauto2=1;
    yrangeauto3=1;
elseif(nargin==4)
    singlestation=0;
    yrangeauto=1;
    yrangeauto2=1;
    yrangeauto3=1;
elseif(nargin>=5) && (~strcmp(extractstation,'-'))
    disp(['Merge-single-station mode is on. Only merge model for: ',extractstation]);
    singlestation=1;
    yrangeauto=1;
    yrangeauto2=1;
    yrangeauto3=1;
    if(nargin>=6)
        yrangeauto1=0;
        if(nargin>=7)
            yrangeauto2=0;
            if(nargin==8)
                yrangeauto3=0;
            end
        end
    end
end

savetofileflag=1;
%use '-' for outfile name to turn off save to file.
if(strcmp(outfilename,'-'));savetofileflag=0;end

if(strcmp(plotornot,'yes') || strcmp(plotornot,'plotyes') || strcmp(plotornot,'yesplot'));close all;end

[model1,header1]=read_layered_vmodel(modelfile1);
[model2,header2]=read_layered_vmodel(modelfile2);
% size(header1)
% size(header2)
headerall=[header1, header2];
if(~singlestation)
    stations_unique(1)={headerall(1).station};

    for i=2:length(headerall)
        if(~strcmp(headerall(i).station,headerall(i-1).station))
            stations_unique=[stations_unique;{headerall(i).station}];
        end
    end
    stations_unique=unique(stations_unique);
else
    stations_unique={extractstation};
end
%nstation=size(stations_unique)
stations_unique=sort(stations_unique);
nstation=length(stations_unique);

if(savetofileflag); fid=fopen(outfilename,'w');end

for k=1:nstation
    
    %if(strcmp(char(stations_unique(k)),
%test for first station
    foundin1=0;
    foundin2=0;
    for i0=1:length(header1)
        if(strcmp(header1(i0).station,char(stations_unique(k))))
            n1=header1(i0).nlayers;
            m1=model1(i0,1:n1,:);
            foundin1=1;
            break;
        end
    end
    for i0=1:length(header2)
        if(strcmp(header2(i0).station,char(stations_unique(k))))
        %if(strcmp(header2(i0).station,'LG19')) %example of merge model1 with all LG19 for model2 
                                                %(which equals to a 1d
                                                %model). The user needs to
                                                %manually modify this line.
            n2=header2(i0).nlayers;
            m2=model2(i0,1:n2,:);
            foundin2=1;
            break;
        end
    end
    sta=char(stations_unique(k));
    if(foundin1) && (foundin2)
        if(strcmp(plotornot,'yes') || strcmp(plotornot,'plotyes') || strcmp(plotornot,'yesplot'))
            figure(k);
            subplot(1,3,1)
            [z1,vmodel1]=plot_layered_vmodel(m1);
            %title([sta,', ',modelfile1]);
            title([sta,', Model 1']);
            if(~yrangeauto1); ylim(yrange1);end
            
            subplot(1,3,2)
            [z2,vmodel2]=plot_layered_vmodel(m2);
            %title([sta,', ',modelfile2]);
            title([sta,', Model 2']);
            if(~yrangeauto2); ylim(yrange2);end
        else
            [z1,vmodel1]=plot_layered_vmodel(m1,'noplot');
            [z2,vmodel2]=plot_layered_vmodel(m2,'noplot');
        end
        zall=z1;
        vmodelall=vmodel1;

        for i=1:length(z2)
            if(z2(i)==zall(length(z1)))
               zall=[zall,z2(i)]; 
               vmodelall=[vmodelall;vmodel2(i,:)];
            elseif(z2(i)>zall(length(z1)))
               if(z1(length(z1)-z1(length(z1)-1)>0.001))  
                   zall=[zall,z1(length(z1))]; 
                   vmodelall=[vmodelall;vmodel2(i,:)];
               else  %1 meter thin layer is ignored here.
                   zall=zall(1:(length(z1)-1)); % trim last member.
                   vmodelall=vmodelall(1:(length(z1)-1),:);
                   %zall(length(zall))=zall(length(zall)-1);
                   vmodelall(length(zall),:)=vmodel2(i,:); %replace last member with model2 value below this point.
                   zall=[zall,z2(i)];
                   vmodelall=[vmodelall;vmodel2(i,:)];
               end
               for j=(i+1):length(z2)
                   zall=[zall,z2(j)]; 
                   vmodelall=[vmodelall;vmodel2(j,:)];
               end
               break;
            end
        end
    elseif(foundin1)
        vmodelall=vmodel1;
        zall=z1;
    elseif(foundin2)
        vmodelall=vmodel2;
        zall=z2;
    else
        error(['ERROR: model information not found for station: ',sta]);
    end
    
    [vmodel_merged,N]=convert_layered_vmodel(vmodelall,zall);

    if(strcmp(plotornot,'yes') || strcmp(plotornot,'plotyes') || strcmp(plotornot,'yesplot'))
        subplot(1,3,3);
        plot_layered_vmodel(vmodel_merged);
        title([sta,', merged']);
        if(~yrangeauto3); ylim(yrange3);end
    end
%     disp(['Saving merged model for station: ',sta,' ... [',num2str(k),'/',num2str(nstation)]);
    if(savetofileflag)
        fprintf(fid,'%s   %d\n',sta,N);
        [c1,c2,c3]=size(vmodel_merged);
        for jj=1:c2
            fprintf(fid,'  %g  %g  %g  %g\n',vmodel_merged(1,jj,:));
        end
    end
end
if(savetofileflag)
    disp(['Saved ',num2str(nstation),' stations to file: ',outfilename]);
    fclose(fid);
end
end