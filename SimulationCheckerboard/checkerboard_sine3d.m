function checkerboard_sine3d(ngridx,ngridy,basevalues,hscales,vscales,outfilebase,chktype)
%Generates sine wave checkerboard used by migsimulation with given grid
%size and the scale of perturbation. Assumes squre grid.
%USAGE: 
%checkerboard_sine3d(ngridx,ngridy,basevalue,hscale,vscale,outfilebase)
%checkerboard_sine3d(ngridx,ngridy,basevalue,hscale,vscale,outfilebase,chktype)
%               chktype: standard or sin3d
%
%The following core code and the concept are originally from Yinzhi Wang: 
%checksinegrid=bsxfun(@(x,y)(sin((x+y)*4)+sin((x-y)*4))*15,-pi:0.1:1.52*pi,(-pi:0.1:3.05*pi)')';
%
%Based on the above core code, I wrote this funciton to make it more
%generic.
%
%Xiaotao Yang@Indiana University
%2016.02.15
%

if(nargin==6)
    chktype='standard';
end
fn=1; %figure number
for i=1:length(basevalues)
    basevalue=basevalues(i);
    for j=1:length(hscales)     
        hscale=round(hscales(j));
        for k=1:length(vscales)
            vscale=vscales(k);
            if(strcmp(chktype,'sine3d'))
                outfilename=[outfilebase,'_SineBase',num2str(basevalue),'H',num2str(hscale),'V',num2str(vscale),'.dat'];
                outfigurename=[outfilebase,'_SineBase',num2str(basevalue),'H',num2str(hscale),'V',num2str(vscale),'.eps'];
                % vscale=30;
                % hscale=5; %number of grids
                % ngridx=101;
                % ngridy=41;
                xmin=-pi;
                xmax=pi;
                ymin=xmin*ngridy/ngridx;
                ymax=xmax*ngridy/ngridx;

                x=linspace(xmin,xmax,ngridx);
                y=linspace(ymin,ymax,ngridy);

                checkgrid=bsxfun(@(x,y)(sin((x+y)*(ngridx/(2*hscale)))+sin((x-y)*(ngridx/(2*hscale))))*vscale/2,x,y');
            elseif(strcmp(chktype,'standard')) %checkerboard
                outfilename=[outfilebase,'_ChkBase',num2str(basevalue),'H',num2str(hscale),'V',num2str(vscale),'.dat'];
                outfigurename=[outfilebase,'_ChkBase',num2str(basevalue),'H',num2str(hscale),'V',num2str(vscale),'.eps'];


                checkgrid0=vscale*((checkerboard(hscale,round(ngridy/2)+1,round(ngridx/2)+1)>0.5)-0.5)*2;

                checkgrid=checkgrid0(1:ngridy,1:ngridx);

            else
                error('ERROR: wrong chktype in input!');
            end
            % size(checksinegrid)
            checkgrid=checkgrid+basevalue;
            
            figure(fn);
            imagesc(checkgrid); colorbar;
            daspect([1 1 1]);
            title(['Base:',num2str(basevalue),',H:',num2str(hscale),',V:',num2str(vscale)]);
            saveas(gca,outfigurename,'epsc');
            %pcolor(checkgrid)
            % figure(2)
            % checkgrid=bsxfun(@(x,y)sin(x+y),-5:0.1:0,(-7:0.1:0)')';
            % imagesc(checkgrid)
            fid=fopen(outfilename,'w');

            for p=1:ngridy
                for q=1:ngridx
                    fprintf(fid,'%d   %d   %g\n',p-1,q-1,checkgrid(p,q));
                end
            end

            fclose(fid);
            fn=fn+1;
        end
    end
end
end
