function checkerboard_sine3d(ngridx,ngridy,basevalue,hscale,vscale,outfilebase)
%Generates sine wave checkerboard used by migsimulation with given grid
%size and the scale of perturbation. Assumes squre grid.
%USAGE: checkerboard_sine3d(ngridx,ngridy,basevalue,hscale,vscale,outfilebase)
%
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

checksinegrid=bsxfun(@(x,y)(sin((x+y)*(ngridx/(2*hscale)))+sin((x-y)*(ngridx/(2*hscale))))*vscale/2,x,y');

% size(checksinegrid)
checksinegrid=checksinegrid+basevalue;
imagesc(checksinegrid); colorbar;
daspect([1 1 1]);
saveas(gca,outfigurename,'epsc');
%pcolor(checksinegrid)
% figure(2)
% checksinegrid=bsxfun(@(x,y)sin(x+y),-5:0.1:0,(-7:0.1:0)')';
% imagesc(checksinegrid)
fid=fopen(outfilename,'w');

for j=1:ngridy
    for i=1:ngridx
        fprintf(fid,'%d   %d   %g\n',j-1,i-1,checksinegrid(j,i));
    end
end

fclose(fid);