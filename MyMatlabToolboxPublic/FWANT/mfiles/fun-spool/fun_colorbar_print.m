function cid=fun_colorbar_print(sid,FIG_ROOT,fnmpart,dfmt,filefmt,xtick)

if ~isdir(FIG_ROOT), mkdir(FIG_ROOT); end

drvfmt=['-d' dfmt];
if exist('filefmt','var')
   fnmfmt=filefmt;
else
   fnmfmt=dfmt;
end

cid=colorbar('vert','location','SouthOutSide');

set(gca,'visible','off');
set(sid,'visible','off');

if exist('xtick','var')
   set(cid,'xtick',xtick);
end
scl_xtick=get(cid,'xtick');

print_south([]);

set(gcf,'InvertHardcopy','off')
print_south('wbg_');
set(gcf,'InvertHardcopy','on')

set(gca,'visible','on');
set(sid,'visible','on');

% ----------- south ------------------
function print_south(wbg)

dnm='South'; set(cid,'location',dnm); pos=get(cid,'position');
   
% whole
set(cid,'XAxisLocation','bottom');
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'B.' fnmfmt]);
set(cid,'XAxisLocation','top');
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'T.' fnmfmt]);
   
% 1/2
set(cid,'position',[pos(1),pos(2),pos(3)/2,pos(4)/2]);
if exist('xtick','var'), set(cid,'xtick',scl_xtick); end

set(cid,'XAxisLocation','bottom');
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'B2.' fnmfmt]);
set(cid,'XAxisLocation','top');
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'T2.' fnmfmt]);
   
% 1/4
set(cid,'position',[pos(1),pos(2),pos(3)/4,pos(4)/4]);
if exist('xtick','var'), set(cid,'xtick',scl_xtick); end

set(cid,'XAxisLocation','bottom');
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'B4.' fnmfmt]);
set(cid,'XAxisLocation','top');
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'T4.' fnmfmt]);
   
% 1/4+0
set(cid,'xtick',[scl_xtick(1),0,scl_xtick(end)]);
set(cid,'XAxisLocation','bottom');
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'B43.' fnmfmt]);
set(cid,'XAxisLocation','top');
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'T43.' fnmfmt]);
   
% 1/4+notick
set(cid,'xtick',[]);
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'E4.' fnmfmt]);
end

end
