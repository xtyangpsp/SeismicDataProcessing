function gmtlines2mat(infilebaselist)
%this program converts gmt formated line segments to matlab *.MAT file. the
%resultant *.MAT file is a cell array containing all of the line segments.
%
%By Xiaotao Yang @ UMass Amherst
%Created 2018.03.02
% 20190128: 
% 1. Improved for mass processing of many files. 
% 2. Changed to function to avoid duplicating codes.
%
% infilebaselist={'CentralMaine','LaurentiaMargin','WestHartfordBasinMargin',...
%     'EastHartfordBasinMargin','WestAvalonMargin'};
for m = 1:length(infilebaselist)
    infilename=[infilebaselist{m}, '.txt'];
    outfilename=[infilebaselist{m}, '.mat'];
    lineinfo={};
    fidin=fopen(infilename,'r');
    dline=fgetl(fidin);
    tag='-';
    nseg=0; %number of segments.
    while ischar(dline)
        if ~strcmp(dline(1),'#')
            clear idxc;
            idxc=strfind(dline,'>');
            if ~isempty(idxc)
               if nseg>0
    %                nseg
                   lineinfo{nseg}.tag=tag;
    %                size(data0)
                   lineinfo{nseg}.data=data0;
    %                size(faults{nseg}.data)
               end
               nseg=nseg+1;
               clear idxs;
               idxs=strfind(dline,'"');
               tag=dline(idxs(1)+1:idxs(2)-1);
               i=1; %start count data
               clear data0;
            else
               data0(i,1:2)=sscanf(dline,'%f %f\n');
    %            disp(dline);
    %            data0(i,1:2)
               i=i+1;
            end
        end
        clear dline;
       dline=fgetl(fidin);
       if ~ischar(dline)
            lineinfo{nseg}.tag=tag;
            %                size(data0)
            lineinfo{nseg}.data=data0;
       end
    end

    save(outfilename,'lineinfo');
end

end