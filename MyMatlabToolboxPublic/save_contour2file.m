function save_contour2file(cdata,outfile)
% Save contour lines to file.
% USAGE: save_contour2file(cdata,outfile)
% Options:
%   cdata: contour line matrix from calling MATLAB contour() function.
%   outfile: output file name. The line segments will be saved in GMT
%   format. Segments are separated by a line starting with '>
%   value,npoints', where value is the contour level, npoints is the number
%   of points for that segment/line.

[~,c]=size(cdata);

fidout=fopen(outfile,'w');
i=1;
while i<=c
    step=int16(cdata(2,i));
    v=cdata(1,i);
    fprintf(fidout,'> %g,%d\n',v,step);
    for j =i+1:i+step
        fprintf(fidout,'%g %g\n',cdata(1,j),cdata(2,j));
    end
    i=i+step+1;
end
fclose(fidout);

disp(['Contour lines saved to: ',outfile])
end