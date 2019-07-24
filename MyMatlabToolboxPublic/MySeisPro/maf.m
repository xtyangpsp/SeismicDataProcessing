function [dout]=maf(din,npts,flag)
% Moving Average Filter (MAF).
% USAGE: dout=maf(din,npts,flag)
%
% SUBROUTINES USED:
%       conv()  --> matlab convolution function.
%
% IN:
%       din     --> input data need to be filtered.
%       npts    --> number of points used to define the filter. The actual
%                   value of npts needs to be odd when flag=2. In this
%                   case, the actual number of points will be (npts + 1).
%       flag    --> -1 for backward oneside averaging; 1 for forward oneside averaging; 2 for symetrical two-side
%                   averaging.
%
% OUT:
%       dout    --> output filtered data.
%
% Author: Xiaotao Yang
% Feb 20, 2013
% 
% History:
%       Feb 20, 2013        First wrote with flag=1 only.

if flag==2
    error('Error: flag=2 for symetrical averaging is not ready yet!');
    return;
end

if flag==-1
    error('Error: flag=-1 for backward oneside averaging is not ready yet!');
    return;
end

dout=conv(ones(npts,1)/npts,din);

return;
end