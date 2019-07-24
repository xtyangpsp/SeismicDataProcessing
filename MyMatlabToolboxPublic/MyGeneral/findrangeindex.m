function [pstart, pend]=findrangeindex(data,range)
% To find the end positions for values fall into the range.
% Very usefull in computing of confidence intervals.
% For example: for a cdf data set, you want to get the 95% confidence
% interval. just call this function like: 
% [pstart pend]=findrangeindex(cdf_data, [0.025 0.975])
% pstart and pend are the start and end positions for the corresponding
% range. Use this two indeces to find the corresponding for the parameter
% you got cdf for.
%
% By Xiaotao Yang @ Indiana University Bloominton
% May 29, 2013

data2=sort(data);
if data2 ~= data
    error('Input data has to be sorted by increasing order!')
end

vmin=range(1);
vmax=range(2);

if vmin>vmax
    error('input range is not correct. range=[min max]');
end

if min(data) >= vmin
    pstart=1;
else
    for i=1:(length(data2)-1)
        if data2(i)< vmin && data2(i+1)>= vmin
            pstart=i+1; 
            break;
        else if data2(i)== vmin
                pstart=i; break;
            end
        end
    end
end

if max(data) <= vmax
    pend=length(data);
else
    for i=1:length(data2)
        if data2(i) <= vmax && data2(i+1) >= vmax
            %disp('test');
            pend=i; 
            break; 
        end
    end
end

return;
end