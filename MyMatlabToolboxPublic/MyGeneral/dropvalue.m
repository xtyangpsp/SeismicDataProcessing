function r=dropvalue(d,v,e)
% Drop elements with the given value.
% INPUT:
%   d -> dadta vector;
%   v -> value to be dropped
%   e -> error, determines the range of the dropped values.
% OUTPUT:
%   r -> remained vector after dropping values inside [v-e v+e].
% Xiaotao Yang
% April 11, 2013
minv=v-e;
maxv=v+e;

j=1;
for i=1:length(d)
   if(d(i)<minv || d(i) >maxv)
       r(j)=d(i);
       j=j+1;
   end
end

return;
end