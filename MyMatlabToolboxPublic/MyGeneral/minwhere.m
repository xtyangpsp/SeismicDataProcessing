function [V,I]=minwhere(A)
% find the minimum value in all of the elements of A and the indices.
% If there are more than one min values, only returns the position of the
% first meet one.
% USAGE: V=minwhere(A), gives you the minimum element in A;
%        [V,I]=minwhere(A), returns the minimum value as well as it
%        position in A.
%Xiaotao Yang Oct 30, 2012 @ Indiana University

[C0,I0]=min(A);

[V,I1]=min(C0);

I(2)=I1;
I(1)=I0(I1);

return
end