function [V,I]=maxwhere(A)
% find the maximum value in all of the elements of A and the indices.
% If there are more than one max values, only returns the position of the first meet one.
% USAGE: V=maxiwhere(A), gives you the maximum element in A;
%        [V,I]=maxwhere(A), returns the maximum value as well as it
%        position in A.
%Xiaotao Yang Oct 30, 2012 @ Indiana University

[C0,I0]=max(A);

[V,I1]=max(C0);

I(2)=I1;
I(1)=I0(I1);

return
end