function I=findclosevalue(A,x)
% To find the closest index for the x in A. A is a vector.

for i= 1:length(A)
    d(i)=abs(A(i)-x);

end

[V I]=min(d);

return;
end