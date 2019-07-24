function f=cumulative_frequency(data, vector)
% calculate the cumulative frequency smaller than the elements given in
% vector.
% data: series of data; vector: values used as index in the calculation.
% e.g. of vector: 1, 2, 3, 4, 7 , 10. The function will calculate number of
% values/frequency in percentage that greater than 1,2,....
% By Xiaotao Yang    Indiana University  Oct 24, 2012
for i=1:length(vector)
    count=0;
    for j=1:length(data)
        if(data(j) >= vector(i))
            count = count + 1;
        end
    end
    
    f(i)=count;
    
end
return;         
end
