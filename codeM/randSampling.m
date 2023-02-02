function output_Index = randSampling(input_Index,N)

% rand('state',sum(100*clock));
if N>length(input_Index)
    disp('The number of sample N should be less than the total input number');
    output_Index = sort(input_Index);
    
else
    rInd = randperm(length(input_Index));
    output_Index = sort(input_Index(rInd(1:N)));    
end




