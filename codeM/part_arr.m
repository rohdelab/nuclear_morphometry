function [out] = part_arr(arr,N)

out=[]; sz=round(length(arr)/N);
for a=1:N-1
    out{a}=randsample(arr,sz);
    arr=setdiff(arr,out{a});
end
out{N}=arr;

end