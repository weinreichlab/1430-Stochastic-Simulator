% Return the number of '1' bits. number is 1's based.

function out = MutationCount(number,NumLoci)
out = 0;
number = number - 1;
for i=1:NumLoci
    if bitand(number,2^(i-1))
        out = out + 1;  
    end
end