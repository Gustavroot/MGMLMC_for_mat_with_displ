% Function that builds a stochastic probing vector (spv) from a normal probing
% vector (pv) filled with 0 and 1. Where pv=0 spv=0 as well, where pv = 1
% spv=1 or spv=-1 with probability 1/2
function spv = build_spv(vec)

spv = zeros(size(vec,1), 1);
indices = find(vec ~= 0);
for l = 1:length(indices)
    num = randi([0,1]);
                
    if num >= 0.5
        spv(indices(l)) = 1;
    else
        spv(indices(l)) = -1;
    end
end
    
end