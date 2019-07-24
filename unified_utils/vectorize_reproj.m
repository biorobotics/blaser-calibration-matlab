function [v_out] = vectorize_reproj(A, world, st, summarize)
%VECTORIZE Summary of this function goes here
%   Detailed explanation goes here

v_out = [];

for i=1:size(A,1)
    if size(A{i,2}, 1) == 0
        continue
    end
    v_out = [v_out; reproj_err(A{i,2}, world, st)];  
end

if exist('summarize','var') && summarize
    v_out = norm(v_out);
end



end
