function t = timeit_reps(f, r, narg)
t = 0;

if ~exist('narg')
    narg = 1;
end

for i=1:r
    t = t + timeit(f, narg)/r;
end
end