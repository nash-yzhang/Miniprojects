function out = my_binvec2dec(vec)
l=length(vec);

out=0;
for i=1:l
    out=out+vec(i)*2^(i-1);
end