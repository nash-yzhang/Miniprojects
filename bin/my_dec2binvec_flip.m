function out = my_dec2binvec_flip(dec, bits)

out=zeros(1,bits);
for m=bits-1:-1:0
    tmp = 2^m;
    if (dec==0)
        break;
    elseif (dec >= tmp)
        out(bits-m)=1;
        dec = dec - tmp;
    end
end