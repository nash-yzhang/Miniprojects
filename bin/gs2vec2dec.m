function out = gs2vec2dec(vec, gs)
% GS2VEC2DEC Convert a 2 level grayscale vector to the appropriate decimal number pair
%
%    This function is only appropriate for values in the range of 0 - 3
 
% Error if G is not defined.
if isempty(vec)
   error('The input vector is empty');
end

% Error if G is not a double.
if (~isa(vec, 'double') | ( any(vec > (2^(gs)-1)) | any(vec < 0) ) )
   error('G must be a gsvec, 0 - 2^(gs)-1');
end

% here we turn each value in the gsvec into a binvec itself
% suppose we are using 3 bits, e.g.
% 2 -> 0 1 1
% 4 -> 1 0 0
for j = 1:length(vec)
    %binvec(j,:) = fliplr(dec2binvec(vec(j),gs));
    binvec(j,:) = my_dec2binvec_flip(vec(j),gs);
end

% then just turn each binvec into a dec value
for j = 1:gs
    %out(j) = binvec2dec(binvec(:,j)');
    out(j) = my_binvec2dec(binvec(:,j)');
end