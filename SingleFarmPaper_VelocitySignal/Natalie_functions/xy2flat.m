function [c] = xy2flat(x,y,datapath,bn,rn)
% convert a given x,y coordinate into the array position of 
% the flat version. 

xg = x;
yg = y;

dataname = sprintf('Block_%d_Recording_%d_Cropped.mat',bn,rn);

% cd(datapath)
% load(dataname,'X','Y')

load(fullfile(datapath, dataname), 'X', 'Y')

if xg < round(min(X(1,:)))
    error('x is less than plane minimum. decrease plane number')
elseif xg <= round(max(X(1,:)))
    [mx,Ix] = min(abs(X(1,:) - xg));
    [my,Iy] = min(abs(Y(:,1) - yg));
    c = Ix*Iy;
else 
    error('x is greater than plane maximum. increase plane number')
end

end