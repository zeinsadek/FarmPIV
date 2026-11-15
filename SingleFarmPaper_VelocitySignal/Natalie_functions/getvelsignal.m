function [ut] = getvelsignal(x,y,bn,rn,datapath)
% with a given coordinate, get the velocity signal at that position for all
% time

xg = x;
yg = y;

% check to make sure value is in the plane 
dataname = sprintf('Block_%d_Recording_%d_Cropped.mat',bn,rn(1));
cd(datapath)
load(dataname,'X','Y','U')
if xg < round(min(X(1,:)))
    error('x is less than plane minimum. decrease plane number')
elseif xg <= round(max(X(1,:)))
    [mx,Ix] = min(abs(X(1,:) - xg));
    [my,Iy] = min(abs(Y(:,1) - yg));
else 
    error('x is greater than plane maximum. increase plane number')
end

% get the signal for all time 
for rr = 1:length(rn)
    dataname = sprintf('Block_%d_Recording_%d_Cropped.mat',bn,rn(rr));
    load(dataname,'U')
    ut(rr,:) = U(Iy,Ix,:);
end

end