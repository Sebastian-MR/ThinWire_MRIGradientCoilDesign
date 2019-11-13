function out = nmz(data)
% function out = nmz(data)
% normalise the input data between 0 and 1

dat = reshape(data,1,[]);
minData = min(dat);
dat = dat - min(dat);
newMaxData = max(dat);
out = (data - minData)/newMaxData;