% contourdata from MathWorks File Exchange to export contours
% Copyright (c) 2007, Duane Hanselman

function s=contourdata(c)
%CONTOURDATA Extract Contour Data from Contour Matrix C.
% CONTOUR, CONTOURF, CONTOUR3, and CONTOURC all produce a contour matrix
% C that is traditionally used by CLABEL for creating contour labels.
%
% S = CONTOURDATA(C) extracts the (x,y) data pairs describing each contour
% line and other data from the contour matrix C. The vector array structure
% S returned has the following fields:
%
% S(k).level contains the contour level height of the k-th line.
% S(k).numel contains the number of points describing the k-th line.
% S(k).isopen is True if the k-th contour is open and False if it is closed.
% S(k).xdata contains the x-axis data for the k-th line as a column vector.
% S(k).ydata contains the y-axis data for the k-th line as a column vector.
%
% For example: PLOT(S(k).xdata,S(k).ydata)) plots just the k-th contour.
%
% See also CONTOUR, CONTOURF, CONTOUR3, CONTOURC.
% From the help text of CONTOURC:
%   The contour matrix C is a two row matrix of contour lines. Each
%   contiguous drawing segment contains the value of the contour, 
%   the number of (x,y) drawing pairs, and the pairs themselves.  
%   The segments are appended end-to-end as
% 
%       C = [level1 x1 x2 x3 ... level2 x2 x2 x3 ...;
%            pairs1 y1 y2 y3 ... pairs2 y2 y2 y3 ...]
% D.C. Hanselman, University of Maine, Orono, ME 04469
% MasteringMatlab@yahoo.com
% Mastering MATLAB 7
% 2007-05-22
if nargin<1 || ~isfloat(c) || size(c,1)~=2 || size(c,2)<4
   error('CONTOURDATA:rhs',...
         'Input Must be the 2-by-N Contour Matrix C.')
end
tol=1e-12;
k=1;     % contour line number
col=1;   % index of column containing contour level and number of points
while col<size(c,2); % while less than total columns in c
   s(k).level = c(1,col); %#ok
   s(k).numel = c(2,col); %#ok
   idx=col+1:col+c(2,col);
   s(k).xdata = c(1,idx).'; %#ok
   s(k).ydata = c(2,idx).'; %#ok
   s(k).isopen = abs(diff(c(1,idx([1 end]))))>tol || ...
                 abs(diff(c(2,idx([1 end]))))>tol; %#ok
   k=k+1;
   col=col+c(2,col)+1;
end
