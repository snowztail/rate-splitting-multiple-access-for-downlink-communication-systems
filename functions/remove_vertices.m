function [validRate] = remove_vertices(rate)
% Function:
%   - remove vertices that do not contribute to the convex hull
%   - remove similar edge points and keep the best point
%
% InputArg(s):
%   - rate (instance * user): all achievable user rate points
%
% OutputArg(s):
%   - validRate (validInstance * user): optimal rate points that contribute to the convex hull
%
% Comment(s):
%   - please adjust threshold if necessary
%
% Author & Date: Yang (i@snowztail.com) - 31 Dec 19


% any rate smaller than this value will be cleared
threshold = 1e-3;
% clear insignificant entries
rate(rate <= threshold) = 0;
% find indexes of all edge points
edgeIndex = find(~all(rate, 2));
% find indexes of optimal edge points
[~, optIndex] = max(rate);
% remove similar edge points and keep the optimal point
rate(setdiff(edgeIndex, optIndex), :) = [];
% add origin for valid convex hull
rate(end + 1, :) = 0;
% find all point indexes that contribute to the convex hull
index = convhull(rate);
% obtain corresponding points
validRate = rate(index, :);

end
