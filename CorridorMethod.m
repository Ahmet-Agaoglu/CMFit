function V = CorridorMethod(x,lb,ub)
% CorridorMethod: Determines the Feasible Parameter Set (FPS) for a lineUb
% fitting problem.


% INPUT:
% x  : Independent variable values
% lb : Lower bounds of the dependent variable values
% ub : Upper bounds of the dependent variable values

% OUTPUT:
% V       : Vertices of the FPS (2D matrix, each row is a vertex [intercept, slope])


% Initialization
% Sort values with respect to x.
Mat = [x(:),lb(:),ub(:)];
Mat = sortrows(Mat,1,"ascend");
x = Mat(:,1); lb = Mat(:,2); ub = Mat(:,3);

%% Main Script
% Determine the vertices of the lower bound (VL)
[VL,a,b] = findvertices_VL(x,lb,ub);
% Determine the vertices of the upper bound (VU)
[VU,c,d] = findvertices_VU(x,lb,ub);
% Find the extreme lines (lmin and lmax) enclosing the feasible parameter set
[lmin,lmax] = find_lminmax(x,lb,ub,VL,VU,a,b,c,d);
% Construct the complete vertex list V for the FPS
V = [VL;lmin;VU;lmax];
end


%% Determination of VL (Lower Bound Vertices)
function [VL,a,b] = findvertices_VL(x,lb,ub)
    VL = []; a = []; b = [];
    errTol = -1e-12; % Tolerance to handle numerical precision issues
    i = 1;
    while i<length(x)
        % Calculate slopes between current and subsequent points
        m = (lb(i+1:end)-lb(i))./(x(i+1:end)-x(i));
        [MaxSlope, ind] = max(m); % Find the maximum slope
        J = ind+i;
        n = lb(i)-MaxSlope*x(i); % Calculate intercept
        % Check if the line satisfies the upper bound condition
        if ub-(MaxSlope*x+n)>errTol
            VL = [VL;n MaxSlope];
	        if isempty(a)
		        a = i; % Record the starting index
	        end
            b = J; % Record the ending index
        end
        i = J;
    end
end

%% Determination of VU (Upper Bound Vertices)
function [VU,c,d] = findvertices_VU(x,lb,ub)
    VU = []; c = []; d = [];
    errTol = -1e-12; % Tolerance to handle numerical precision issues
    i = 1;
    while i<length(x)
        % Calculate slopes between current and subsequent points
        m = (ub(i+1:end)-ub(i))./(x(i+1:end)-x(i));
        [MinSlope, ind] = min(m); % Find the minimum slope
        J = ind+i;
        n = ub(i)-MinSlope*x(i); % Calculate intercept
        % Check if the line satisfies the lower bound condition
        if (MinSlope*x+n)-lb>errTol
            VU = [VU;n MinSlope];
	        if isempty(c)
		        c = i; % Record the starting index
	        end
            d = J; % Record the ending index
        end
        i = J;
    end
end

%% Determination of l_min and l_max (Extreme Lines)
function [lmin,lmax] = find_lminmax(x,lb,ub,VL,VU,a,b,c,d)
    % Determine extreme lines (lmin and lmax) depending on the availability of VL and VU
    if ~isempty(VL) && ~isempty(VU)
        P1 = [x(b) lb(b)]; P2 = [x(c) ub(c)];
        P3 = [x(a) lb(a)]; P4 = [x(d) ub(d)];
    elseif isempty(VU)
        m = (ub(a+1:b-1)-lb(a))./(x(a+1:b-1)-x(a));
        [~, ind] = min(m);
        w = ind+a;
        P1 = [x(b) lb(b)]; P2 = [x(w) ub(w)];
        P3 = [x(a) lb(a)]; P4 = [x(w) ub(w)];
    elseif isempty(VL)
        m = (ub(c)-lb(c+1:d-1))./(x(c)-x(c+1:d-1));
        [~, ind] = max(m);
        z = ind+c;
        P1 = [x(z) lb(z)]; P2 = [x(c) ub(c)];
        P3 = [x(z) lb(z)]; P4 = [x(d) ub(d)];
    end
    % Calculate slopes and intercepts of lmin and lmax
    lmin_slope = (P2(2)-P1(2))/(P2(1)-P1(1));
    lmin_intercept = P1(2)-lmin_slope*P1(1);

    lmax_slope = (P4(2)-P3(2))/(P4(1)-P3(1));
    lmax_intercept = P3(2)-lmax_slope*P3(1);

    lmin = [lmin_intercept lmin_slope];
    lmax = [lmax_intercept lmax_slope];
end
