function x = Pol2Cart( Bearing, Range )
%CART2POL Converts 2D polar bearing and range to cartesian x, y

% Bearing measured anti-clockwise from east, in the interval (-pi, pi]

x = zeros(1,2);

x(1) = Range*cos(Bearing);
x(2) = Range*sin(Bearing);

end

