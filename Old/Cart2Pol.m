function [Bearing, Range] = Cart2Pol( x )
%CART2POL Converts 2D cartesian coordinate vector to bearing and range

% Bearing measured anti-clockwise from east, in the interval (-pi, pi]

Bearing = atan(x(2)/x(1));
if x(1)<0
    Bearing = Bearing + pi;
end
if Bearing > pi
    Bearing = Bearing - 2*pi;
end

assert( (Bearing<pi)&&(Bearing>-pi), 'Invalid bearing');

Range = sqrt(sum(x.^2));

end

