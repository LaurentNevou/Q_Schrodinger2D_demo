a=8e-9;     %% side length of the hexagon [m]

idx1=  (abs(X)<a*sqrt(3)/2);
idx2=(tan(pi/6)*X+a>Y) .* (tan(pi/6)*X-a<Y) .* (-tan(pi/6)*X-a<Y) .* (-tan(pi/6)*X+a>Y);
idx=idx1.*idx2;

