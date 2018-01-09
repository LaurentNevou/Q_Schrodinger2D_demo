Rx=4e-9;       %% radius in the the x-direction of the ellipse [m]
Ry=8e-9;       %% radius in the the y-direction of the ellipse [m]
x0=0;y0=0;     %% center positions of the ellipse [m]

idx=((X-x0)/Rx).^2 + ((Y-y0)/Ry).^2 < 1;
