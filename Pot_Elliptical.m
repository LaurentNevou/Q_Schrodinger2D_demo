a=4e-9;       %% radius in the the x-direction of the ellipse [m]
b=8e-9;       %% radius in the the y-direction of the ellipse [m]
x0=0;y0=0;    %% center positions of the ellipse [m]

idx=((X-x0)/a).^2 + ((Y-y0)/b).^2 < 1;