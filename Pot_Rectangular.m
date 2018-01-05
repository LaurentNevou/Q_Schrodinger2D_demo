Lx=5e-9;        %% well width X [m]
Ly=7e-9;        %% well width Y [m]
x0=0;y0=0;      %% center positions of the rectangle [m]

idx= ( abs(X-x0)<Lx ) .* ( abs(Y-y0)<Ly ) ;