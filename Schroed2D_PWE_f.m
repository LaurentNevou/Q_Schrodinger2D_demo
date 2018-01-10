function[E,psi]=Schroed2D_PWE_f(x,y,V0,Mass,n,Nx,Ny,NGx,NGy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant [J.s]
hbar=h/(2*pi);
e=1.602176487E-19;              %% electron charge [C]
me=9.10938188E-31;              %% electron mass [kg]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Interpolation on a grid that have 2^N points %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NGx = 2*floor(NGx/2);           %% round to lower even number
NGy = 2*floor(NGy/2);           %% round to lower even number

[X,Y] = meshgrid(x,y);
xx=linspace(x(1),x(end),Nx);
yy=linspace(y(1),y(end),Ny);
[XX,YY] = meshgrid(xx,yy);
%size(V0)
V=interp2(X,Y,V0,XX,YY);
%size(V)

dx=x(2)-x(1);
dxx=xx(2)-xx(1);
dy=y(2)-y(1);
dyy=yy(2)-yy(1);
Ltotx=xx(end)-xx(1);
Ltoty=yy(end)-yy(1);

[XX,YY] = meshgrid(xx,yy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Building of the potential in Fourier space %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vk = fftshift(fft2(V))*dxx*dyy/Ltotx/Ltoty;
Vk =Vk(Ny/2-NGy+1:Ny/2+NGy+1 , Nx/2-NGx+1:Nx/2+NGx+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Reciprocal lattice vectors %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gx = (-NGx/2:NGx/2)'*2*pi/Ltotx;
Gy = (-NGy/2:NGy/2)'*2*pi/Ltoty;

NGx=length(Gx);
NGy=length(Gy);
NG=NGx*NGy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Building Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx_x = repmat((1:NGx), [NGy 1 ]);
idx_x = idx_x(:);

idx_y = repmat((1:NGy)', [1 NGx]);
idx_y = idx_y(:);

idx_X = (idx_x-idx_x') + NGx;      %% work only in Octave
idx_Y = (idx_y-idx_y') + NGy;      %% work only in Octave
%idx_X = (repmat(idx_x,[1 NG])-repmat(idx_x',[NG 1])) + NGx;     %% work in Octave and Matlab
%idx_Y = (repmat(idx_y,[1 NG])-repmat(idx_y',[NG 1])) + NGy;     %% work in Octave and Matlab

idx = sub2ind(size(Vk), idx_Y(:), idx_X(:));
idx = reshape(idx, [NG NG]);

GX = diag(Gx(idx_x));
GY = diag(Gy(idx_y));

D2 = GX.^2 + GY.^2 ;

H =  hbar^2/(2*me*Mass)*D2  +  Vk(idx)*e ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Solving Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = sparse(H);
[psik, Ek] = eigs(H,n,'SM');
E = diag(Ek)  / e;
%E=abs(E);
E=real(E);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Transforming & Scaling the waves functions %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:n
    PSI = reshape(psik(:,i),[NGy,NGx]);
    PSI = invFFT2D(PSI,Ny,Nx)/(dxx*dyy) ;
    psi_temp = interp2(XX,YY,PSI,X,Y);
    psi(:,:,i) = psi_temp / sqrt( trapz( y' , trapz(x,abs(psi_temp).^2 ,2) , 1 )  );  % normalisation of the wave function psi
end

% in Octave, the order of the eigen values are reversed...
%psi=psi(:,:,end:-1:1);
%E=E(end:-1:1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Vxy] = invFFT2D(Vk2D,Ny,Nx)

Nkx=length(Vk2D(1,:));
Nky=length(Vk2D(:,1));

Nx1=Nx/2-floor(Nkx/2);
Nx2=Nx/2+ceil(Nkx/2);
Ny1=Ny/2-floor(Nky/2);
Ny2=Ny/2+ceil(Nky/2);

Vk2D00=zeros(Ny,Nx);
Vk2D00( Ny1+1:Ny2 , Nx1+1:Nx2)=Vk2D;
Vxy=ifft2(ifftshift(Vk2D00));

end
