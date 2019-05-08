%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% last update 5Jan2018, lne %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% model activation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 for turn off
% 1 for turn on

FE_Method=1;            % Diagonalization of the Hamiltonian (FEM)
PWE_Method=1;           % Plane Wave Expansion (PWE)

saveXY=0;
saveV=0;
savePSI=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=15;                   %% number of solution asked 
Mass = 0.067;           %% effective mass, constant over all the structure...
Fx=0;%5e7;              %% Electric field [V/m] in the x-direction
Fy=0;%5e7;              %% Electric field [V/m] in the y-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Potential definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two vectors (x and y) and one matrix V0 must be defined with homogeneous grid
% x and y [meter]
% V0 [eV]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx=40;                  %% Meshing point in x-direction
Ny=50;                  %% Meshing point in y-direction
Mx=20e-9;               %% map X [m]
My=20e-9;               %% map Y [m]

x=linspace(-Mx/2,Mx/2,Nx);
y=linspace(-My/2,My/2,Ny);

[X,Y]=meshgrid(x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose your between the next 3 potentials or build your own!

% Pot_Rectangular
% Pot_Elliptical
 Pot_Hexagonal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vb=1.5;                 %% Potential barrier height[eV]
V0=(idx)*0 + (1-idx)*Vb ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% NOTHING TO CHANGE ANYMORE!!! %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V0=(Fx*X)+V0;        % adding the electric field Fx to the potential in the x-direction
V0=(Fy*Y)+V0;        % adding the electric field Fx to the potential in the y-direction
V0=V0-min(min(V0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Selection of the model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E1=[];E2=[];
display('=======================================')

if FE_Method==1
    tic
    if length(x)*length(y)>1e4
      N=length(x)*length(y);
      display(strcat('Warning: Take care, H=',num2str(N),'x',num2str(N),'elements'))
    end
    [E1,psi1] = Schroed2D_FEM_f(x,y,V0,Mass,n);  % m=cste
    display(strcat('-> Finite Elements method =',num2str(toc),'sec'))
end

if PWE_Method==1
    Nx = 64 ;        % number of points on the x grid % has to be a power of 2 (32,64,128,256,512,...) (smaller => faster)
    Ny = 32 ;        % number of points on the y grid % has to be a power of 2 (32,64,128,256,512,...) (smaller => faster)
    NGx = 11;%Nx/2-1  ;    % number of harmonics % must be at least 2 times -1 smaller than Nz (smaller => faster)
    NGy = 10;%Ny/2-1  ;    % number of harmonics % must be at least 2 times -1 smaller than Nz (smaller => faster)
    
    tic
    [E2,psi2] = Schroed2D_PWE_f(x,y,V0,Mass,n,Nx,Ny,NGx,NGy);
    display(strcat('-> PWE method =',num2str(toc),'sec'))
    if Fx~=0 || Fy~=0
      display('Warning: The PWE method is not the best for non-periodic potential')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E=nan(n,2);
E(1:length(E1),1)=E1;
E(1:length(E2),2)=E2;

display('=======================================')
display('Results:')
display('=======================================')
display(strcat('E(eV)='))
display(strcat(num2str(E)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','Potential','position',[100 100 1200 400])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,1)
hold on;grid on;
surf(x*1e9,y*1e9,V0)

colormap(jet)
colorbar
view(30,30)
%shading flat

xlabel('x (nm)')
ylabel('y (nm)')
zlabel('Energy (eV)')
title(strcat('Potential, Fx=',num2str(Fx,'%.1e'),'[V/m], Fy=',num2str(Fy,'%.1e'),'[V/m]'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,2)
hold on;grid on;
pcolor(x*1e9,y*1e9,V0)

colormap(jet)
colorbar
%shading flat

xlabel('x (nm)')
ylabel('y (nm)')
zlabel('Energy (eV)')
title(strcat('Potential, Fx=',num2str(Fx,'%.1e'),'[V/m], Fy=',num2str(Fy,'%.1e'),'[V/m]'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FE_Method==1
c=0;
ii=0;
for i=1:n
    if i>45
      break
    end
    if i==1 || i==16 || i==31
      figure('Name','FEM method','position',[100 100 1600 900])
      c=c+1;
      ii=0;
    end
    ii=ii+1;
    
    subplot(3,5,ii,'fontsize',10)
    hold on
    
    pcolor(x*1e9,y*1e9,(psi1(:,:,i)) )  
    contour(x*1e9,y*1e9,V0,1,'linewidth',3,'linecolor','w')
    
    xlabel('x (nm)')
    ylabel('y (nm)')
    title(strcat('E',num2str(i),'=',num2str(E(i,1)*1000,'%.1f'),'meV'))
    %axis equal
    shading flat
    colormap(jet)
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PWE_Method==1
c=0;
ii=0;
for i=1:n
    if i>45
      break
    end
    if i==1 || i==16 || i==31
      figure('Name','PWE method','position',[100 100 1600 900])
      c=c+1;
      ii=0;
    end
    ii=ii+1;
    
    subplot(3,5,ii,'fontsize',10)
    hold on
    
    pcolor(x*1e9,y*1e9,real(psi2(:,:,i)) )
    %pcolor(x*1e9,y*1e9,abs(psi2(:,:,i)) )  
    contour(x*1e9,y*1e9,V0,1,'linewidth',3,'linecolor','w')
    
    xlabel('x (nm)')
    ylabel('y (nm)')
    title(strcat('E',num2str(i),'=',num2str(E(i,2)*1000,'%.1f'),'meV'))
    %axis equal
    shading flat
    colormap(jet)
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data save %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveXY==1;
    x=x';y=y';
    save('data_x.txt','x','-ascii')
    save('data_y.txt','y','-ascii')
end

if saveV==1;
    save('data_V.txt','V0','-ascii')
end

if savePSI==1;
    
  if FE_Method==1
    for i=1:n
      M1 = psi1(:,:,i);
      save(strcat('data_psi',num2str(i),'_FEM.txt'),'M1','-ascii')
    end
  end
  if PWE_Method==1
    for i=1:n
      M2 = real(psi2(:,:,i));
      save(strcat('data_psi',num2str(i),'_PWE.txt'),'M2','-ascii')
    end
  end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%