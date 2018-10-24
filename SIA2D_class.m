clear, clf
tic

% physics
Lx      = 250000;
Ly      = 200000;
% Lx      =   173514.571;
% Ly      =   123888.341;

B0      = 3500;
ro      = 910.0;
g       = 9.81; 
yr      = 31556926.0;
A1      = (1.9*1e-24)*((ro*g)^3.0)*yr;
A2      = (5.7*1e-20)*((ro*g)^3.0)*yr;

% numerics
nx      = 100;
ny      = 100;
nt      = 8e4;
dt      = 1;
dx      = Lx/(nx-1);
dy      = Ly/(ny-1);
x       = -Lx/2:dx:Lx/2;
y       = -Ly/2:dy:Ly/2;

% Bedrock
[X, Y]  = meshgrid(x,y); 
% B       = B0.*exp(-X.*X/1e10 -Y.*Y/1e9);
% B       = B0.*exp(-X.*X/1e10 - 2.*Y.*Y/1e9) + B0.*exp(- 2.*X.*X/1e9 -Y.*Y/1e10);
B       = B0.*exp(-X.*X/1e10 - Y.*Y/1e9) + B0.*exp(- X.*X/1e9 -(Y-Ly/8).*(Y-Ly/8)/1e10);
S       = B;
% load('B_SN.mat');
% B       = imresize(B,[nx,ny]);
% B([1,end],:) = B([2,end-1],:);               
% B(:,[1,end]) = B(:,[2,end-1]); 
% B = medfilt2(B,[2,2]);
% mass balance
beta    = 0.01*ones(nx,ny);
% Ela     = 2000.*ones(nx,ny);
% Ela     = 3050.*ones(nx,ny) + 1500*atan(Y/Ly);  2000.*ones(nx,ny) + 3000.*tan(Y/Ly); %3350.*ones(nx,ny) + 1500*atan(Y/Ly); %
Ela     = 2150.*ones(nx,ny) + 900*atan(Y/Ly);
c       = 2.0.*ones(nx,ny);

a       = min(beta.*(S-Ela),c);

% initialization 
H       = zeros(nx,ny);
% dq      = zeros(nx-2,ny-2);
tic();
for it=1:nt
    H0     = H;
    Havg   = 0.25*(H(1:end-1,1:end-1) + H(2:end,1:end-1) + ...
                   H(2:end,2:end)     + H(1:end-1,2:end));
    Sx     = diff(S,1,1)/dx;             %% same as (S(2:end,:) - S(1:end-1,:))/dx;
    Sy     = diff(S,1,2)/dy;             %% same as (S(:,2:end) - S(:,1:end-1))/dy;
    Sx_avy = 0.5*(Sx(:,1:end-1) +  Sx(:,2:end));
    Sy_avx = 0.5*(Sy(1:end-1,:) +  Sy(2:end,:));
    SNorm  =    ((Sx_avy.^2.0)  + (Sy_avx.^2.0)).^0.5;
    D      = ((A1.*Havg.^5.0)   + (A2.*Havg.^3.0)).*(SNorm).^2.0;
    Dxy    = 0.5.*(D(:,1:end-1) + D(:,2:end));
    Dyx    = 0.5.*(D(1:end-1,:) + D(2:end,:));
    qx     = Dxy.*diff(S(:,2:end-1),1,1)/dx;
    qy     = Dyx.*diff(S(2:end-1,:),1,2)/dy;
    dqx    = diff(qx,1,1)/dx;            %% same as (qx(2:end,:) - qx(1:end-1,:))/dx;
    dqy    = diff(qy,1,2)/dy;            %% same as (qy(:,2:end) - qy(:,1:end-1))/dy;
   
%     Davg   = 0.25*(D(1:end-1,1:end-1) + D(2:end,1:end-1) + ... 
%                    D(2:end,2:end)     + D(1:end-1,2:end));
%     dt     = min(1/4.1*(min(dx*dx,dy*dy))./(1e-3 + Davg),1.0);
    dt     = min(1/8.1*(min(dx*dx,dy*dy))./(max(max(D))),1.0);
    
    H(2:end-1, 2:end-1)= H(2:end-1, 2:end-1) + dt.*(dqx + dqy + a(2:end-1,2:end-1));
    H      = max(H,0.0);
    S      = B + H;
    a      = min(beta.*(S-Ela),c);
  
    %add plot
    
    if max(max((abs(H-H0))))<1e-2   
        break;
    end
end
toc();
figure(2); imagesc(H)
Ela_data         = Ela; 
H_data           = H;
% H_data_mask      = H; 
% H_data_mask(H>0) = 1.0;
% H_data_mask(H<0) = 0.0;
H_data_mask                = H>0;

save SIA_example.mat 'B' 'H_data' 'H_data_mask' 'Ela_data'