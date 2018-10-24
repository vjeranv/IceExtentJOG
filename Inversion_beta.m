clear, figure(1), clf
tic
load('SIA_example.mat')

%% physics
Lx      = 250000;
Ly      = 200000;
B0      = 3500;
ro      = 910.0;
g       = 9.81; 
yr      = 31556926.0;
A1      = (1.9*1e-24)*((ro*g)^3.0)*yr;
A2      = (5.7*1e-20)*((ro*g)^3.0)*yr;
%% numerics
nx      = 100;
ny      = 100;
nt      = 8e4;
dt      = 1.0;
dx      = Lx/(nx-1);
dy      = Ly/(ny-1);
x       = -Lx/2:dx:Lx/2;
y       = -Ly/2:dy:Ly/2;
%% plot numerics
x0     = 10;
y0     = 10;
width  = 1024;
height = 1024;
%% Bedrock
[X, Y]  = meshgrid(x,y); 
% B       = B0.*exp(-X.*X/1e9 - Y.*Y/1e9);
% B       = B0.*exp(-(X + 0.5.*Lx).*(X + 0.5.*Lx)/1e10 -Y.*Y/1e10)+ B0.*exp(-(X-0.6.*Lx).*(X-0.6.*Lx)/1e10 -Y.*Y/1e10);
B       = B0.*exp(-X.*X/1e10 - Y.*Y/1e9) + B0.*exp(- X.*X/1e9 -(Y-Ly/8).*(Y-Ly/8)/1e10);
S       = B;
%% mass balance guess
Ela     = 2150.*ones(nx,ny) + 900*atan(Y/Ly);
beta    = 0.0153 + 0.007.*rand(nx,ny);
c       = 2.0.*ones(nx,ny);
a       = min(beta.*(S-Ela),c);
%% inversion parameters 
niter   = 150;
tau1    = 1.5e-2; %1.5e-2 %440 %510s 180 iter 75m
nsmooth = 30; %40
tau2    = 0.25*min(dx*dx,dy*dy);
%% initialization 
H       = zeros(nx,ny);
res     = zeros(1,niter);
for iter  = 1:niter
  H       = zeros(nx,ny);
  S       = B;
  a       = min(beta.*(S-Ela),c);
 for it=1:nt
    H0     = H;
    Havg   = 0.25*(H(1:end-1,1:end-1) + H(2:end,1:end-1) + ...
                   H(2:end,2:end)     + H(1:end-1,2:end));
    Sx     = diff(S,1,1)/dx;          %% (S(2:end,:) - S(1:end-1,:))/dx;
    Sy     = diff(S,1,2)/dy;          %% (S(:,2:end) - S(:,1:end-1))/dy;
    Sx_avy = 0.5*(Sx(:,1:end-1) +  Sx(:,2:end));
    Sy_avx = 0.5*(Sy(1:end-1,:) +  Sy(2:end,:));
    SNorm  =    ((Sx_avy.^2.0)  + (Sy_avx.^2.0)).^0.5;
    D      = ((A1.*Havg.^5.0)   + (A2.*Havg.^3.0)).*(SNorm).^2.0;
    Dxy    = 0.5.*(D(:,1:end-1) + D(:,2:end));
    Dyx    = 0.5.*(D(1:end-1,:) + D(2:end,:));
    qx     = Dxy.*diff(S(:,2:end-1),1,1)/dx;
    qy     = Dyx.*diff(S(2:end-1,:),1,2)/dy;
    dqx    = diff(qx,1,1)/dx;         %% (qx(2:end,:) - qx(1:end-1,:))/dx;
    dqy    = diff(qy,1,2)/dy;         %% (qy(:,2:end) - qy(:,1:end-1))/dy;
    dt     = min(1/4.1*(min(dx*dx,dy*dy))./(max(max(D))),1.0);
    
    H(2:end-1, 2:end-1) = H(2:end-1, 2:end-1)+dt.*(dqx + dqy + a(2:end-1,2:end-1));
    H                   = max(H,0.0);
    S                   = B + H;
    a                   = min(beta.*(S-Ela),2.0);
    if max(max(abs(H-H0)))<1e-2
        break;
    end
    if it == nt
        disp('end')
    end
end
    H_in_mask                = H>0;

    gamma                    = H_data_mask - H_in_mask; 
    beta0                    = beta;
    beta                     = beta         - tau1*(gamma);
 for ismooth=1:nsmooth
     beta(2:end-1,2:end-1)   = beta(2:end-1,2:end-1)     ...
                             + tau2*(diff(beta( : ,2:end-1),2,1)/(dx*dx) ...
                             +       diff(beta(2:end-1, :),2,2)/(dy*dy)); 
     beta([1,end],:)       = beta([2,end-1],:);               
     beta(:,[1,end])       = beta(:,[2,end-1]);  
 end
        res(iter)    = sum(sum(abs(gamma)));
     if mod(iter,10) == 0
        figure(1);
        subplot(4,1,1)
        imagesc(H_data_mask.*(beta_data-beta));colorbar,colormap(parula); set(gca,'YDir','normal'); xlabel('Lx','FontSize', 24); ylabel('Ly','FontSize', 24); title(iter); %caxis([min(min(Ela_in)),max(max(Ela_in))]);
        subplot(4,1,2)
        imagesc(gamma); colorbar; set(gca,'YDir','normal'); colorbar; colormap(parula);
        subplot(4,1,3)
        imagesc(H_in_mask.*(beta-beta0));colorbar,colormap(jet); set(gca,'YDir','normal'); xlabel('Lx','FontSize', 24); ylabel('Ly','FontSize', 24); title(iter); %caxis([min(min(Ela_in)),max(max(Ela_in))]);
        subplot(4,1,4)
        plot(iter,sum(sum(abs(gamma))),'r*'); set(gca,'fontsize',20); title(sum(sum(abs(gamma)))); hold on
        drawnow
     end
     
%      if max(max(abs((beta_data-beta)).*H_data_mask)) < 25 || sum(sum(abs(gamma))) < 10
%          break
%      end
     if sum(sum(abs(gamma))) < 10
         break
     end


end
%      figure(2); imagesc(Hmbin.*Ela_in-Hmbin.*Ela_data); colorbar; set(gca,'YDir','normal'); colorbar;
 figure(1);
        subplot(4,1,1)
        imagesc(H_data_mask.*(beta_data-beta));colorbar,colormap(parula); set(gca,'YDir','normal'); xlabel('Lx','FontSize', 24); ylabel('Ly','FontSize', 24); title(iter); %caxis([min(min(Ela_in)),max(max(Ela_in))]);
        subplot(4,1,2)
        imagesc(-gamma); colorbar; set(gca,'YDir','normal'); colorbar; colormap(parula);
        subplot(4,1,3)
        imagesc(H_in_mask.*(beta-beta0));colorbar,colormap(jet); set(gca,'YDir','normal'); xlabel('Lx','FontSize', 24); ylabel('Ly','FontSize', 24); title(iter); %caxis([min(min(Ela_in)),max(max(Ela_in))]);
        subplot(4,1,4)
        plot(iter,sum(sum(abs(gamma))),'r*'); set(gca,'fontsize',20); title(sum(sum(abs(gamma)))); hold on
       
toc()
