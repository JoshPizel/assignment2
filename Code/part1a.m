
global C

C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665;
    
    
L = 75;
W = 50;

G = sparse(L*W,L*W);
R = zeros(L*W,1);

for i =1:1:L
    for j =1:1:W
        n = j+(i-1)*W;
        nxm = j+(i-2)*W;
        nxp = j+i*W;
        nyp = j+1+ (i-1)*W;
        nym = j-1+ (i-1)*W;
        
        if(i==1)
            G(n,:) = 0;
            G(n,n) = 1;
            R(n) = 1;
        elseif(i==L)
            G(n,:) = 0;
            G(n,n) = 1;
            R(n) = 0;
        elseif(j==1)
            G(n,:) = 0;
            G(n,n) = -3;
            G(n,nyp) = 1;
            G(n,nxp) = 1;
            G(n,nxm) = 1;
        elseif(j==W)
            G(n,:) = 0;
            G(n,n) = -3;
            G(n,nym) = 1;
            G(n,nxp) = 1;
            G(n,nxm) = 1;
        else          
            G(n,:)=0;
            G(n,nxm) = 1; 
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1;
            G(n,n) = -4;
        end
    end
end

V = G\R;

%remap
K = zeros(L,W);
for i =1:1:L
    for j =1:1:W
        n=j+(i-1)*W;
        K(i,j) =V(n);
    end
end

figure(1);
surf(K)
colorbar
title('Voltage map of V(x=0)=1, V(x=L)=0'),xlabel('Width'),ylabel('Length'),zlabel('Voltage');
