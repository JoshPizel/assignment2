clearvars
clearvars -GLOBAL
close all 


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
        elseif(i==L)
            G(n,:) = 0;
            G(n,n) = 0;
        elseif(j==1)
            G(n,:) = 0;
            G(n,n) = 1;
        elseif(j==W)
            G(n,:) = 0;
            G(n,n) = 1;
        else          
            G(n,:)=0;
            G(n,nxm) = 1; 
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1;
            G(n,n) = -(G(n,nxm)+G(n,nxp)+G(n,nym)+G(n,nyp));
        end
    end
end

[V,D] = eigs(G);

%remap
Vs = zeros(L,W);
for i =1:1:L
    for j =1:1:W
        n=j+(i-1)*W;
        
        Vs(i,j) =V(n,1);
        %+V(n,2)+V(n,3)+V(n,4)+V(n,5)+V(n,6)
    end
end

%http://www.southampton.ac.uk/~feeg6002/lecturenotes/feeg6002_numerical_methods05.pdf

F = zeros(L,W);
F(1,:) =1;

Vs2 =F./Vs;

surf(Vs2);

