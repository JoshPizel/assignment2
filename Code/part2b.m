function I = part2b(pointsPerSpace)

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
    

%working area
W = 10;
L = W*3/2;

%matrices
G = sparse(L*W,L*W);
B = zeros(L*W,1);

%conductivity
s1 = 1;
s2 = 0.01;

%for mesh density
dr = 1/pointsPerSpace;

pointsL = pointsPerSpace*L;
pointsW = pointsPerSpace*W;

%resistive regions size
rL = L*1/4;
rW = W*2/5;
centreX = L/2;
centreY = W/2;

Smap = zeros(L,W);
for i =1:1:pointsL
    for j =1:1:pointsW
        n = j+(i-1)*pointsW;
        nxm = j+(i-2)*pointsW;
        nxp = j+i*pointsW;
        nyp = j+1+ (i-1)*pointsW;
        nym = j-1+ (i-1)*pointsW;
        
        if(i==1)
            G(n,:) = 0;
            G(n,n) = 1/(dr^2);
            B(n) = 1;
            Smap(i,j) = s1;
        elseif(i==pointsL)
            G(n,:) = 0;
            G(n,n) = 1/(dr^2);
            B(n) = 1;
            Smap(i,j) = s1;
        elseif(j==1)
            G(n,:) = 0;
            G(n,n) = -3/(dr^2);
            if(i/pointsPerSpace > centreX-(rL/2) && i/pointsPerSpace < centreX+(rL/2))%x boundary
                G(n,nxm) = s2/(dr^2);
                G(n,nxp) = s2/(dr^2);
                G(n,nyp) = s2/(dr^2);
                Smap(i,j) = s2;
            else
                G(n,nxm) = s1/(dr^2);
                G(n,nxp) = s1/(dr^2);
                G(n,nyp) = s1/(dr^2);
                Smap(i,j) = s1;
            end
        elseif(j==pointsW)
            G(n,:) = 0;
            G(n,n) = -3/(dr^2);
            if(i/pointsPerSpace > centreX-(rL/2) && i/pointsPerSpace < centreX+(rL/2))%x boundary
                G(n,nxm) = s2/(dr^2);
                G(n,nxp) = s2/(dr^2);
                G(n,nym) = s2/(dr^2);
                Smap(i,j) = s2;
            else
                G(n,nxm) = s1/(dr^2);
                G(n,nxp) = s1/(dr^2);
                G(n,nym) = s1/(dr^2);
                Smap(i,j) = s1;
            end
        else          
            G(n,:) = 0;
            G(n,n) = -4/(dr^2);
            %  |                         X Boundaries                                | 
            if((i/pointsPerSpace > centreX-(rL/2) && i/pointsPerSpace < centreX+(rL/2)) && ...
                    (j/pointsPerSpace > centreY+(rW/2) || j/pointsPerSpace < centreY-(rW/2)))
            %       |                                Y Boundaries                          |
                G(n,nxp) = s2/(dr^2);
                G(n,nxm) = s2/(dr^2);
                G(n,nyp) = s2/(dr^2);
                G(n,nym) = s2/(dr^2);
                Smap(i,j) = s2;
            else
                G(n,nxp) = s1/(dr^2);
                G(n,nxm) = s1/(dr^2);
                G(n,nyp) = s1/(dr^2);
                G(n,nym) = s1/(dr^2);
                Smap(i,j) = s1;
            end     
        end
    end
end

V = G\B;

%remap
Vmap = zeros(pointsL,pointsW);
for i =1:1:pointsL
    for j =1:1:pointsW
        n=j+(i-1)*pointsW;
        Vmap(i,j) =V(n);
    end
end

E = gradient(Vmap);

J = Smap.*E;

area = L*W;
Javg = sum(sum(J))/(pointsL*pointsW);
I = Javg/area;

end
