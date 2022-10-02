%Code for finding ply failure loads using complete degradation

clc;
clear all;

%Properties of plies in the laminate
n = 8;       %No.of layers in the laminate
E1 = [38.6e9, 38.6e9, 38.6e9, 38.6e9, 38.6e9, 38.6e9, 38.6e9, 38.6e9];    %Longitudinal Youngs modulus in Pa. for ply 1 to ply n
E2 = [8.27e9, 8.27e9, 8.27e9, 8.27e9, 8.27e9, 8.27e9, 8.27e9, 8.27e9];    %Transverse Youngs modulus in Pa. for ply 1 to ply n
nu12 = [0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28];                  %Major poissons ratio for ply 1 to ply n
G12 = [4.14e9, 4.14e9, 4.14e9, 4.14e9, 4.14e9, 4.14e9, 4.14e9, 4.14e9];   %Shear modulus in Pa. for ply 1 to ply n

%Calculating minor poissons ratio
for i = 1:n
    nu21(i) = (E2(i)/E1(i))*nu12(i);
end

S1Tu = [1062e6, 1062e6, 1062e6, 1062e6, 1062e6, 1062e6, 1062e6, 1062e6];  %Longitudinal tensile strength in Pa. for ply 1 to ply n
S1Cu = [610e6, 610e6, 610e6, 610e6, 610e6, 610e6, 610e6, 610e6];          %Longitudinal compressive strength in Pa. for ply 1 to ply n
S2Tu = [31e6, 31e6, 31e6, 31e6, 31e6, 31e6, 31e6, 31e6];                  %Transverse tensile strength in Pa. for ply 1 to ply n
S2Cu = [118e6, 118e6, 118e6, 118e6, 118e6, 118e6, 118e6, 118e6];          %Transverse compressive strength in Pa. for ply 1 to ply n
T12u = [72e6, 72e6, 72e6, 72e6, 72e6, 72e6, 72e6, 72e6];                  %Shear strength in Pa. for ply 1 to ply n

a1 = [8.6e-6, 8.6e-6, 8.6e-6, 8.6e-6, 8.6e-6, 8.6e-6, 8.6e-6, 8.6e-6];    %Longitudinal coefficient of thermal expansion in m/m/deg C for ply 1 to ply n
a2 = [22.1e-6, 22.1e-6, 22.1e-6, 22.1e-6, 22.1e-6, 22.1e-6, 22.1e-6, 22.1e-6];   %Transverse coefficient of thermal expansion in m/m/deg C for ply 1 to ply n
b1 = [0, 0, 0, 0, 0, 0, 0, 0];                                            %Longitudinal coefficient of moisture expansion in m/m/deg C for ply 1 to ply n
b2 = [0, 0, 0, 0, 0, 0, 0, 0];                                            %Transverse coefficient of moisture expansion in m/m/deg C for ply 1 to ply n

%CTE & CME vectors for plies in material axes(Column i = for ply i) 
for i= 1:n
    alpham(1,i) = a1(i);
    alpham(2,i) = a2(i);
    alpham(3,i) = 0;
    betam(1,i) = b1(i);
    betam(2,i) = b2(i);
    betam(3,i) = 0;
end

dT = 50;         %Change in temperature
dC = 0;          %Change in moisture concentration

theta = [0; 45; -45; 90; 90; -45; 45; 0];                   %Orientation w.r.t to global axes for ply 1 to ply n
t = [0.125e-3; 0.125e-3; 0.125e-3; 0.125e-3; 0.125e-3; 0.125e-3; 0.125e-3; 0.125e-3];   %thickness of plies from ply 1 to ply n

%Reduced stiffness matrix[Q](Q(:,:,i) = for ply i)
for i = 1:n
    Q11 = E1(i)/(1-nu12(i)*nu21(i));
    Q12 = E2(i)*nu12(i)/(1-nu12(i)*nu21(i));
    Q22 = E2(i)/(1-nu12(i)*nu21(i));
    Q66 = G12(i);
    
    Q(:,:,i) = [Q11, Q12, 0;
                Q12, Q22, 0;
                0, 0, Q66];
end

%Transformed reduced stiffness matrix, CTE & CME vectors in global axes
for i = 1:n
    th = theta(i);
    c = cosd(th);
    s = sind(th);
    Q11b = Q(1,1,i)*c^4 + 2*s^2*c^2*(Q(1,2,i) + 2*Q(3,3,i)) + Q(2,2,i)*s^4;
    Q12b = (Q(1,1,i) + Q(2,2,i) - 4*Q(3,3,i))*s^2*c^2 + Q(1,2,i)*(s^4 + c^4);
    Q22b = Q(1,1,i)*s^4 + 2*s^2*c^2*(Q(1,2,i) + 2*Q(3,3,i)) + Q(2,2,i)*c^4;
	Q16b = (Q(1,1,i) - Q(1,2,i) - 2*Q(3,3,i))*s*c^3 + (Q(1,2,i) - Q(2,2,i) + 2*Q(3,3,i))*s^3*c;
    Q26b = (Q(1,1,i) - Q(1,2,i) - 2*Q(3,3,i))*s^3*c + (Q(1,2,i) - Q(2,2,i) + 2*Q(3,3,i))*s*c^3;
    Q66b = (Q(1,1,i) + Q(2,2,i) - 2*Q(1,2,i) - 2*Q(3,3,i))*s^2*c^2 + Q(3,3,i)*(s^4 + c^4);
    Qb(:,:,i) = [Q11b, Q12b, Q16b;
                 Q12b, Q22b, Q26b;
                 Q16b, Q26b, Q66b];
    T = [c^2, s^2, 2*s*c;                   %Transformation matrix
         s^2, c^2, -2*s*c;
        -s*c, s*c, c^2 - s^2];
    alphax(:,i) = T\alpham(:,i);            %[alphax, alphay, alphaxy/2]
    betax(:,i) = T\betam(:,i);              %[betax, betay, betaxy/2]
end

%Calculation of Z values(distance to a ply from reference plane)
Z = zeros(n+1,1);
for i = 1:n/2            %For ply 1
    Z(1) = Z(1) - t(i);
end
for i = 2:n+1
    Z(i) = Z(i-1) + t(i-1);
end

%Applied force and moment resultants
N = zeros(6,1);
N(1,1) = 100e3;              %Nx = 100 kN/m

PFL = [];        %Ply failure loads
DL = [];         %Degraded layers
FM = [];         %Failure modes
ii = 0;

Output = fopen('Complete Degradation.txt', 'w');

while (size(DL,1) < n)
    ii = ii + 1;
    fprintf (Output, 'Iteration %d\n', ii);
    fprintf (Output, '--------------\n');
    
    %Complete degradation of failed plies
    for i=1:size(DL,1)
        Q(:,:,DL(i)) = zeros(3,3);
        Qb(:,:,DL(i)) = zeros(3,3);
    end

    %Finding A, B and D matrices
    A = zeros(3,3);
    B = zeros(3,3);
    D = zeros(3,3);
    for i = 1:n
        A = A + Qb(:,:,i)*(Z(i+1)-Z(i));
        B = B + 0.5*Qb(:,:,i)*((Z(i+1))^2 - (Z(i))^2);
        D = D + (Qb(:,:,i)/3)*((Z(i+1))^3 - (Z(i))^3);
    end

    ABBD = [A, B; B, D];
    
    %Considering only mechanical loading
    EK0 = ABBD\N;                 %Mid surface strains and curvatures 
    E0 = EK0(1:3);                %Mid surface strains
    K0 = EK0(4:6);                %Mid surface curvatures

    for i=1:n
        if (i<=n/2)                    %Strains in global axes
            Ex(:,i) = E0 + Z(i)*K0;
        else
            Ex(:,i) = E0 + Z(i+1)*K0;
        end
        th = theta(i);
        c = cosd(th);
        s = sind(th);
        T = [c^2, s^2, 2*s*c;           %Transformation matrix
            s^2, c^2, -2*s*c;
            -s*c, s*c, c^2 - s^2];
        Em(:,i) = T*[Ex(1,i); Ex(2,i); Ex(3,i)/2];      %Strains in material axes - [Epsl1, Epsl2, Gamma12/2]
        Sm(:,i) = Q(:,:,i)*[Em(1,i); Em(2,i); 2*Em(3,i)];      %Stresses in material axes
    end
    fprintf (Output, 'Stresses(MPa) due to mechanical load\n');
    fprintf (Output, 'Ply 1\tPly 2\tPly 3\tPly 4\tPly 5\tPly 6\tPly 7\tPly 8\n');
    for i=1:3
        for j = 1:n
            fprintf (Output, '%.2f\t', Sm(i, j)/1e6);
        end
        fprintf (Output, '\n');
    end
    fprintf (Output, '\n');

    %Finding strength ratios using Max stress criteria
    for j=1:n
        if (Sm(1,j)> 0)
            SR(1,j) = Sm(1,j)/S1Tu(j);
        else
            SR(1,j) = -Sm(1,j)/S1Cu(j);
        end
        if (Sm(2,j)> 0)
            SR(2,j) = Sm(2,j)/S2Tu(j);
        else
            SR(2,j) = -Sm(2,j)/S2Cu(j);
        end
        SR(3,j) = abs(Sm(3,j))/T12u(j);
    end
    fprintf (Output, 'Strength ratios for above stress values\n');
    fprintf (Output, 'Ply 1\tPly 2\tPly 3\tPly 4\tPly 5\tPly 6\tPly 7\tPly 8\n');
    for i=1:3
        for j = 1:n
            fprintf (Output, '%.2f\t', SR(i, j));
        end
        fprintf (Output, '\n');
    end
    fprintf (Output, '\n');
    
    M = max(max(SR));        %Maximum strength ratio
    fprintf (Output, 'Maximum strength ratio = %.2f\n', M);
    fprintf (Output, '\n');
    
    %Finding which plies will fail and in what failure mode(FM)
    for j = 1:n
        if (SR(1,j) == M && Sm(1,j) > 0) 
            DL = [DL;j];
            FM = [FM;1];             %Longitudinal tensile
        elseif (SR(1,j) == M && Sm(1,j) < 0)
            DL = [DL;j];
            FM = [FM;2];             %Longitudinal compressive
        end
    
        if (SR(2,j) == M && Sm(2,j) > 0) 
            DL = [DL;j];
            FM = [FM;3];             %Transverse tensile
        elseif (SR(2,j) == M && Sm(2,j) < 0)
            DL = [DL;j];
            FM = [FM;4];             %Transverse compressive
        end
    
        if (SR(3,j) == M && Sm(3,j) > 0) 
            DL = [DL;j];
            FM = [FM;5];             %Positive shear
        elseif (SR(3,j) == M && Sm(3,j) < 0)
            DL = [DL;j];
            FM = [FM;6];             %Negative shear
        end
    end
    fprintf (Output, 'Degraded layers after iteration %d \n', ii);
    for i = 1:size(DL,1)
        fprintf (Output, '%d\t', DL(i));
    end
    fprintf (Output, '\n');
    fprintf (Output, 'Failure modes of degraded layers \n');
    for i = 1:size(FM,1)
        fprintf (Output, '%d\t', FM(i));
    end
    fprintf (Output, '\n\n');
    
    NT = zeros(6,1);
    NH = zeros(6,1);
    %Calculating equivalent thermal and hygroscopic loads
    for i = 1:n
        NT(1:3) = NT(1:3) + dT*Qb(:,:,i)*[alphax(1,i); alphax(2,i); 2*alphax(3,i)]*(Z(i+1)-Z(i));
        NT(4:6) = NT(4:6) + 0.5*dT*Qb(:,:,i)*[alphax(1,i); alphax(2,i); 2*alphax(3,i)]*((Z(i+1))^2 - (Z(i))^2);
        NH(1:3) = NH(1:3) + dC*Qb(:,:,i)*[betax(1,i); betax(2,i); 2*betax(3,i)]*(Z(i+1)-Z(i));
        NH(4:6) = NH(4:6) + 0.5*dC*Qb(:,:,i)*[betax(1,i); betax(2,i); 2*betax(3,i)]*((Z(i+1))^2 - (Z(i))^2);
    end

    %Finding stresses considering only Hygro-thermal loading
    EK0 = ABBD\(NT+NH);           %Mid surface strains and curvatures 
    E0 = EK0(1:3);                %Mid surface strains
    K0 = EK0(4:6);                %Mid surface curvatures
    for i=1:n
        if (i<=n/2)                    %Strains in global axes
            Ex(:,i) = E0 + Z(i)*K0;
        else
            Ex(:,i) = E0 + Z(i+1)*K0;
        end
        Ef(:,i) = dT*[alphax(1,i); alphax(2,i); 2*alphax(3,i)] + dC*[betax(1,i); betax(2,i); 2*betax(3,i)];    %Free strains due to Hygro-thermal loading
        Erx(:,i) = Ex(:,i) - Ef(:,i);      %Residual strains in global axes
        Srx(:,i) = Qb(:,:,i)*Erx(:,i);     %Residual stresses in global axes
        th = theta(i);
        c = cosd(th);
        s = sind(th);
        T = [c^2, s^2, 2*s*c;              %Transformation matrix
            s^2, c^2, -2*s*c;
            -s*c, s*c, c^2 - s^2];
        Srm(:,i) = T*Srx(:,i);            %Residual stresses in global axes
    end
    fprintf (Output, 'Residual stresses(MPa) due to Hygro-thermal load\n');
    fprintf (Output, 'Ply 1\tPly 2\tPly 3\tPly 4\tPly 5\tPly 6\tPly 7\tPly 8\n');
    for i=1:3
        for j = 1:n
            fprintf (Output, '%.2f\t', Srm(i, j)/1e6);
        end
        fprintf (Output, '\n');
    end
    fprintf (Output, '\n');

    dl = DL(size(DL,1));       %recently added degraded layer
    fm = FM(size(FM,1));       %recently added degraded layer's failure mode
    if (fm == 1)   %Longitudinal tensile
        a = S1Tu(dl) - Srm(1,dl);    %Modified longitudinal tensile strength pertaining to degraded layer
        M = Sm(1,dl)/a;              %Modified strength ratio pertaining to failure mode
    elseif (fm == 2)   %Longitudinal compressive
        a = S1Cu(dl) + Srm(1,dl);    %Modified longitudinal compressive strength pertaining to degraded layer
        M = -Sm(1,dl)/a;              %Modified strength ratio pertaining to failure mode
    elseif (fm == 3)   %Transverse tensile
        a = S2Tu(dl) - Srm(2,dl);    %Modified transverse tensile strength pertaining to degraded layer
        M = Sm(2,dl)/a;              %Modified strength ratio pertaining to failure mode
    elseif (fm == 4)   %Transverse compressive
        a = S2Cu(dl) + Srm(2,dl);    %Modified transverse compressive strength pertaining to degraded layer
        M = -Sm(2,dl)/a;              %Modified strength ratio pertaining to failure mode
    elseif (fm == 5)   %Positive shear
        a = T12u(dl) - Srm(3,dl);    %Modified shear strength pertaining to degraded layer
        M = Sm(3,dl)/a;              %Modified strength ratio pertaining to failure mode
    elseif (fm == 6)   %Negative shear
        a = T12u(dl) + Srm(3,dl);    %Modified shear strength pertaining to degraded layer
        M = Sm(3,dl)/a;              %Modified strength ratio pertaining to failure mode
    end
    
    fprintf (Output, 'Corrected strength ratio considering Hygro-thermal effects = %.2f\n', M);
    fprintf (Output, '\n');
    
    N(1,1) = N(1,1)/M;
    PFL = [PFL; N(1,1)];
    fprintf (Output, 'Ply Failure Load values(kN/m):\n', M);
    for i=1:size(PFL,1)
        fprintf (Output, '%.3f\t', PFL(i)/1e3);
    end
    fprintf (Output, '\n\n');
end
fclose(Output);