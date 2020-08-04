clear
nloop = 101;
%flistx = linspace(-0.13,0.13,nloop);
flistx = linspace(0, 0, 1);
flisty = linspace(0.37, 0.63, nloop);

%flistx = linspace(-0.04,0.04,nloop);flisty = linspace(0.46,0.54,nloop);
%flistx = linspace(0,1,nloop);flisty = linspace(0,1,nloop);
EC = 27;%1/(6450*5.7E-15)/1E9; %EC=(2e)^2/CJ in GHz
%EC = 1/(6450*10E-15)/1E9;
EJ = 91; % Josephson energy of the standard junction in GHz
alphaC = 1.0; % Ratio between standard and alpha-junction
alphaJ = 1.0;
Nch = 3; % Number of5 charge states used for calculations. Generally (EC*((Nch-1)/2)^2)/EJ > 1
C01 = 1/EC; 
C12 = C01; 
C23 = C01; 
C03 = C12;
C02 = C01 * alphaC;
EJa = alphaJ * EJ;

b=2;
dm=Nch*Nch*Nch;
H=spalloc(dm,dm,3*dm);

vq=[1 Nch Nch^2];

% Off-diagonal elements (except transitions through alpha-junction)
for n1 = 1:Nch-1; for n2 = 1:Nch-1; for n3 = 1:Nch
           j0 = vq * [n1-1; n2-1; n3-1] + 1; 
           j = j0 + vq * [1; 0; 0]; 
           k = j0 + vq * [0; 1; 0]; 

           fprintf('(%d, %d, %d) --> %d --> j=%d, k=%d\n', n1, n2, n3, j0, j, k);
end; end; end