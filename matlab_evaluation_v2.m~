clear
nloop = 101; %flistx = linspace(0.45,0.55,nloop); flisty = linspace(0.55,0.45,nloop);
flistx = linspace(0.37,0.63,nloop);flisty = linspace(0.37,0.63,nloop);
EC = 27; % 1/(6450*5.7E-15)/1E9; %EC=(2e)^2/CJ in GHz
EJ = 91; % Josephson energy of the standard junction in GHz
% Ratio between standard and alpha-junction
alphaJ = 1.023; alphaC = 1.023;
Nch = 3; % Number of5 charge states used for calculations. Generally (EC*((Nch-1)/2)^2)/EJ > 1
C01=1/EC;C12=C01;C23=C01;C03=C12;
C02=C01*alphaC;
EJa=alphaJ*EJ;
b=2;

dm=Nch*Nch*Nch; H=spalloc(dm,dm,12*dm); H1=spalloc(dm,dm,12*dm);
P01a=spalloc(dm,dm,3*dm); P10a=spalloc(dm,dm,3*dm); P03a=spalloc(dm,dm,3*dm); P30a=spalloc(dm,dm,3*dm);
P01=spalloc(dm,dm,3*dm); P10=spalloc(dm,dm,3*dm); P12=spalloc(dm,dm,3*dm); P21=spalloc(dm,dm,3*dm);
P02=spalloc(dm,dm,3*dm); P20=spalloc(dm,dm,3*dm); P03=spalloc(dm,dm,3*dm); P30=spalloc(dm,dm,3*dm);
P32=spalloc(dm,dm,3*dm); P23=spalloc(dm,dm,3*dm); U=spalloc(dm,dm,3*dm); VT2=spalloc(dm,dm,3*dm);

%Cp = [C01+C12 -C12 0;-C12 C12+C23+C02 -C23;0 -C23 C23+C03]; Cinv = inv(Cp); % Definition of capacitance matrix
vq=[1 Nch Nch^2];

q = [-(Nch-1)/2:(Nch-1)/2];
q1 = linspace(1,1,Nch);

C02=C01*alphaC;
Cp = [C01+C12 -C12 0;-C12 C12+C23+C02 -C23;0 -C23 C23+C03]; Cinv = inv(Cp); % Definition of capacitance matrix
EJa=alphaJ*EJ;
% Diagonal elements:
for n1 = 1:Nch for n2 = 1:Nch for n3 = 1:Nch 
   qv = [q(n1) q(n2) q(n3)]; j=vq*[n1-1 n2-1 n3-1]'+1; U(j,j)=0.5*qv*Cinv*qv';
   qv1 = [0 1 0]; VT2(j,j)=qv1*Cinv*qv';
end;end;end; 

% Phi operators: P01, P02, P03
for n1 = 1:Nch-1 
    for n2 = 1:Nch 
        for n3 = 1:Nch 
            j0 = vq * [n1-1;n2-1;n3-1] + 1; 
            j = j0 + vq * [1;0;0]; 
            k = j0 + vq * [0;0;0]; 
            %% Swapped so that it matches mine
            P03(j,k)=1;  
        end; 
    end; 
end;  
for n1 = 1:Nch for n2 = 1:Nch-1 for n3 = 1:Nch j0=vq*[n1-1;n2-1;n3-1]+1; j=j0+vq*[0;1;0]; k=j0+vq*[0;0;0]; P02(j,k)=1; end; end; end;  
for n1 = 1:Nch 
    for n2 = 1:Nch 
        for n3 = 1:Nch-1 
            j0=vq*[n1-1;n2-1;n3-1]+1; 
            j=j0+vq*[0;0;1]; 
            k=j0+vq*[0;0;0]; 
            
            %% Swapped so that it matches mine
            P01(j,k)=1; 
        end; 
    end; 
end;  

%% Continue
P10=P01'; P20=P02'; P30=P03';
P12=P10*P02; P23=P20*P03; P21=P12'; P32=P23';

%The Hamiltonian except varying part of P12, P21 and P32, P23 
% H1=U-0.5*EJ*(P01+P10) - 0.5*EJ*(P30+P03) - 0.5*EJa*(P02+P20);
H1= U -0.5*EJ*(P01+P10) -0.5*EJ*(P30+P03) - 0.5*EJa*(P02+P20);

[rf, cf, vf] = find(VT2);
dlmwrite('matlab_dumps/voltage_matrix.txt',[rf-1 cf-1 vf], 'delimiter', ',')
[rf, cf, vf] = find(P01+P10);
dlmwrite('matlab_dumps/phi01_rows_and_cols.txt',[rf-1 cf-1 vf], 'delimiter', ',')
[rf, cf, vf] = find(P03+P30);
dlmwrite('matlab_dumps/phi03_rows_and_cols.txt',[rf-1 cf-1 vf], 'delimiter', ',')
[rf, cf, vf] = find(P12);
dlmwrite('matlab_dumps/phi12_rows_and_cols.txt',[rf-1 cf-1 vf], 'delimiter', ',')
[rf, cf, vf] = find(P21);
dlmwrite('matlab_dumps/phi21_rows_and_cols.txt',[rf-1 cf-1 vf], 'delimiter', ',')
[rf, cf, vf] = find(P23);
dlmwrite('matlab_dumps/phi23_rows_and_cols.txt',[rf-1 cf-1 vf], 'delimiter', ',')
[rf, cf, vf] = find(P32);
dlmwrite('matlab_dumps/phi32_rows_and_cols.txt',[rf-1 cf-1 vf], 'delimiter', ',')

for nfx = 1:nloop
%  for nfy = 1:nloop
      nfy=nfx;
      f=flistx(nfx); f1=flisty(nfy); 
      f = 0.4;
      f1 = 0.4;
      %display(f)
      %display(f1)

      for n1 = 1:Nch for n2 = 1:Nch for n3 = 1:Nch 
        qv = [q(n1) q(n2) q(n3)]; 
        j=vq*[n1-1 n2-1 n3-1]'+1; 
        H(j,j)=0.5*qv*Cinv*qv'; 
        %qv1 = [1 0 0]; P1(j,j)=qv1*Cinv*qv';
      end;end;end; 

      % Trasitions trough phase-droped-junction 
      % f = 0.422;
      % f1 = 0.1;
      H = H1 - 0.5*EJ* (P12*exp(-i*2*pi*f) + P21*exp(i*2*pi*f)) - 0.5*EJ*(P23*exp(-i*2*pi*f1) + P32*exp(+i*2*pi*f1));
      % H = H1 - 0.5*EJ*(P12*exp(-i*2*pi*f)) - 0.5*EJ*P21*exp(i*2*pi*f) + 10003*P23 + 10004*P32;      % H = H1 + 10001*P12 + 10002*P21 + 10003*P23 + 10004*P32;

      Nstat = 6;
      opts.disp=0;
      [vd,D] = eigs(H,Nstat,'sr',opts);
      ds=diag(D); svec=vd;
        
      [rf, cf, vf] = find(H);
      dlmwrite('matlab_dumps/H_new.txt',[rf-1 cf-1 vf], 'delimiter', ',')
      [rf, cf, vf] = find(ds);
      dlmwrite('matlab_dumps/eigvals_new.txt',[rf, vf], 'delimiter', ',')
      for vecidx = 1:Nstat
          [vecrow, veccol, vecval] = find(svec(:, vecidx));
          dlmwrite(sprintf('matlab_dumps/eigenvec_new%d.txt',vecidx-1), [vecrow - 1, vecval], 'delimiter', ',')
      end
 
      
      ds1(nfx) = ds(1); 
      ds2(nfx) = ds(2); ds3(nfx) = ds(3); ds4(nfx) = ds(4); ds5(nfx) = ds(5);
      ds21(nfx) = ds(2)-ds(1);
      dds1(nfx)=ds(2)-ds(1);
      
      v = svec(:,1); v2 = svec(:,2); v3 = svec(:,3); v4 = svec(:,4); v5 = svec(:,5);
      
      ph01(nfx) = v'*P01*v; ph12(nfx) = v'*P12*v; ph02(nfx) = v'*P02*v; ph012(nfx) = v'*(P01+P12+P20)*v;           
      ph03(nfx) = v'*P03*v; ph32(nfx) = v'*P32*v; 
      
      a = v'*VT2*v2; t12(nfx)=abs(a); 
      display(abs(a));
      dlmwrite(sprintf('matlab_dumps/transition.txt',vecidx-1), [vecrow - 1, vecval], 'delimiter', ',')
      
      a = v'*VT2*v3; t13(nfx)=abs(a); 
      a = v'*VT2*v4; t14(nfx)=abs(a); 
      a = v'*VT2*v5; t15(nfx)=abs(a); 
      a = v2'*VT2*v3; t23(nfx)=abs(a); 
      a = v2'*VT2*v4; t24(nfx)=abs(a); 
      a = v2'*VT2*v5; t25(nfx)=abs(a); 
      a = v3'*VT2*v4; t34(nfx)=abs(a); 
      a = v3'*VT2*v5; t35(nfx)=abs(a); 
      return
%  end;
  fprintf('%d ',nfx);
end;

p2=2*pi;

% f1 = figure(1); plot(flistx,ds1,flistx,ds2,flistx,ds3,flistx,ds4,'LineWidth',2); set(gca,'LineWidth',1.5,'FontSize',36); 
% xlim([flistx(1) flistx(nloop)]);
% f11 = figure(11); plot(flistx,ds2-ds1,flistx,ds3-ds1,flistx,ds4-ds1,'LineWidth',2); set(gca,'LineWidth',1.5,'FontSize',36); 
% xlim([flistx(1) flistx(nloop)]);

%f2 = figure(2); plot(flistx,b*t12,flistx,b*t23,flistx,b*t34,flistx,b*t14,flistx,b*t13,flistx,b*t24,'LineWidth',2); set(gca,'LineWidth',1.5,'FontSize',36);
f2 = figure(2); plot(flistx,b*t12,flistx,b*t23,'LineWidth',2); set(gca,'LineWidth',1.5,'FontSize',36);
xlim([flistx(1) flistx(nloop)]);