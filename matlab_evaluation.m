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
alphaC = 1.023; % Ratio between standard and alpha-junction
alphaJ = 1.023;
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
dH1 = spalloc(dm,dm,dm);
dH01 = spalloc(dm,dm,2*dm);
dH10 = spalloc(dm,dm,2*dm);
dH12 = spalloc(dm,dm,2*dm);
dH21 = spalloc(dm,dm,2*dm);
dH02 = spalloc(dm,dm,2*dm);
dH20 = spalloc(dm,dm,2*dm);
dH03 = spalloc(dm,dm,2*dm);
dH30 = spalloc(dm,dm,2*dm);
dH32 = spalloc(dm,dm,2*dm);
dH23 = spalloc(dm,dm,2*dm);
dH20p = spalloc(dm,dm,2*dm);
%Cp = [C01+C12 -C12 0;-C12 C12+C23+C02 -C23;0 -C23 C23+C03]; Cinv = inv(Cp); % Definition of capacitance matrix
vq=[1 Nch Nch^2];

q = [-(Nch-1)/2:(Nch-1)/2];
q1 = linspace(1,1,Nch);

% Off-diagonal elements (except transitions through alpha-junction)
for n1 = 1:Nch-1 
    for n2 = 1:Nch-1 
        for n3 = 1:Nch
           j0 = vq * [n1-1; n2-1; n3-1] + 1; 
           j = j0 + vq * [1; 0; 0]; 
           k = j0 + vq * [0; 1; 0]; 

           H(j,k) = 10000 % -0.5 * EJ; 
           H(k,j) = 10000 % -0.5 * EJ; 
           
           dH12(j,k) = 1;
           dH21(k,j) = 1;
end; end; end
for n1 = 1:Nch 
    for n2 = 1:Nch-1 
        for n3 = 1:Nch-1
           j0= vq * [n1-1; n2-1; n3-1] + 1; 
           j = j0 + vq * [0; 0; 1]; 
           k = j0 + vq * [0; 1; 0]; 
           
           H(j,k) = 10000%3 % -0.5 * EJ; 
           H(k,j) = 10000%3 % -0.5 * EJ; 
           
           dH32(j,k) = 1; 
           dH23(k,j) = 1; 
end; end; end
for n3 = 1:Nch 
    for n2 = 1:Nch-1 
        for n1 = 1:Nch
           j0= vq * [n1-1; n2-1; n3-1] + 1; 
           j = j0; 
           k= j0 + vq * [0; 1; 0]; 
           
           H(j,k) = -0.5 * EJa; 
           H(k,j) = -0.5 * EJa; 
           
           dH02(j,k) = 1; 
           dH20p(j,k) = exp(-i * pi); 
           dH20(k,j) = 1; 
end; end; end

C02 = C01 * alphaC;
Cp = [C01+C12 -C12 0;-C12 C12+C23+C02 -C23;0 -C23 C23+C03]; 
Cinv = inv(Cp); % Definition of capacitance matrix
EJa= alphaJ * EJ;
% Diagonal elements:
for n1 = 1:Nch 
    for n2 = 1:Nch 
        for n3 = 1:Nch
           qv = [q(n1) q(n2) q(n3)]; 
           j = vq * [n1-1 n2-1 n3-1]'+1; 
           % H(j,j) = 10000;
           H(j,j) = 0.5 * qv * Cinv * qv';
           qv1 = [0 1 0]; 
           dH1(j,j) = qv1 * Cinv*qv';
end;end;end

for nfx = 1:nloop
  for nfy = 1:nloop
      %nfy=nfx;

      %f = 0.5; f1 =0.5; qq=flistx(nfx);
      %%%% f=flistx(nfx); f1=flisty(nfy);
      % display(f * 2 * pi)
      % display(f1 * 2 * pi)
      qq=0;
      %%%% f=flistx(nfx)+flisty(nfy); f1=flistx(nfx)-flisty(nfy);
      % display(f * 2 * pi)
      % display(f1 * 2 * pi)
      %f=flistx(nfx); f1=flisty(nfy);
      f = 0.3;
      f1 = 0.1;
      %alpha=flist(nf);

      %for n1 = 1:Nch for n2 = 1:Nch for n3 = 1:Nch
      %  qv = [q(n1) q(n2)+qq q(n3)];
      %  j=vq*[n1-1 n2-1 n3-1]'+1;
      %  H(j,j) = 0.5*qv*Cinv*qv';
      %  %qv1 = [1 0 0]; dH1(j,j)=qv1*Cinv*qv';
      %end;end;end;

      % Trasitions trough phase-droped-junction

      for n1 = 1:Nch-1 for n2 = 1:Nch for n3 = 1:Nch
        j0=vq*[n1-1;n2-1;n3-1]+1; j=j0; k=j0+vq*[1;0;0]; 
        H(j,k) = 10000%8%-0.5*EJ*exp(-i*2*pi*f); 
        H(k,j) = 10000%8%-0.5*EJ*exp(i*2*pi*f);
        
        dH01(j,k)=exp(i*2*pi*f);
        dH10(k,j)=exp(-i*2*pi*f);
      end;end;end;
      for n1 = 1:Nch for n2 = 1:Nch for n3 = 1:Nch-1
        j0=vq*[n1-1;n2-1;n3-1]+1; j=j0; k=j0+vq*[0;0;1]; 
        H(j,k) = 10000%-0.5*EJ*exp(i*2*pi*f1);
        H(k,j) = 10000%-0.5*EJ*exp(-i*2*pi*f1);
        
        dH03(j,k)=exp(-i*2*pi*f1); dH30(k,j)=exp(i*2*pi*f1);
      end;end;end;

      Nstat = 5;

      opts.disp=0;
      [vd,D] = eigs(H, Nstat,'sr',opts);

      % Find three lowest eigen energies
%       for k = 1:Nstat vecDre(:,k) = [D(k,k) real(vd(:,k)')]; vecDim(:,k) = [D(k,k) imag(vd(:,k)')]; end;
%       svecDre = sortrows(vecDre')';  svecDim = sortrows(vecDim')';
%       svecD = svecDre + i*svecDim;
%       ds=real(svecD(1,:)); svec = svecD(2:size(vd(:,1))+1,:);
        ds=diag(D); svec=vd;
        
      %opts.disp=0;
      %[vp,Dp] = eigs(dH01+dH10+dH12+dH21+dH23+dH32+dH30+dH03,Nstat,'sr',opts);

      %vg = vd(:,1); ve = vd(:,2);
      %Nm1(nfx) = svec(:,1)'*No1*svec(:,2); Nm2(nfx) = svec(:,1)'*No2*svec(:,2); Nm3(nfx) = svec(:,1)'*No3*svec(:,2);
      %vv = Cinv*[Nm1(nfx) Nm2(nfx) Nm3(nfx)]';
      %v2(nfx)=vv(2);

      %ds=sort(diag(D)); %ds21(nfx,nfy)=ds(2)-ds(1);

      ds1(nfx,nfy) = ds(1);
      ds2(nfx,nfy) = ds(2); 
      ds3(nfx,nfy) = ds(3); 
      ds4(nfx,nfy) = ds(4); 
      ds5(nfx,nfy) = ds(5);
      %ds4(nfx) = ds(4); ds5(nfx) = ds(5); ds6(nfx) = ds(6); ds7(nfx) = ds(7);
      %ds21(nfx,nfy) = ds(2)-ds(1);

      %dds1(nfx,nfy)=ds(2)-ds(1);

      v = svec(:,1); v2 = svec(:,2); v3 = svec(:,3); v4 = svec(:,4); v5 = svec(:,5);
      
      [Hrow, Hcol, Hval] = find(H);
      dlmwrite('matlab_dumps/H_old.txt',[Hrow-1 Hcol-1 Hval], 'delimiter', ',')
      for vecidx = 1:Nstat
          [vecrow, veccol, vecval] = find(svec(:, vecidx));
          dlmwrite(sprintf('matlab_dumps/eigenvec_old%d.txt',vecidx-1), [vecrow - 1, vecval], 'delimiter', ',')
      end
      display(ds)
      return
      %display(svec(:,1))
      %display(svec(:,2))
      %display(a)
      
      %display(v)
      % return
      %v = (vp1+vp2)/sqrt(2); v2 =(vp1-vp2)/sqrt(2);
      ph01(nfx,nfy) = v'*dH01*v;
      ph12(nfx,nfy) = v'*dH12*v;
      ph02(nfx,nfy) = v'*dH02*v;
      ph012(nfx,nfy) = v'*(dH01+dH12+dH20)*v;

%       a20_1(nfx,nfy) = vp3'*v;
%       thet=angle(a20_1(nfx,nfy));
%       a20_1(nfx,nfy)=a20_1(nfx,nfy)*exp(-i*thet);
%       a20_2(nfx,nfy) = vp4'*v;
%       a20_2(nfx,nfy)=a20_2(nfx,nfy)*exp(-i*thet);

      ph03(nfx,nfy) = v'*dH03*v;
      ph32(nfx,nfy) = v'*dH32*v;

      a = v'*dH1*v2; t12(nfx,nfy)=abs(a);
      display(abs(a))
      display(a)

      return

      a = v'*dH1*v3; t13(nfx,nfy)=abs(a);
      a = v'*dH1*v4; t14(nfx,nfy)=abs(a);   
      a = v'*dH1*v5; t15(nfx,nfy)=abs(a);
      a = v2'*dH1*v3; t23(nfx,nfy)=abs(a);
      a = v2'*dH1*v4; t24(nfx,nfy)=abs(a);
      a = v2'*dH1*v5; t25(nfx,nfy)=abs(a);
      a = v3'*dH1*v4; t34(nfx,nfy)=abs(a);
      a = v3'*dH1*v5; t35(nfx,nfy)=abs(a);
      
      


      %if (nfx==nloop)&&(nfy==1) vp1 = svec(:,1); vp2 = svec(:,2); end;
      %if (nfx==26)&&(nfy==26) vp1 = svec(:,1); end;
      %if (nfx==26)&&(nfy==26) vp2 = svec(:,2); end;
      %if (nfx==36)&&(nfy==16) vp3 = (svec(:,1)-svec(:,2))/sqrt(2); end;
      %if (nfx==16)&&(nfy==36) vp4 = (svec(:,1)+svec(:,2))/sqrt(2); end;
  end
  fprintf('%d ',nfx);
end;
%f11 = figure(11); plot(flistx,ds1,flistx,ds2,flistx,ds3,flistx,ds4,flistx,ds5);
%f21 = figure(21); plot(flistx,ds2-ds1,flistx,ds3-ds1,flistx,ds4-ds1,flistx,ds5-ds1,flistx,ds6-ds1,flistx,ds7-ds1);

%f3 = figure(3); plot(flistx,abs(t12),flistx,abs(t13),flistx,abs(t23));

%f31 = figure(31); plot(flistx,abs(t12));
%f4 = figure(4); plot(flistx,abs(v2));

f13 = figure(13); surf(flistx,flisty,dds1); colormap 'Hot'; colorbar; shading flat; view(2); zlim([0 max(dds1,[],'all')]);

p2=2*pi;

f20 = figure(20); surf(flistx,flisty,angle(ph01)/p2); colorbar; shading flat; set(gca,'FontSize',20);  view(2);
f21 = figure(21); surf(flistx,flisty,angle(ph12)/p2); colorbar; shading flat; set(gca,'FontSize',20);  view(2);
f22 = figure(22); surf(flistx,flisty,angle(ph02)/p2); colorbar; shading flat; set(gca,'FontSize',20);  view(2);

%f26 = figure(26); surf(flistx,flisty,rem(angle(ph01)/p2+angle(ph12)/p2-angle(ph03)/p2-angle(ph32)/p2+3.0,1)); colorbar; shading flat; set(gca,'FontSize',20);  view(2);

% f30 = figure(30); surf(flistx,flisty,real(t12)); colorbar; shading flat; set(gca,'FontSize',20);  view(2);
% f31 = figure(31); surf(flistx,flisty,real(t13)); colorbar; shading flat; set(gca,'FontSize',20);  view(2);
% f32 = figure(32); surf(flistx,flisty,real(t23)); colorbar; shading flat; set(gca,'FontSize',20);  view(2);
% f33 = figure(33); surf(flistx,flisty,real(t14)); colorbar; shading flat; set(gca,'FontSize',20);  view(2);
% f34 = figure(34); surf(flistx,flisty,real(t24)); colorbar; shading flat; set(gca,'FontSize',20);  view(2);

f40 = figure(40); surf(flistx,flisty,ds21); colorbar; colormap('hot'); shading flat; set(gca,'LineWidth',2,'FontSize',30); zlim([0 40]);

nl=51; f41 = figure(41); plot(flistx,ds1(nl,:),flistx,ds2(nl,:),flistx,ds3(nl,:),flistx,ds4(nl,:),'LineWidth',2); set(gca,'LineWidth',1.5,'FontSize',36);
xlim([flistx(1) flistx(nloop)]);

nl=51; f42 = figure(42); plot(flistx,ds2(nl,:)-ds1(nl,:),flistx,ds3(nl,:)-ds1(nl,:),flistx,ds4(nl,:)-ds1(nl,:),'LineWidth',2); set(gca,'LineWidth',1.5,'FontSize',20); ylim([0 50]);
nl=26; f44 = figure(44); plot(flisty,ds21(:,51),'LineWidth',2); ylim([0 13]); set(gca,'LineWidth',2,'FontSize',40);
nl=51; f410 = figure(410); plot(flisty,ds1(:,nl),flisty,ds2(:,nl),flisty,ds3(:,nl),flisty,ds4(:,nl),'LineWidth',2); set(gca,'LineWidth',1.5,'FontSize',36);
xlim([flisty(1) flisty(nloop)]);
nl=51; f420 = figure(420); plot(flisty,ds2(:,nl)-ds1(:,nl),flisty,ds3(:,nl)-ds1(:,nl),flisty,ds4(:,nl)-ds1(:,nl),'LineWidth',2); set(gca,'LineWidth',1.5,'FontSize',20);

nl=51; f43 = figure(43); plot(flistx,b*t12(nl,:),flistx,b*t23(nl,:),flistx,b*t34(nl,:),flistx,b*t14(nl,:),flistx,b*t13(nl,:),flistx,b*t24(nl,:),'LineWidth',2); set(gca,'LineWidth',1.5,'FontSize',36);
xlim([flistx(1) flistx(nloop)]);
%nl=51; f44 = figure(44); plot(flisty,t13(nl,:),flisty,t24(nl,:),'LineWidth',2); set(gca,'FontSize',20);
nl=51; f430 = figure(430); plot(flisty,b*t12(:,nl),flisty,b*t23(:,nl),flisty,b*t34(:,nl),flisty,b*t14(:,nl),flisty,b*t13(:,nl),flisty,b*t24(:,nl),'LineWidth',2); set(gca,'LineWidth',1.5,'FontSize',36);
xlim([flisty(1) flisty(nloop)]);
%nl=51; f44 = figure(44); plot(flisty,t13(:,nl),flisty,t24(:,nl),'LineWidth',2); set(gca,'FontSize',20);

dd=flistx(2)-flistx(1); f45 = figure(45); plot(flisty(2:nloop-1),diff(diff(ds21(:,51)))/dd^2,'LineWidth',2); set(gca,'LineWidth',2,'FontSize',36);
