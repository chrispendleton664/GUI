function out = GUIFailureFaceSheetDeflection(FailABD, FailureCriteria, GUIInput)
 
t1 = GUIInput.t1LF;
t2 = GUIInput.t2LF;

pt = FailureCriteria.pt;

Xt = t2(1,1);
Xc = t2(1,2);
Yt = t2(1,3);
Yc = t2(1,4);
S = t2(1,5);
  
ext = t2(1,6); 
exc = t2(1,7);
eyt = t2(1,8);
eyc = t2(1,9);
es = t2(1,10);

  A11 = FailABD.A(1,1);
  A12 = FailABD.A(1,2);
  A16 = FailABD.A(1,3);
  A22 = FailABD.A(2,2);
  A26 = FailABD.A(2,3);
  A66 = FailABD.A(3,3);

  B11 = 0;
  B12 = 0; 
  B16 = 0; 
  B22 = 0; 
  B26 = 0;
  B66 = 0; 
  
  D11 = FailABD.D(1,1);
  D12 = FailABD.D(1,2);
  D16 = FailABD.D(1,3);
  D22 = FailABD.D(2,2);
  D26 = FailABD.D(2,3);
  D66 = FailABD.D(3,3);
  
% Wave parameters
  alf=0.35;
  tp=0.0018;

  tt=0.01;
  dt=0.000002;
  
  Pi = pi; 
            
  RefStresses = FailABD.Q;
  
%   time difference
  NM=tt/dt+1;
  NM=int32(NM); 

  U11 = 1;
  V11 = 1;
  W11 = 1;
   
%   Laminate physical dimensions

  N = FailABD.numPlies;
  a = t1(1,5);
  b = t1(1,6);
  h = t1(1,7);
  rho = t1(1,8);

  M = rho*h;
  z = FailABD.z(1,N+1);
  x = a/2;
  y = b/2;

%   N = FailABD.numPlies;
%   a = 0.22;
%   b = 0.22;
%   rho = 1800;
%   h = 1.96e-3;
%   M = rho*h;
%   z = FailABD.z(1,N+1);
%   x = a/2;
%   y = b/2;

  m = FailABD.m(1,1);
  n = FailABD.n(1,1);
  trans = [m^2 n^2 m*n;
           n^2 m^2 -m*n;
          -2*m*n 2*m*n m^2 - n^2];
      
%   Xt = 2207e6
%   Xc = -1531e6
%   Yt = 80.7e6
%   Yc = -199.8e6
%   S = 114.5e6
%   
%   ext = 4e-3 
%   exc = -3e-3
%   eyt = 4e-3
%   eyc = -3e-3
%   es = 6e-3
  
      
% Coefficients in the time-dependent nonlinear differential equations
  a0=(a*b^9*M)/1260;
  b0=(a^9*b*M)/1260;
  c0=(a*b*M)/4;

  a1=((3*a^2*A66*b^7 + A11*b^9*Pi^2)/(315*a))/a0;
  a2=((9*a^4*(A12 + A66)*b^4)/Pi^6)/a0;
  a3=(16*b^3*(b^2*B11 + a^2*(B12 + 2*B66))*(-12 + Pi^2))/(3*a^2*Pi^3)/a0;
  a4=((b^3*(A11*b^2*(45 + Pi^4) + a^2*(90*A66 - A12*(-45 + Pi^4))))/(240*a^2*Pi))/a0;
  a5=0/a0;

  b1=(((3*a^7*A66*b^2 + a^9*A22*Pi^2))/(315*b))/b0;
  b2=((9*a^4*(A12 + A66)*b^4)/Pi^6)/b0;
  b3=(16*a^3*(a^2*B22 + b^2*(B12 + 2*B66))*(-12 + Pi^2))/(3*b^2*Pi^3)/b0;
  b4=((a^3*(a^2*A22*(45 + Pi^4) + b^2*(90*A66 - A12*(-45 + Pi^4))))/(240*b^2*Pi))/b0;
  b5=0/b0;

  c1=(((b^4*D11 + a^4*D22 + 2*a^2*b^2*(D12 + 2*D66))*Pi^4)/(4*a^3*b^3))/c0;
  c2=((32*(B12 - B66)*Pi^2)/(9*a*b) +(8*(-B12 + B66)*Pi^2)/(9*a*b))/c0;
  c3=((9*a*A22*Pi^4)/(128*b^3) - (A66*Pi^4)/(32*a*b)+(3*(A12+2*A66)*Pi^4)/(64*a*b)+(9*A11*b*Pi^4)/(128*a^3))/c0;
  c4=((16*b^3*(4*b^2*B11 + a^2*(B12+2*B66))*(-12+Pi^2))/(3*a^2*Pi^3))/c0;
  c5=((16*a^3*(4*a^2*B22 + b^2*(B12+2*B66))*(-12+Pi^2))/(3*b^2*Pi^3))/c0;
  c6=((3*(A12+A66)*b^3)/(4*Pi)+(A11*(b^5/60+(3*b^5)/(4*Pi^4))*Pi^3)/(2*a^2)-(A12*(b^5/60+(3*b^5)/(4*Pi^4))*Pi^3)/(2*b^2))/c0;
  c7=((3*a^3*(A12+A66))/(4*Pi)-(A12*(a^5/60+(3*a^5)/(4*Pi^4))*Pi^3)/(2*a^2)+(A22*(a^5/60+(3*a^5)/(4*Pi^4))*Pi^3)/(2*b^2))/c0;
  c8=((-4*a*b)/Pi^2)/c0;
  
%   IC's
  un=0;
  vn=0;
  wn=0;
  
  u = zeros(1,(NM-1));
  v = zeros(1,(NM-1));
  w = zeros(1,(NM-1));
  timesteps = zeros(1,(NM-1));

  dt1=1/dt;
  dt2=1/dt^2;

  t=0;
  
%   Iteration of w, v and u for each time step
  for i = 1:(NM)-1

  u(1)=0;
  v(1)=0;
  w(1)=0;    
      
  t= t+dt;
  timesteps(i) = t;

  www=w(i);
  
  eps = 1;
  
%   Coefficients in the finite difference equations
    while (eps > 1.0e-12)
  AA1=a1+dt2+a0*dt1;
  AA2=a2;
  AA3=a3+a4*www;
  AA4=dt1*un+(dt2+a0*dt1)*u(i)-a5;

  BB1=b2;
  BB2=dt2+b1+b0*dt1;
  BB3=b3+b4*www;
  BB4=dt1*vn+(dt2+b0*dt1)*v(i)-b5;

  CC1=c4+c6*www;
  CC2=c5+c7*www;
  CC3=dt2+c1+c2*www+c3*www^2+c0*dt1;
  CC4=dt1*wn+(dt2+c0*dt1)*w(i)+c8*pt;
 
  w(i+1) = (AA1*BB2*CC4 - AA1*BB4*CC2 - AA2*BB1*CC4 + AA2*BB4*CC1 + AA4*BB1*CC2 - AA4*BB2*CC1)/(AA1*BB2*CC3 - AA1*BB3*CC2 - AA2*BB1*CC3 + AA2*BB3*CC1 + AA3*BB1*CC2 - AA3*BB2*CC1);
  v(i+1) = -(AA1*BB3*CC4 - AA1*BB4*CC3 - AA3*BB1*CC4 + AA3*BB4*CC1 + AA4*BB1*CC3 - AA4*BB3*CC1)/(AA1*BB2*CC3 - AA1*BB3*CC2 - AA2*BB1*CC3 + AA2*BB3*CC1 + AA3*BB1*CC2 - AA3*BB2*CC1);
  u(i+1) = (AA2*BB3*CC4 - AA2*BB4*CC3 - AA3*BB2*CC4 + AA3*BB4*CC2 + AA4*BB2*CC3 - AA4*BB3*CC2)/(AA1*BB2*CC3 - AA1*BB3*CC2 - AA2*BB1*CC3 + AA2*BB3*CC1 + AA3*BB1*CC2 - AA3*BB2*CC1);
    
  eps=abs(w(i+1)-www);

  www=w(i+1);

    end

  un=(u(i+1)-u(i))/dt;
  vn=(v(i+1)-v(i))/dt;
  wn=(w(i+1)-w(i))/dt;    
    
   
  end

   epsx(1, :) = ((w(1,:)).^2*pi()^2*cos((pi()*x)/a)^2*sin((pi()*y)/b)^2)/(2*a^2) + ((w(1,:))*z*(pi())^2*sin((pi()*x)/a)*sin((pi()*y)/b))/a^2 + (2*(u(1,:))*y^2*pi()*cos((2*pi()*x)/a)*(b - y)^2)/a;
   epsy(1, :)  = ((w(1,:)).^2*pi^2*cos((pi*y)/b)^2*sin((pi*x)/a)^2)/(2*b^2) + ((w(1,:))*z*pi^2*sin((pi*x)/a)*sin((pi*y)/b))/b^2 + (2*(v(1,:))*x^2*pi*cos((2*pi*y)/b)*(a - x)^2)/b;
   epsxy(1, :)  = 2*(u(1,:))*y*sin((2*pi*x)/a)*(b - y)^2 + 2*(v(1,:))*x*sin((2*pi*y)/b)*(a - x)^2 - (u(1,:))*y^2*sin((2*pi*x)/a)*(2*b - 2*y) - (v(1,:))*x^2*sin((2*pi*y)/b)*(2*a - 2*x) - (2*(w(1,:))*z*pi^2*cos((pi*x)/a)*cos((pi*y)/b))/(a*b) + ((w(1,:)).^2*pi^2*cos((pi*x)/a)*cos((pi*y)/b)*sin((pi*x)/a)*sin((pi*y)/b))/(a*b);

   RefStrains = [epsx;
                 epsy;
                 epsxy];
            
% Conversion of Reference axes Strains to Material axes Strains
    MatStrains(:,:) = trans * RefStrains;
    
% Calc of Material Principal Stresses
    MatStresses(:,:) = RefStresses * MatStrains(:,:);  

% Hoffman Criterion

F1 = 1/Xt - 1/(abs(Xc));
F2 = 1/Yt - 1/(abs(Yc));
F11 = 1/(Xt*(abs(Xc)));
F22 = 1/(Yt*(abs(Yc)));
F33 = 1/S^2;
F12 = -1/(2*Xt*(abs(Xc)));

Hoffman = (MatStresses(1,:) * F1)  + (MatStresses(2,:) * F2) + (MatStresses(1,:).^2 * F11)  + (MatStresses(2,:).^2 * F22) + (MatStresses(3,:).^2 * F33) + (MatStresses(1,:).* MatStresses(2,:) * 2 * F12);

Hoffman = Hoffman';

% Tsai-Hill Criterion

TsaiHill = (MatStresses(1,:) / abs(Xt)).^2 + (MatStresses(2,:) / abs(Yt)).^2  + (MatStresses(3,:) / S).^2  - ((MatStresses(1,:) / abs(Xt)).* (MatStresses(2,:) / abs(Xt)));

TsaiHill = TsaiHill';

% Tsai-Wu Criterion

TsaiWu = ((MatStresses(1,:).^2 ) / (Xt*abs(Xc))) - ((MatStresses(1,:).*MatStresses(2,:)) / (sqrt(Xt*Xc*Yt*Yc))) + ((MatStresses(2,:).^2 ) / (Yt*abs(Yc))) + ((1/Xt + 1/Xc)*(MatStresses(1,:))) + ((MatStresses(3,:).^2) / S);

TsaiWu = TsaiWu';

    
  out = [Hoffman, TsaiHill, TsaiWu];

%   out = TsaiHill;
    
%   out = TsaiWu;    
end