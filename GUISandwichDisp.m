function out = GUISandwichDisp(GUIInput)
    
t1s = GUIInput.t1s;
t2s = GUIInput.t2s;

Xt = GUIInput.t3s(1,1);
Xc = GUIInput.t3s(1,2);
Yt = GUIInput.t3s(1,3);
Yc = GUIInput.t3s(1,4);
S = GUIInput.t3s(1,5);
  
ext = GUIInput.t3s(1,6); 
exc = GUIInput.t3s(1,7);
eyt = GUIInput.t3s(1,8);
eyc = GUIInput.t3s(1,9);
es = GUIInput.t3s(1,10);

% ABD Calc input
San = GUIABD_and_Strain_Sandwich(GUIInput);

  A11 = San.Aij(1,1);
  A12 = San.Aij(1,2);
  A16 = San.Aij(1,6);
  A22 = San.Aij(2,2);
  A26 = San.Aij(2,6);
  A44 = San.Aij(4,4);
  A55 = San.Aij(5,5);
  A66 = San.Aij(6,6);

  B11 = 0;
  B12 = 0; 
  B16 = 0; 
  B22 = 0; 
  B26 = 0;
  B66 = 0; 
  
  D11 = San.Dij(1,1);
  D12 = San.Dij(1,2);
  D16 = San.Dij(1,6);
  D22 = San.Dij(2,2);
  D26 = San.Dij(2,6);
  D66 = San.Dij(6,6);

  K44=0.833;
  K55=0.833;
  
  
% Wave parameters
  alf=0.35;
  tp=0.0018;

  tt=0.025;
  dt=0.0000021;
  
  Pi = pi; 

  RefStresses = San.Q;
  
% time difference
  NM=tt/dt+1;
  NM=int32(NM); 

  U11 = 1;
  V11 = 1;
  WB1 = 1;
  WS1 = 1;
  
  
% Laminate physical dimensions
  a = t1s(1,5);
  b = t1s(1,6);
  tf = t1s(1,7);
  tc = t2s(1,7);
  rhof = t1s(1,8);
  rhoc = t2s(1,8);
  
  M = (rhof*tf*2)+(rhoc*tc);
  z = San.tf1;
  x = a/2;
  y = b/2;
  m = San.m;
  n = San.n;
  trans = [m^2 n^2 m*n 0 0 0;
           n^2 m^2 -m*n 0 0 0;
          -2*m*n 2*m*n m^2 - n^2 0 0 0;
           0 0 0 0 0 0;
           0 0 0 0 0 0;
           0 0 0 0 0 0];
% Coefficients in the time-dependent nonlinear differential equations      
      
%   a1=((3*a^2*A66*b^7 + A11*b^9*Pi^2)/(315*a))/a0;
%   a2=((9*a^4*(A12 + A66)*b^4)/Pi^6)/a0;
%   a3=((16*b^3*(b^2*B11 + a^2*(B12 + 2*B66))*(-12 + Pi^2))/(3*a^2*Pi^3))/a0;
%   a4=0;
%   a5=((b^3*(A11*b^2*(45 + Pi^4) + a^2*(90*A66 - A12*(-45 + Pi^4))))/(240*a^2*Pi))/a0;
%   a6=((b^3*(A11*b^2*(45 + Pi^4) + a^2*(90*A66 - A12*(-45 + Pi^4))))/(240*a^2*Pi))/a0;
%   a7=((b^3*(A11*b^2*(45 + Pi^4) + a^2*(90*A66 - A12*(-45 + Pi^4))))/(120*a^2*Pi))/a0;
%   a8=0;  
%   
%   b1=((9*a^4*(A12 + A66)*b^4)/Pi^6)/b0;
%   b2=((3*a^2*A66*b^7 + A11*b^9*Pi^2)/(315*a))/a0;
%   b3=((16*a^3*(a^2*B22 + b^2*(B12 + 2*B66))*(-12 + Pi^2))/(3*b^2*Pi^3))/b0;
%   b4=0;
%   b5=((a^3*(a^2*A22*(45 + Pi^4) + b^2*(90*A66 - A12*(-45 + Pi^4))))/(240*b^2*Pi))/b0;
%   b6=((a^3*(a^2*A22*(45 + Pi^4) + b^2*(90*A66 - A12*(-45 + Pi^4))))/(240*b^2*Pi))/b0;
%   b7=((a^3*(a^2*A22*(45 + Pi^4) + b^2*(90*A66 - A12*(-45 + Pi^4))))/(120*b^2*Pi))/b0;
%   b8=0;
% 
%   c1= ((16*b^3*(4*b^2*B11 + a^2*(B12+2*B66))*(-12+Pi^2))/(3*a^2*Pi^3))/c0;
%   c2= ((16*a^3*(4*a^2*B22 + b^2*(B12+2*B66))*(-12+Pi^2))/(3*b^2*Pi^3))/c0;
%   c3= (((b^4*D11 + a^4*D22 + 2*a^2*b^2*(D12 + 2*D66))*Pi^4)/(4*a^3*b^3))/c0;
%   c4= ((a*D22*pi^4)/(4*b^3))/c0;
%   c5= ((8*pi^2*(B12 - B66))/ (3*a*b))/c0;
%   c6= ((-8*(b^4*B11 + a^4*B22 + a^2*b^2*(-B12 + B66))*pi^2)/(9*a^3*b^3))/c0;
%   c7= ((-8*(b^4*B11 + a^4*B22 + 4*a^2*b^2*(-B12 + B66))*pi^2)/(9*a^3*b^3))/c0;
%   c8= ((b^3*(A11*b^2*(45 + Pi^4) + a^2*(90*A66 - A12*(-45 + Pi^4))))/(120*a^2*Pi))/c0;
%   c9= ((b^3*(A11*b^2*(45 + Pi^4) + a^2*(90*A66 - A12*(-45 + Pi^4))))/(120*a^2*Pi))/c0;
%   c10= ((a^3*(2*a^2*A22*(45 + Pi^4) + b^2*(90*A66 - A12*(-45 + Pi^4))))/(120*b^2*Pi))/c0;
%   c11= ((a^3*(2*a^2*A22*(45 + Pi^4) + b^2*(90*A66 - A12*(-45 + Pi^4))))/(120*b^2*Pi))/c0;
%   c12= (((9*a^4*A22 + 2*a^2*(3*A12 + 4*A66)*b^2 + 9*A11*b^4)*Pi^4)/(128*a^3*b^3))/c0; 
%   c13= (((9*a^4*A22 + 2*a^2*(3*A12 + 4*A66)*b^2 + 9*A11*b^4)*Pi^4)/(128*a^3*b^3))/c0;
%   c14= (3*((9*a^4*A22 + 2*a^2*(3*A12 + 4*A66)*b^2 + 9*A11*b^4)*Pi^4)/(128*a^3*b^3))/c0;
%   c15= (3*((9*a^4*A22 + 2*a^2*(3*A12 + 4*A66)*b^2 + 9*A11*b^4)*Pi^4)/(128*a^3*b^3))/c0;
% %   c16= ((-4*a*b)/Pi^2)/c0;
%   
%   d1= 0;
%   d2= 0;
%   d3= 0;
%   d4= ((-(a^2*A44*K44 + b^2*A55*K55)*pi^2)/4*a*b)/d0;
%   d5= ((-16*(b^4*B11 + a^4*B22 + a^2*b^2*(2*B12 + B66))*pi^2)/(9*a^3*b^3))/d0;
%   d6= 0; 
%   d7= ((-16*(b^4*B11 + a^4*B22 + a^2*b^2*(2*B12 + B66))*pi^2)/(9*a^3*b^3))/d0; 
%   d8= ((b^3*(A11*b^2*(45 + Pi^4) + a^2*(90*A66 + A12*(45 + Pi^4))))/(120*a^2*Pi))/d0;
%   d9= ((b^3*(A11*b^2*(45 + Pi^4) + a^2*(90*A66 + A12*(45 + Pi^4))))/(120*a^2*Pi))/d0;
%   d10= ((a^3*(2*a^2*A22*(45 + Pi^4) + b^2*(90*A66 + A12*(45 + Pi^4))))/(120*b^2*Pi))/d0;
%   d11= ((a^3*(2*a^2*A22*(45 + Pi^4) + b^2*(90*A66 + A12*(45 + Pi^4))))/(120*b^2*Pi))/d0;
%   d12= (((-3*a^4*A22 + 2*a^2*(3*A12 - 2*A66)*b^2 + 3*A11*b^4)*Pi^4)/(128*a^3*b^3))/d0;
%   d13= (((-3*a^4*A22 + 2*a^2*(3*A12 - 2*A66)*b^2 + 3*A11*b^4)*Pi^4)/(128*a^3*b^3))/d0;
%   d14= ((3*(-3*a^4*A22 + 2*a^2*(3*A12 - 2*A66)*b^2 + 3*A11*b^4)*Pi^4)/(128*a^3*b^3))/d0;
%   d15= ((3*(-3*a^4*A22 + 2*a^2*(3*A12 - 2*A66)*b^2 + 3*A11*b^4)*Pi^4)/(128*a^3*b^3))/d0;
%   d16= ((-4*a*b)/Pi^2)/d0;  
  
  a0=(9*a* b^9 *M)/1260;
  b0=(9* a^9 *b*M)/1260;
  c0=(a*b*M*WB1^2)/4;
  d0=(a*b*M*WS1^2)/4;

  a1=((3*a^2*A66*b^7 + A11*b^9*Pi^2)*U11^2)/(315*a);
  a2=(9*a^4*(A12 + A66)*b^4*U11*V11)/Pi^6;
  a3=(16*b^3*(b^2*B11 + a^2*(B12 + 2*B66))*(-12 + Pi^2)*U11*WB1)/(3*a^2*Pi^3);
  a4=0;
  a5=(b^3*(A11*b^2*(45 + Pi^4) + a^2*(90*A66 - A12*(-45 + Pi^4)))*U11*WB1^2)/(240*a^2*Pi);
  a6=(b^3*(A11*b^2*(45 + Pi^4) + a^2*(90*A66 - A12*(-45 + Pi^4)))*U11*WS1^2)/(240*a^2*Pi);
  a7=(b^3*(A11*b^2*(45 + Pi^4) + a^2*(90*A66 - A12*(-45 + Pi^4)))*U11*WB1*WS1)/(120*a^2*Pi);
  a8=0;

  b1=(9*a^4*(A12 + A66)*b^4*U11*V11)/Pi^6;
  b2=((3*a^7*A66*b^2 + a^9*A22*Pi^2)*V11^2)/(315*b);
  b3=(16*a^3*(a^2*B22 + b^2*(B12 + 2*B66))*(-12 + Pi^2)*V11*WB1)/(3*b^2*Pi^3);
  b4=0;
  b5=(a^3*(a^2*A22*(45 + Pi^4) + b^2*(90*A66 - A12*(-45 + Pi^4)))*V11*WB1^2)/(240*b^2*Pi);
  b6=(a^3*(a^2*A22*(45 + Pi^4) + b^2*(90*A66 - A12*(-45 + Pi^4)))*V11*WS1^2)/(240*b^2*Pi);
  b7=(a^3*(a^2*A22*(45 + Pi^4) + b^2*(90*A66 - A12*(-45 + Pi^4)))*V11*WB1*WS1)/(120*b^2*Pi);
  b8=0;

  c1=(16*b^3*(4*b^2*B11 + a^2*(B12 + 2*B66))*(-12 + Pi^2)*U11*WB1)/(3*a^2*Pi^3);
  c2=(16*a^3*(4*a^2*B22 + b^2*(B12 + 2*B66))*(-12 + Pi^2)*V11*WB1)/(3*b^2*Pi^3);
  c3=((b^4*D11 + a^4*D22 + 2*a^2*b^2*(D12 + 2*D66))*Pi^4*WB1^2)/(4*a^3*b^3);
  c4=(a*D22*Pi^4*WB1*WS1)/(4*b^3);
  c5=(40*B12*Pi^2*WB1^2)/(9*a*b) + (8*B66*Pi^2*WB1^2)/(9*a*b)-(16*(B12 + 2*B66)*Pi^2*WB1^2)/(9*a*b);
  c6=(-8*b*B11*Pi^2*WS1^2)/(9*a^3) + (8*B12*Pi^2*WS1^2)/(9*a*b) - (8*a*B22*Pi^2*WS1^2)/(9*b^3)-(8*B66*Pi^2*WS1^2)/(9*a*b);
  c7=(-8*b*B11*Pi^2*WB1*WS1)/(9*a^3) + (16*B12*Pi^2*WB1*WS1)/(3*a*b)-(8*a*B22*Pi^2*WB1*WS1)/(9*b^3)-(16*(B12+2*B66)*Pi^2*WB1*WS1)/(9*a*b);
  c8=(3*(A12 + A66)*b^3*U11*WB1)/(4*Pi)+(A11*(b^5/60+(3*b^5)/(4*Pi^4))*Pi^3*U11*WB1)/(2*a^2)-(A12*(b^5/60+(3*b^5)/(4*Pi^4))*Pi^3*U11*WB1)/(2*b^2);
  c9=(3*(A12 + A66)*b^3*U11*WS1)/(4*Pi) + (A11*(b^5/60+(3*b^5)/(4*Pi^4))*Pi^3*U11*WS1)/(2*a^2)-(A12*(b^5/60+(3*b^5)/(4*Pi^4))*Pi^3*U11*WS1)/(2*b^2);
  c10=(3*a^3*(A12 + A66)*V11*WB1)/(4*Pi)-(A12*(a^5/60+(3*a^5)/(4*Pi^4))*Pi^3*V11*WB1)/(2*a^2) + (A22*(a^5/60+(3*a^5)/(4*Pi^4))*Pi^3*V11*WB1)/b^2;
  c11=(3*a^3*(A12 + A66)*V11*WS1)/(4*Pi)-(A12*(a^5/60+(3*a^5)/(4*Pi^4))*Pi^3*V11*WS1)/(2*a^2) + (A22*(a^5/60+(3*a^5)/(4*Pi^4))*Pi^3*V11*WS1)/b^2;
  c12=(9*a*A22*Pi^4*WB1^3)/(128*b^3) + (3*A12*Pi^4*WB1^3)/(64*a*b) + (A66*Pi^4*WB1^3)/(16*a*b) +(9*A11*b*Pi^4*WB1^3)/(128*a^3);
  c13=(9*a*A22*Pi^4*WS1^3)/(128*b^3) + (3*A12*Pi^4*WS1^3)/(64*a*b)+(A66*Pi^4*WS1^3)/(16*a*b)+(9*A11*b*Pi^4*WS1^3)/(128*a^3);
  c14=(27*a*A22*Pi^4*WB1^2*WS1)/(128*b^3)+(9*A12*Pi^4*WB1^2*WS1)/(64*a*b)+(3*A66*Pi^4*WB1^2*WS1)/(16*a*b) + (27*A11*b*Pi^4*WB1^2*WS1)/(128*a^3);
  c15=(27*a*A22*Pi^4*WB1*WS1^2)/(128*b^3)+(9*A12*Pi^4*WB1*WS1^2)/(64*a*b) + (3*A66*Pi^4*WB1*WS1^2)/(16*a*b) + (27*A11*b*Pi^4*WB1*WS1^2)/(128*a^3);

  d1=0;
  d2=0;
  d3=0;
  d4=-((a^2*A44*K44 + A55*b^2*K55)*Pi^2*WS1^2)/(4*a*b);
  d5=(-16*b*B11*Pi^2*WB1^2)/(9*a^3)-(32*B12*Pi^2*WB1^2)/(9*a*b)-(16*a*B22*Pi^2*WB1^2)/(9*b^3)-(16*B66*Pi^2*WB1^2)/(9*a*b);
  d6=0;
  d7=(-16*b*B11*Pi^2*WB1*WS1)/(9*a^3) - (32*B12*Pi^2*WB1*WS1)/(9*a*b)-(16*a*B22*Pi^2*WB1*WS1)/(9*b^3)-(16*B66*Pi^2*WB1*WS1)/(9*a*b);
  d8=(3*A66*b^3*U11*WB1)/(4*Pi) + (A11*(b^5/60 + (3*b^5)/(4*Pi^4))*Pi^3*U11*WB1)/(2*a^2) + (A12*(b^5/60+(3*b^5)/(4*Pi^4))*Pi^3*U11*WB1)/(2*b^2);
  d9=(3*A66*b^3*U11*WS1)/(4*Pi) + (A11*(b^5/60 + (3*b^5)/(4*Pi^4))*Pi^3*U11*WS1)/(2*a^2) + (A12*(b^5/60+(3*b^5)/(4*Pi^4))*Pi^3*U11*WS1)/(2*b^2);
  d10=(3*a^3*A66*V11*WB1)/(4*Pi) + (A12*(a^5/60 + (3*a^5)/(4*Pi^4))*Pi^3*V11*WB1)/(2*a^2) + (A22*(a^5/60+(3*a^5)/(4*Pi^4))*Pi^3*V11*WB1)/(2*b^2);
  d11=(3*a^3*A66*V11*WS1)/(4*Pi) + (A12*(a^5/60 + (3*a^5)/(4*Pi^4))*Pi^3*V11*WS1)/(2*a^2) + (A22*(a^5/60+(3*a^5)/(4*Pi^4))*Pi^3*V11*WS1)/(2*b^2);
  d12=(-3*a*A22*Pi^4*WB1^3)/(128*b^3) - (3*A12*Pi^4*WB1^3)/(64*a*b)+(A66*Pi^4*WB1^3)/(32*a*b)-(3*A11*b*Pi^4*WB1^3)/(128*a^3);
  d13=(-3*a*A22*Pi^4*WS1^3)/(128*b^3) - (3*A12*Pi^4*WS1^3)/(64*a*b)+(A66*Pi^4*WS1^3)/(32*a*b)-(3*A11*b*Pi^4*WS1^3)/(128*a^3);
  d14=(-9*a*A22*Pi^4*WB1^2*WS1)/(128*b^3)-(9*A12*Pi^4*WB1^2*WS1)/(64*a*b)+(3*A66*Pi^4*WB1^2*WS1)/(32*a*b)-(9*A11*b*Pi^4*WB1^2*WS1)/(128*a^3);
  d15=(-9*a*A22*Pi^4*WB1*WS1^2)/(128*b^3)-(9*A12*Pi^4*WB1*WS1^2)/(64*a*b)+(3*A66*Pi^4*WB1*WS1^2)/(32*a*b) - (9*A11*b*Pi^4*WB1*WS1^2)/(128*a^3);

  
%   IC's
  un=0;
  vn=0;
  wbbn=0;
  wssn=0;

  aS=0;
  bS=0;
  cS=0;
  dS=0;
  
  
  u = zeros(1,(NM-1));
  v = zeros(1,(NM-1));
  wbb = zeros(1,(NM-1));
  wss = zeros(1,(NM-1));
  
  zerostrain = zeros(1,(NM));

  dt1=1/dt;
  dt2=1/dt^2;

  t=zeros((NM-1),1);
  
%   Iteration of w, v and u for each time step  
  for i = 1:(NM)-1

  u(1)=0;
  v(1)=0;
  wbb(1)=0;    
  wss(1)=0;    

  t(i+1,:)= t(i,:)+dt;  
  
%   Blast load eqn  
  pm= 29.8e3;  
  pt(i+1,:) = -pm*(1-t(i+1,:)./tp).*exp(-alf*t(i+1,:)./tp);
  
  c16=(-4*a*b*WB1*pt(i,:))/Pi^2;
  d16=(-4*a*b*WS1*pt(i,:))/Pi^2;

  wb=wbb(i);
  ws=wss(i);
  
  ep = 1;
  
%   Coefficients in the finite difference equations  
    while (ep > 1.0e-18)
  AA1=(dt2*a0+dt1*aS+a1);
  AA2=(a2);
  AA3=(a3+a5*wb);
  AA4=(a4+a6*ws+a7*wb);
  AA5=dt1*a0*un+(dt2*a0+dt1*aS)*u(i)-a8;
  
  BB1=(b1);
  BB2=(dt2*b0+dt1*bS+b2);
  BB3=(b3+b5*wb);
  BB4=(b4+b6*ws+b7*wb);
  BB5=dt1*b0*vn+(dt2*b0+dt1*bS)*v(i)-b8;
  
  CC1=(c1+c8*wb+c9*ws);
  CC2=(c2+c10*wb+c11*ws);
  CC3=(dt2*c0+dt1*cS+c3+c5*wb+c7*ws+c12*wb^2+c14*wb*ws);
  CC4=(c4+c6*ws+c13*ws^2+c15*wb*ws);
  CC5=dt1*c0*wbbn+(dt2*c0+dt1*cS)*wbb(i)-c16;
  
  DD1=(d1+d8*wb+d9*ws);
  DD2=(d2+d10*wb+d11*ws);
  DD3=(d3+d5*wb+d7*ws+d12*wb^2+d14*wb*ws);
  DD4=(dt2*d0+dt1*dS+d4+d6*ws+d13*ws^2+d15*wb*ws);
  DD5=dt1*d0*wssn+(dt2*d0+dt1*dS)*wss(i)-d16;

 
  wss(i+1) = -(AA1*BB2*CC4*DD5 - AA1*BB2*CC5*DD4 - AA1*BB4*CC2*DD5 + AA1*BB4*CC5*DD2 + AA1*BB5*CC2*DD4 - AA1*BB5*CC4*DD2 - AA2*BB1*CC4*DD5 + AA2*BB1*CC5*DD4 + AA2*BB4*CC1*DD5 - AA2*BB4*CC5*DD1 - AA2*BB5*CC1*DD4 + AA2*BB5*CC4*DD1 + AA4*BB1*CC2*DD5 - AA4*BB1*CC5*DD2 - AA4*BB2*CC1*DD5 + AA4*BB2*CC5*DD1 + AA4*BB5*CC1*DD2 - AA4*BB5*CC2*DD1 - AA5*BB1*CC2*DD4 + AA5*BB1*CC4*DD2 + AA5*BB2*CC1*DD4 - AA5*BB2*CC4*DD1 - AA5*BB4*CC1*DD2 + AA5*BB4*CC2*DD1)/(AA1*BB2*CC3*DD4 - AA1*BB2*CC4*DD3 - AA1*BB3*CC2*DD4 + AA1*BB3*CC4*DD2 + AA1*BB4*CC2*DD3 - AA1*BB4*CC3*DD2 - AA2*BB1*CC3*DD4 + AA2*BB1*CC4*DD3 + AA2*BB3*CC1*DD4 - AA2*BB3*CC4*DD1 - AA2*BB4*CC1*DD3 + AA2*BB4*CC3*DD1 + AA3*BB1*CC2*DD4 - AA3*BB1*CC4*DD2 - AA3*BB2*CC1*DD4 + AA3*BB2*CC4*DD1 + AA3*BB4*CC1*DD2 - AA3*BB4*CC2*DD1 - AA4*BB1*CC2*DD3 + AA4*BB1*CC3*DD2 + AA4*BB2*CC1*DD3 - AA4*BB2*CC3*DD1 - AA4*BB3*CC1*DD2 + AA4*BB3*CC2*DD1);
  wbb(i+1) = -(AA1*BB2*CC4*DD5 - AA1*BB2*CC5*DD4 - AA1*BB4*CC2*DD5 + AA1*BB4*CC5*DD2 + AA1*BB5*CC2*DD4 - AA1*BB5*CC4*DD2 - AA2*BB1*CC4*DD5 + AA2*BB1*CC5*DD4 + AA2*BB4*CC1*DD5 - AA2*BB4*CC5*DD1 - AA2*BB5*CC1*DD4 + AA2*BB5*CC4*DD1 + AA4*BB1*CC2*DD5 - AA4*BB1*CC5*DD2 - AA4*BB2*CC1*DD5 + AA4*BB2*CC5*DD1 + AA4*BB5*CC1*DD2 - AA4*BB5*CC2*DD1 - AA5*BB1*CC2*DD4 + AA5*BB1*CC4*DD2 + AA5*BB2*CC1*DD4 - AA5*BB2*CC4*DD1 - AA5*BB4*CC1*DD2 + AA5*BB4*CC2*DD1)/(AA1*BB2*CC3*DD4 - AA1*BB2*CC4*DD3 - AA1*BB3*CC2*DD4 + AA1*BB3*CC4*DD2 + AA1*BB4*CC2*DD3 - AA1*BB4*CC3*DD2 - AA2*BB1*CC3*DD4 + AA2*BB1*CC4*DD3 + AA2*BB3*CC1*DD4 - AA2*BB3*CC4*DD1 - AA2*BB4*CC1*DD3 + AA2*BB4*CC3*DD1 + AA3*BB1*CC2*DD4 - AA3*BB1*CC4*DD2 - AA3*BB2*CC1*DD4 + AA3*BB2*CC4*DD1 + AA3*BB4*CC1*DD2 - AA3*BB4*CC2*DD1 - AA4*BB1*CC2*DD3 + AA4*BB1*CC3*DD2 + AA4*BB2*CC1*DD3 - AA4*BB2*CC3*DD1 - AA4*BB3*CC1*DD2 + AA4*BB3*CC2*DD1);
  v(i+1) = (AA1*BB3*CC4*DD5 - AA1*BB3*CC5*DD4 - AA1*BB4*CC3*DD5 + AA1*BB4*CC5*DD3 + AA1*BB5*CC3*DD4 - AA1*BB5*CC4*DD3 - AA3*BB1*CC4*DD5 + AA3*BB1*CC5*DD4 + AA3*BB4*CC1*DD5 - AA3*BB4*CC5*DD1 - AA3*BB5*CC1*DD4 + AA3*BB5*CC4*DD1 + AA4*BB1*CC3*DD5 - AA4*BB1*CC5*DD3 - AA4*BB3*CC1*DD5 + AA4*BB3*CC5*DD1 + AA4*BB5*CC1*DD3 - AA4*BB5*CC3*DD1 - AA5*BB1*CC3*DD4 + AA5*BB1*CC4*DD3 + AA5*BB3*CC1*DD4 - AA5*BB3*CC4*DD1 - AA5*BB4*CC1*DD3 + AA5*BB4*CC3*DD1)/(AA1*BB2*CC3*DD4 - AA1*BB2*CC4*DD3 - AA1*BB3*CC2*DD4 + AA1*BB3*CC4*DD2 + AA1*BB4*CC2*DD3 - AA1*BB4*CC3*DD2 - AA2*BB1*CC3*DD4 + AA2*BB1*CC4*DD3 + AA2*BB3*CC1*DD4 - AA2*BB3*CC4*DD1 - AA2*BB4*CC1*DD3 + AA2*BB4*CC3*DD1 + AA3*BB1*CC2*DD4 - AA3*BB1*CC4*DD2 - AA3*BB2*CC1*DD4 + AA3*BB2*CC4*DD1 + AA3*BB4*CC1*DD2 - AA3*BB4*CC2*DD1 - AA4*BB1*CC2*DD3 + AA4*BB1*CC3*DD2 + AA4*BB2*CC1*DD3 - AA4*BB2*CC3*DD1 - AA4*BB3*CC1*DD2 + AA4*BB3*CC2*DD1);
  u(i+1) = -(AA2*BB3*CC4*DD5 - AA2*BB3*CC5*DD4 - AA2*BB4*CC3*DD5 + AA2*BB4*CC5*DD3 + AA2*BB5*CC3*DD4 - AA2*BB5*CC4*DD3 - AA3*BB2*CC4*DD5 + AA3*BB2*CC5*DD4 + AA3*BB4*CC2*DD5 - AA3*BB4*CC5*DD2 - AA3*BB5*CC2*DD4 + AA3*BB5*CC4*DD2 + AA4*BB2*CC3*DD5 - AA4*BB2*CC5*DD3 - AA4*BB3*CC2*DD5 + AA4*BB3*CC5*DD2 + AA4*BB5*CC2*DD3 - AA4*BB5*CC3*DD2 - AA5*BB2*CC3*DD4 + AA5*BB2*CC4*DD3 + AA5*BB3*CC2*DD4 - AA5*BB3*CC4*DD2 - AA5*BB4*CC2*DD3 + AA5*BB4*CC3*DD2)/(AA1*BB2*CC3*DD4 - AA1*BB2*CC4*DD3 - AA1*BB3*CC2*DD4 + AA1*BB3*CC4*DD2 + AA1*BB4*CC2*DD3 - AA1*BB4*CC3*DD2 - AA2*BB1*CC3*DD4 + AA2*BB1*CC4*DD3 + AA2*BB3*CC1*DD4 - AA2*BB3*CC4*DD1 - AA2*BB4*CC1*DD3 + AA2*BB4*CC3*DD1 + AA3*BB1*CC2*DD4 - AA3*BB1*CC4*DD2 - AA3*BB2*CC1*DD4 + AA3*BB2*CC4*DD1 + AA3*BB4*CC1*DD2 - AA3*BB4*CC2*DD1 - AA4*BB1*CC2*DD3 + AA4*BB1*CC3*DD2 + AA4*BB2*CC1*DD3 - AA4*BB2*CC3*DD1 - AA4*BB3*CC1*DD2 + AA4*BB3*CC2*DD1);
 
    
  ep=abs(wbb(i+1)-wb)+abs(wss(i+1)-ws);

  wb=wbb(i+1);
  ws=wss(i+1);

    end

  un=(u(i+1)-u(i))/dt;
  vn=(v(i+1)-v(i))/dt;
  wbbn=(wbb(i+1)-wbb(i))/dt;
  wssn=(wss(i+1)-wss(i))/dt;

%   epx = 31.0280755910103*(3.141592653589793*(0.00554*wbb(i+1)));
  
%   if(mod(i+1,5) == 0)
%       fileID = fopen('w(i+1).txt','a');
%       fprintf(fileID,'%.24f\n',(wbb(i+1)*1e3) + (wss(i+1)*1e3));
%       fclose(fileID);
%       type w(i+1).txt
%   end

%   if(mod(i+1,5) == 0)
%       fileID = fopen('w(i+1).txt','a');
%       fprintf(fileID,'%.24f\n',t);
%       fclose(fileID);
%       type w(i+1).txt
%   end  

  end

%  out = ((wbb(1,:)) + (wss(1,:)))*1e3;
%  out = out';
%  out = out(1:5:end);

   epsx(1, :) = ((wbb(1,:).*pi*cos((pi*x)/a)*sin((pi*y)/b))/a + (wss(1,:).*pi*cos((pi*x)/a)*sin((pi*y)/b))/a).^2/2 + (wbb(1,:).*z*pi^2*sin((pi*x)/a)*sin((pi*y)/b))/a^2 + (2*un(1,:).*y^2*pi*cos((2*pi*x)/a)*(b - y)^2)/a;   
   epsy(1, :) = ((wbb(1,:).*pi*cos((pi*y)/b)*sin((pi*x)/a))/b + (wss(1,:).*pi*cos((pi*y)/b)*sin((pi*x)/a))/b).^2/2 + (wbb(1,:).*z*pi^2*sin((pi*x)/a)*sin((pi*y)/b))/b^2 + (2*vn(1,:).*x^2*pi*cos((2*pi*y)/b)*(a - x)^2)/b;
   epsxy(1, :) = ((wbb(1,:).*pi*cos((pi*x)/a)*sin((pi*y)/b))/a + (wss(1,:).*pi*cos((pi*x)/a)*sin((pi*y)/b))/a).*((wbb(1,:).*pi*cos((pi*y)/b)*sin((pi*x)/a))/b + (wss(1,:).*pi*cos((pi*y)/b)*sin((pi*x)/a))/b) + 2*un(1,:).*y*sin((2*pi*x)/a)*(b - y)^2 + 2*vn(1,:).*x*sin((2*pi*y)/b)*(a - x)^2 - un(1,:).*y^2*sin((2*pi*x)/a)*(2*b - 2*y) - vn(1,:).*x^2*sin((2*pi*y)/b)*(2*a - 2*x) - (2*wbb(1,:).*z*pi^2*cos((pi*x)/a)*cos((pi*y)/b))/(a*b);
   epsxz(1, :) = (wss(1,:).*pi*cos((pi*x)/a)*sin((pi*y)/b))/a;
   epsyz(1, :) = (wss(1,:).*pi*cos((pi*y)/b)*sin((pi*x)/a))/b;

   
   RefStrains = [epsx;
                 epsy;
                 epsxy;
                 zerostrain;
                 epsxz;
                 epsyz];
            
% Conversion of Reference axes Strains to Material axes Strains
    MatStrains(:,:) = trans * RefStrains;
    MatStrains1 = MatStrains(1,:)';
    MatStrains1 = MatStrains1(1:5:end);
    MatStrains2 = MatStrains(2,:)';
    MatStrains2 = MatStrains2(1:5:end);
    MatStrains3 = MatStrains(3,:)';
    MatStrains3 = MatStrains3(1:5:end);
    
% Calc of Material Principal Stresses
    MatStresses(:,:) = RefStresses * MatStrains(:,:);
    MatStresses1 = MatStresses(1,:)';
    MatStresses1 = MatStresses1(1:5:end);
    MatStresses2 = MatStresses(2,:)';
    MatStresses2 = MatStresses2(1:5:end);
    MatStresses3 = MatStresses(3,:)';
    MatStresses3 = MatStresses3(1:5:end);
    
    
% Max Stress & StrainCriterion
    MaxStressF1Xt = MatStresses(1,:)/Xt;
    MaxStressF1Xt = MaxStressF1Xt(1:5:end);
    
    MaxStressF1Xc = MatStresses(1,:)/Xc;
    MaxStressF1Xc = MaxStressF1Xc(1:5:end);
    
    MaxStressF2Yt = MatStresses(2,:)/Yt;
    MaxStressF2Yt = MaxStressF2Yt(1:5:end);
    
    MaxStressF2Yc = MatStresses(2,:)/Yc;
    MaxStressF2Yc = MaxStressF2Yc(1:5:end);
    
    MaxStressF12S = MatStresses(3,:)/S;
    MaxStressF12S = MaxStressF12S(1:5:end);
   
    
    MaxStrainsF1Xt = MatStrains(1,:)/ext;
    MaxStrainsF1Xt = MaxStrainsF1Xt(1:5:end);
    
    MaxStrainsF1Xc = MatStrains(1,:)/exc;
    MaxStrainsF1Xc = MaxStrainsF1Xc(1:5:end);
    
    MaxStrainsF2Yt = MatStrains(2,:)/eyt;
    MaxStrainsF2Yt = MaxStrainsF2Yt(1:5:end);
    
    MaxStrainsF2Yc = MatStrains(2,:)/eyc;
    MaxStrainsF2Yc = MaxStrainsF2Yc(1:5:end);
    
    MaxStrainsF12S = MatStrains(3,:)/es;
    MaxStrainsF12S = MaxStrainsF12S(1:5:end);
    
    MaxStressF1Xt = MaxStressF1Xt';
    MaxStressF1Xc = MaxStressF1Xc';
    MaxStressF2Yt = MaxStressF2Yt';
    MaxStressF2Yc = MaxStressF2Yc';
    MaxStressF12S = MaxStressF12S';
    
    MaxStrainsF1Xt = MaxStrainsF1Xt';
    MaxStrainsF1Xc = MaxStrainsF1Xc';
    MaxStrainsF2Yt = MaxStrainsF2Yt';
    MaxStrainsF2Yc = MaxStrainsF2Yc';
    MaxStrainsF12S = MaxStrainsF12S';
    
    % Hoffman Criterion

F1 = 1/Xt - 1/(abs(Xc));
F2 = 1/Yt - 1/(abs(Yc));
F11 = 1/(Xt*(abs(Xc)));
F22 = 1/(Yt*(abs(Yc)));
F33 = 1/S^2;
F12 = -1/(2*Xt*(abs(Xc)));

Hoffman = (MatStresses(1,:) * F1)  + (MatStresses(2,:) * F2) + (MatStresses(1,:).^2 * F11)  + (MatStresses(2,:).^2 * F22) + (MatStresses(3,:).^2 * F33) + (MatStresses(1,:).* MatStresses(2,:) * 2 * F12);

Hoffman = Hoffman';

Hoffman = Hoffman(1:5:end);

% Tsai-Hill Criterion

TsaiHill = (MatStresses(1,:) / abs(Xt)).^2 + (MatStresses(2,:) / abs(Yt)).^2  + (MatStresses(3,:) / S).^2  - ((MatStresses(1,:) / abs(Xt)).* (MatStresses(2,:) / abs(Xt)));

TsaiHill = TsaiHill';

TsaiHill = TsaiHill(1:5:end);

% Tsai-Wu Criterion

TsaiWu = ((MatStresses(1,:).^2 )./ (Xt*abs(Xc))) - ((MatStresses(1,:).*MatStresses(2,:))./ (sqrt(Xt*Xc*Yt*Yc))) + ((MatStresses(2,:).^2 )./ (Yt*abs(Yc))) + ((MatStresses(1,:)).*(1/Xt + 1/Xc)) + ((MatStresses(2,:)).*(1/Yt + 1/Yc)) + ((MatStresses(3,:).^2)./ S^2);

TsaiWu = TsaiWu';

TsaiWu = TsaiWu(1:5:end);

pt = pt(1:5:end);

t = t.*1e3;
t = t(1:5:end);

w = (((wbb(1,:)) + (wss(1,:))).*1e3)';
w = w(1:5:end);

    StressAndStrainsArray = [t pt w MaxStressF1Xt MaxStressF1Xc MaxStressF2Yt MaxStressF2Yc MaxStressF12S MaxStrainsF1Xt MaxStrainsF1Xc MaxStrainsF2Yt MaxStrainsF2Yc MaxStrainsF12S MatStresses1 MatStresses2 MatStresses3 MatStrains1 MatStrains2 MatStrains3 Hoffman TsaiHill TsaiWu];
    
    failure_criteria = {'t/ms', 'pt', 'w/mm', 'MaxStressF1Xt', 'MaxStressF1Xc', 'MaxStressF2Yt', 'MaxStressF2Yc', 'MaxStressF12S', 'MaxStrainsF1Xt', 'MaxStrainsF1Xc', 'MaxStrainsF2Yt', 'MaxStrainsF2Yc', 'MaxStrainsF12S', 'MatStresses1', 'MatStresses2', 'MatStresses3','MatStrains1',' MatStrains2','MatStrains3','Hoffman','TsaiHill','TsaiWu'};
  
    SSData = array2table(StressAndStrainsArray, 'VariableNames', failure_criteria);
    
    out = SSData;
    
    filename = 'LamApp Sandwich Sheet Laminate Results.xlsx';
    writetable(SSData,filename,'Sheet',1);

end  