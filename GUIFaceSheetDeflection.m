function out = GUIFaceSheetDeflection(GUIInput)

t1 = GUIInput.t1;
t2 = GUIInput.t2;
NumPlies = GUIInput.NumPlies;

Xt = GUIInput.t4(1,1);
Xc = GUIInput.t4(1,2);
Yt = GUIInput.t4(1,3);
Yc = GUIInput.t4(1,4);
S = GUIInput.t4(1,5);
  
ext = GUIInput.t4(1,6); 
exc = GUIInput.t4(1,7);
eyt = GUIInput.t4(1,8);
eyc = GUIInput.t4(1,9);
es = GUIInput.t4(1,10);

% ABD Calc input
ABD = GUIABD_and_Strain(GUIInput);

  A11 = ABD.A(1,1);
  A12 = ABD.A(1,2);
  A16 = ABD.A(1,3);
  A22 = ABD.A(2,2);
  A26 = ABD.A(2,3);
  A66 = ABD.A(3,3);

  B11 = 0;
  B12 = 0; 
  B16 = 0; 
  B22 = 0; 
  B26 = 0;
  B66 = 0; 
  
  D11 = ABD.D(1,1);
  D12 = ABD.D(1,2);
  D16 = ABD.D(1,3);
  D22 = ABD.D(2,2);
  D26 = ABD.D(2,3);
  D66 = ABD.D(3,3);
  
% Wave parameters
  alf=0.35;
  tp=0.0018;

  tt=0.1;
  dt=0.000002;
  
  Pi = pi; 
            
  RefStresses = ABD.Q;
  
%   time difference
  NM=tt/dt+1;
  NM=int32(NM); 

  U11 = 1;
  V11 = 1;
  W11 = 1;
   
%   Laminate physical dimensions
  N = ABD.numPlies;
  a = t2(1,5);
  b = t2(1,6);
  h = t2(1,7);
  rho = t2(1,8);

  M = rho*h;
  z = ABD.z(1,N+1);
  x = a/2;
  y = b/2;
  m = ABD.m(1,N);
  n = ABD.n(1,N);
  trans = [m^2 n^2 m*n;
           n^2 m^2 -m*n;
          -2*m*n 2*m*n m^2 - n^2];    
      
%   Xt = 700.11e6;
%   Xc = 570.37e6;
%   Yt = 69.67e6;
%   Yc = 122.12e6;
%   S = 68.89e6;
%   
%   ext = 4e-3; 
%   exc = -3e-3;
%   eyt = 4e-3;
%   eyc = -3e-3;
%   es = 6e-3;
  
      
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
  
  dt1=1/dt;
  dt2=1/dt^2;

  t=zeros((NM-1),1);
  pt=zeros((NM-1),1);
  
%   Iteration of w, v and u for each time step
  for i = 1:(NM)-1

  u(1)=0;
  v(1)=0;
  w(1)=0;    
      
  t(i+1,:)= t(i,:)+dt;
  
%   Blast load eqn  

  pm= 29.8e3;
  
  pt(i+1,:) = -pm*(1-t(i+1,:)./tp)*exp(-alf*t(i+1,:)/tp);

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
  CC4=dt1*wn+(dt2+c0*dt1)*w(i)+c8*pt(i,:);
 
  w(i+1) = (AA1*BB2*CC4 - AA1*BB4*CC2 - AA2*BB1*CC4 + AA2*BB4*CC1 + AA4*BB1*CC2 - AA4*BB2*CC1)/(AA1*BB2*CC3 - AA1*BB3*CC2 - AA2*BB1*CC3 + AA2*BB3*CC1 + AA3*BB1*CC2 - AA3*BB2*CC1);
  v(i+1) = -(AA1*BB3*CC4 - AA1*BB4*CC3 - AA3*BB1*CC4 + AA3*BB4*CC1 + AA4*BB1*CC3 - AA4*BB3*CC1)/(AA1*BB2*CC3 - AA1*BB3*CC2 - AA2*BB1*CC3 + AA2*BB3*CC1 + AA3*BB1*CC2 - AA3*BB2*CC1);
  u(i+1) = (AA2*BB3*CC4 - AA2*BB4*CC3 - AA3*BB2*CC4 + AA3*BB4*CC2 + AA4*BB2*CC3 - AA4*BB3*CC2)/(AA1*BB2*CC3 - AA1*BB3*CC2 - AA2*BB1*CC3 + AA2*BB3*CC1 + AA3*BB1*CC2 - AA3*BB2*CC1);
    
  eps=abs(w(i+1)-www);

  www=w(i+1);

    end

  un=(u(i+1)-u(i))/dt;
  vn=(v(i+1)-v(i))/dt;
  wn=(w(i+1)-w(i))/dt;    

%   if(mod(i+49,50) == 0)
%       fileID = fopen('w(i+1).txt','a');
%       fprintf(fileID,'%.24f\n',w(i+1)*1e3);
%       fclose(fileID);
%       type w(i+1).txt
%   end
   
  end
  
%   Reference Strains Calc
%    epsx(1, :) = ((w(i+1))^2*pi()^2*cos((pi()*x)/a)^2*sin((pi()*y)/b)^2)/(2*a^2) + ((w(i+1))*z*(pi())^2*sin((pi()*x)/a)*sin((pi()*y)/b))/a^2 + (2*(u(i+1))*y^2*pi()*cos((2*pi()*x)/a)*(b - y)^2)/a;
%    epsy(1, :)  = ((w(i+1))^2*pi^2*cos((pi*y)/b)^2*sin((pi*x)/a)^2)/(2*b^2) + ((w(i+1))*z*pi^2*sin((pi*x)/a)*sin((pi*y)/b))/b^2 + (2*(v(i+1))*x^2*pi*cos((2*pi*y)/b)*(a - x)^2)/b;
%    epsxy(1, :)  = 2*(u(i+1))*y*sin((2*pi*x)/a)*(b - y)^2 + 2*(v(i+1))*x*sin((2*pi*y)/b)*(a - x)^2 - (u(i+1))*y^2*sin((2*pi*x)/a)*(2*b - 2*y) - (v(i+1))*x^2*sin((2*pi*y)/b)*(2*a - 2*x) - (2*(w(i+1))*z*pi^2*cos((pi*x)/a)*cos((pi*y)/b))/(a*b) + ((w(i+1))^2*pi^2*cos((pi*x)/a)*cos((pi*y)/b)*sin((pi*x)/a)*sin((pi*y)/b))/(a*b);

   epsx(1, :) = ((w(1,:)).^2*pi^2*cos((pi*x)/a)^2*sin((pi*y)/b)^2)/(2*a^2) + ((w(1,:))*z*pi^2*sin((pi*x)/a)*sin((pi*y)/b))/a^2 + (2*(u(1,:))*y^2*pi*cos((2*pi*x)/a)*(b - y)^2)/a;
   epsy(1, :) = ((w(1,:)).^2*pi^2*cos((pi*y)/b)^2*sin((pi*x)/a)^2)/(2*b^2) + ((w(1,:))*z*pi^2*sin((pi*x)/a)*sin((pi*y)/b))/b^2 + (2*(v(1,:))*x^2*pi*cos((2*pi*y)/b)*(a - x)^2)/b;
   epsxy(1, :) = 2*(u(1,:))*y*sin((2*pi*x)/a)*(b - y)^2 + 2*(v(1,:))*x*sin((2*pi*y)/b)*(a - x)^2 - (u(1,:))*y^2*sin((2*pi*x)/a)*(2*b - 2*y) - (v(1,:))*x^2*sin((2*pi*y)/b)*(2*a - 2*x) - (2*(w(1,:))*z*pi^2*cos((pi*x)/a)*cos((pi*y)/b))/(a*b) + ((w(1,:)).^2*pi^2*cos((pi*x)/a)*cos((pi*y)/b)*sin((pi*x)/a)*sin((pi*y)/b))/(a*b);

   RefStrains = [epsx;
                 epsy;
                 epsxy];
            
% Conversion of Reference axes Strains to Material axes Strains
    MatStrains(:,:) = trans * RefStrains;
    MatStrains1 = MatStrains(1,:)';
%     MatStrains1 = MatStrains1(1:50:end);
    MatStrains2 = MatStrains(2,:)';
%     MatStrains2 = MatStrains2(1:50:end);
    MatStrains3 = MatStrains(3,:)';
%     MatStrains3 = MatStrains3(1:50:end);
    
% Calc of Material Principal Stresses
    MatStresses(:,:) = RefStresses * MatStrains(:,:);
    MatStresses1 = MatStresses(1,:)';
%     MatStresses1 = MatStresses1(1:50:end);
    MatStresses2 = MatStresses(2,:)';
%     MatStresses2 = MatStresses2(1:50:end);
    MatStresses3 = MatStresses(3,:)';
%     MatStresses3 = MatStresses3(1:50:end);
    
    
% Max Stress & StrainCriterion
    MaxStressF1Xt = MatStresses(1,:)/Xt;
%     MaxStressF1Xt = MaxStressF1Xt(1:50:end);
    
    MaxStressF1Xc = MatStresses(1,:)/Xc;
%     MaxStressF1Xc = MaxStressF1Xc(1:50:end);
    
    MaxStressF2Yt = MatStresses(2,:)/Yt;
%     MaxStressF2Yt = MaxStressF2Yt(1:50:end);
    
    MaxStressF2Yc = MatStresses(2,:)/Yc;
%     MaxStressF2Yc = MaxStressF2Yc(1:50:end);
    
    MaxStressF12S = MatStresses(3,:)/S;
%     MaxStressF12S = MaxStressF12S(1:50:end);
   
    
    MaxStrainsF1Xt = MatStrains(1,:)/ext;
%     MaxStrainsF1Xt = MaxStrainsF1Xt(1:50:end);
    
    MaxStrainsF1Xc = MatStrains(1,:)/exc;
%     MaxStrainsF1Xc = MaxStrainsF1Xc(1:50:end);
    
    MaxStrainsF2Yt = MatStrains(2,:)/eyt;
%     MaxStrainsF2Yt = MaxStrainsF2Yt(1:50:end);
    
    MaxStrainsF2Yc = MatStrains(2,:)/eyc;
%     MaxStrainsF2Yc = MaxStrainsF2Yc(1:50:end);
    
    MaxStrainsF12S = MatStrains(3,:)/es;
%     MaxStrainsF12S = MaxStrainsF12S(1:50:end);
    
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

% Hoffman = Hoffman(1:50:end);

% Tsai-Hill Criterion

TsaiHill = (MatStresses(1,:) / abs(Xt)).^2 + (MatStresses(2,:) / abs(Yt)).^2  + (MatStresses(3,:) / S).^2  - ((MatStresses(1,:) / abs(Xt)).* (MatStresses(2,:) / abs(Xt)));

TsaiHill = TsaiHill';

% TsaiHill = TsaiHill(1:50:end);

% Tsai-Wu Criterion

TsaiWu = ((MatStresses(1,:).^2 )./ (Xt*abs(Xc))) - ((MatStresses(1,:).*MatStresses(2,:))./ (sqrt(Xt*Xc*Yt*Yc))) + ((MatStresses(2,:).^2 )./ (Yt*abs(Yc))) + ((MatStresses(1,:)).*(1/Xt + 1/Xc)) + ((MatStresses(2,:)).*(1/Yt + 1/Yc)) + ((MatStresses(3,:).^2)./ S^2);

TsaiWu = TsaiWu';

% TsaiWu = TsaiWu(1:50:end);

% pt = pt(1:50:end);

t = t.*1e3;
% t = t(1:50:end);

w = (w.*1e3)';
% w = w(1:50:end);

    StressAndStrainsArray = [t pt w MaxStressF1Xt MaxStressF1Xc MaxStressF2Yt MaxStressF2Yc MaxStressF12S MaxStrainsF1Xt MaxStrainsF1Xc MaxStrainsF2Yt MaxStrainsF2Yc MaxStrainsF12S MatStresses1 MatStresses2 MatStresses3 MatStrains1 MatStrains2 MatStrains3 Hoffman TsaiHill TsaiWu];
    
    failure_criteria = {'t/ms', 'pt', 'w/mm', 'MaxStressF1Xt', 'MaxStressF1Xc', 'MaxStressF2Yt', 'MaxStressF2Yc', 'MaxStressF12S', 'MaxStrainsF1Xt', 'MaxStrainsF1Xc', 'MaxStrainsF2Yt', 'MaxStrainsF2Yc', 'MaxStrainsF12S', 'MatStresses1', 'MatStresses2', 'MatStresses3','MatStrains1',' MatStrains2','MatStrains3','Hoffman','TsaiHill','TsaiWu'};
  
    SSData = array2table(StressAndStrainsArray, 'VariableNames', failure_criteria);
    
    out = SSData;
    
    filename = 'LamApp Face Sheet Laminate Results.xlsx';
    writetable(SSData,filename,'Sheet',1);
    
    
%     failed = StressAndStrainsArray > 1;
%     failed_times_mask = any(failed,2);
%     failed_times = timesteps(failed_times_mask);
%     failure_types = cell(numel(failed_times),1);
%     failure_values = cell(numel(failed_times),1);
%     
%     SSArray_sub = StressAndStrainsArray(failed_times_mask,:);
    
%     for i =1:numel(failed_times)
%         
%         idxs = find(SSArray_sub(i,:)>1);
%         failure_types{i} = failure_criteria(idxs);
%         failure_values{i} = SSArray_sub(i,idxs);
%         
%     end
    
    
    
%     SSData = mat2dataset(StressAndStrainsArray);
%     
%     SSData.Properties.VarNames = {'MaxStressF1Xt', 'MaxStressF1Xc', 'MaxStressF2Yt', 'MaxStressF2Yc', 'MaxStressF12S', 'MaxStrainsF1Xt', 'MaxStrainsF1Xc', 'MaxStrainsF2Yt', 'MaxStrainsF2Yc', 'MaxStrainsF12S'};
%     
%     SSData.Properties.Description = 'Max Stress and Strain Criterion Indices Values';
%     
%     SSData(1:50,:)
    
end