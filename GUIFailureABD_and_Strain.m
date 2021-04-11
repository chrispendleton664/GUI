function FailABD = GUIFailureABD_and_Strain(theta, GUIInput)

% Laminate Engineering Material Constants
t1 = GUIInput.t1LF;

FailABD.numPlies = GUIInput.NumPlies;

E1 = t1(1,1);
E2 = t1(1,2);
G12 = t1(1,3);
v12 = t1(1,4);
v21 = (E2*v12)/E1;

% E1 = 24.14E9;
% E2 = 24.14E9;
% v12 = 0.11;
% v21 = (E2*v12)/E1;
% G12 = 3.79E9;
% fx = 0;
% fy = 0;
% fxy = 0;
% 
% FailABD.numPlies = 7;

theta = zeros(1, FailABD.numPlies) + theta;

FailABD.z = zeros(1,FailABD.numPlies+1);
h = t1(1,7);

for i = 0:1:FailABD.numPlies
    FailABD.z(1,i+1) = -h/2 + ((h/FailABD.numPlies)*i);
end

% FailABD.z = [-0.98e-3, -0.7e-3, -0.42e-3, -0.14e-3, 0.14e-3, 0.42e-3, 0.7e-3, 0.98E-3];

FailABD.m = cosd(theta);  
FailABD.n = sind(theta);

% Q Values
Q11 = E1/(1-(v12*v21));
Q12 = (v21*E1)/(1-(v12*v21));
Q21 = (v12*E2)/((1-v12*v21));
Q22 = E2/(1-(v12*v21));
Q66 = G12;

% Q Matrix
FailABD.Q =  [Q11 Q12 0;
          Q12 Q22 0;
          0   0   Q66];

% Q Tranformation Matrix, using 3 x 3 x i matrix, with i being 1 to
% numPlies
Qb = zeros(3,3,FailABD.numPlies);
Aij = zeros(3,3,FailABD.numPlies);
Bij = zeros(3,3,FailABD.numPlies);
Dij = zeros(3,3,FailABD.numPlies);

for i = 1:FailABD.numPlies
    Qt1 = [FailABD.m(i)^2 FailABD.n(i)^2 -2*FailABD.m(i)*FailABD.n(i);
          FailABD.n(i)^2 FailABD.m(i)^2 2*FailABD.m(i)*FailABD.n(i);
          FailABD.m(i)*FailABD.n(i) -FailABD.m(i)*FailABD.n(i) FailABD.m(i)^2 - FailABD.n(i)^2];
    Qt2 = [FailABD.m(i)^2 FailABD.n(i)^2 FailABD.m(i)*FailABD.n(i);
          FailABD.n(i)^2 FailABD.m(i)^2 -FailABD.m(i)*FailABD.n(i);
          -2*FailABD.m(i)*FailABD.n(i) 2*FailABD.m(i)*FailABD.n(i) FailABD.m(i)^2 - FailABD.n(i)^2];
    
% Qbar Matrix for each Ply
Qb(:,:,i) = Qt1*FailABD.Q*Qt2;

% ABD Matricies for each Ply
Aij(:,:,i) = Qb(:,:,i)*(FailABD.z(i+1)-FailABD.z(i));
Bij(:,:,i) = Qb(:,:,i)*-0.5*((FailABD.z(i+1))^2 - (FailABD.z(i))^2);
Dij(:,:,i) = Qb(:,:,i)*(1/3)*((FailABD.z(i+1))^3 - (FailABD.z(i))^3);
end

% Qbar Matricies sum for all plies
AQbarSum = sum(Qb,3);

% ABD Matricies sum for all plies
Aijsum = sum(Aij,3);
Bijsum = sum(Bij,3);
Dijsum = sum(Dij,3);
% ABD Matrix 
ABDMatrix = [Aijsum Bijsum;
             Bijsum Dijsum];

FailABD.A = Aijsum;
FailABD.B = Bijsum;
FailABD.D = Dijsum;
end
