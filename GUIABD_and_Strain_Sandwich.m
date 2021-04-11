function San = GUIABD_and_Strain_Sandwich(GUIInput)


t1s = GUIInput.t1s;
t2s = GUIInput.t2s;


theta = GUIInput.Orientation;


% Laminate Engineering Material Constants

E1 = t1s(1,1);
E2 = t1s(1,2);
G12 = t1s(1,3);
v12 = t1s(1,4);
v21 = (E2*v12)/E1;

tf = t1s(1,7);
tc = t2s(1,7);

San.tf1 = -(tf + tc/2);
San.tf2 = -(tc/2);

San.tf3 = tc/2;
San.tf4 = tc/2 + tf;

numPlies = 1;

for j = 1:2

z = [-(tf + tc/2), tc/2;    
     -(tc/2), (tc/2 + tf)];
m = cosd(theta);  
n = sind(theta);
San.m = m;
San.n = n;


% Q Values
Q11 = E1/(1-(v12*v21));
Q12 = (v21*E1)/(1-(v12*v21));
Q21 = (v12*E2)/((1-v12*v21));
Q22 = E2/(1-(v12*v21));
Q66 = G12;

% Q Matrix
Q =  [Q11   Q12   0   0   0   0;
      Q12   Q22   0   0   0   0;
      0   0   0   0   0   0;
      0   0   0   0   0   0;
      0   0   0   0   0   0;
      0   0   0   0   0   Q66];

% Q =  [Q11 Q12 0;
%       Q12 Q22 0;
%       0   0   Q66];

% Q Tranformation Matrix, using 3 x 3 x i matrix, with i being 1 to
% numPlies
Qb = zeros(6,6,numPlies);
Aij = zeros(6,6,numPlies);
Bij = zeros(6,6,numPlies);
Dij = zeros(6,6,numPlies);
for i = 1:numPlies
    
    Qt1 =  [m(i)^2       n(i)^2      0   0   0   -2*m(i)*n(i);
            n(i)^2       m(i)^2      0   0   0   2*m(i)*n(i);
            0            0           0   0   0   0;
            0            0           0   0   0   0;
            0            0           0   0   0   0;
            m(i)*n(i)   -m(i)*n(i)   0   0   0   m(i)^2 - n(i)^2];
        
    Qt2 =  [m(i)^2       n(i)^2      0   0   0   m(i)*n(i);
            n(i)^2       m(i)^2      0   0   0   -m(i)*n(i);
            0            0           0   0   0   0;
            0            0           0   0   0   0;
            0            0           0   0   0   0;
           -2*m(i)*n(i)  2*m(i)*n(i) 0   0   0   m(i)^2 - n(i)^2];    
    
% Qbar Matrix for each Ply
Qb(:,:,i) = Qt1*Q*Qt2;

% ABD Matricies for each Ply
Aij(:,:,i) = Qb(:,:,i)*(z(i+1,j)-z(i,j));
Bij(:,:,i) = Qb(:,:,i)*-0.5*((z(i+1,j))^2 - (z(i,j))^2);
Dij(:,:,i) = Qb(:,:,i)*(1/3)*((z(i+1,j))^3 - (z(i,j))^3);
end

% Qbar Matricies sum for all plies
AQbarSum(:,:,j) = sum(Qb,3);

% ABD Matricies sum for all plies
Aijsum(:,:,j) = sum(Aij,3);
Bijsum(:,:,j) = sum(Bij,3);
Dijsum(:,:,j) = sum(Dij,3);
% ABD Matrix 
ABDMatrix(:,:,j) = [Aijsum(:,:,j) Bijsum(:,:,j);
             Bijsum(:,:,j) Dijsum(:,:,j)];

ABDSandwich.A(:,:,j) = Aijsum(:,:,j);
ABDSandwich.B(:,:,j) = Bijsum(:,:,j);
ABDSandwich.D(:,:,j) = Dijsum(:,:,j);

Aijface1(:,:,j) = Aijsum(:,:,j);
Aijface2(:,:,j) = Aijsum(:,:,j);


ABDSanQ44 = GUIABD_and_Strain_Sandwich_Q44_Q55(GUIInput);
Aijcore(:,:,j) = ABDSanQ44.Aij;

Q44 = ABDSanQ44.Q; 
San.Q = Q + Q44;

San.Aij(:,:,j) = Aijface1(:,:,j) + Aijface2(:,:,j) + Aijcore(:,:,j);
San.Dij(:,:,j) = Aijface1(:,:,j)*((tc+tf)*0.5)^2 +  Aijface2(:,:,j)*((tc+tf)*0.5)^2;

end
San.Aij = San.Aij(:,:,1) + San.Aij(:,:,2); 
San.Dij = San.Dij(:,:,1) + San.Dij(:,:,2);

end
