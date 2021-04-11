function ABDSanQ44 = GUIABD_and_Strain_Sandwich_Q44_Q55(GUIInput)

t2s = GUIInput.t2s;

% Laminate Engineering Material Constants

E1c = t2s(1,1);
E2c = t2s(1,2);
G12c = t2s(1,3);
v12c = t2s(1,4);
v21c = (E2c*v12c)/E1c;

tc = t2s(1,7);

% E1c = 29.6E6;
% E2c = 14.5E6;
% v12c = 0.3;
% v21c = (E2c*v12c)/E1c;
% G12c = 14E6;
% tc = 0.00508;

fx = 0;
fy = 0;
fxy = 0;
theta1 = 0;



numPlies = 1;
theta = theta1;
z = [-tc/2,tc/2];
m = cosd(theta);  
n = sind(theta);


% Laminate Force and Moment Instensities
Nx = 0;
Ny = 0;
Nxy = 0;
Mx = 0;
My = 0;
Mxy = 0;

% Q Values for Sandwich plate assuming core is isentropic material
Q44 = G12c;
Q55 = Q44;

% Q Matrix
Q =  [0   0   0   0   0   0;
      0   0   0   0   0   0;
      0   0   0   0   0   0;
      0   0   0   Q44 0   0;
      0   0   0   0   Q55 0;
      0   0   0   0   0   0];

% Q =  [Q44 0
%       0   Q55];

% Q Tranformation Matrix, using 3 x 3 x i matrix, with i being 1 to
% numPlies
Qb = zeros(6,6,numPlies);
Aij = zeros(6,6,numPlies);

for i = 1:numPlies
     
    Qt1 =  [0   0   0   0   0   0;
            0   0   0   0   0   0;
            0   0   0   0   0   0;
            0   0   0   m(i) -n(i) 0;
            0   0   0   n(i)  m(i) 0;
            0   0   0   0   0   0];
        
    Qt2 =  [0   0   0   0   0   0;
            0   0   0   0   0   0;
            0   0   0   0   0   0;
            0   0   0   m(i) -n(i) 0;
            0   0   0   n(i)  m(i) 0;
            0   0   0   0   0   0];
%     Qt1 = [m(i) -n(i);
%            n(i)  m(i)];
%     Qt2 = [m(i) -n(i);
%            n(i)  m(i)];
    
% Qbar Matrix for each Ply
Qb(:,:,i) = Qt1*Q*Qt2;

% ABD Matricies for each Ply
Aij(:,:,i) = Qb(:,:,i)*(z(i+1)-z(i));

end

% Qbar Matricies sum for all plies
AQbarSum = sum(Qb,3);

% A Matrices sum for all plies
Aijsum = sum(Aij,3);

% A Matrix 
ABDMatrix = Aijsum;

ABDSanQ44.Aij = Aijsum;
ABDSanQ44.Q = Q;
end
