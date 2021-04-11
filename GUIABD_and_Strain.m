function ABD = GUIABD_and_Strain(GUIInput)

% Laminate Engineering Material Constants
t1 = GUIInput.t1;
t2 = GUIInput.t2;


ABD.numPlies = GUIInput.NumPlies;

E1 = t2(1,1);
E2 = t2(1,2);
G12 = t2(1,3);
v12 = t2(1,4);
v21 = (E2*v12)/E1;

theta = zeros(1,ABD.numPlies);

for i = 1:ABD.numPlies
    theta(1,i) = t1(i,1);
end

ABD.z = zeros(1,ABD.numPlies+1);
h = t2(1,7);

for i = 0:1:ABD.numPlies
    ABD.z(1,i+1) = -h/2 + ((h/ABD.numPlies)*i);
end

ABD.m = cosd(theta);  
ABD.n = sind(theta);


% Laminate Force and Moment Instensities
Nx = 0;
Ny = 0;
Nxy = 0;
Mx = 0;
My = 0;
Mxy = 0;

% Q Values
Q11 = E1/(1-(v12*v21));
Q12 = (v21*E1)/(1-(v12*v21));
Q21 = (v12*E2)/((1-v12*v21));
Q22 = E2/(1-(v12*v21));
Q66 = G12;

% Q Matrix
ABD.Q =  [Q11 Q12 0;
          Q12 Q22 0;
          0   0   Q66];

% Q Tranformation Matrix, using 3 x 3 x i matrix, with i being 1 to
% numPlies
Qb = zeros(3,3,ABD.numPlies);
Aij = zeros(3,3,ABD.numPlies);
Bij = zeros(3,3,ABD.numPlies);
Dij = zeros(3,3,ABD.numPlies);
for i = 1:ABD.numPlies
% for i = 1
    Qt1 = [ABD.m(i)^2 ABD.n(i)^2 -2*ABD.m(i)*ABD.n(i);
          ABD.n(i)^2 ABD.m(i)^2 2*ABD.m(i)*ABD.n(i);
          ABD.m(i)*ABD.n(i) -ABD.m(i)*ABD.n(i) ABD.m(i)^2 - ABD.n(i)^2];
    Qt2 = [ABD.m(i)^2 ABD.n(i)^2 ABD.m(i)*ABD.n(i);
          ABD.n(i)^2 ABD.m(i)^2 -ABD.m(i)*ABD.n(i);
          -2*ABD.m(i)*ABD.n(i) 2*ABD.m(i)*ABD.n(i) ABD.m(i)^2 - ABD.n(i)^2];
    
% Qbar Matrix for each Ply
Qb(:,:,i) = Qt1*ABD.Q*Qt2;

% ABD Matricies for each Ply
Aij(:,:,i) = Qb(:,:,i)*(ABD.z(i+1)-ABD.z(i));
Bij(:,:,i) = Qb(:,:,i)*-0.5*((ABD.z(i+1))^2 - (ABD.z(i))^2);
Dij(:,:,i) = Qb(:,:,i)*(1/3)*((ABD.z(i+1))^3 - (ABD.z(i))^3);
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

ABD.A = Aijsum;
ABD.B = Bijsum;
ABD.D = Dijsum;
end
