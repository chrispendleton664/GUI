function FailureLoad = GUILaminate_Orientation_Hoffman(GUIInput)

theta = zeros(90,1);
FailureLoad = zeros(90,2);
Mmax = zeros(90,1);

for j = 0:1:90
    
theta(j+1,:) = j;

theta = theta(j+1,1);

Hoffman = GUIFailLoad(theta, GUIInput);

FailureLoad(j+1,:) = [theta,Hoffman.pt];


end

% Mmax = FailureLoad(:,:).* (3/8 * 48e-3 * (96e-3)^2);

% Mmax = FailureLoad(:,:); 

    failure_criteria = {'Theta', 'Failure Pressure'};
  
    FailureLoad = array2table(FailureLoad, 'VariableNames', failure_criteria);

    filename = 'LamApp Hoffman Failure Results.xlsx';
    writetable(FailureLoad,filename,'Sheet',1);

end
