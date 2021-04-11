function FailureCriteria = GUIFailLoadTsaiWu(theta, GUIInput)


dif = 0.1;
FailureCriteria.pt = 10000;

ABD =  GUIFailureABD_and_Strain(theta, GUIInput);

while (abs(dif)>1e-3)

% FailureFaceSheetDeflection Inputs
FailureCriteria.Inp = GUIFailureFaceSheetDeflection(ABD, FailureCriteria, GUIInput);

inc = 10000 * dif;

FailureCriteria.pt = FailureCriteria.pt + inc;

dif = 1 - max(FailureCriteria.Inp(:,3));

end

end