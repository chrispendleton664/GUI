function Pressure = GUIPressurePreview(GUIInput)

t1p = GUIInput.t1p;

%   time difference
  tt=0.1;
  dt=0.000002;

  NM=tt/dt+1;
  NM=int32(NM); 
  
  t=zeros((NM-1),1);
  pt=zeros((NM-1),1);
  
for i = 1:(NM)-1
     
  t(i+1,:)= t(i,:)+dt;


% Wave parameters

  pm = t1p(1,1);
  alf = t1p(1,2);
  tp = t1p(1,3);

% Blast load eqn    
  pt(i+1,:) = -pm*(1-t(i+1,:)./tp)*exp(-alf*t(i+1,:)/tp);
  
  
  Pressure.pt = pt;
  Pressure.t = t.*1e3;
end

end
