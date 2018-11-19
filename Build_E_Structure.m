  function E = Build_E_Structure(T, D, W, N_xy)
  
  global  THERMAL A_GLEN n_GLEN A_0 ETA_0 R_gas Q_creep YRSEC n_XI
  
  if THERMAL
    A      = A_0*exp(-Q_creep./(R_gas*T));
  else
    A      = A_GLEN*ones(n_XI, N_xy);
  end
  
  % Convert strain rates to 1/s from 1/yr for viscosity calculation

  E1     = squeeze(D(1,1,:,:) + D(2,2,:,:) + D(3,3,:,:));
  
  if any(abs(E1)>0)
    fprintf(1, 'BuildSpecialStructure(): Non-zero E1 invariant - max(E1)=%f\n', max(max(E1)));
  end  
  
  E2     = 0.5*(D(1,1,:,:).^2 + D(2,2,:,:).^2 + D(3,3,:,:).^2) + D(1,2,:,:).^2 + D(1,3,:,:).^2 + D(2,3,:,:).^2;
  E2     = squeeze(E2);
    
  s2     = (E2./A.^2).^(1/n_GLEN);         % A has units of Pa^{-3} yr^{-1} and E2 units of 1/yr^2. Thus s2 has Pa^2;
  TAU    = sqrt(s2);                       % Units of TAU are Pa
  ETA    = 1./(1/ETA_0 + 2*A.*TAU.^(n_GLEN-1)/YRSEC);   % Convert Glen coefficient from 1/yr to 1/sec
  
  % Only now convert the strain rate to 1/sec so that stress deviator is Pa
  % Note that at this point ETA (Pa s) and D(i,j)/YRSEC (1/s)

  s      = NaN(3,3,n_XI, N_xy);
  
  s(1,1,:,:) = 2*ETA.*squeeze(D(1,1,:,:))/YRSEC;
  s(1,2,:,:) = 2*ETA.*squeeze(D(1,2,:,:))/YRSEC;
  s(1,3,:,:) = 2*ETA.*squeeze(D(1,3,:,:))/YRSEC;
  
  s(2,1,:,:) = 2*ETA.*squeeze(D(2,1,:,:))/YRSEC;
  s(2,2,:,:) = 2*ETA.*squeeze(D(2,2,:,:))/YRSEC;
  s(2,3,:,:) = 2*ETA.*squeeze(D(2,3,:,:))/YRSEC;

  s(3,1,:,:) = 2*ETA.*squeeze(D(3,1,:,:))/YRSEC;
  s(3,2,:,:) = 2*ETA.*squeeze(D(3,2,:,:))/YRSEC;
  s(3,3,:,:) = 2*ETA.*squeeze(D(3,3,:,:))/YRSEC;
  
  E.A     = A;
  E.D     = D;
  E.W     = W;
  E.ETA   = ETA;
  E.TAU   = TAU;
  E.s     = s;
  
