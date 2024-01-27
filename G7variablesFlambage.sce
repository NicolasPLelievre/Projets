n_inputs = 7;
n_outputs = 1;

function out=G(X) 
  // Cas test du flambage de coque
  e=X(:,1);
  Ls=X(:,2);
  hw=X(:,3);
  ew=X(:,4);
  wf=X(:,5);
  ef=X(:,6);
  E=X(:,7);
  //p0=6.35;
  p0=7;
  n=2
  Lc=15e3
  R=2.488e3
  Le = R*(1.56*sqrt(e/R)).*(sqrt(sqrt(1+((n^4)/2)*(e/R)^2)+((n^2)/sqrt(3))*(e/R))).^(-1)
  Ass = hw.*ew + wf.*ef
  Rs = ((R-(e/2+hw/2)).*hw.*ew + (R-(e/2+hw+ef/2)).*wf.*ef ).*Ass^(-1)
  Iss = ew.*hw^3/12 + ew.*hw.*(R-(e/2+hw/2)-Rs)^2 + wf.*ef^3/12 + wf.*ef.*(R-(e/2+hw+e/2)-Rs)^2
  Ac = Ass + e.*Le
  Xc = ( Le.*(e^2)/2 + Ass.*(e/2+(R-Rs))).*Ac^(-1)
  Ic = (e^3).*Le/3 + Iss + Ass.*(e/2+(R-Rs))^2 - Ac.*Xc^2
  be = (n^2-1+0.5*(%pi*R.*Lc^(-1))^2)^(-1).*(n^2*(Lc/%pi/R)^2+1)^(-2)
  pn = E.*e.*be.*R^(-1) + (n^2-1).*E.*Ic.*(R^-3).*Ls^(-1);
  out=pn-p0
endfunction
