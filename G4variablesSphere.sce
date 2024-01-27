n_inputs = 4;
n_outputs = 1;

function out=G(X) 
  // Cas test de la sphere
  // Variable 1 : fy
  // Variable 2 : P0
  // Variable 3 : r1
  // Variabe 4 : r0
  out=X(:,1)-3/2*X(:,2).*(X(:,3)^3).*(((X(:,3)^3-X(:,4)^3)).^(-1));
endfunction


