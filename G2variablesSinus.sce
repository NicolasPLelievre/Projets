n_inputs = 2;
n_outputs = 1;

function out=G(X) 
  Nc=2*10^5;
  out=Nc-(0.1*(X(:,1)-X(:,2)).^2+(X(:,1)+X(:,2))/sqrt(2)).*(1+0.05*cos((X(:,1)+X(:,2))*5))*10^5;
endfunction


