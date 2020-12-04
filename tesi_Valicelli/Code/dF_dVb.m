function der_M = dF_dVb(k1,M,Ci,t0,C0,t,alpha)
          
der_M = -[1,1]*concentration(k1,M,Ci,t0,C0,t)+Ci(t);
%der_M = der_M;

end 