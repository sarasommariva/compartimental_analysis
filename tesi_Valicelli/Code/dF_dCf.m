function der_M = dF_dCf(k1,M,Ci,t0,C0,t,alpha)
          
for i = 1:length(t) 
z(:,i) = expm((t(i)-t0)*M)*[1,0]';
end
for i = 1:length(t)
der_M(i) = alpha*z(:,i);
end
   
end  