function der_M = dF_dCp(k1,M,Ci,t0,C0,t,alpha)
          
for i = 1:length(t)
z(:,i) = expm((t(i)-t0)*M)*[0,1]';
end
for i = 1:length(t)
der_M(i) = alpha*z(:,i);
end
   
end 