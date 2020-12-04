function [Cxdata, Cx] = forward_Skf(Ca, Vb, Vi, t, t0, C0, k1, k2, k3, k4)

alpha = [1-Vb, 1-Vb-Vi];

Mx = [[-(k2+k3);k3],[k4;-k4]];
Cx = concentration(k1,Mx,Ca,t0,C0,t);
Cxdata = (alpha*Cx + Vb*Ca(t))';


end