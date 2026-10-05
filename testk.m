function r = testk(rn1,vn1,vn,an1,an,w,T,Rn)





%Rn = utils.q2R(q);

r = rn1 + (vn1+vn)/2*T;%-Rn*(a-alpha)*T^2/2;
v = vn1 + Rn*(an+an1)/2*T;

%q = utils.Gamma(utils.qwT(w - 0,T))*q;



end