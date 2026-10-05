function qc = correctq(delta_theta,qx)


            R = utils.q2R(qx);
       
        R_delta = eye(3) + utils.Omega(delta_theta);
        R_new = R*R_delta;
        qc = utils.R2q(R_new);
        %qc = utils.qMult(qx,delta_q);
        qc = qc/norm(qc);