

classdef utils
    methods(Static)
        function Omegav = Omega(v)

            Omegav = [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0 ];
        end
        function Gammaq = Gamma(q)
            q0 = q(1);
            v = q(2:4);
            Gammaq = q0*eye(4) + [0 -v.';v -utils.Omega(v)];
        end
       
        function Xiq = Xi(q) 
            w = q(1);x=q(2);y=q(3);z=q(4);
            Xiq = [-x -y -z;w -z y;z w -x;-y x w];
        end
            
            
        function r = FF(p,R,t)
            n = size(p,2);
            r = R*p + t*ones(1,n);
        end
        function p = TF(r,R,t)
            n = size(r,2);
            p = R.'*(r-t*ones(1,n));
        end

         function y = PIy(p,f) 
            y(1,1) = f*p(1)/p(2);
            y(2,1) = -f*p(3)/p(2);
         end
         function y = PIz(p,f) 
            y(1,1) = f*p(1)/p(3);
            y(2,1) = f*p(2)/p(3);
         end
         function p = invPIz(y,L,z)
             p = z*L\[y;1];
             p = p(1:2);
         end
    

    function q = qwT(w,T)
        if sum(abs(w))==0
            q = [1;0;0;0];
        else
            s = norm(w);
            u = w/s;

            q = [cos(s*T/2); u*sin(s*T/2)];
        end
    end

    function pq = qMult(p,q)
        u = p(2:4);
        v = q(2:4);
        p0 = p(1);
        q0 = q(1);
        pq = [p0*q0-u.'*v;q0*u + p0*v + cross(u,v)];
    end

    function q = R2q(R)
        T = 1 + trace(R);

        r11 = R(1,1);r12=R(1,2);r13 = R(1,3);
        r21 = R(2,1);r22=R(2,2);r23 = R(2,3);
        r31 = R(3,1);r32=R(3,2);r33 = R(3,3);
        if (T >1e-8)
            S = 2*sqrt(T);
            a= S/4;
            b= (r32 -r23)/S;
            c= (r13- r31)/S;
            d= (r21- r12)/S;
        else
            if (r11 >r22) && (r11 >r33)
                S = 2*sqrt(1 + r11- r22- r33);
                a= (r23- r32)/S;
                b=-S/4;
                c= (r21- r12)/S;
                d= (r13- r31)/S;
            elseif (r22 >r33)
                S = 2*sqrt(1- r11 + r22- r33);
                a= (r31- r13)/S;
                b= (r21- r12)/S;
                c=-S/4;
                d= (r32- r23)/S;
            else
                S = 2*sqrt(1- r11- r22 + r33);
                a= (r12- r21)/S;
                b= (r13- r31)/S;
                c= (r32- r23)/S;
                d=-S/4;
            end
        end
        q = [a;b;c;d];


    end

    function Rq = q2R(q)
        v = q(2:4);
        Rq = (q(1)^2 - (q(2)^2 + q(3)^2 + q(4)^2))*eye(3) + 2*(v*v.') + 2*q(1)*utils.Omega(v);
    end
   function R = v2R(v)
        t = norm(v);
        A = utils.Omega(v);
        R = eye(3) + (sin(t)/t)*A + ((1-cos(t))/t^2)*A^2;
    end

    function R = e2R(phi,theta,psi)%x,y,z
        c1 = cos(phi);
        s1 = sin(phi);
        R1 = [1 0 0;0 c1 s1;0 -s1 c1];
        c2 = cos(theta);
        s2 = sin(theta);
        R2 = [c2 0 -s2;0 1 0;s2 0 c2];
        c3 = cos(psi);
        s3 = sin(psi);
        R3 = [c3 s3 0;-s3 c3 0;0 0 1];
        R = R1*R2*R3;
    end
    function e = q2e(q)
        q = q/norm(q);
        R = utils.q2R(q);

        % Matches utils.e2R(phi, theta, psi)
        theta = asin(max(-1, min(1, -R(1,3))));

        phi = atan2(R(2,3), R(3,3));
        psi = atan2(R(1,2), R(1,1));

        e = [phi; theta; psi];   % radians
    end

    function dtheta = rightAttitudeError(qHat, qTrue)
    % Returns dtheta such that:
    % RTrue = RHat * Exp(Omega(dtheta))

    qHat  = qHat/norm(qHat);
    qTrue = qTrue/norm(qTrue);

    % qErr = qHat^{-1} ⊗ qTrue
    qHatInv = [qHat(1); -qHat(2:4)];
    qErr = utils.qMult(qHatInv, qTrue);
    qErr = qErr/norm(qErr);

    % q and -q represent the same rotation.
    % Select the shortest rotation: angle in [0, pi].
    if qErr(1) < 0
        qErr = -qErr;
    end

    s = norm(qErr(2:4));

    if s < 1e-8
        % qErr ≈ [1; dtheta/2]
        dtheta = 2*qErr(2:4);
    else
        angle = 2*atan2(s, qErr(1));
        dtheta = angle*(qErr(2:4)/s);
    end
end
    function r = TriangulateLS(p,P)
    %estimate landmark positions from camera projections and 
    %known vehicle position and attitude.
        K= size(p,2);
        H = [];d=[];
        for k=1:K
            if sum(abs(p(:,k)))>0
                uk = p(1,k);
                vk = p(2,k);

                Pk = P(:,:,k);
                a1 = Pk(1,1:3);
                a2 = Pk(2,1:3);
                a3 = Pk(3,1:3);



                Hk = [a3*uk-a1; a3*vk-a2];
                dk = [Pk(1,4)-uk*Pk(3,4); Pk(2,4)-vk*Pk(3,4)];

                H = [H; Hk];
                d = [d; dk];
            end
        end

        %r = inv(H.'*H)*H.'*d;
        r = H  \ d;
    end
     function r = TriangulateLSy(p,P)

        K= size(p,2);
        H = [];d=[];
        for k=1:K
            if sum(abs(p(:,k)))>0
                uk = p(1,k);
                vk = p(2,k);

                Pk = P(:,:,k);
                a1 = Pk(1,1:3);
                a2 = Pk(2,1:3);
                a3 = Pk(3,1:3);



                Hk = [a2*uk-a1; a2*vk-a3];
                dk = [Pk(1,4)-uk*Pk(2,4); Pk(3,4)-vk*Pk(2,4)];

                H = [H; Hk];
                d = [d; dk];
            end
        end

        r = inv(H.'*H)*H.'*d;
     end

     function qc = correctq(delta_theta, qx)
         dq = utils.qwT(delta_theta, 1);  % Exp(delta_theta)
         qc = utils.qMult(qx, dq);        % qx ⊗ dq: right/body-local update
         qc = qc / norm(qc);
     end
    end

    
    
end