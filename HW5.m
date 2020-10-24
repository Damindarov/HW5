syms q1 q2 s l m1 m2 I1 I2 g dq1 dq2 ddq1 ddq2 real

T00 = eye(4);
T01 = Rz(q1) * Tx(q2)
T02 = Rz(q1) * Tx(q2) * Tx(2*l)

Z0 = [0;0;1];
Z1 = [cos(q1);sin(q1);0];%first column

O0 = [0;0;0];
Oc1 = [0;0;0];
O1 = [q2*cos(q1); q2*sin(q1); 0];
Oc2 =[(q2 + l)*cos(q1);(q2 + l)*sin(q1);0];
O2 = [(q2 + 2*l)*cos(q1);(q2 + 2*l)*sin(q1);0];
R1 = Rz(q1);
R1 = R1(1:3,1:3);

q = [q1,q2];

Zero = [0;0;0];

Jv1 = [cross(Z0,Oc1-O0),Zero];
Jv2 = [cross(Z0,Oc2-O0),Z1];%Z1 tk prismatic
Jw1 = [Z0,Zero]
Jw2 = [Z0,Zero]

D1 = m1*Jv1'*Jv1 + Jw1'*R1*I1*R1'*Jw1;
D2 = m2*Jv2'*Jv2 + Jw2'*R1*I2*R1'*Jw2;
D = D1 + D2;
D = simplify(D)

%potential energy
P1 = 0;
P2 = m2*(q2 + l)*g*sin(q1);
P = P1 + P2
G1 = diff(P,q1);
G2 = diff(P,q2);
G = [G1;G2];
dq = [dq1;dq2];
ddq = [ddq1; ddq2];
C = Coriolis(D,q,dq,2);
C = simplify(C);
tor = D*ddq + C*dq + G;%ure dvish invdinam

%re
D(q1,q2) = subs(D,{m1, m2, I1, I2, l},{2 2 1 2 0.2});
C(q1,q2,dq1,dq2) = subs(C*dq ,{m1, m2, I1, I2, l},{2 2 1 2 0.2});
G(q1,q2) = subs(G ,{m1, m2, I1, I2, l,g},{2 2 1 2 0.2 9.81});

q1_0 = -pi/2;
q2_0 = 1;

dq1_0 = 0;
dq2_0 = 0;

delta_t = 0.01;

U = [0,0];
% U = [0;0];
n=100;
ddq_1 = [];
ddq_2 = [];
tor_1p = [];
tor_2p = [];
t = 0;
U_1p = [];
U_2p = [];
for i = 1:n
 U = [sin(t^2),t];
 U_1p = [U_1p U(1)];
 U_2p = [U_2p U(2)];
 q1p(i)=q1_0; 
 q2p(i)=q2_0;
 dq1p(i)=dq1_0;
 dq2p(i)=dq2_0;
 ddq = inv(D(q1_0, q2_0))*(U-C(q1_0, q2_0,dq1_0,dq2_0)-G(q1_0,q2_0));
 ddq_1 = [ddq_1 ddq(1)]; 
 ddq_2 = [ddq_2 ddq(2)]; 
    
 dq1_0=dq1p(i) + double(ddq(1)*delta_t);
 dq2_0=dq2p(i) + double(ddq(2)*delta_t);
   
 q1_0 = q1p(i) + dq1_0*delta_t;
 q2_0 = q2p(i) + dq2_0*delta_t;
 t = t+delta_t;
 if q2_0 < 0 
     q2_0 = 0;
 end
end
t = 0:delta_t:(delta_t*(n-1));
plot(t,q1p);
title('q1 vs time');
figure;
plot(t,q2p);
title('q2 vs time');
figure;
plot(t,ddq_1);
title('ddq1 vs time');
figure;
plot(t,ddq_2);
title('ddq2 vs time');

figure;
plot(t,dq1p);
title('dq1 vs time');
figure;
plot(t,dq2p);
title('dq2 vs time');

figure;
plot(t,U_1p);
title('U1 vs time');

figure;
plot(t,U_2p);
title('U2 vs time');



function C = Coriolis(D,q,dq,n)
sym C;
for k = 1:n
     for j =1:n
     C(k,j) = sym(0);
         for i=1:n
            c_ijk = (1/2)*(diff(D(k,j),q(i)) + diff(D(k,i),q(j))-diff(D(i,j),q(k)));
            C(k,j) = C(k,j) + c_ijk*dq(i);
         end
     end
end
end