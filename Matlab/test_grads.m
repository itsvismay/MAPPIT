syms q [1 9] real
syms K real
v1 = q(4:6) - q(1:3);
v2 = q(7:9) - q(4:6);
% syms v1 [1 3] real
% syms v2 [1 3] real
% syms K real
v1xv2 = cross(v1,v2);
v1dv2 = dot(v1, v2);
v1xv2norm = norm(v1xv2);

%%
%Dave's energy
E = 0.5*K*atan2(v1xv2norm, v1dv2)^2
v2m = [0 v2(3) -v2(2); -v2(3) 0 v2(1); v2(2) -v2(1) 0];
v1m = [0 -v1(3) v1(2); v1(3) 0 -v1(1); -v1(2) v1(1) 0]
lodhi1 = v1dv2*(1/v1xv2norm)*(v1xv2*v2m);
hidlo1 = v1xv2norm*v2;
lolo = v1dv2^2;
lodhi2 = v1dv2*(1/v1xv2norm)*(v1xv2*v1m);
hidlo2 = v1xv2norm*v1;
Gp = -K*atan2(v1xv2norm, v1dv2)*(1/(1+(v1xv2norm/v1dv2)^2))*(lodhi1-hidlo1)/lolo;
Gn = K*atan2(v1xv2norm, v1dv2)*(1/(1+(v1xv2norm/v1dv2)^2))*(lodhi2-hidlo2)/lolo;
G = [Gp -Gp-Gn Gn];

q1 = 0; q2 = 0; q3 = 0;
q4 = 0.5; q5 = 0.3; q6 = 0.5;
q7 = 1.2; q8 = 1.7; q9 = 1.1;
K = 1;
v1 = [ q4 - q1, q5 - q2, q6 - q3]; v11 = v1(1); v12 = v1(2); v13 = v1(3);
v2 = [ q7 - q4, q8 - q5, q9 - q6]; v21 = v2(1); v22 = v2(2); v23 = v2(3);
v1xv2 = cross(v1,v2);
v1dv2 = dot(v1, v2);
v1xv2norm = norm(v1xv2);
symGrad = [ (K*atan2((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2), v11*v21 + v12*v22 + v13*v23)*(v21*(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2) - ((v22*(v11*v22 - v12*v21) + v23*(v11*v23 - v13*v21))*(v11*v21 + v12*v22 + v13*v23))/(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2)))/(((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)/(v11*v21 + v12*v22 + v13*v23)^2 + 1)*(v11*v21 + v12*v22 + v13*v23)^2), (K*atan2((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2), v11*v21 + v12*v22 + v13*v23)*(v22*(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2) + ((v21*(v11*v22 - v12*v21) - v23*(v12*v23 - v13*v22))*(v11*v21 + v12*v22 + v13*v23))/(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2)))/(((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)/(v11*v21 + v12*v22 + v13*v23)^2 + 1)*(v11*v21 + v12*v22 + v13*v23)^2), (K*atan2((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2), v11*v21 + v12*v22 + v13*v23)*(v23*(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2) + ((v21*(v11*v23 - v13*v21) + v22*(v12*v23 - v13*v22))*(v11*v21 + v12*v22 + v13*v23))/(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2)))/(((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)/(v11*v21 + v12*v22 + v13*v23)^2 + 1)*(v11*v21 + v12*v22 + v13*v23)^2), (K*atan2((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2), v11*v21 + v12*v22 + v13*v23)*(v11*(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2) + ((v12*(v11*v22 - v12*v21) + v13*(v11*v23 - v13*v21))*(v11*v21 + v12*v22 + v13*v23))/(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2)))/(((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)/(v11*v21 + v12*v22 + v13*v23)^2 + 1)*(v11*v21 + v12*v22 + v13*v23)^2) - (K*atan2((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2), v11*v21 + v12*v22 + v13*v23)*(v21*(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2) - ((v22*(v11*v22 - v12*v21) + v23*(v11*v23 - v13*v21))*(v11*v21 + v12*v22 + v13*v23))/(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2)))/(((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)/(v11*v21 + v12*v22 + v13*v23)^2 + 1)*(v11*v21 + v12*v22 + v13*v23)^2), (K*atan2((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2), v11*v21 + v12*v22 + v13*v23)*(v12*(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2) - ((v11*(v11*v22 - v12*v21) - v13*(v12*v23 - v13*v22))*(v11*v21 + v12*v22 + v13*v23))/(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2)))/(((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)/(v11*v21 + v12*v22 + v13*v23)^2 + 1)*(v11*v21 + v12*v22 + v13*v23)^2) - (K*atan2((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2), v11*v21 + v12*v22 + v13*v23)*(v22*(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2) + ((v21*(v11*v22 - v12*v21) - v23*(v12*v23 - v13*v22))*(v11*v21 + v12*v22 + v13*v23))/(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2)))/(((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)/(v11*v21 + v12*v22 + v13*v23)^2 + 1)*(v11*v21 + v12*v22 + v13*v23)^2), (K*atan2((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2), v11*v21 + v12*v22 + v13*v23)*(v13*(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2) - ((v11*(v11*v23 - v13*v21) + v12*(v12*v23 - v13*v22))*(v11*v21 + v12*v22 + v13*v23))/(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2)))/(((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)/(v11*v21 + v12*v22 + v13*v23)^2 + 1)*(v11*v21 + v12*v22 + v13*v23)^2) - (K*atan2((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2), v11*v21 + v12*v22 + v13*v23)*(v23*(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2) + ((v21*(v11*v23 - v13*v21) + v22*(v12*v23 - v13*v22))*(v11*v21 + v12*v22 + v13*v23))/(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2)))/(((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)/(v11*v21 + v12*v22 + v13*v23)^2 + 1)*(v11*v21 + v12*v22 + v13*v23)^2), -(K*atan2((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2), v11*v21 + v12*v22 + v13*v23)*(v11*(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2) + ((v12*(v11*v22 - v12*v21) + v13*(v11*v23 - v13*v21))*(v11*v21 + v12*v22 + v13*v23))/(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2)))/(((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)/(v11*v21 + v12*v22 + v13*v23)^2 + 1)*(v11*v21 + v12*v22 + v13*v23)^2), -(K*atan2((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2), v11*v21 + v12*v22 + v13*v23)*(v12*(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2) - ((v11*(v11*v22 - v12*v21) - v13*(v12*v23 - v13*v22))*(v11*v21 + v12*v22 + v13*v23))/(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2)))/(((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)/(v11*v21 + v12*v22 + v13*v23)^2 + 1)*(v11*v21 + v12*v22 + v13*v23)^2), -(K*atan2((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2), v11*v21 + v12*v22 + v13*v23)*(v13*(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2) - ((v11*(v11*v23 - v13*v21) + v12*(v12*v23 - v13*v22))*(v11*v21 + v12*v22 + v13*v23))/(abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)^(1/2)))/(((abs(v11*v22 - v12*v21)^2 + abs(v11*v23 - v13*v21)^2 + abs(v12*v23 - v13*v22)^2)/(v11*v21 + v12*v22 + v13*v23)^2 + 1)*(v11*v21 + v12*v22 + v13*v23)^2)];   
symGrad

%%
%Etienne's energy
syms z [1 3] real
%assume(v1xv2norm > 0)
%z = v1xv2/v1xv2norm;
X = norm(v1)*norm(v2) + dot(v1, v2);
Y = dot(cross(v1, v2),z);
angle = 2*atan2(Y, X);
E = 0.5*K*angle^2;
% G1 = [
%       (4*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))*((z3*(q5 - q8) - z2*(q6 - q9))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9)) - ((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))*(q4 - q7 + (abs(q1 - q4)*sign(q1 - q4)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2))/(abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)*((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)/((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))^2 + ((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)
%      -(4*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))*((z3*(q4 - q7) - z1*(q6 - q9))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9)) + ((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))*(q5 - q8 + (abs(q2 - q5)*sign(q2 - q5)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2))/(abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)*((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)/((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))^2 + ((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)
%       (4*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))*((z2*(q4 - q7) - z1*(q5 - q8))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9)) - ((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))*(q6 - q9 + (abs(q3 - q6)*sign(q3 - q6)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2))/(abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)*((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)/((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))^2 + ((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)
%      ];
% G3 = [
%       (4*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))*((z3*(q2 - q5) - z2*(q3 - q6))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9)) + ((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))*(q1 - q4 + (abs(q4 - q7)*sign(q4 - q7)*(abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2))/(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2)))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)*((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)/((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))^2 + ((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)
%      -(4*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))*((z3*(q1 - q4) - z1*(q3 - q6))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9)) - ((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))*(q2 - q5 + (abs(q5 - q8)*sign(q5 - q8)*(abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2))/(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2)))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)*((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)/((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))^2 + ((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)
%       (4*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))*((z2*(q1 - q4) - z1*(q2 - q5))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9)) + ((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))*(q3 - q6 + (abs(q6 - q9)*sign(q6 - q9)*(abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2))/(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2)))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)*((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)/((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))^2 + ((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)
%      ];
%  G2 = -G1 - G3;
%  G = [G1' G2' G3'];
Gjp = -K*(angle)*(cross(v1, z)/(norm(v1)^2));
Gjn = -K*(angle)*(cross(v2, z)/(norm(v2)^2));
Gj = -Gjp - Gjn;
size(Gj)
G = [Gjp Gj Gjn]

K = 1;

q1 = 0; q2 = 0; q3 = 0;
q4 = 0.5; q5 = 0.3; q6 = 0.5;
q7 = 1; q8 = 1; q9 = 1;
v1 = [ q4 - q1, q5 - q2, q6 - q3];
v2 = [ q7 - q4, q8 - q5, q9 - q6];
v1xv2 = cross(v1,v2);
v1dv2 = dot(v1, v2);
v1xv2norm = norm(v1xv2);
z = v1xv2/v1xv2norm;
z1 = z(1); z2 = z(2); z3 = z(3);
e = 2*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2;
 
g1 = [
      (4*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))*((z3*(q5 - q8) - z2*(q6 - q9))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9)) - ((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))*(q4 - q7 + (abs(q1 - q4)*sign(q1 - q4)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2))/(abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)*((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)/((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))^2 + ((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)
     -(4*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))*((z3*(q4 - q7) - z1*(q6 - q9))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9)) + ((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))*(q5 - q8 + (abs(q2 - q5)*sign(q2 - q5)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2))/(abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)*((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)/((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))^2 + ((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)
      (4*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))*((z2*(q4 - q7) - z1*(q5 - q8))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9)) - ((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))*(q6 - q9 + (abs(q3 - q6)*sign(q3 - q6)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2))/(abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)*((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)/((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))^2 + ((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)
     ];
g3 = [
      (4*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))*((z3*(q2 - q5) - z2*(q3 - q6))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9)) + ((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))*(q1 - q4 + (abs(q4 - q7)*sign(q4 - q7)*(abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2))/(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2)))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)*((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)/((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))^2 + ((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)
     -(4*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))*((z3*(q1 - q4) - z1*(q3 - q6))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9)) - ((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))*(q2 - q5 + (abs(q5 - q8)*sign(q5 - q8)*(abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2))/(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2)))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)*((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)/((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))^2 + ((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)
      (4*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))*((z2*(q1 - q4) - z1*(q2 - q5))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9)) + ((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))*(q3 - q6 + (abs(q6 - q9)*sign(q6 - q9)*(abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2))/(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2)))/((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)*((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)/((z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)))^2 + ((abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))^2)
     ];
 g2 = -g1 - g3;
 g = [g1' g2' g3'];
 [ (2*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))*(z3*(q2 - q5) - z2*(q3 - q6)))/(abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2), -(2*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))*(z3*(q1 - q4) - z1*(q3 - q6)))/(abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2), (2*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))*(z2*(q1 - q4) - z1*(q2 - q5)))/(abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)]
  [ (2*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))*(z3*(q5 - q8) - z2*(q6 - q9)))/(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2), -(2*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))*(z3*(q4 - q7) - z1*(q6 - q9)))/(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2), (2*K*atan2(z3*((q1 - q4)*(q5 - q8) - (q2 - q5)*(q4 - q7)) - z2*((q1 - q4)*(q6 - q9) - (q3 - q6)*(q4 - q7)) + z1*((q2 - q5)*(q6 - q9) - (q3 - q6)*(q5 - q8)), (abs(q1 - q4)^2 + abs(q2 - q5)^2 + abs(q3 - q6)^2)^(1/2)*(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)^(1/2) + (q1 - q4)*(q4 - q7) + (q2 - q5)*(q5 - q8) + (q3 - q6)*(q6 - q9))*(z2*(q4 - q7) - z1*(q5 - q8)))/(abs(q4 - q7)^2 + abs(q5 - q8)^2 + abs(q6 - q9)^2)]