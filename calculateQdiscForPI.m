% Input ad,bd shall be the discretized PI controller.
% Input Kp and Ki shall be positive PI parameters.
% Input R is normally one.
function [Q, P] = calculateQdiscForPI(ad,bd,Ki,Kp,R)

a=ad;
b=bd;

% State feedback is u=-Kx=-[-ki,-kp]x=ki*x1+kp*x2 = PI
%Ki = -Ki;
%Kp = -Kp;

P12 = (Ki*R*a)/(b^2*(Ki-Kp)+b*a);
P22 = R*(Kp-Ki)/(b^2*(Ki-Kp)+b*a);
P11 = P12*(1-a) + b^2*(P12^2 + a*P12*P22)/(R+b^2*P22);


Q1 = (P12^2)/((R/b^2)+P22);
Q2 = P22*(1-a^2) - P11 - 2*a*P12 + ((P12+a*P22)^2)/(R/(b^2) + P22);

Q=[Q1 0;0 Q2];
P=[P11 P12;P12 P22];


%[Kri, Pri,] = dlqr([1 1;0 a], [0;b], Q, R, 0);

%disp(['Q1: ', num2str(Q1), '. Q2: ', num2str(Q2), ...
%      '. If forward Riccati is solved with Q, we get Kp: ', ...
%      num2str(Kri(2)), ' and Ki: ', num2str(Kri(1))]);

end