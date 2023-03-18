clc; clear variables; close all;

%Distances
d1 = 800; d2 = 200;

% %Power allocation coefficients
% a1 = 0.75; a2 = 0.25;

N = 10^5;
eta = 4;%Path loss exponent

%Rayleigh fading channels
h11 = sqrt(d1^-eta)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h12 = sqrt(d1^-eta)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h13 = sqrt(d1^-eta)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h21 = sqrt(d2^-eta)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h22 = sqrt(d2^-eta)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h23 = sqrt(d2^-eta)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);

h1_2 = h11+h12;
h2_2 = h21+h22;

h1_3 = h11+h12+h13;
h2_3 = h21+h22+h23;

%Channel gains
g1_2 = (abs(h1_2)).^2;
g2_2 = (abs(h2_2)).^2;

g1_3 = (abs(h1_3)).^2;
g2_3 = (abs(h2_3)).^2;

%Transmit power
Pt = 10; %in dBm
pt = (10^-3)*db2pow(Pt); %linear scale

BW = 10^6;  %bandwidth
%Noise power
No = -174 + 10*log10(BW);   %in dBm
no = (10^-3)*db2pow(No);    %in linear scale

%Target rates
R1 = 0.5:0.5:6;
R2 = R1;
% R2 = R1+3;

p1n_2 = zeros(1,length(R1));
p2n_2 = zeros(1,length(R1));
p1n_3 = zeros(1,length(R1));
p2n_3 = zeros(1,length(R1));
p1f_2 = zeros(1,length(R1));
p2f_2 = zeros(1,length(R1));
p1f_3 = zeros(1,length(R1));
p2f_3 = zeros(1,length(R1));

epsilon = (2.^R1)-1;

aa1 = 0.75; aa2 = 0.25;
for u = 1:length(R1)
% % %     %BASIC FAIR PA%
%      a1 = min(1,epsilon(u)*(no + pt*g1)./(pt*g1*(1+epsilon(u))));
%      a2 = 1 - a1;
% %    
   %IMPROVED FAIR PA MIMO 2*2
   a1_2 = epsilon(u)*(no + pt*g1_2)./(pt*g1_2*(1+epsilon(u)));
   a1_2(a1_2>1) = 0;
   a2_2 = 1 - a1_2;
   
   %IMPROVED FAIR PA MIMO 3*3
   a1_3 = epsilon(u)*(no + pt*g1_3)./(pt*g1_3*(1+epsilon(u)));
   a1_3(a1_3>1) = 0;
   a2_3 = 1 - a1_3;

   %Achievable rates for MIMO-NOMA 2*2 Fair PA
   R1n_2 = log2(1 + pt*a1_2.*g1_2./(pt*a2_2.*g1_2 + no));
   R12n_2 = log2(1 + pt*a1_2.*g2_2./(pt*a2_2.*g2_2 + no));
   R2n_2 = log2(1 + pt*a2_2.*g2_2/no);
   
   %Achievable rates for MIMO-NOMA 3*3 Fair PA
   R1n_3 = log2(1 + pt*a1_3.*g1_3./(pt*a2_3.*g1_3 + no));
   R12n_3 = log2(1 + pt*a1_3.*g2_3./(pt*a2_3.*g2_3 + no));
   R2n_3 = log2(1 + pt*a2_3.*g2_3/no);
   
   %Achievable rates for MIMO-NOMA 2*2 fixed PA
   R1f_2 = log2(1 + pt*aa1.*g1_2./(pt*aa2.*g1_2 + no));
   R12f_2 = log2(1 + pt*aa1.*g2_2./(pt*aa2.*g2_2 + no));
   R2f_2 = log2(1 + pt*aa2.*g2_2/no);
   
   %Achievable rates for MIMO-NOMA 3*3 fixed PA
   R1f_3 = log2(1 + pt*aa1.*g1_3./(pt*aa2.*g1_3 + no));
   R12f_3 = log2(1 + pt*aa1.*g2_3./(pt*aa2.*g2_3 + no));
   R2f_3 = log2(1 + pt*aa2.*g2_3/no);
   
   %Outage calculation
   for k = 1:N
       %MIMO-NOMA 2*2 FAIR POWER USER 1 (FAR)
       if R1n_2(k) < R1(u)
           p1n_2(u) = p1n_2(u)+1;
       end
       %MIMO-NOMA 2*2 FAIR POWER USER 2 (NEAR)
       if a1_2(k) ~= 0
           if (R12n_2(k)<R1(u))||(R2n_2(k) < R2(u))
               p2n_2(u) = p2n_2(u)+1;
           end
       else
           if R2n_2(k) < R2(u)
               p2n_2(u) = p2n_2(u)+1;
           end
       end
       
       %MIMO-NOMA 3*3 FAIR POWER USER 1 (FAR)
       if R1n_3(k) < R1(u)
           p1n_3(u) = p1n_3(u)+1;
       end
       %MIMO-NOMA 3*3 FAIR POWER USER 2 (NEAR)
       if a1_3(k) ~= 0
           if (R12n_3(k)<R1(u))||(R2n_3(k) < R2(u))
               p2n_3(u) = p2n_3(u)+1;
           end
       else
           if R2n_3(k) < R2(u)
               p2n_3(u) = p2n_3(u)+1;
           end
       end
       
       %MIMO-NOMA 2*2 FIXED POWER USER 1 (FAR)
       if R1f_2(k)< R1(u)
           p1f_2(u) = p1f_2(u)+1;
       end
       %MIMO-NOMA 2*2 FIXED POWER USER 2 (NEAR)
       if (R12f_2(k)<R1(u))||(R2f_2(k) < R2(u))
           p2f_2(u) = p2f_2(u)+1;
       end
       
       %MIMO-NOMA 3*3 FIXED POWER USER 1 (FAR)
       if R1f_3(k)< R1(u)
           p1f_3(u) = p1f_3(u)+1;
       end
       %MIMO-NOMA 3*3 FIXED POWER USER 2 (NEAR)
       if (R12f_3(k)<R1(u))||(R2f_3(k) < R2(u))
           p2f_3(u) = p2f_3(u)+1;
       end
   end
end

figure;
semilogy(R1,p1n_2/N,'-*b','linewidth',1.5); hold on; grid on;
semilogy(R1,p2n_2/N,'-ob','linewidth',1.5);
semilogy(R1,p1n_3/N,'-*g','linewidth',1.5); hold on; grid on;
semilogy(R1,p2n_3/N,'-og','linewidth',1.5);
semilogy(R1,p1f_2/N,'-*r','linewidth',1.5); hold on; grid on;
semilogy(R1,p2f_2/N,'-or','linewidth',1.5);
semilogy(R1,p1f_3/N,'-*c','linewidth',1.5); hold on; grid on;
semilogy(R1,p2f_3/N,'-oc','linewidth',1.5);
legend({'MIMO-NOMA 2\times2 far - FAIR PA','MIMO-NOMA 2\times2 near FAIR PA','MIMO-NOMA 3\times3 far - FAIR PA','MIMO-NOMA 3\times3 near FAIR PA','MIMO-NOMA 2\times2 far FIXED PA', 'MIMO-NOMA 2\times2 near FIXED PA','MIMO-NOMA 3\times3 far FIXED PA', 'MIMO-NOMA 3\times3 near FIXED PA'},'location','southeast')
xlabel('Target rate (bps/Hz)');
ylabel('Outage probability');
str = '$$ Bandwidth = 1MHz ,  (Far User = 800m),(Near User = 200m) $$';
text(1.1,0.0003,str,'Interpreter','latex')
title('Outage comparison');
xlim([R1(1) R1(end)])

% figure;
% plot(Pt,R1n_av,'-*b','linewidth',1.5); hold on; grid on;
% plot(Pt,R2n_av,'-ob','linewidth',1.5); 
% plot(Pt,R1o_av,'-*r','linewidth',1.5); 
% plot(Pt,R2o_av,'-or','linewidth',1.5); 
% legend('MIMO-NOMA weak','MIMO-NOMA strong','MIMO-OMA weak', 'MIMO-OMA strong')
% xlabel('Transmit power (dBm)');
% ylabel('Achievable rates (bps/Hz)');
% title('Individual user rates')