function [ ps ] = successProbV3( d , pl_model, PLOT_COMMAND )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%close all;

[pl, beta, one_meter_loss] = pathloss( d, pl_model, PLOT_COMMAND );

lamda = 100;
B = 2e9; %mmWave bandwidth


Pb = 30;%dBm
Pb_lin = 10^(Pb/10);

Gmax = 18;%dB
Gmax_lin = 10^(Gmax/10);

C = 0.11;
D = d;
pn = -174 + 10*log10(B) + 10; % noise power
pn_lin = 10^(pn/10);

%%%%%%%%%%%%%
%{
if strcmp(pl_model, '28') || strcmp(pl_model, '73')
      pl_compute = pl;
elseif strcmp(pl_model, '3gpp')
      pl_compute = pl;
end
%}
%%%%%%%%%%%%%

SNR0 = Pb + Gmax - pn;
SNR0_lin=10^(SNR0/10);
SNR = SNR0 - pl;
SNR_lin=10.^(SNR/10);

if strcmp(PLOT_COMMAND, 'yes')
    figure;
    plot(d,SNR)
end

pl_lin = 10.^(pl/10);

xsi_l_lin = 5.2;
%xsi_l_lin = 10^(xsi_l/10);

xsi_n_lin = 7.6;
%xsi_n_lin = 10^(xsi_n/10);

% beta_l_n = 22;%dB
% beta_l_n_lin = 10^(beta_l_n/10);
% 
% beta_l_lin = beta_l_n_lin;
% beta_n_lin = beta_l_n_lin;
% 
% ml = -0.1*beta_l_lin * log(10);
% sigma_l = 0.1* xsi_l_lin * log(10);
% 
% mn = -0.1*beta_n_lin * log(10);
% sigma_n = 0.1* xsi_n_lin * log(10);

beta_l_n = one_meter_loss;%dB
beta_l_n_lin = 10^(beta_l_n/10);

%beta_l_lin = beta_l_n;
%beta_n_lin = beta_l_n;

%ml_db = -0.1*beta_l_lin * log(10);
ml = -log(beta_l_n_lin);  %10^(ml_db/10);
sigma_l = 0.1* xsi_l_lin * log(10);

%mn_db = -0.1*beta_n_lin * log(10);
mn = -log(beta_l_n_lin); %10^(mn_db/10);
sigma_n = 0.1* xsi_n_lin * log(10);

tau = 3;%pl(i)
tau_lin = 10^(tau/10);
    

%for i = 1:length(d)
%factor(i)= SNR_lin(i)/(tau_lin);
% end

for i = 1:length(d)

    %tau = 8;%pl(i)
    %tau_lin = 10^(tau/10);
    
    factor(i)= SNR_lin(i)/(tau_lin);%((Pb_lin*Gmax_lin) / (pl_lin(i)*pn_lin)); 
    

    pc1(i) = D(i)^2*(   qfunc_( (log(D(i)^beta/factor(i))- ml)/sigma_l) - qfunc_( (log(D(i)^beta/factor(i))- mn)/sigma_n)    );
    
    pc2(i) = factor(i)^(2/beta) * exp(2*(sigma_l^2/beta^2) + 2*(ml/beta)) *...
            qfunc_( (sigma_l^2*(2/beta) - log(D(i)^beta/factor(i)) + ml ) / sigma_l   );
        
    pc3(i) = factor(i)^(2/beta) * exp(2*(sigma_n^2/beta^2) + 2*(mn/beta)) *...
        (1/C - qfunc_( (sigma_n^2*(2/beta) - log(D(i)^beta/factor(i)) + mn ) / sigma_n   ));
    
    Lamda_a(i) = lamda*pi*C*(pc1(i) +  pc2(i) + pc3(i));
    
    
    Ma(i) = Lamda_a(i)/lamda;
    
    
    ps(i) = 1 - exp(-lamda*Ma(i)); 
    
    %ps(i) = 1 - exp(-lamda*Ma(i)*factor(i)); 

end

if strcmp(PLOT_COMMAND, 'yes')
    figure 
    plot(d,ps)
end


end

