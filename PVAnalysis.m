clc;
close all;



wavelength = 1550e-9;
n = 32; %number of wavelength channels
spac = 0.63e-9; %spacing between the wavelength channels
Rb = 5e9:1e9:30e9; %bitrate
MAOP = 20; %maximum allowable optical power
Sen = -20; %detector sensitivity
c = 3*(10^8); %speed of light
q0 = 0.04;




%% This snippet of code aids in importing PV maps data into the workspace

% for z = 0:79
%   myfilename = sprintf('InterDie_Q_%d_block100_100.csv', z);
%   mydata_Q{z+1} = importdata(myfilename);
% end
% 
% for z = 0:79
%   myfilename = sprintf('InterDie_ER_%d_block100_100.csv', z);
%   mydata_ER{z+1} = importdata(myfilename);
% end
% 
% for z = 1:80
%     mydata_QyFr{z}.Q = transpose(mydata_Q{1,z});
% end
% 
% for z = 1:80
%     mydata_ExRo{z}.ER = transpose(mydata_ER{1,z});
% end




%% Crosstalk Penalty for different source and destination pairs

for s = 1:16
for d = 1:16
if s ~= d
for br = 1:26
for z = 1:80
for j = 0:n-1
        wl(j+1)=wavelength+j*spac; 
end

for j = 1:n
    freq(j) = (c/wl(j));
end


for i = 1:n
    
    for j = 1:n
        wli(j) = wl(j)+(spac/2);
    end
    
    for j = 1:n
        freql(j)=(c/wli(j));
    end
    
    for j = 1:n
       FWHM_m(z,j) = (freq(j)/mydata_QyFr{1,z}.Q(r(s),j)); 
    end
    
    for j = 1:n
        FWHM(z,j) = (freq(j)/mydata_QyFr{1,z}.Q(r(s),cl(d)+j));
    end
     
    for j = 1:n
        Q_intrinsic(j) = 2*pi*4.2/(100*wl(j));
    end
    
    for j = 1:n
        ring_drop_IL(z,j) = -10*log10(1-(mydata_QyFr{1,z}.Q(r(s),cl(d)+j)/Q_intrinsic(j)));
    end
     
    for j = 1:n
        fd(i,j) = (freq(i)-freq(j));
    end
    
    for j = 1:n
        Fd(i,j) = (freq(i)-freq(j))/((Rb(br)));
    end
    
    for j = 1:n
        Fd1(i,j) = (freq(j)-freql(i))/(1*(Rb(br)));
    end
    
    for j = 1:n
        F_normal(j) = freq(j)/Rb(br);
    end
    
    for j = 1:n
        F_normali(j) = freql(j)/Rb(br);
    end
    
    for j = 1:n
        Ei(z,j) = FWHM(z,j)/(2*(Rb(br)));
    end
    
    for j = 1:n
        Eii(z,j) = FWHM(z,j)/(2*(Rb(br)));
    end
    
    for j = 1:n
        eps(i,j) = fd(i,j);
    end
    
    for j = 1:n
        pp(i,j) = -5*log10((((2*eps(i,j)/FWHM_m(z,i))^2)+q0)/(((2*eps(i,j)/FWHM_m(z,i))^2)+1));
    end
    
    sum = 0;
    for j = 1:n
        if(i ~= j)
            pp_mux(z,i) = sum+pp(i,j);
            sum = pp_mux(z,i);
        end
    end

    for j = 1:n
        beta(i,j) = (2*(fd(i,j)))/(FWHM(z,i));
    end
    
    for j = 1:n
        beta_square(i,j) = beta(i,j)^2;
    end
    
    for j = 1:n
        gamma(i,j) = (1/(1+((beta_square(i,j)))))-((1/(2*pi*Ei(i))).*real((1-exp(-2*pi*Ei(i)*(1-(1j*beta(i,j)))))./((1-beta_square(i,j)-1j*2*beta(i,j)))));
    end
    
    gamma_sqrt = sqrt(gamma);
    
    sum1 = 0;
    for j = 1:n
        p1(i) = sum1+(gamma(i,j));
        sum1 = p1(i);
    end
    
    sum_beta_square = 1./(1+beta_square);
    
    sum2 = 0;
    for j = 1:n
        p2(i) = sum2+(sum_beta_square(i,j));
        sum2 = p2(i);
    end
    
    pp_demux1(z,i) = -5*log10((0.8^2)*(p1(i))*(p2(i)));
    
    sum3 = 0;
    for j = 1:n
        if(i ~= j)
            p3(i) = sum3+(gamma(i,j));
            sum3 = p3(i);
        end
    end
    
    pp_demux2(z,i) = -10*log10(1-(0.5*p3(i)));
    
    
    for j = 1:n 
        for k = 1:i
            fun_1 = @(F) ((((F)+(F_normal(i)-F_normal(j))).^2)/((Ei(z,k).^2)+(((F-(F_normal(i)-F_normal(j))).^2))));
            fun_2 = @(F) ((sinc(F).^2)./(1+(((F-(Fd(i,j)))./Ei(z,i)).^2))).*(fun_1(F).^k);
        end
        fr1(i,j) = integral(fun_2,-Inf,Inf);
    end
    
    for j = 1:n 
        for k = 1:i
            fun_3 = @(F) ((((F)+(F_normali(i)-F_normali(j))).^2)/((Eii(z,k).^2)+(((F+(F_normali(i)-F_normali(j))).^2))));
            fun_4 = @(F) ((sinc(F).^2)./(1+(((F-(Fd1(i,j)))./Eii(z,i)).^2))).*(fun_1(F).^k);
        end
        fr2(i,j) = integral(fun_4,-Inf,Inf);
    end
    
end

for i = 1:n
    pp_demux(z,i) = (pp_demux1(z,i))+(pp_demux2(z,i))+(ring_drop_IL(z,i));
end
    
%% Active Loss
for j = 1:n
    sum = 0;
        for i = 1:n
            if(j ~= i)
                p4(j) = sum+(fr1(i,j));
                sum = p4(j);
            end
        end
end

Active_Loss_dB(z) = abs(10*log10(1-max(p4)));

%% Inactive Loss
for j = 1:n
    sum = 0;
        for i = 1:n
            if(j ~= i)
                p5(j) = sum+(fr2(i,j));
                sum = p5(j); 
            end
        end
end

Inactive_Loss_dB(z) = abs(10*log10(1-max(p5)));


%% Power Budget
Per_waveguide_Budget = MAOP - Sen; 

%% Other Losses
Splitter_Coupling_Loss_dB = 2.5;
WG_Propogation_Loss_dB = WG_Loss_dB(s+d,1);  

%% Power Budget Violation Check
if(imag(pp_demux(z,1:n)) == 0)
    Total_Loss(z) = WG_Propogation_Loss_dB + Splitter_Coupling_Loss_dB + (2*Active_Loss_dB(z)) + ((d-s-1)*Inactive_Loss_dB(z)) + max(pp_mux(z,1:n)) + max(pp_demux(z,1:n));
    Demux_Penalty(z) = max(pp_demux(z,1:n));
    mux_Penalty(z) = max(pp_mux(z,1:n));
    error_function(z) = Per_waveguide_Budget - Total_Loss(z) - 10*log10(n);
    Per_Wavelength_IP(z) = Total_Loss(z) + Sen;

if(error_function(z) > 0)
    disp('Y');
    DC(z,1) = char('N');
else
    disp('Power Budget Violation');
    DC(z,1) = char('V');
end

else
    disp('Power Budget Violation');
    DC(z,1) = char('V');
end

%% BER Evaluation
    for i = 1:n       
        sum = 0;
        for j = 1:n
            if(j ~= i)
                p11(z,i) = sum+(fr1(i,j));
                sum = p11(z,i);
            end
        end
    end
    
Attenuation(z) = WG_Propogation_Loss_dB + Splitter_Coupling_Loss_dB + (2*Active_Loss_dB(z)) + ((d-s-1)*Inactive_Loss_dB(z));
Psignal(z) = (10^(-Attenuation(z)/10))*1e-3;
Pnoise(z) = Psignal(z)*max(p11(z,1:n));
SNR(z) = Psignal(z)/Pnoise(z);
BER(z) = 0.5*erfc(sqrt(SNR(z))/(2*sqrt(2)));


end
%%
brate = Rb(br)/1e9;
filename = sprintf('AnalysisResults_%d_%d_%dGbps.xlsx',s,d,brate);
labels = {'PVMaps','Total Loss (dB)','Per-Wavelength Input Optical Power (dBm)','Bit-Error Rate','Power Budget Violation Check'};
xlswrite(filename,labels,1,'A1:D1');

xlswrite(filename,transpose(pvmaps),1,'A2');
xlswrite(filename,transpose(Total_Loss),1,'B2');
xlswrite(filename,transpose(Per_Wavelength_IP),1,'C2');
xlswrite(filename,transpose(BER),1,'D2');
xlswrite(filename,DC,1,'E2');

end
end
end
end