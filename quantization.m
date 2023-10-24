clear
close all
clc

%% mi kompanding kvantizator

[x, fs] = audioread('sekvenca2_1.wav');

b = [4 8 12]; % broj bita
mi = [100 500]; 

Xmax = max(abs(x));
M = 2.^b; % broj kvantizacionih nivoa

d = 2*Xmax./M; % korak kvantizacije

t = 1/fs:1/fs:length(x)/fs;

for i = 1:length(mi)
    for j = 1:length(b)
        x_comp = Xmax*log10(1 + mi(i)*abs(x)/Xmax)/(log10(1+mi(i))).*sign(x);
        xq_mi = round(x_comp/d(j))*d(j); % uniformna kvantizacija kompresovanog signala
        xq_mi(x_comp > (M(j) - 1)/2*d(j)) = (M(j)/2 - 1)*d(j);
    
        x_decomp = 1/mi(i)*sign(xq_mi).*((1+mi(i)).^(abs(xq_mi)/Xmax)-1)*Xmax;

        figure()
        plot(x, xq_mi, '*');
        xlabel('x','Interpreter','Latex');
        ylabel('$\hat{x}$','Interpreter','Latex')
        title('Kvantizaciona karakteristika')

        figure()

        hold all
        plot(t,x, 'b');
        plot(t,x_decomp, 'r--')
        xlabel('t[s]')
        title('Vremenski oblik signala')
        legend('original','kvantizovan')

        fprintf('SNR = %f\n', 10*log10(var(x)/var(x-x_decomp)));
    end
end

%% SNR karakteristika mi kompanding kvantizatora

for i = 1:length(mi)
    for j = 1:length(b)
        aten = 0.1:0.01:1.5; 
        SNR1 = [];
        xvar1 = [];

        for k=1:length(aten)
            x1 = aten(k)*x; % skaliranje (pojacavanje/slabljenje)
            xvar1 = [xvar1 var(x1)];
            x_comp = Xmax*log10(1+mi(i)*abs(x1)/Xmax)/(log10(1+mi(i))).*(sign(x1));
            xq_mi = round(x_comp/d(j))*d(j);
            x_mi_decomp =1/mi(i)*sign(xq_mi).*((1+mi(i)).^(abs(xq_mi)/Xmax)-1)*Xmax;
            SNR1 = [SNR1 10*log10(var(x1)/var(x1-x_mi_decomp))]; % eksperimentalni SNR
        end

        figure(12+j);
        semilogx(Xmax./(sqrt(xvar1)),SNR1,'b');
        hold on
        % teorijski SNR
        semilogx(Xmax./(sqrt(xvar1)),4.77+6*b(j)-20*log10(log(1+mi(i)))-10*log10(1+(Xmax./mi(i))^2./xvar1+sqrt(2)*Xmax./mi(i)./sqrt(xvar1)),'r--');
        legend('Eksperimentalno','Teorijski')
        title('SNR karakteristika')
        xlabel('$X_{max}/\sigma_x$','Interpreter','Latex')
        ylabel('SNR[dB]','Interpreter','Latex')
    end
end

figure(13)
text(10,15.2,'mi = 100')
text(10,13,'mi = 500')

figure(14)
text(10,39.5,'mi = 100')
text(10,37.25,'mi = 500')

figure(15)
text(10,63.5,'mi = 100')
text(10,61.25,'mi = 500')

%% Delta kvantizator

N = length(x);

Q = 0.01; % korak kvantizacije

d = zeros(1,N); % pamtimo prirastaj da bismo crtali histogram
d(1) = x(1);

x_delta = zeros(1,N);
x_delta(1) = Q;

for i = 2:N
    d(i) = x(i)-x_delta(i-1); % prirastaj
    x_delta(i) = x_delta(i-1)+Q*(2*(d(i) > 0)-1);
end

t = 1/fs:1/fs:length(x)/fs; 
figure();
plot(t,x)
hold on
plot(t,x_delta)
hold off
legend('Originalni', 'Kvantizovani');
xlabel('t[s]')
title('Vremenski oblik signala')

figure();
histogram(d);
xlim([-0.15,0.15])
title('Histogram prirastaja')
