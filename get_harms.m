function [a0,cn,phi_n,freq] = get_harms(t,y,fund_freq,signal_name)
%GET_HARMS Computes fourier coefficients by fft method
%   INPUTS:
%       - t = time values
%       - y = sample values
%       - fs = sampling frequency
%       - fund_freq = frequency of fundamental component
%       - signal_name = name of variable for for plots
%       - power = average output power
%
%   OUTPUTS:
%       - a0 = 2 times mean value of signal
%       - cn = amplitude coefficients of sine series
%       - phi_n = phase of sine series
%       - freq = freq. vector only for harmonics (excludes mean)
% f(t) = a0/2 + sum(n=1 to N_trunc) [ cn sin(nw0t + phi_n)]
% Output plots always present 50 harmonics

% Gets sampling frequency
fs = 1/(t(2) - t(1));

% Checks if signal_name was passed, if not defaults to 'Sinal'
if ~exist('signal_name','var')
    signal_name = 'Sinal [%]';
end

N_fft = length(y);
df = fs/N_fft;
F = fft(y,N_fft)/N_fft;
a0 = 2*F(1);

cn = 2*abs(F(2:floor(N_fft/2)));
% phi_n = phase(F(2:floor(N_fft/2))) + pi/2;
phi_n = unwrap(angle(F(2:floor(N_fft/2)))) + pi/2;

% Reconstroi sinal
y_teste = a0/2 * ones(length(t),1);

for i=1:length(cn)
    y_teste = y_teste + cn(i)*sin(i*2*pi*df*t + phi_n(i));
end

figure
plot(t,y,'-b')
xlabel('Tempo')
ylabel('Sinal')
hold on
plot(t,y_teste,'--r')
hold off
legend('Original', 'Reconstruído')

freq = (1:N_fft-1)*fs/N_fft;
ix = zeros(1,length(fund_freq:fund_freq:50*fund_freq));
i = 1;
for ref=fund_freq:fund_freq:50*fund_freq
    ix1 = find(freq >= ref,1);
    ix2 = find(freq <= ref, 1, 'last');
    if abs(ref - freq(ix1)) <= abs(ref - freq(ix2))
        ix(i) = ix1;
    else
        ix(i) = ix2;
    end
    i=i+1;
end

% Possibly remove line below
%[~,ix] = findpeaks(cn);

figure
b = bar(cn(ix)/cn(ix(1)) * 100); grid on;
b.FaceColor = [0 0.45 0.85];
ylabel(signal_name,'Fontsize',14,'Fontname','Times New Roman')
xlabel("Índice da harmônica [n]",'Fontsize',14,'Fontname','Times New Roman')
set(gca,'Fontsize',12)


% Plots abs value of specter and selected components
figure
plot(freq(1:ix(end)),cn(1:ix(end))); hold on

cn = cn(ix);
phi_n = phi_n(ix);
freq = freq(ix);

plot(freq,cn,'rs'); hold off

end

