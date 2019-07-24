function [q , C, PSD] = psd_2D(z , PixelWidth)
% Calculates radially averaged 2D power spectrum. 
% It is better to remove mean and tilt in your topography before applying
% this function
% Here, same pixelwidth in x and y direction is assumed which is typical
% of measurement instruments.

% ===
% inputs:
% z : height topography, a matrix of n*m size (SI units, i.e. meters)

% PixelWidth: size of each Pixel in topography/image (SI units, i.e.
% meters). If you don't know your pixelwidth just devide topography length by
% number of pixels in length.

% ===
% outputs:
% q: wavevectors, which is 2pi/lambda. lambda is wavelength of your
% roughness components.

% C: Radially averaged power spectrum (2D PSD)

% PSD: this structure contain any other interesting information you may
% need.

% in order to plot final result just use:
% loglog(q,C)

% 2D FFT of surface topography could be seen by:
% imagesc(1+log10(abs(PSD.Hm)))

% =========================================================================
% make the size an even number (For PSD)
[n,m] = size(z);
if mod(n,2)
    z = z(2:end,:);
    n = n -1;
end
if mod(m,2)
    z = z(:,2:end);
    m = m -1;
end

Ln = n * PixelWidth; % width of image
Lm = m * PixelWidth; % length of image

a = PixelWidth; % lattice spacing in meter

% =========================================================================
% Window function (up to user)
% win = tukeywin(n,0.25)*tukeywin(m,0.25)';
% win = hann(n)*hann(m)';
% win = kaiser(n,10)*kaiser(m,10)';
% win = bartlett(n) * bartlett(m)';
win = (1 - (((0:n-1)-((n-1)/2))/((n+1)/2)).^2)'*(1 - (((0:m-1)-((m-1)/2))/((m+1)/2)).^2); % Welch window

recInt = trapz(0:n-1,trapz(0:m-1,((rectwin(n)* rectwin(m)').^2),2),1); % integral of squared rectangular window
winInt = trapz(0:n-1,trapz(0:m-1,(win.^2),2),1); % integral of square of selected window function
U = (winInt/recInt); % Normalization constant

z_win = z.*win;
% =========================================================================
% Calculate 2D PSD
Hm = fftshift(fft2(z_win,n,m));
Cq = (1/U)*(a^2/((n*m)*((2*pi)^2)).*((abs((Hm))).^2));
Cq(n/2+1,m/2+1) = 0; % remove the mean 

% =========================================================================
% corresponding wavevectors to Cq values after fftshift has been applied
qx_1=zeros(m,1);
for k=0:m-1
    qx_1(k+1)=(2*pi/m)*(k);
end
qx_2 = fftshift(qx_1);
qx_3 = unwrap(qx_2-2*pi);
qx = qx_3/a;


qy_1=zeros(n,1);
for k=0:n-1
    qy_1(k+1)=(2*pi/n)*(k);
end
qy_2 = fftshift(qy_1);
qy_3 = unwrap(qy_2-2*pi);
qy = qy_3/a;


% =========================================================================
% Radial Averaging
[qxx,qyy] = meshgrid(qx,qy);
[~,rho] = cart2pol(qxx,qyy);
rho = floor(rho);
J = 100; % resolution in q space (increase if you want)
qrmin = log10(sqrt((((2*pi)/Lm)^2+((2*pi)/Ln)^2)));
qrmax = log10(sqrt(qx(end).^2 + qy(end).^2)); % Nyquist
q = floor(10.^linspace(qrmin,qrmax,J));

% =========================================================================
% Averaging Cq values
C_AVE = zeros(1,length(q));
ind = cell(length(q)-1 , 1);
for j = 1:length(q)-1
    ind{j} = find( rho > q(j) & rho <=(q(j+1)));
    C_AVE(j) = nanmean(Cq(ind{j}));
end
ind = ~isnan(C_AVE);
C = C_AVE(ind);
q = q(ind);

% =========================================================================
PSD.Hm = Hm; % 2D FFT
PSD.C = C;
PSD.q = q;
PSD.Cq = Cq;
PSD.qx = qx;
PSD.qy = qy;
PSD.z_win = z_win; % z profile after window function applied

