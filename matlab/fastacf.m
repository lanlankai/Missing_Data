function [c,lag,c2,c1]=fastacf(x,c2,c1)
% [c,lag]=fastacf(x);
% This function is to estimate the autocorrelation function via xcorr



xi=find(isinf(x)==1);
x(xi)=nan;

xi=find(isnan(x)==0);
x=x(xi(1):xi(end));
xi=find(isnan(x)==1);




M=length(x);

x=x-nanmean(x);
x(xi)=0;
y=ones(size(x));
L=2^nextpow2(2*M-1);
X = fft(x,L); % this is the trick to  have the right correlation function
c = ifft(abs(X).^2);

c(M-1:end-M+2)=0;

X = fft(y,L); % this is the trick to  have the right correlation function
c2 = ifft(abs(X).^2);
c2(c2<1)=1;


y(xi)=0;
X = fft(y,L); % this is the trick to  have the right correlation function
c1 = ifft(abs(X).^2);
c1(c1<1)=1;

c=c./c1.*c2; % correct the missing data effect.
lag=-L/2:L/2-1;

c=fftshift(c);



