
format long

x = 2
%h = [1e-100 1e-50 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1]
h=[-100:0];
h=10.^h;

f = exp(x)* (sin(x) + cos (x)) + 3*cos(x)*x^2 - sin(x)*x^3

% Taylorentwicklung

f_1 =(exp(x+h).*sin(x+h)+cos(x+h).*(x+h).^3 ...
      -(exp(x).*sin(x)+cos(x).*(x)^3))./h

f_2 = imag(exp(x+i*h).*sin(x+i*h)+cos(x+i*h).*((x+i*h).^3))./h

figure
loglog(h,abs(f-f_1))
figure
loglog(h,abs(f-f_2))

