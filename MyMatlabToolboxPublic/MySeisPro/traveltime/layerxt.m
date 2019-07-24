function [dx, dt, irtr] = layerxt(p, h, utop, ubot)

if (p >= utop)
   dx = 0.0;
   dt = 0.0;
   irtr = 0;
   return
elseif (h == 0)
  dx = 0.0;
  dt = 0.0;
  irtr = -1;
  return
end

u1 = utop;
u2 = ubot;
v1 = 1.0/u1;
v2 = 1.0/u2;
b = (v2 - v1)/h;

eta1 = sqrt(u1^2 - p^2);

if (b == 0.0)
   dx = h*p/eta1;
   dt = h*u1^2/eta1;
   irtr = 1;
   return
end

x1 = eta1/(u1*b*p);
tau1 = (log((u1+eta1)/p)-eta1/u1)/b;

if (p >= ubot)
   dx = x1;
   dtau = tau1;
   dt = dtau+p*dx;
   irtr = 2;
   return
end

irtr = 1;

eta2 = sqrt(u2^2-p^2);
x2 = eta2/(u2*b*p);
tau2 = (log((u2+eta2)/p) - eta2/u2)/b;

dx = x1-x2;
dtau = tau1-tau2;

dt = dtau+p*dx;

return;
