.gumb.123.density = function(x, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2];u3 = u[3]}else{u1 = u[,1];u2 = u[,2];u3 = u[,3]}
  p123 = (-log(u1))^x + (-log(u2))^x + (-log(u3))^x
  p123.x = p123^(1/x)
  ((1 + 2*(x^2) + 3*x*(-1 + p123.x) - 3*p123.x + p123.x^2)*p123^(-3)*((-(log(u1)*log(u2)*log(u3)))^(-1 + x)))/(exp(p123.x)*u1*u2*u3)
}

.fran.123.density = function(x, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2];u3 = u[3]}else{u1 = u[,1];u2 = u[,2];u3 = u[,3]}
  
  au1 = -1 + exp(-x * u1)
  au2 = -1 + exp(-x * u2)
  au3 = -1 + exp(-x * u3)
  au123 = au1*au2*au3
  au  = -1 + exp(-x)
  a = (au^2)*(1 + au123/(au^2))
  e123 = exp(-x*(u1+u2+u3))
  
  (2 * e123 * (au123^2) * x^2) / (a^3) - (3 * e123 * au123 * x^2) / (a^2) + (e123 * x^2) / a
}

.clay.123.density = function(x, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2];u3 = u[3]}else{u1 = u[,1];u2 = u[,2];u3 = u[,3]}
  u1pt = u1^(-x)
  u2pt = u2^(-x)
  u3pt = u3^(-x)
  u123 = u1*u2*u3
  ((1 + x)*(1 + 2 * x)*((u123)^(2 * x - 1)))/(((u1pt + u2pt + u3pt-2)^(1/x))*((u1*u2)^x + (u1*u3)^x + (u2*u3)^x - 2*(u123)^x)^3)
}

.l.gumb.123.density = function(x, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2];u3 = u[3]}else{u1 = u[,1];u2 = u[,2];u3 = u[,3]}
  plu1 = (-log(u1))^x
  plu2 = (-log(u2))^x
  plu3 = (-log(u3))^x
  p123 = plu1 + plu2 + plu3
  p123.x = p123^(1/x)
  sum(log(((1 + 2*(x^2) + 3*x*(-1 + p123.x) - 3*p123.x + p123.x^2)*((-(log(u1)*log(u2)*log(u3)))^(-1 + x)))/(p123^(3) * exp(p123.x)*u1*u2*u3)))
}

.l.gumb.123.density.t= function(x, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2];u3 = u[3]}else{u1 = u[,1];u2 = u[,2];u3 = u[,3]}
  x = 1 / (1 - x)
  p123 = (-log(u1))^x + (-log(u2))^x + (-log(u3))^x
  p123.x = p123^(1/x)
  sum(log(1 + 2*(x^2) + 3*x*(p123.x-1) - 3*p123.x + p123.x^2)+(x-1)*log(-log(u1)*log(u2)*log(u3)) - 3*log(p123) - p123.x - log(u1*u2*u3))
}

.l.fran.123.density = function(x, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2];u3 = u[3]}else{u1 = u[,1];u2 = u[,2];u3 = u[,3]}
  
  au123 = (-1 + exp(-x * u1))*(-1 + exp(-x * u2))*(-1 + exp(-x * u3))
  au  = -1 + exp(-x)
  a = (au^2)*(1 + au123/(au^2))
  e123 = exp(-x*(u1+u2+u3))
  
  sum(log((2 * e123 * (au123^2) * x^2) / (a^3) - (3 * e123 * au123 * x^2) / (a^2) + (e123 * x^2) / a))
}

.l.clay.123.density = function(x, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2];u3 = u[3]}else{u1 = u[,1];u2 = u[,2];u3 = u[,3]}
  u1pt = u1^(-x)
  u2pt = u2^(-x)
  u3pt = u3^(-x)
  sum(log(((1 + x)*(1 + 2 * x)*((u1*u2*u3)^(2 * x - 1)))/(((u1pt + u2pt + u3pt-2)^(1/x))*((u1*u2)^x + (u1*u3)^x + (u2*u3)^x - 2*(u1*u2*u3)^x)^3)))
}

.l.clay.123.density.t.ls = function(x, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2];u3 = u[3]}else{u1 = u[,1];u2 = u[,2];u3 = u[,3]}
  x = 2*x/(1-x)
  u1pt = u1^(-x)
  u2pt = u2^(-x)
  u3pt = u3^(-x)
  u123 = u1*u2*u3
  sum(log(1 + x) + log(1 + 2 * x) + (2 * x - 1)*log(u123) - log(u1pt + u2pt + u3pt-2) / x + 3*log((u1*u2)^x + (u1*u3)^x + (u2*u3)^x - 2*(u123)^x))
}

.gumb.123.e2 = function(theta, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2]}else{u1 = u[,1];u2 = u[,2]}
  lu1 = -log(u1)
  plu12t = lu1^theta + (-log(u2))^theta
  ((lu1^(theta-1))*(plu12t^(1/theta-1)))/(exp((plu12t)^(1/theta))*u1)
}

.gumb.123.e3 = function(theta, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2];u3 = u[3]}else{u1 = u[,1];u2 = u[,2];u3 = u[,3]}
  plu1 = (-log(u1))^theta
  plu2 = (-log(u2))^theta
  plu3 = (-log(u3))^theta
  plu12 = plu1 + plu2
  plu123 = plu12 + plu3
  plu12.t = plu12^(1/theta)
  plu123.t = plu123^(1/theta)
  (exp(plu12.t - plu123.t)*(-1 + theta + plu123.t)*((plu12/plu123)^(2 - 1/theta)))/(-1 + theta + plu12.t)
}

.fran.123.e2 = function(theta, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2]}else{u1 = u[,1];u2 = u[,2]}
  (exp(theta)*(-1 + exp(theta*u2)))/(-exp(theta) + exp(theta*(1 + u1)) - exp(theta*(u1 + u2)) + exp(theta*(1 + u2)))
}

.fran.123.e3 = function(theta, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2];u3 = u[3]}else{u1 = u[,1];u2 = u[,2];u3 = u[,3]}
  
  (exp(theta*(1+u3))*(-1 + exp(theta))*((exp(theta) - exp(theta*(1+u1)) + exp(theta*(u1 + u2)) - exp(theta*(1+u2)))^2)*(-1 + exp(theta*u3)))/((exp(2*theta) - exp(theta*(2 + u1)) - exp(theta*(2 + u2)) + exp(theta*(2 + u1 + u2)) - exp(theta*(2 + u3)) + exp(theta*(2 + u1 + u3)) + exp(theta*(2 + u2 + u3)) + exp(theta*(u1 + u2 + u3)) - 2*exp(theta*(1 + u1 + u2 + u3)))^2)
}


.clay.123.e2 = function(theta, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2]}else{u1 = u[,1];u2 = u[,2]}
  u1pt = u1^theta
  u2pt = u2^theta
  -(u2pt/(u1*((-1 + 1/u1pt + 1/u2pt)^(1/theta))*(-u2pt + u1pt*(-1 + u2pt))))
}

.clay.123.e3 = function(theta, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2];u3 = u[3]}else{u1 = u[,1];u2 = u[,2];u3 = u[,3]}
  u1pt = u1^(-theta)
  u2pt = u2^(-theta)
  ((-1 + u1pt + u2pt)^(2 + 1/theta))*((-2 + u1pt + u2pt + u3^(-theta))^(-2 - 1/theta))
}

.clay.S.123 = function(theta, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2];u3 = u[3]}else{u1 = u[,1];u2 = u[,2];u3 = u[,3]}
  
  t2 = theta^2
  t3 = theta^3
  t4 = theta^4
  t5 = theta^5
  t6 = theta^6
  u12 = u1*u2
  u13 = u1*u3
  u23 = u2*u3
  u123 = u1*u2*u3
  
  u1t = .Power(u1,-theta)
  u2t = .Power(u2,-theta)
  u3t = .Power(u3,-theta)
  u1tt = .Power(u1,theta)
  u2tt = .Power(u2,theta)
  u3tt = .Power(u3,theta)
  u1tt2 = u1tt^2
  u2tt2 = u2tt^2
  u3tt2 = u3tt^2
  l3 = log(u3)
  l3p2 = .Power(l3,2)
  lu123t = log(-2 + u1t + u2t + u3t)
  
  a12t = .Power(u12,theta)
  a13t = .Power(u13,theta)
  a23t = .Power(u23,theta)
  
  a12t2 = a12t^2
  a13t2 = a13t^2
  a23t2 = a23t^2
  
  -((5*t3*a12t2 + 12*t4*a12t2 + 8*t5*a12t2 - 20*t3*a12t2*u3tt - 48*t4*a12t2*u3tt - 32*t5*a12t2*u3tt + 10*t3*a12t*u3tt2 + 24*t4*a12t*u3tt2 + 16*t5*a12t*u3tt2 + 10*t3*u2tt2*a13t + 24*t4*u2tt2*a13t + 16*t5*u2tt2*a13t + 5*t3*a13t2 + 12*t4*a13t2 + 8*t5*a13t2 - 20*t3*u2tt*a13t2 - 48*t4*u2tt*a13t2 - 32*t5*u2tt*a13t2 + 10*t3*u1tt2*a23t + 24*t4*u1tt2*a23t + 16*t5*u1tt2*a23t + 5*t3*a23t2 + 12*t4*a23t2 + 8*t5*a23t2 - 20*t3*u1tt*a23t2 - 48*t4*u1tt*a23t2 - 32*t5*u1tt*a23t2 + 20*t3*.Power(u123,2*theta) + 48*t4*.Power(u123,2*theta) + 32*t5*.Power(u123,2*theta) + t2*(1 + 3*theta)*.Power(1 + 3*theta + 2*t2,2)*(a12t*u3tt2 + u2tt2*a13t - 2*u1tt*a23t2)*.Power(log(u1),2) + t2*(1 + 3*theta)*.Power(1 + 3*theta + 2*t2,2)*(a12t*u3tt2 - 2*u2tt*a13t2 + u1tt2*a23t)*.Power(log(u2),2) + 2*theta*a12t2*l3 + 12*t2*a12t2*l3 + 26*t3*a12t2*l3 + 24*t4*a12t2*l3 + 8*t5*a12t2*l3 - 4*theta*a12t2*u3tt*l3 - 24*t2*a12t2*u3tt*l3 - 52*t3*a12t2*u3tt*l3 - 48*t4*a12t2*u3tt*l3 - 16*t5*a12t2*u3tt*l3 + 2*theta*u2tt2*a13t*l3 + 12*t2*u2tt2*a13t*l3 + 26*t3*u2tt2*a13t*l3 + 24*t4*u2tt2*a13t*l3 + 8*t5*u2tt2*a13t*l3 + 2*theta*u1tt2*a23t*l3 + 12*t2*u1tt2*a23t*l3 + 26*t3*u1tt2*a23t*l3 + 24*t4*u1tt2*a23t*l3 + 8*t5*u1tt2*a23t*l3 - 2*t2*a12t2*u3tt*l3p2 - 18*t3*a12t2*u3tt*l3p2 - 62*t4*a12t2*u3tt*l3p2 - 102*t5*a12t2*u3tt*l3p2 - 80*t6*a12t2*u3tt*l3p2 - 24*.Power(theta,7)*a12t2*u3tt*l3p2 + t2*u2tt2*a13t*l3p2 + 9*t3*u2tt2*a13t*l3p2 + 31*t4*u2tt2*a13t*l3p2 + 51*t5*u2tt2*a13t*l3p2 + 40*t6*u2tt2*a13t*l3p2 + 12*.Power(theta,7)*u2tt2*a13t*l3p2 + t2*u1tt2*a23t*l3p2 + 9*t3*u1tt2*a23t*l3p2 + 31*t4*u1tt2*a23t*l3p2 + 51*t5*u1tt2*a23t*l3p2 + 40*t6*u1tt2*a23t*l3p2 + 12*.Power(theta,7)*u1tt2*a23t*l3p2 - 2*theta*.Power(1 + 3*theta + 2*t2,2)*log(u1)*(-(a12t*u3tt2) - u2tt2*a13t - a23t2 + 2*u1tt*a23t2 +   theta*(1 + 3*theta)*a12t*u3tt2*log(u2) + theta*(1 + 3*theta)*u2tt2*a13t*l3) + 2*theta*.Power(1 + 3*theta + 2*t2,2)*log(u2)*(a12t*u3tt2 + a13t2 - 2*u2tt*a13t2 + u1tt2*a23t -   theta*(1 + 3*theta)*u1tt2*a23t*l3) + 2*a12t2*lu123t + 12*theta*a12t2*lu123t + 26*t2*a12t2*lu123t + 24*t3*a12t2*lu123t + 8*t4*a12t2*lu123t - 8*a12t2*u3tt*lu123t - 48*theta*a12t2*u3tt*lu123t - 104*t2*a12t2*u3tt*lu123t - 96*t3*a12t2*u3tt*lu123t - 32*t4*a12t2*u3tt*lu123t + 4*a12t*u3tt2*lu123t + 24*theta*a12t*u3tt2*lu123t + 52*t2*a12t*u3tt2*lu123t + 48*t3*a12t*u3tt2*lu123t + 16*t4*a12t*u3tt2*lu123t + 4*u2tt2*a13t*lu123t + 24*theta*u2tt2*a13t*lu123t + 52*t2*u2tt2*a13t*lu123t + 48*t3*u2tt2*a13t*lu123t + 16*t4*u2tt2*a13t*lu123t + 2*a13t2*lu123t + 12*theta*a13t2*lu123t + 26*t2*a13t2*lu123t + 24*t3*a13t2*lu123t + 8*t4*a13t2*lu123t - 8*u2tt*a13t2*lu123t - 48*theta*u2tt*a13t2*lu123t - 104*t2*u2tt*a13t2*lu123t - 96*t3*u2tt*a13t2*lu123t - 32*t4*u2tt*a13t2*lu123t + 4*u1tt2*a23t*lu123t + 24*theta*u1tt2*a23t*lu123t + 52*t2*u1tt2*a23t*lu123t + 48*t3*u1tt2*a23t*lu123t + 16*t4*u1tt2*a23t*lu123t + 2*a23t2*lu123t + 12*theta*a23t2*lu123t + 26*t2*a23t2*lu123t + 24*t3*a23t2*lu123t + 8*t4*a23t2*lu123t - 8*u1tt*a23t2*lu123t - 48*theta*u1tt*a23t2*lu123t - 104*t2*u1tt*a23t2*lu123t - 96*t3*u1tt*a23t2*lu123t - 32*t4*u1tt*a23t2*lu123t + 8*.Power(u123,2*theta)*lu123t + 48*theta*.Power(u123,2*theta)*lu123t + 104*t2*.Power(u123,2*theta)*lu123t + 96*t3*.Power(u123,2*theta)*lu123t + 32*t4*.Power(u123,2*theta)*lu123t)/(t3*.Power(1 + theta,2)*.Power(1 + 2*theta,2)*.Power(a12t + a13t + a23t - 2*.Power(u123,theta),2)))
}

.clay.V.123 = function(theta, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2];u3 = u[3]}else{u1 = u[,1];u2 = u[,2];u3 = u[,3]}
  
  t2 = theta^2
  t3 = theta^3
  t4 = theta^4
  t5 = theta^5
  t6 = theta^6
  u12 = u1*u2
  u13 = u1*u3
  u23 = u2*u3
  u123 = u1*u2*u3
  l3 = log(u3)
  u12t = u12^theta
  u13t = u13^theta
  u23t = u23^theta
  u1t = u1^(-theta)
  u2t = u2^(-theta)
  u3t = u3^(-theta)
  a = log(-2 + u1t + u2t + u3t)
  p123 = u123^theta
  
  (3*t2*u12t + 4*t3*u12t + 3*t2*u13t + 4*t3*u13t + 3*t2*u23t + 4*t3*u23t - 6*t2*p123 - 8*t3*p123 + theta*(1 + 3*theta + 2*t2)*(u23t - theta*(u12t + u13t - 2*u23t - 2*p123))*log(u1) + theta*(1 + 3*theta + 2*t2)*(u13t - theta*(u12t - 2*u13t + u23t - 2*p123))*log(u2) + theta*u12t*l3 + 5*t2*u12t*l3 + 8*t3*u12t*l3 + 4*t4*u12t*l3 - t2*u13t*l3 - 3*t3*u13t*l3 - 2*t4*u13t*l3 - t2*u23t*l3 - 3*t3*u23t*l3 - 2*t4*u23t*l3 + 2*t2*p123*l3 + 6*t3*p123*l3 + 4*t4*p123*l3 + u12t*a + 3*theta*u12t*a + 2*t2*u12t*a + u13t*a + 3*theta*u13t*a + 2*t2*u13t*a + u23t*a + 3*theta*u23t*a + 2*t2*u23t*a - 2*p123*a - 6*theta*p123*a - 4*t2*p123*a)/(t2*(1 + theta)*(1 + 2*theta)*(u12t + u13t + u23t - 2*p123))
}

.gumb.S.123 = function(theta, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2];u3 = u[3]}else{u1 = u[,1];u2 = u[,2];u3 = u[,3]}
  
  l1 = log(u1)
  l2 = log(u2)
  l3 = log(u3)
  ll1 = log(-l1)
  ll2 = log(-l2)
  ll3 = log(-l3)
  
  ti = 1/theta
  
  pu1t = (-l1)^theta
  pu2t = (-l2)^theta
  pu3t = (-l3)^theta
  
  pl13t = (l1*l3)^theta
  pl12t = (l1*l2)^theta
  pl23t = (l2*l3)^theta
  
  pu123t = pu1t + pu2t + pu3t
  lt = log(pu123t)
  b = pu123t^ti
  b2 = b^2
  b3 = b^3
  b4 = b^4
  ll32 = ll3^2
  a = b*ll32
  lpu123t = lt^2
  
  
  t2 = theta^2
  t3 = theta^3
  t4 = theta^4
  t5 = theta^5
  t6 = theta^6
  
  (-.Power(-3*t2*pu1t + 4*t3*pu1t - 3*t2*pu2t + 4*t3*pu2t + 3*t2*pu1t*b + 3*t2*pu2t*b - 3*t2*pu3t + 4*t3*pu3t + 3*t2*b*pu3t - theta*(pu1t*(-1 + 7*b - 6*b2 + b3) + t2*(-1 + b)*(8*pu1t - 3*pu2t - 3*pu3t) + t3*(4*pu1t - 2*pu2t - 2*pu3t) + theta*(1 - 3*b + b2)*(5*pu1t - pu2t - pu3t))*ll1 - theta*(pu2t*(-1 + 7*b - 6*b2 + b3) - theta*(1 - 3*b + b2)*(pu1t - 5*pu2t + pu3t) - 2*t3*(pu1t - 2*pu2t + pu3t) - t2*(-1 + b)*(3*pu1t - 8*pu2t + 3*pu3t))*ll2 - pu1t*lt + 3*theta*pu1t*lt - 2*t2*pu1t*lt - pu2t*lt + 3*theta*pu2t*lt - 2*t2*pu2t*lt + 7*pu1t*b*lt - 9*theta*pu1t*b*lt + 2*t2*pu1t*b*lt + 7*pu2t*b*lt - 9*theta*pu2t*b*lt + 2*t2*pu2t*b*lt - 6*pu1t*b2*lt + 3*theta*pu1t*b2*lt - 6*pu2t*b2*lt + 3*theta*pu2t*b2*lt + pu1t*b3*lt + pu2t*b3*lt - pu3t*lt + 3*theta*pu3t*lt - 2*t2*pu3t*lt + 7*b*pu3t*lt - 9*theta*b*pu3t*lt + 2*t2*b*pu3t*lt - 6*b2*pu3t*lt + 3*theta*b2*pu3t*lt + b3*pu3t*lt + t2*pu1t*ll3 - 3*t3*pu1t*ll3 + 2*t4*pu1t*ll3 + t2*pu2t*ll3 - 3*t3*pu2t*ll3 + 2*t4*pu2t*ll3 - 3*t2*pu1t*b*ll3 + 3*t3*pu1t*b*ll3 - 3*t2*pu2t*b*ll3 + 3*t3*pu2t*b*ll3 + t2*pu1t*b2*ll3 + t2*pu2t*b2*ll3 + theta*pu3t*ll3 - 5*t2*pu3t*ll3 + 8*t3*pu3t*ll3 - 4*t4*pu3t*ll3 - 7*theta*b*pu3t*ll3 + 15*t2*b*pu3t*ll3 - 8*t3*b*pu3t*ll3 + 6*theta*b2*pu3t*ll3 - 5*t2*b2*pu3t*ll3 - theta*b3*pu3t*ll3,2) + (1 + 2*t2 + 3*theta*(-1 + b) - 3*b + b2)*(4*t4*(pu1t^2) + 4*t4*(pu2t^2) + 8*t4*pl12t + 4*t4*(pu3t^2) + 8*t4*pl13t + 8*t4*pl23t + t2*((pu1t^2)*(1 - 15*b + 25*b2 - 10*b3 + b4) + theta*(-1 + 7*b - 6*b2 + b3)*(7*(pu1t^2) - 3*pl12t - 3*pl13t) + t2*(1 - 3*b + b2)*(18*(pu1t^2) + (pu2t^2) - 16*pl12t + (pu3t^2) - 16*pl13t + 2*pl23t) + 2*t4*(4*(pu1t^2) + (pu2t^2) - 7*pl12t + (pu3t^2) - 7*pl13t + 2*pl23t) + t3*(-1 + b)*(20*(pu1t^2) + 3*(pu2t^2) - 27*pl12t + 3*(pu3t^2) - 27*pl13t + 6*pl23t))*(ll1^2) + t2*((pu2t^2)*(1 - 15*b + 25*b2 - 10*b3 + b4) + t3*(-1 + b)*(3*(pu1t^2) + 20*(pu2t^2) - 27*pl12t + 3*(pu3t^2) + 6*pl13t - 27*pl23t) + t2*(1 - 3*b + b2)*((pu1t^2) + 18*(pu2t^2) - 16*pl12t + (pu3t^2) + 2*pl13t - 16*pl23t) + 2*t4*((pu1t^2) + 4*(pu2t^2) - 7*pl12t + (pu3t^2) + 2*pl13t - 7*pl23t) + theta*(-1 + 7*b - 6*b2 + b3)*(7*(pu2t^2) - 3*pl12t - 3*pl23t))*(ll2^2) + 2*theta*(pu1t^2)*lt - 4*t3*(pu1t^2)*lt + 2*theta*(pu2t^2)*lt - 4*t3*(pu2t^2)*lt + 4*theta*pl12t*lt - 8*t3*pl12t*lt - 14*theta*(pu1t^2)*b*lt + 4*t3*(pu1t^2)*b*lt - 14*theta*(pu2t^2)*b*lt + 4*t3*(pu2t^2)*b*lt - 28*theta*pl12t*b*lt + 8*t3*pl12t*b*lt + 12*theta*(pu1t^2)*b2*lt + 12*theta*(pu2t^2)*b2*lt + 24*theta*pl12t*b2*lt - 2*theta*(pu1t^2)*b3*lt - 2*theta*(pu2t^2)*b3*lt - 4*theta*pl12t*b3*lt + 2*theta*(pu3t^2)*lt - 4*t3*(pu3t^2)*lt - 14*theta*b*(pu3t^2)*lt + 4*t3*b*(pu3t^2)*lt + 12*theta*b2*(pu3t^2)*lt - 2*theta*b3*(pu3t^2)*lt + 4*theta*pl13t*lt - 8*t3*pl13t*lt - 28*theta*b*pl13t*lt + 8*t3*b*pl13t*lt + 24*theta*b2*pl13t*lt - 4*theta*b3*pl13t*lt + 4*theta*pl23t*lt - 8*t3*pl23t*lt - 28*theta*b*pl23t*lt + 8*t3*b*pl23t*lt + 24*theta*b2*pl23t*lt - 4*theta*b3*pl23t*lt + (pu1t^2)*lpu123t - 3*theta*(pu1t^2)*lpu123t + 2*t2*(pu1t^2)*lpu123t + (pu2t^2)*lpu123t - 3*theta*(pu2t^2)*lpu123t + 2*t2*(pu2t^2)*lpu123t + 2*pl12t*lpu123t - 6*theta*pl12t*lpu123t + 4*t2*pl12t*lpu123t - 15*(pu1t^2)*b*lpu123t + 21*theta*(pu1t^2)*b*lpu123t - 6*t2*(pu1t^2)*b*lpu123t - 15*(pu2t^2)*b*lpu123t + 21*theta*(pu2t^2)*b*lpu123t - 6*t2*(pu2t^2)*b*lpu123t - 30*pl12t*b*lpu123t + 42*theta*pl12t*b*lpu123t - 12*t2*pl12t*b*lpu123t + 25*(pu1t^2)*b2*lpu123t - 18*theta*(pu1t^2)*b2*lpu123t + 2*t2*(pu1t^2)*b2*lpu123t + 25*(pu2t^2)*b2*lpu123t - 18*theta*(pu2t^2)*b2*lpu123t + 2*t2*(pu2t^2)*b2*lpu123t + 50*pl12t*b2*lpu123t - 36*theta*pl12t*b2*lpu123t + 4*t2*pl12t*b2*lpu123t - 10*(pu1t^2)*b3*lpu123t + 3*theta*(pu1t^2)*b3*lpu123t - 10*(pu2t^2)*b3*lpu123t + 3*theta*(pu2t^2)*b3*lpu123t - 20*pl12t*b3*lpu123t + 6*theta*pl12t*b3*lpu123t + (pu1t^2)*b4*lpu123t + (pu2t^2)*b4*lpu123t + 2*pl12t*b4*lpu123t + (pu3t^2)*lpu123t - 3*theta*(pu3t^2)*lpu123t + 2*t2*(pu3t^2)*lpu123t - 15*b*(pu3t^2)*lpu123t + 21*theta*b*(pu3t^2)*lpu123t - 6*t2*b*(pu3t^2)*lpu123t + 25*b2*(pu3t^2)*lpu123t - 18*theta*b2*(pu3t^2)*lpu123t + 2*t2*b2*(pu3t^2)*lpu123t - 10*b3*(pu3t^2)*lpu123t + 3*theta*b3*(pu3t^2)*lpu123t + b4*(pu3t^2)*lpu123t + 2*pl13t*lpu123t - 6*theta*pl13t*lpu123t + 4*t2*pl13t*lpu123t - 30*b*pl13t*lpu123t + 42*theta*b*pl13t*lpu123t - 12*t2*b*pl13t*lpu123t + 50*b2*pl13t*lpu123t - 36*theta*b2*pl13t*lpu123t + 4*t2*b2*pl13t*lpu123t - 20*b3*pl13t*lpu123t + 6*theta*b3*pl13t*lpu123t + 2*b4*pl13t*lpu123t + 2*pl23t*lpu123t - 6*theta*pl23t*lpu123t + 4*t2*pl23t*lpu123t - 30*b*pl23t*lpu123t + 42*theta*b*pl23t*lpu123t - 12*t2*b*pl23t*lpu123t + 50*b2*pl23t*lpu123t - 36*theta*b2*pl23t*lpu123t + 4*t2*b2*pl23t*lpu123t - 20*b3*pl23t*lpu123t + 6*theta*b3*pl23t*lpu123t + 2*b4*pl23t*lpu123t - 6*t4*(pu1t^2)*ll3 + 8*t5*(pu1t^2)*ll3 - 6*t4*(pu2t^2)*ll3 + 8*t5*(pu2t^2)*ll3 - 12*t4*pl12t*ll3 + 16*t5*pl12t*ll3 + 6*t4*(pu1t^2)*b*ll3 + 6*t4*(pu2t^2)*b*ll3 + 12*t4*pl12t*b*ll3 - 2*t2*(pu3t^2)*ll3 + 16*t4*(pu3t^2)*ll3 - 16*t5*(pu3t^2)*ll3 + 14*t2*b*(pu3t^2)*ll3 - 16*t4*b*(pu3t^2)*ll3 - 12*t2*b2*(pu3t^2)*ll3 + 2*t2*b3*(pu3t^2)*ll3 - 2*t2*pl13t*ll3 + 10*t4*pl13t*ll3 - 8*t5*pl13t*ll3 + 14*t2*b*pl13t*ll3 - 10*t4*b*pl13t*ll3 - 12*t2*b2*pl13t*ll3 + 2*t2*b3*pl13t*ll3 - 2*t2*pl23t*ll3 + 10*t4*pl23t*ll3 - 8*t5*pl23t*ll3 + 14*t2*b*pl23t*ll3 - 10*t4*b*pl23t*ll3 - 12*t2*b2*pl23t*ll3 + 2*t2*b3*pl23t*ll3 - 2*t2*(pu1t^2)*lt*ll3 + 6*t3*(pu1t^2)*lt*ll3 - 4*t4*(pu1t^2)*lt*ll3 - 2*t2*(pu2t^2)*lt*ll3 + 6*t3*(pu2t^2)*lt*ll3 - 4*t4*(pu2t^2)*lt*ll3 - 4*t2*pl12t*lt*ll3 + 12*t3*pl12t*lt*ll3 - 8*t4*pl12t*lt*ll3 + 14*t2*(pu1t^2)*b*lt*ll3 - 18*t3*(pu1t^2)*b*lt*ll3 + 4*t4*(pu1t^2)*b*lt*ll3 + 14*t2*(pu2t^2)*b*lt*ll3 - 18*t3*(pu2t^2)*b*lt*ll3 + 4*t4*(pu2t^2)*b*lt*ll3 + 28*t2*pl12t*b*lt*ll3 - 36*t3*pl12t*b*lt*ll3 + 8*t4*pl12t*b*lt*ll3 - 12*t2*(pu1t^2)*b2*lt*ll3 + 6*t3*(pu1t^2)*b2*lt*ll3 - 12*t2*(pu2t^2)*b2*lt*ll3 + 6*t3*(pu2t^2)*b2*lt*ll3 - 24*t2*pl12t*b2*lt*ll3 + 12*t3*pl12t*b2*lt*ll3 + 2*t2*(pu1t^2)*b3*lt*ll3 + 2*t2*(pu2t^2)*b3*lt*ll3 + 4*t2*pl12t*b3*lt*ll3 - 2*theta*(pu3t^2)*lt*ll3 + 10*t2*(pu3t^2)*lt*ll3 - 16*t3*(pu3t^2)*lt*ll3 + 8*t4*(pu3t^2)*lt*ll3 + 30*theta*b*(pu3t^2)*lt*ll3 - 70*t2*b*(pu3t^2)*lt*ll3 + 48*t3*b*(pu3t^2)*lt*ll3 - 8*t4*b*(pu3t^2)*lt*ll3 - 50*theta*b2*(pu3t^2)*lt*ll3 + 60*t2*b2*(pu3t^2)*lt*ll3 - 16*t3*b2*(pu3t^2)*lt*ll3 + 20*theta*b3*(pu3t^2)*lt*ll3 - 10*t2*b3*(pu3t^2)*lt*ll3 - 2*theta*b4*(pu3t^2)*lt*ll3 - 2*theta*pl13t*lt*ll3 + 8*t2*pl13t*lt*ll3 - 10*t3*pl13t*lt*ll3 + 4*t4*pl13t*lt*ll3 + 30*theta*b*pl13t*lt*ll3 - 56*t2*b*pl13t*lt*ll3 + 30*t3*b*pl13t*lt*ll3 - 4*t4*b*pl13t*lt*ll3 - 50*theta*b2*pl13t*lt*ll3 + 48*t2*b2*pl13t*lt*ll3 - 10*t3*b2*pl13t*lt*ll3 + 20*theta*b3*pl13t*lt*ll3 - 8*t2*b3*pl13t*lt*ll3 - 2*theta*b4*pl13t*lt*ll3 - 2*theta*pl23t*lt*ll3 + 8*t2*pl23t*lt*ll3 - 10*t3*pl23t*lt*ll3 + 4*t4*pl23t*lt*ll3 + 30*theta*b*pl23t*lt*ll3 - 56*t2*b*pl23t*lt*ll3 + 30*t3*b*pl23t*lt*ll3 - 4*t4*b*pl23t*lt*ll3 - 50*theta*b2*pl23t*lt*ll3 + 48*t2*b2*pl23t*lt*ll3 - 10*t3*b2*pl23t*lt*ll3 + 20*theta*b3*pl23t*lt*ll3 - 8*t2*b3*pl23t*lt*ll3 - 2*theta*b4*pl23t*lt*ll3 + t4*(pu1t^2)*ll32 - 3*t5*(pu1t^2)*ll32 + 2*t6*(pu1t^2)*ll32 + t4*(pu2t^2)*ll32 - 3*t5*(pu2t^2)*ll32 + 2*t6*(pu2t^2)*ll32 + 2*t4*pl12t*ll32 - 6*t5*pl12t*ll32 + 4*t6*pl12t*ll32 - 3*t4*(pu1t^2)*a + 3*t5*(pu1t^2)*a - 3*t4*(pu2t^2)*a + 3*t5*(pu2t^2)*a - 6*t4*pl12t*a + 6*t5*pl12t*a + t4*(pu1t^2)*b2*ll32 + t4*(pu2t^2)*b2*ll32 + 2*t4*pl12t*b2*ll32 + t2*(pu3t^2)*ll32 - 7*t3*(pu3t^2)*ll32 + 18*t4*(pu3t^2)*ll32 - 20*t5*(pu3t^2)*ll32 + 8*t6*(pu3t^2)*ll32 - 15*t2*b*(pu3t^2)*ll32 + 49*t3*b*(pu3t^2)*ll32 - 54*t4*b*(pu3t^2)*ll32 + 20*t5*b*(pu3t^2)*ll32 + 25*t2*b2*(pu3t^2)*ll32 - 42*t3*b2*(pu3t^2)*ll32 + 18*t4*b2*(pu3t^2)*ll32 - 10*t2*b3*(pu3t^2)*ll32 + 7*t3*b3*(pu3t^2)*ll32 + t2*b4*(pu3t^2)*ll32 + 3*t3*pl13t*ll32 - 16*t4*pl13t*ll32 + 27*t5*pl13t*ll32 - 14*t6*pl13t*ll32 - 21*t3*b*pl13t*ll32 + 48*t4*b*pl13t*ll32 - 27*t5*b*pl13t*ll32 + 18*t3*b2*pl13t*ll32 - 16*t4*b2*pl13t*ll32 - 3*t3*b3*pl13t*ll32 + 3*t3*pl23t*ll32 - 16*t4*pl23t*ll32 + 27*t5*pl23t*ll32 - 14*t6*pl23t*ll32 - 21*t3*b*pl23t*ll32 + 48*t4*b*pl23t*ll32 - 27*t5*b*pl23t*ll32 + 18*t3*b2*pl23t*ll32 - 16*t4*b2*pl23t*ll32 - 3*t3*b3*pl23t*ll32 - 2*theta*ll1*(-(theta*(pl12t*(1 - 15*b + 25*b2 - 10*b3 + b4) - theta*(-1 + 7*b - 6*b2 + b3)*((pu1t^2) + (pu2t^2) - 8*pl12t + pl13t + pl23t)- 2*t4*(2*(pu1t^2) + 2*(pu2t^2) - 8*pl12t - (pu3t^2) + pl13t + pl23t) - t2*(1 - 3*b + b2)*(5*(pu1t^2) + 5*(pu2t^2) - 25*pl12t - (pu3t^2) + 4*pl13t + 4*pl23t) - t3*(-1 + b)*(8*(pu1t^2) + 8*(pu2t^2) - 34*pl12t - 3*(pu3t^2) + 5*pl13t + 5*pl23t))*ll2) + ((1 - 15*b + 25*b2 - 10*b3 + b4)*((pu1t^2) + pl12t + pl13t) + t2*(1 - 3*b + b2)*(8*(pu1t^2) - 3*(pu2t^2) + 5*pl12t - 3*(pu3t^2) + 5*pl13t - 6*pl23t) + 2*t3*(-1 + b)*(2*(pu1t^2) - (pu2t^2) + pl12t - (pu3t^2) + pl13t - 2*pl23t) + theta*(-1 + 7*b - 6*b2 + b3)*(5*(pu1t^2) - (pu2t^2) + 4*pl12t - (pu3t^2) + 4*pl13t - 2*pl23t))*lt + theta*((pu1t^2) - 8*t2*(pu1t^2) + 8*t3*(pu1t^2) + 3*t2*(pu2t^2) - 4*t3*(pu2t^2) + pl12t - 5*t2*pl12t + 4*t3*pl12t - 7*(pu1t^2)*b + 8*t2*(pu1t^2)*b - 3*t2*(pu2t^2)*b - 7*pl12t*b + 5*t2*pl12t*b + 6*(pu1t^2)*b2 + 6*pl12t*b2 - (pu1t^2)*b3 - pl12t*b3 + 3*t2*(pu3t^2) - 4*t3*(pu3t^2) - 3*t2*b*(pu3t^2) + pl13t - 5*t2*pl13t + 4*t3*pl13t - 7*b*pl13t + 5*t2*b*pl13t + 6*b2*pl13t - b3*pl13t + 6*t2*pl23t - 8*t3*pl23t - 6*t2*b*pl23t - ((1 - 15*b + 25*b2 - 10*b3 + b4)*pl13t - theta*(-1 + 7*b - 6*b2 + b3)*((pu1t^2) + pl12t + (pu3t^2) - 8*pl13t + pl23t)- 2*t4*(2*(pu1t^2) - (pu2t^2) + pl12t + 2*(pu3t^2) - 8*pl13t + pl23t) - t2*(1 - 3*b + b2)*(5*(pu1t^2) - (pu2t^2) + 4*pl12t + 5*(pu3t^2) - 25*pl13t + 4*pl23t) - t3*(-1 + b)*(8*(pu1t^2) - 3*(pu2t^2) + 5*pl12t + 8*(pu3t^2) - 34*pl13t + 5*pl23t))*ll3)) - 2*theta*ll2*((-(t2*(1 - 3*b + b2)*(3*(pu1t^2) - 8*(pu2t^2) - 5*pl12t + 3*(pu3t^2) + 6*pl13t - 5*pl23t)) - theta*(-1 + 7*b - 6*b2 + b3)*((pu1t^2) - 5*(pu2t^2) - 4*pl12t + (pu3t^2) + 2*pl13t - 4*pl23t) - 2*t3*(-1 + b)*((pu1t^2) - 2*(pu2t^2) - pl12t + (pu3t^2) + 2*pl13t - pl23t) + (1 - 15*b + 25*b2 - 10*b3 + b4)*((pu2t^2) + pl12t + pl23t))*lt + theta*(3*t2*(pu1t^2) - 4*t3*(pu1t^2) + (pu2t^2) - 8*t2*(pu2t^2) + 8*t3*(pu2t^2) + pl12t - 5*t2*pl12t + 4*t3*pl12t - 3*t2*(pu1t^2)*b - 7*(pu2t^2)*b + 8*t2*(pu2t^2)*b - 7*pl12t*b + 5*t2*pl12t*b + 6*(pu2t^2)*b2 + 6*pl12t*b2 - (pu2t^2)*b3 - pl12t*b3 + 3*t2*(pu3t^2) - 4*t3*(pu3t^2) - 3*t2*b*(pu3t^2) + 6*t2*pl13t - 8*t3*pl13t - 6*t2*b*pl13t + pl23t - 5*t2*pl23t + 4*t3*pl23t - 7*b*pl23t + 5*t2*b*pl23t + 6*b2*pl23t - b3*pl23t - ((1 - 15*b + 25*b2 - 10*b3 + b4)*pl23t - theta*(-1 + 7*b - 6*b2 + b3)*((pu2t^2) + pl12t + (pu3t^2) + pl13t - 8*pl23t)+ 2*t4*((pu1t^2) - 2*(pu2t^2) - pl12t - 2*(pu3t^2) - pl13t + 8*pl23t) + t2*(1 - 3*b + b2)*((pu1t^2) - 5*(pu2t^2) - 4*pl12t - 5*(pu3t^2) - 4*pl13t + 25*pl23t) + t3*(-1 + b)*(3*(pu1t^2) - 8*(pu2t^2) - 5*pl12t - 8*(pu3t^2) - 5*pl13t + 34*pl23t))*ll3))))/(t4*((1 + 2*t2 + 3*theta*(-1 + b) - 3*b + b2)^(2))*((pu123t)^(2)))
}

.gumb.V.123 = function(theta, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2];u3 = u[3]}else{u1 = u[,1];u2 = u[,2];u3 = u[,3]}
  
  l1 = log(u1)
  l2 = log(u2)
  l3 = log(u3)
  ll1 = log(-l1)
  ll2 = log(-l2)
  ll3 = log(-l3)
  ti = 1/theta
  
  pu1t = (-l1)^theta
  pu2t = (-l2)^theta
  pu3t = (-l3)^theta
  pu123t = pu1t + pu2t + pu3t
  lt = log(pu123t)
  b = .Power(pu123t,ti)
  b2 = b^2
  b3 = b^3
  
  t2 = theta^2
  t3 = theta^3
  t4 = theta^4
  t5 = theta^5
  
  (-3*t2*pu1t + 4*t3*pu1t - 3*t2*pu2t + 4*t3*pu2t + 3*t2*pu1t*b + 3*t2*pu2t*b - 3*t2*pu3t + 4*t3*pu3t + 3*t2*b*pu3t - theta*(pu1t*(-1 + 7*b - 6*b2 + b3) + t2*(-1 + b)*(8*pu1t - 3*pu2t - 3*pu3t) + t3*(4*pu1t - 2*pu2t - 2*pu3t) + theta*(1 - 3*b + b2)*(5*pu1t - pu2t - pu3t))*ll1 - theta*(pu2t*(-1 + 7*b - 6*b2 + b3) - theta*(1 - 3*b + b2)*(pu1t - 5*pu2t + pu3t) - 2*t3*(pu1t - 2*pu2t + pu3t) - t2*(-1 + b)*(3*pu1t - 8*pu2t + 3*pu3t))*ll2 - pu1t*lt + 3*theta*pu1t*lt - 2*t2*pu1t*lt - pu2t*lt + 3*theta*pu2t*lt - 2*t2*pu2t*lt + 7*pu1t*b*lt - 9*theta*pu1t*b*lt + 2*t2*pu1t*b*lt + 7*pu2t*b*lt - 9*theta*pu2t*b*lt + 2*t2*pu2t*b*lt - 6*pu1t*b2*lt + 3*theta*pu1t*b2*lt - 6*pu2t*b2*lt + 3*theta*pu2t*b2*lt + pu1t*b3*lt + pu2t*b3*lt - pu3t*lt + 3*theta*pu3t*lt - 2*t2*pu3t*lt + 7*b*pu3t*lt - 9*theta*b*pu3t*lt + 2*t2*b*pu3t*lt - 6*b2*pu3t*lt + 3*theta*b2*pu3t*lt + b3*pu3t*lt + t2*pu1t*ll3 - 3*t3*pu1t*ll3 + 2*t4*pu1t*ll3 + t2*pu2t*ll3 - 3*t3*pu2t*ll3 + 2*t4*pu2t*ll3 - 3*t2*pu1t*b*ll3 + 3*t3*pu1t*b*ll3 - 3*t2*pu2t*b*ll3 + 3*t3*pu2t*b*ll3 + t2*pu1t*b2*ll3 + t2*pu2t*b2*ll3 + theta*pu3t*ll3 - 5*t2*pu3t*ll3 + 8*t3*pu3t*ll3 - 4*t4*pu3t*ll3 - 7*theta*b*pu3t*ll3 + 15*t2*b*pu3t*ll3 - 8*t3*b*pu3t*ll3 + 6*theta*b2*pu3t*ll3 - 5*t2*b2*pu3t*ll3 - theta*b3*pu3t*ll3)/(t2*(1 + 2*t2 + 3*theta*(-1 + b) - 3*b + b2)*(pu123t))
}

.fran.S.123 = function(theta, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2];u3 = u[3]}else{u1 = u[,1];u2 = u[,2];u3 = u[,3]}
  et = exp(theta)
  et2 = exp(2*theta)
  t2 = theta^2
  u12 = u1 + u2
  u13 = u1 + u3
  u23 = u2 + u3
  et12 = exp(theta*(2 + u12))
  et13 = exp(theta*(2 + u13))
  et23 = exp(theta*(2 + u23))
  
  u123 = u12 + u3
  t2u1 = exp(theta*(2 + u1))
  t2u2 = exp(theta*(2 + u2))
  t2u3 = exp(theta*(2 + u3))
  e123 = exp(theta*(u123))
  e1123 = exp(theta*(1 + u123))
  e2123 = exp(theta*(2 + u123))
  a = -et2 + t2u1 + t2u2 - et12 + t2u3 - et13 - et23 + e123 - 2*e1123 + 2*e2123
  c1 =  et2 - t2u1 - t2u2 + et12 - t2u3 + et13 + et23 + e123 - 2*e1123
  b = c1^2
  
  ((-2*(-1 + et)*c1*(2*(-1 + et)*c1*a + 2*et*c1*a*theta + (-1 + et)*c1*a*theta*(2 + u123)- 3*(-1 + et)*a*theta*(2*et2 - t2u1*(2 + u1) - t2u2*(2 + u2) + et12*(2 + u12) - t2u3*(2 + u3) + et13*(2 + u13) + et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123)) + (-1 + et)*c1*theta*(-2*et2 + t2u1*(2 + u1) + t2u2*(2 + u2) - et12*(2 + u12) + t2u3*(2 + u3) - et13*(2 + u13) - et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123) + 2*e2123*(2 + u123))))/a - (2*et*c1*theta*(2*(-1 + et)*c1*a + 2*et*c1*a*theta + (-1 + et)*c1*a*theta*(2 + u123)- 3*(-1 + et)*a*theta*(2*et2 - t2u1*(2 + u1) - t2u2*(2 + u2) + et12*(2 + u12) - t2u3*(2 + u3) + et13*(2 + u13) + et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123)) + (-1 + et)*c1*theta*(-2*et2 + t2u1*(2 + u1) + t2u2*(2 + u2) - et12*(2 + u12) + t2u3*(2 + u3) - et13*(2 + u13) - et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123) + 2*e2123*(2 + u123))))/a + ((-1 + et)*c1*theta*(-2 - u1 - u2 - u3)*(2*(-1 + et)*c1*a + 2*et*c1*a*theta + (-1 + et)*c1*a*theta*(2 + u123)- 3*(-1 + et)*a*theta*(2*et2 - t2u1*(2 + u1) - t2u2*(2 + u2) + et12*(2 + u12) - t2u3*(2 + u3) + et13*(2 + u13) + et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123)) + (-1 + et)*c1*theta*(-2*et2 + t2u1*(2 + u1) + t2u2*(2 + u2) - et12*(2 + u12) + t2u3*(2 + u3) - et13*(2 + u13) - et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123) + 2*e2123*(2 + u123))))/a + (3*(-1 + et)*theta*(2*et2 - t2u1*(2 + u1) - t2u2*(2 + u2) + et12*(2 + u12) - t2u3*(2 + u3) + et13*(2 + u13) + et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123))*(2*(-1 + et)*c1*a + 2*et*c1*a*theta + (-1 + et)*c1*a*theta*(2 + u123)- 3*(-1 + et)*a*theta*(2*et2 - t2u1*(2 + u1) - t2u2*(2 + u2) + et12*(2 + u12) - t2u3*(2 + u3) + et13*(2 + u13) + et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123)) + (-1 + et)*c1*theta*(-2*et2 + t2u1*(2 + u1) + t2u2*(2 + u2) - et12*(2 + u12) + t2u3*(2 + u3) - et13*(2 + u13) - et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123) + 2*e2123*(2 + u123))))/a - ((-1 + et)*c1*theta*(-2*et2 + t2u1*(2 + u1) + t2u2*(2 + u2) - et12*(2 + u12) + t2u3*(2 + u3) - et13*(2 + u13) - et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123) + 2*e2123*(2 + u123))*(2*(-1 + et)*c1*a + 2*et*c1*a*theta + (-1 + et)*c1*a*theta*(2 + u123)- 3*(-1 + et)*a*theta*(2*et2 - t2u1*(2 + u1) - t2u2*(2 + u2) + et12*(2 + u12) - t2u3*(2 + u3) + et13*(2 + u13) + et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123)) + (-1 + et)*c1*theta*(-2*et2 + t2u1*(2 + u1) + t2u2*(2 + u2) - et12*(2 + u12) + t2u3*(2 + u3) - et13*(2 + u13) - et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123) + 2*e2123*(2 + u123))))/.Power(et2 - t2u1 - t2u2 + et12 - t2u3 + et13 + et23 - e123 + 2*e1123 - 2*e2123,2) + (2*.Power(-1 + et,2)*b*a + 8*et*(-1 + et)*b*a*theta + 2*et2*b*a*t2 + 4*.Power(-1 + et,2)*b*a*theta*(2 + u123) + 2*et*(-1 + et)*b*a*t2*(2 + u123) + .Power(-1 + et,2)*b*a*t2*.Power(2 + u123,2) + 2*et*(-1 + et)*b*a*t2*(3 + u123) - 12*.Power(-1 + et,2)*c1*a*theta*(2*et2 - t2u1*(2 + u1) - t2u2*(2 + u2) + et12*(2 + u12) - t2u3*(2 + u3) + et13*(2 + u13) + et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123)) - 12*et*(-1 + et)*c1*a*t2*(2*et2 - t2u1*(2 + u1) - t2u2*(2 + u2) + et12*(2 + u12) - t2u3*(2 + u3) + et13*(2 + u13) + et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123)) - 6*.Power(-1 + et,2)*c1*a*t2*(2 + u123)*(2*et2 - t2u1*(2 + u1) - t2u2*(2 + u2) + et12*(2 + u12) - t2u3*(2 + u3) + et13*(2 + u13) + et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123)) + 12*.Power(-1 + et,2)*a*t2*.Power(2*et2 - t2u1*(2 + u1) - t2u2*(2 + u2) + et12*(2 + u12) - t2u3*(2 + u3) + et13*(2 + u13) + et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123),2) - 3*.Power(-1 + et,2)*c1*a*t2*(4*et2 - t2u1*.Power(2 + u1,2) - t2u2*.Power(2 + u2,2) + et12*.Power(2 + u12,2) - t2u3*.Power(2 + u3,2) + et13*.Power(2 + u13,2) + et23*.Power(2 + u23,2) + e123*.Power(u123,2) - 2*e1123*.Power(1 + u123,2)) + 4*.Power(-1 + et,2)*b*theta*(-2*et2 + t2u1*(2 + u1) + t2u2*(2 + u2) - et12*(2 + u12) + t2u3*(2 + u3) - et13*(2 + u13) - et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123) + 2*e2123*(2 + u123)) + 4*et*(-1 + et)*b*t2*(-2*et2 + t2u1*(2 + u1) + t2u2*(2 + u2) - et12*(2 + u12) + t2u3*(2 + u3) - et13*(2 + u13) - et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123) + 2*e2123*(2 + u123)) + 2*.Power(-1 + et,2)*b*t2*(2 + u123)*(-2*et2 + t2u1*(2 + u1) + t2u2*(2 + u2) - et12*(2 + u12) + t2u3*(2 + u3) - et13*(2 + u13) - et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123) + 2*e2123*(2 + u123)) - 6*.Power(-1 + et,2)*c1*t2*(2*et2 - t2u1*(2 + u1) - t2u2*(2 + u2) + et12*(2 + u12) - t2u3*(2 + u3) + et13*(2 + u13) + et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123))*(-2*et2 + t2u1*(2 + u1) + t2u2*(2 + u2) - et12*(2 + u12) + t2u3*(2 + u3) - et13*(2 + u13) - et23*(2 + u23) + e123*(u123) - 2*e1123*(1 + u123) + 2*e2123*(2 + u123)) + .Power(-1 + et,2)*b*t2*(-4*et2 + t2u1*.Power(2 + u1,2) + t2u2*.Power(2 + u2,2) - et12*.Power(2 + u12,2) + t2u3*.Power(2 + u3,2) - et13*.Power(2 + u13,2) - et23*.Power(2 + u23,2) + e123*.Power(u123,2) - 2*e1123*.Power(1 + u123,2) + 2*e2123*.Power(2 + u123,2)))/a)/(.Power(-1 + et,2)*b*t2)
}

.fran.V.123 = function(theta, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2];u3 = u[3]}else{u1 = u[,1];u2 = u[,2];u3 = u[,3]}
  et = .Power(E,theta)
  et2 = .Power(E,2*theta)
  et12 = .Power(E,theta*(2 + u1 + u2))
  et13 = .Power(E,theta*(2 + u1 + u3))
  et23 = .Power(E,theta*(2 + u2 + u3))
  
  u123 = u1 + u2 + u3
  t2u1 = .Power(E,theta*(2 + u1))
  t2u2 = .Power(E,theta*(2 + u2))
  t2u3 = .Power(E,theta*(2 + u3))
  e123 = .Power(E,theta*(u123))
  e1123 = .Power(E,theta*(1 + u123))
  e2123 = .Power(E,theta*(2 + u123))
  a =  et2 - t2u1 - t2u2 + et12 - t2u3 + et13 + et23 + e123 - 2*e1123
  b = -et2 + t2u1 + t2u2 - et12 + t2u3 - et13 - et23 + e123 - 2*e1123 + 2*e2123
  
  (2*(-1 + et)*a*b + 2*et*a*b*theta + (-1 + et)*a*b*theta*(2 + u123) - 3*(-1 + et)*b*theta*(2*et2 - t2u1*(2 + u1) - t2u2*(2 + u2) + et12*(2 + u1 + u2) - t2u3*(2 + u3) + et13*(2 + u1 + u3) + et23*(2 + u2 + u3) + e123*u123 - 2*e1123*(1 + u123))+ (-1 + et)*a*theta*(-2*et2 + t2u1*(2 + u1) + t2u2*(2 + u2) - et12*(2 + u1 + u2) + t2u3*(2 + u3) - et13*(2 + u1 + u3) - et23*(2 + u2 + u3) + e123*u123 - 2*e1123*(1 + u123) + 2*e2123*(2 + u123)))/((-1 + et)*a*b*theta)
}

.S.ga.12.12 = function(u, sig){(-((-1 + (sig[1,2]^2) + (sig[1,3]^2) - 2*sig[1,2]*sig[1,3]*sig[2,3] + (sig[2,3]^2))*(-1 + 3*(sig[1,2]^2) - (u[2]^2)*(-1 + (sig[1,3]^2)) - 2*u[2]*u[3]*sig[2,3] + (sig[2,3]^2) + (u[3]^2)*(sig[2,3]^2) + 2*sig[1,2]*sig[1,3]*(u[2]*u[3] - (3 + (u[3]^2))*sig[2,3]) - (u[1]^2)*(-1 + (sig[2,3]^2)) + (sig[1,3]^2)*(1 + (u[3]^2) + 2*(sig[2,3]^2)) - 2*u[1]*(u[3]*(sig[1,3] - sig[1,2]*sig[2,3]) + u[2]*(sig[1,2] - sig[1,3]*sig[2,3])))) + 4*(sig[1,2] - sig[1,3]*sig[2,3])*((sig[1,2]^3) + (sig[1,2]^2)*sig[1,3]*(u[2]*u[3] - (3 + (u[3]^2))*sig[2,3]) - (u[1]^2)*(sig[1,2] - sig[1,3]*sig[2,3])*(-1 + (sig[2,3]^2)) + sig[1,3]*((u[2]^2)*(-1 + (sig[1,3]^2))*sig[2,3] + u[2]*u[3]*(1 - (sig[1,3]^2) + (sig[2,3]^2)) - sig[2,3]*(-1 + (u[3]^2) + (sig[1,3]^2) + (sig[2,3]^2))) + sig[1,2]*(-1 - (u[2]^2)*(-1 + (sig[1,3]^2)) - 2*u[2]*u[3]*sig[2,3] + (sig[2,3]^2) + (u[3]^2)*(sig[2,3]^2) + (sig[1,3]^2)*(1 + (u[3]^2) + 2*(sig[2,3]^2))) + u[1]*(u[2]*(-1 - (sig[1,2]^2) + (sig[1,3]^2) + 2*sig[1,2]*sig[1,3]*sig[2,3] + (sig[2,3]^2) - 2*(sig[1,3]^2)*(sig[2,3]^2)) + u[3]*(-2*sig[1,2]*sig[1,3] + sig[2,3] + (sig[1,2]^2)*sig[2,3] + (sig[1,3]^2)*sig[2,3] - (sig[2,3]^3)))))/.Power(-1 + (sig[1,2]^2) + (sig[1,3]^2) - 2*sig[1,2]*sig[1,3]*sig[2,3] + (sig[2,3]^2),3)}
.S.ga.12.13 = function(u, sig){(-((-1 + (sig[1,2]^2) + (sig[1,3]^2) - 2*sig[1,2]*sig[1,3]*sig[2,3] + (sig[2,3]^2))*(2*sig[1,2]*sig[1,3] + 2*(u[3]^2)*sig[1,2]*sig[1,3] + sig[2,3] - (u[3]^2)*sig[2,3] - 3*(sig[1,2]^2)*sig[2,3] - (u[3]^2)*(sig[1,2]^2)*sig[2,3] - 3*(sig[1,3]^2)*sig[2,3] + 4*sig[1,2]*sig[1,3]*(sig[2,3]^2) - (sig[2,3]^3) - 2*u[1]*u[3]*(sig[1,2] - sig[1,3]*sig[2,3]) - (u[2]^2)*(2*sig[1,2]*sig[1,3] + sig[2,3] - 3*(sig[1,3]^2)*sig[2,3]) + (u[1]^2)*sig[2,3]*(-1 + (sig[2,3]^2)) + u[2]*(u[3]*(1 + (sig[1,2]^2) - 3*(sig[1,3]^2) + (sig[2,3]^2)) + 2*u[1]*(sig[1,3] + sig[1,2]*sig[2,3] - 2*sig[1,3]*(sig[2,3]^2))))) + 4*(sig[1,3] - sig[1,2]*sig[2,3])*((sig[1,2]^3) + (sig[1,2]^2)*sig[1,3]*(u[2]*u[3] - (3 + (u[3]^2))*sig[2,3]) - (u[1]^2)*(sig[1,2] - sig[1,3]*sig[2,3])*(-1 + (sig[2,3]^2)) + sig[1,3]*((u[2]^2)*(-1 + (sig[1,3]^2))*sig[2,3] + u[2]*u[3]*(1 - (sig[1,3]^2) + (sig[2,3]^2)) - sig[2,3]*(-1 + (u[3]^2) + (sig[1,3]^2) + (sig[2,3]^2))) + sig[1,2]*(-1 - (u[2]^2)*(-1 + (sig[1,3]^2)) - 2*u[2]*u[3]*sig[2,3] + (sig[2,3]^2) + (u[3]^2)*(sig[2,3]^2) + (sig[1,3]^2)*(1 + (u[3]^2) + 2*(sig[2,3]^2))) + u[1]*(u[2]*(-1 - (sig[1,2]^2) + (sig[1,3]^2) + 2*sig[1,2]*sig[1,3]*sig[2,3] + (sig[2,3]^2) - 2*(sig[1,3]^2)*(sig[2,3]^2)) + u[3]*(-2*sig[1,2]*sig[1,3] + sig[2,3] + (sig[1,2]^2)*sig[2,3] + (sig[1,3]^2)*sig[2,3] - (sig[2,3]^3)))))/.Power(-1 + (sig[1,2]^2) + (sig[1,3]^2) - 2*sig[1,2]*sig[1,3]*sig[2,3] + (sig[2,3]^2),3)}
.S.ga.12.23 = function(u, sig){(-((-1 + (sig[1,2]^2) + (sig[1,3]^2) - 2*sig[1,2]*sig[1,3]*sig[2,3] + (sig[2,3]^2))*(-((3 + (u[3]^2))*(sig[1,2]^2)*sig[1,3]) + 2*(u[1]^2)*sig[2,3]*(-sig[1,2] + sig[1,3]*sig[2,3]) + 2*sig[1,2]*(-(u[2]*u[3]) + (1 + (u[3]^2) + 2*(sig[1,3]^2))*sig[2,3]) + sig[1,3]*(1 - (u[3]^2) - (sig[1,3]^2) + (u[2]^2)*(-1 + (sig[1,3]^2)) + 2*u[2]*u[3]*sig[2,3] - 3*(sig[2,3]^2)) + (u[1]^2)*sig[1,3]*(-1 + (sig[2,3]^2)) + u[1]*(2*u[2]*(sig[1,2]*sig[1,3] + sig[2,3] - 2*(sig[1,3]^2)*sig[2,3]) + u[3]*(1 + (sig[1,2]^2) + (sig[1,3]^2) - 3*(sig[2,3]^2))))) + 2*(-2*sig[1,2]*sig[1,3] + 2*sig[2,3])*((sig[1,2]^3) + (sig[1,2]^2)*sig[1,3]*(u[2]*u[3] - (3 + (u[3]^2))*sig[2,3]) - (u[1]^2)*(sig[1,2] - sig[1,3]*sig[2,3])*(-1 + (sig[2,3]^2)) + sig[1,3]*((u[2]^2)*(-1 + (sig[1,3]^2))*sig[2,3] + u[2]*u[3]*(1 - (sig[1,3]^2) + (sig[2,3]^2)) - sig[2,3]*(-1 + (u[3]^2) + (sig[1,3]^2) + (sig[2,3]^2))) + sig[1,2]*(-1 - (u[2]^2)*(-1 + (sig[1,3]^2)) - 2*u[2]*u[3]*sig[2,3] + (sig[2,3]^2) + (u[3]^2)*(sig[2,3]^2) + (sig[1,3]^2)*(1 + (u[3]^2) + 2*(sig[2,3]^2))) + u[1]*(u[2]*(-1 - (sig[1,2]^2) + (sig[1,3]^2) + 2*sig[1,2]*sig[1,3]*sig[2,3] + (sig[2,3]^2) - 2*(sig[1,3]^2)*(sig[2,3]^2)) + u[3]*(-2*sig[1,2]*sig[1,3] + sig[2,3] + (sig[1,2]^2)*sig[2,3] + (sig[1,3]^2)*sig[2,3] - (sig[2,3]^3)))))/.Power(-1 + (sig[1,2]^2) + (sig[1,3]^2) - 2*sig[1,2]*sig[1,3]*sig[2,3] + (sig[2,3]^2),3)}
.S.ga.13.13 = function(u, sig){-(((-1 + (sig[1,2]^2) + (sig[1,3]^2) - 2*sig[1,2]*sig[1,3]*sig[2,3] + (sig[2,3]^2))*(-1 + (sig[1,2]^2) + (u[2]^2)*(sig[1,2]^2) - (u[3]^2)*(-1 + (sig[1,2]^2)) + 3*(sig[1,3]^2) + 2*u[2]*u[3]*(sig[1,2]*sig[1,3] - sig[2,3]) - 6*sig[1,2]*sig[1,3]*sig[2,3] - 2*(u[2]^2)*sig[1,2]*sig[1,3]*sig[2,3] + (sig[2,3]^2) + (u[2]^2)*(sig[2,3]^2) + 2*(sig[1,2]^2)*(sig[2,3]^2) - (u[1]^2)*(-1 + (sig[2,3]^2)) - 2*u[1]*(u[3]*(sig[1,3] - sig[1,2]*sig[2,3]) + u[2]*(sig[1,2] - sig[1,3]*sig[2,3]))) + 4*(sig[1,3] - sig[1,2]*sig[2,3])*((u[1]^2)*(sig[1,3] - sig[1,2]*sig[2,3])*(-1 + (sig[2,3]^2)) - (sig[1,3] - sig[1,2]*sig[2,3])*(-1 + (sig[1,2]^2) - (u[3]^2)*(-1 + (sig[1,2]^2)) + (sig[1,3]^2) - 2*sig[1,2]*sig[1,3]*sig[2,3] + (sig[2,3]^2)) - (u[2]^2)*((sig[1,2]^2)*sig[1,3] - sig[1,2]*(1 + (sig[1,3]^2))*sig[2,3] + sig[1,3]*(sig[2,3]^2)) - u[2]*u[3]*(-(sig[1,2]^3) - 2*sig[1,3]*sig[2,3] + sig[1,2]*(1 + (sig[1,3]^2) + (sig[2,3]^2))) + u[1]*(-(u[2]*(-2*sig[1,2]*sig[1,3] + sig[2,3] + (sig[1,2]^2)*sig[2,3] + (sig[1,3]^2)*sig[2,3] - (sig[2,3]^3))) + u[3]*(1 + (sig[1,3]^2) - 2*sig[1,2]*sig[1,3]*sig[2,3] - (sig[2,3]^2) + (sig[1,2]^2)*(-1 + 2*(sig[2,3]^2))))))/.Power(-1 + (sig[1,2]^2) + (sig[1,3]^2) - 2*sig[1,2]*sig[1,3]*sig[2,3] + (sig[2,3]^2),3))}
.S.ga.13.23 = function(u, sig){((-1 + (sig[1,2]^2) + (sig[1,3]^2) - 2*sig[1,2]*sig[1,3]*sig[2,3] + (sig[2,3]^2))*(-((-1 + (u[3]^2))*(sig[1,2]^3)) - 4*(sig[1,2]^2)*sig[1,3]*sig[2,3] - 2*sig[1,3]*(-(u[2]*u[3]) + sig[2,3] + (u[2]^2)*sig[2,3]) + sig[1,2]*(-1 + (u[3]^2) + 3*(sig[1,3]^2) + (u[2]^2)*(1 + (sig[1,3]^2)) - 2*u[2]*u[3]*sig[2,3] + 3*(sig[2,3]^2)) + (u[1]^2)*(sig[1,2] + 2*sig[1,3]*sig[2,3] - 3*sig[1,2]*(sig[2,3]^2)) - u[1]*(2*u[3]*(sig[1,2]*sig[1,3] + sig[2,3] - 2*(sig[1,2]^2)*sig[2,3]) + u[2]*(1 + (sig[1,2]^2) + (sig[1,3]^2) - 3*(sig[2,3]^2)))) - 2*(-2*sig[1,2]*sig[1,3] + 2*sig[2,3])*((u[1]^2)*(sig[1,3] - sig[1,2]*sig[2,3])*(-1 + (sig[2,3]^2)) - (sig[1,3] - sig[1,2]*sig[2,3])*(-1 + (sig[1,2]^2) - (u[3]^2)*(-1 + (sig[1,2]^2)) + (sig[1,3]^2) - 2*sig[1,2]*sig[1,3]*sig[2,3] + (sig[2,3]^2)) - (u[2]^2)*((sig[1,2]^2)*sig[1,3] - sig[1,2]*(1 + (sig[1,3]^2))*sig[2,3] + sig[1,3]*(sig[2,3]^2)) - u[2]*u[3]*(-(sig[1,2]^3) - 2*sig[1,3]*sig[2,3] + sig[1,2]*(1 + (sig[1,3]^2) + (sig[2,3]^2))) + u[1]*(-(u[2]*(-2*sig[1,2]*sig[1,3] + sig[2,3] + (sig[1,2]^2)*sig[2,3] + (sig[1,3]^2)*sig[2,3] - (sig[2,3]^3))) + u[3]*(1 + (sig[1,3]^2) - 2*sig[1,2]*sig[1,3]*sig[2,3] - (sig[2,3]^2) + (sig[1,2]^2)*(-1 + 2*(sig[2,3]^2))))))/.Power(-1 + (sig[1,2]^2) + (sig[1,3]^2) - 2*sig[1,2]*sig[1,3]*sig[2,3] + (sig[2,3]^2),3)}
.S.ga.23.23 = function(u, sig){(-((-1 + (sig[1,2]^2) + (sig[1,3]^2) - 2*sig[1,2]*sig[1,3]*sig[2,3] + (sig[2,3]^2))*(-1 + (sig[1,2]^2) + (u[1]^2)*(sig[1,2]^2) - (u[3]^2)*(-1 + (sig[1,2]^2)) + (sig[1,3]^2) + (u[1]^2)*(sig[1,3]^2) + 2*(sig[1,2]^2)*(sig[1,3]^2) - (u[2]^2)*(-1 + (sig[1,3]^2)) - 6*sig[1,2]*sig[1,3]*sig[2,3] - 2*(u[1]^2)*sig[1,2]*sig[1,3]*sig[2,3] + 3*(sig[2,3]^2) - 2*u[1]*u[3]*(sig[1,3] - sig[1,2]*sig[2,3]) + 2*u[2]*(-(u[1]*sig[1,2]) + u[3]*sig[1,2]*sig[1,3] - u[3]*sig[2,3] + u[1]*sig[1,3]*sig[2,3]))) + 2*(-2*sig[1,2]*sig[1,3] + 2*sig[2,3])*((u[2]^2)*(-1 + (sig[1,3]^2))*(sig[1,2]*sig[1,3] - sig[2,3]) + (sig[1,2]*sig[1,3] - sig[2,3])*(1 - (sig[1,2]^2) + (u[3]^2)*(-1 + (sig[1,2]^2)) - (sig[1,3]^2) + 2*sig[1,2]*sig[1,3]*sig[2,3] - (sig[2,3]^2)) + (u[1]^2)*((sig[1,2]^2)*sig[2,3] + (sig[1,3]^2)*sig[2,3] - sig[1,2]*sig[1,3]*(1 + (sig[2,3]^2))) + u[1]*u[3]*(-(sig[1,2]^3) - 2*sig[1,3]*sig[2,3] + sig[1,2]*(1 + (sig[1,3]^2) + (sig[2,3]^2))) + u[2]*(u[3]*(-1 + (sig[1,3]^2) + (sig[1,2]^2)*(1 - 2*(sig[1,3]^2)) + 2*sig[1,2]*sig[1,3]*sig[2,3] - (sig[2,3]^2)) + u[1]*(-(sig[1,3]^3) - 2*sig[1,2]*sig[2,3] + sig[1,3]*(1 + (sig[1,2]^2) + (sig[2,3]^2))))))/.Power(-1 + (sig[1,2]^2) + (sig[1,3]^2) - 2*sig[1,2]*sig[1,3]*sig[2,3] + (sig[2,3]^2),3)}

