
rm(list = ls())

a = 0.1786
b = 0.85
c = 0.47
d = 1

T = a + d
D = a*d - b*c


L1 = T/2 + ((T^2/4) - D)^(0.5)
L2 = T/2 - ((T^2/4) - D)^(0.5)

lam_e = 0.2506842


se1 = lam_e/L1
se2 = lam_e/L2

me = 1 - se1
