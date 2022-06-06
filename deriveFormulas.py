import sympy as smp
import dill
dill.settings['recurse'] = True

t, g = smp.symbols("t g")
m1, m2 = smp.symbols("m1 m2")
l1, l2 =smp.symbols("L1, L2")

theta1, theta2 = smp.symbols(r"\theta_1 \theta_2", cls=smp.Function)

theta1 = theta1(t)
theta2 = theta2(t)

theta1_d = smp.diff(theta1,t)
theta1_dd = smp.diff(theta1_d,t)
theta2_d = smp.diff(theta2,t)
theta2_dd = smp.diff(theta2_d,t)

x1 = smp.sin(theta1)*l1
y1 = -1*smp.cos(theta1)*l1
x2 = smp.sin(theta1)*l1 + smp.sin(theta2)*l2
y2 = -1*smp.cos(theta1)*l1-smp.cos(theta2)*l2

ekin_1 = 1/2 * m1 * (smp.diff(x1,t)**2+ smp.diff(y1,t)**2 )
ekin_2 = 1/2 * m2 * (smp.diff(x2,t)**2+ smp.diff(y2,t)**2 )
ekin = ekin_1 + ekin_2

epot_1 = m1*g*y1
epot_2 = m2*g*y2
epot = epot_1 + epot_2
L = ekin - epot

LG1 = smp.diff(L,theta1)-smp.diff(smp.diff(L,theta1_d),t).simplify()
LG2 = smp.diff(L,theta2)-smp.diff(smp.diff(L,theta2_d),t).simplify()

sols = smp.solve([LG1,LG2],(theta1_dd,theta2_dd),simplify=False,rational=False) #DGL 2. Ordnung

dz1dt_f = smp.lambdify((t,g,m1,m2,l1,l2,theta1,theta2,theta1_d,theta2_d),sols[theta1_dd])
dz2dt_f = smp.lambdify((t,g,m1,m2,l1,l2,theta1,theta2,theta1_d,theta2_d),sols[theta2_dd])
dtheta1dt_f = smp.lambdify(theta1_d,theta1_d)
dtheta2dt_f = smp.lambdify(theta2_d,theta2_d)

formulas = [dz1dt_f,dz2dt_f,dtheta1dt_f,dtheta2dt_f]
dill.dump(formulas, open("derivedFormulas", "wb"))