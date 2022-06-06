import dearpygui.dearpygui as dpg
import threading
import time
import numpy as np
import sympy as smp
from scipy.integrate import odeint
from math import sin, cos


t, g = smp.symbols("t g")
m1, m2 = smp.symbols("m1 m2")
l1, l2 =smp.symbols("L1, L2")

theta1, theta2 = smp.symbols("theta_1 theta_2", cls=smp.Function)

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

def dSdt(S,t,g,m1,m2,l1,l2):
    theta1, z1, theta2, z2 = S
    return [dtheta1dt_f(z1), 
            dz1dt_f(t,g,m1,m2,l1,l2,theta1,theta2,z1,z2),
            dtheta2dt_f(z2),
            dz2dt_f(t,g,m1,m2,l1,l2,theta1,theta2,z1,z2)
           ]

currTime = 15
t = np.linspace(0,15,500)
g = 9.81
m1 = 2
m2 = 2
l1 = 1.5
l2 = 1.5
theta1_0 = 0
theta2_0 = 0

ans = []
theta1 = []
theta2 = []

def get_pos(theta1,theta2,l1,l2):
    l1 *= 100
    l2 *= 100
    return(l1*sin(theta1),
            l1*cos(theta1),
            l1*sin(theta1)+l2*sin(theta2),
            l1*cos(theta1)+l2*cos(theta2))

def start():
    global running,theta1, theta2, t, ans
    t = np.linspace(0,15,500)
    ans = odeint(dSdt,y0=[theta1_0,0,theta2_0,0],t=t,args=(g,m1,m2,l1,l2))
    print(1)
    theta1 = ans.T[0].tolist()
    theta2 = ans.T[2].tolist()
    dpg.configure_item(startBtn,enabled=False)
    dpg.configure_item(the1Slider,enabled=False)
    dpg.configure_item(the2Slider,enabled=False)
    dpg.configure_item(m1Slider,enabled=False)
    dpg.configure_item(m2Slider,enabled=False)
    dpg.configure_item(l1Slider,enabled=False)
    dpg.configure_item(l2Slider,enabled=False)
    running = True

def reset():
    global running, plot1Data, plot2Data
    running = False
    dpg.configure_item(startBtn,enabled=True)
    dpg.configure_item(the1Slider,enabled=True)
    dpg.configure_item(the2Slider,enabled=True)
    dpg.configure_item(m1Slider,enabled=True)
    dpg.configure_item(m2Slider,enabled=True)
    dpg.configure_item(l1Slider,enabled=True)
    dpg.configure_item(l2Slider,enabled=True)
    plot1Data = [0 for i in range(300)]
    plot2Data = [0 for i in range(300)]
    update()

def update():
    global theta1_0, theta2_0, l1, l2, m1, m2
    theta1_0 = dpg.get_value(the1Slider)
    theta2_0 = dpg.get_value(the2Slider)
    l1 = dpg.get_value(l1Slider)
    l2 = dpg.get_value(l2Slider)
    m1 = dpg.get_value(m1Slider)
    m2 = dpg.get_value(m2Slider)

    pos = get_pos(theta1_0,theta2_0,l1,l2)
    dpg.configure_item(pendulum1,p1=[350,350],p2=[350+pos[0],350+pos[1]])
    dpg.configure_item(pendulum2,p1=[350+pos[0],350+pos[1]],p2=[350+pos[2],350+pos[3]])
    dpg.configure_item(circle1,center=[350+pos[0],350+pos[1]],radius=m1*7)
    dpg.configure_item(circle2,center=[350+pos[2],350+pos[3]],radius=m2*7)

dpg.create_context()
header1 = dpg.add_window(label="header",width=350,height=105,no_close=True,no_title_bar=True,no_resize=True,no_move=True)
header2 = dpg.add_window(label="header",width=350,height=105,pos=(350,0),no_close=True,no_title_bar=True,no_resize=True,no_move=True)
body = dpg.add_window(label="main",width=700,height=700,pos=(0,105),no_close=True,no_title_bar=True,no_resize=True,no_move=True)
sidebar = dpg.add_window(label="sidebar",width=300,height=820,pos=(700,0),no_close=True,no_title_bar=True,no_resize=True,no_move=True)
plot1 = dpg.add_simple_plot(overlay="Angle1",parent=sidebar,width=285,height=393)
plot2 = dpg.add_simple_plot(overlay="Angle2",parent=sidebar,width=285,height=393)

running = False
startBtn = dpg.add_button(label="Start",parent=header1,callback=start)
resetBtn = dpg.add_button(label="Reset",parent=header1,callback=reset)

the1Slider = dpg.add_slider_float(label="Angle1",min_value=0,max_value=2*3.1415926,default_value=0,parent=header1,callback=update)
the2Slider = dpg.add_slider_float(label="Angle2",min_value=0,max_value=2*3.1415926,default_value=0,parent=header1,callback=update)
m1Slider = dpg.add_slider_float(label="Mass1",min_value=0.3,max_value=10,default_value=2,parent=header2,callback=update)
m2Slider = dpg.add_slider_float(label="Mass2",min_value=0.3,max_value=10,default_value=2,parent=header2,callback=update)
l1Slider = dpg.add_slider_float(label="Length1",min_value=0.5,max_value=3,default_value=1.5,parent=header2,callback=update)
l2Slider = dpg.add_slider_float(label="Length2",min_value=0.5,max_value=3,default_value=1.5,parent=header2,callback=update)

pos = get_pos(theta1_0,theta2_0,l1,l2)

pendulum1 = dpg.draw_line([350,350],[350+pos[0],350+pos[1]],parent=body)
pendulum2 = dpg.draw_line([350+pos[0],350+pos[1]],[350+pos[2],350+pos[3]],parent=body)
circle1 = dpg.draw_circle([350+pos[0],350+pos[1]],m1*7,fill=(255,105,97),parent=body) 
circle2 = dpg.draw_circle([350+pos[2],350+pos[3]],m2*7,fill=(255,105,97),parent=body)
plot1Data = [0 for i in range(300)]
plot2Data = [0 for i in range(300)]

def computingLoop():
    global ans, currTime
    t = np.linspace(currTime,currTime+15,500)
    ans = odeint(dSdt,y0=[ans[-1][0],ans[-1][1],ans[-1][2],ans[-1][3]],t=t,args=(g,m1,m2,l1,l2))

def updateLoop():
    global theta1, theta2
    while 1:
        if running:
            computingThread = threading.Thread(target=computingLoop)
            computingThread.start()
            for the1, the2 in zip(theta1,theta2):
                if not running:
                    break
                plot1Data.pop(0)
                plot1Data.append(the1)
                plot2Data.pop(0)
                plot2Data.append(the2)
                dpg.set_value(plot1,plot1Data)
                dpg.set_value(plot2,plot2Data)
                pos = get_pos(the1,the2,l1,l2)
                dpg.configure_item(pendulum1,p1=[350,350],p2=[350+pos[0],350+pos[1]])
                dpg.configure_item(pendulum2,p1=[350+pos[0],350+pos[1]],p2=[350+pos[2],350+pos[3]])
                dpg.configure_item(circle1,center=[350+pos[0],350+pos[1]])
                dpg.configure_item(circle2,center=[350+pos[2],350+pos[3]])
                time.sleep(0.03)
            computingThread.join()
            theta1 = ans.T[0].tolist()
            theta2 = ans.T[2].tolist()


mainLoop = threading.Thread(target=updateLoop)
mainLoop.start()


dpg.create_viewport(title="Double Pendulum",width=1000,height=805,resizable=False)
dpg.setup_dearpygui()
dpg.show_viewport()
dpg.start_dearpygui()
dpg.destroy_context()

