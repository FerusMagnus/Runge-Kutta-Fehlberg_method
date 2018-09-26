# -*- coding: utf-8 -*-
from numpy import*
import matplotlib.pyplot as plt

def rungeKutta(f, to, yo, tEnd, tau):
         def increment(f, t, y, tau):  # поиск приближённого решения методом Рунге—Кутта—Фельберга.
                  k1 = tau*f(t, y)
                  k2 = tau*f(t+(1/4)*tau, y+(1/4)*k1)
                  k3 = tau*f(t+(3/8)*tau, y+(3/32)*k1+(9/32)*k2)
                  k4 = tau*f(t+(12/13)*tau, y+(1932/2197)*k1-(7200/2197)*k2+(7296/2197)*k3)
                  k5 = tau*f(t+tau, y+(439/216)*k1-8*k2+(3680/513)*k3 -(845/4104)*k4)
                  k6 = tau*f(t+(1/2)*tau, y-(8/27)*k1+2*k2-(3544/2565)*k3 +(1859/4104)*k4-(11/40)*k5)
                  return (16/135)*k1+(6656/12825)*k3+(28561/56430)*k4-(9/50)*k5+(2/55)*k6

         t = []  # подготовка пустого списка t
         y = []  # подготовка пустого списка y
         t.append(to)  # внесение в список t начального значения to
         y.append(yo)  # внесение в список y начального значения yo
         while to < tEnd:  # внесение результатов расчёта в массивы t,y
                  tau = min(tau, tEnd - to)  # определение минимального шага tau
                  yo = yo + increment(f, to, yo, tau)  # расчёт значения в точке t0,y0 для задачи Коши
                  to = to + tau  # приращение времени
                  t.append(to)  # заполнение массива t
                  y.append(yo)  # заполнение массива y
         return array(t), array(y)

def f(t, y): # функция правых частей системы ОДУ
         f = zeros([6])
         f[0]= y[2]
         f[1]= y[3]
         f[2]= -t * y[1]
         f[3]= y[1] - 2 * y[2]
         f[4]= y[3] - 1.5 * y[4]
         f[5] = y[3] + 2 * y[4] - y[5]
         return f

to = 0  # начальный момент отсчёта времени
tEnd = 2  # конечный момент отсчёта времени
#yo = array([0,1,1,0])  # начальные условия
#yo = array([0,1,1,0,1])
yo = array([0.2, 1, 1, 0, 1, 1.5])
tau = 0.2  # шаг

t, y = rungeKutta(f, to, yo, tEnd, tau)
y1 = array([i[0] for i in y])
y2 = array([i[1] for i in y])
y3 = array([i[2] for i in y])
y4 = array([i[3] for i in y])
y5 = array([i[4] for i in y])
y6 = array([i[5] for i in y])

# визуализация
plt.plot(t, y1, label = 'y1')
plt.plot(t, y2, label = 'y2')
plt.plot(t, y3, label = 'y3')
plt.plot(t, y4, label = 'y4')
plt.plot(t, y5, label = 'y5')
plt.plot(t, y6, label = 'y6')

print(y1)
print(y2)
print(y3)
print(y4)
print(y5)
print(y6)

plt.title("Результаты численного решения системы ОДУ \n с использованием \n метода Рунге—Кутта—Фельберга ")
plt.xlabel('Значение t')
plt.ylabel('Значение Yn')
plt.legend(loc='best')
# plt.xlim(0, 8)
# plt.ylim(-0.1, 2)
plt.grid(True)
plt.show()