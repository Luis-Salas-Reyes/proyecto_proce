from sumatoria import *

def promedio(x,y):
    return (x/len(lisx),y/len(lisy))

prox,proy=promedio(x,y)
print("Promedio de x:",prox)
print("Promedio de y:",proy)

def pendiente(x,y,cuadx,cuady,sumxy):
    a1 = (len(lisx)*sumxy-x*y)/(len(lisx)*cuadx-x**2)
    return a1

def intersaccion(prox,proy,a1):
    a0=proy-a1*prox
    return a0


