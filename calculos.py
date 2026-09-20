from sumatoria import *
from math import sqrt

def promedio(x,y):
    return (x/len(lisx),y/len(lisy))

prox,proy=promedio(x,y)

def pendiente(x,y,cuadx,sumxy):
    a1 = (len(lisx)*sumxy-x*y)/(len(lisx)*cuadx-x**2)
    return a1

def intersaccion(prox,proy,a1):
    a0=proy-a1*prox
    return a0

def desviacion_estandar():
    desviacion_sum=0
    for i in range(len(lisx)):
        desviaciony= (lisy[i]-proy)**2
        desviacion_sum+=desviaciony
    deviacion_est= sqrt(desviacion_sum/(len(lisy)-1))
    return deviacion_est

def error_estandar():
    error_sum=0
    for i in range(len(lisx)):
        error= (lisy[i]-a0-a1*lisx[i])**2
        error_sum+=error
    error_est= sqrt(error_sum/(len(lisy)-2))
    return error_est

def coeficiente_correlacion():
    r=(des-err)/des
    return r

a1=pendiente(x,y,cuadx,sumxy)
a0=intersaccion(prox,proy,a1)
des=desviacion_estandar()
err=error_estandar()
r=coeficiente_correlacion()

print(f"Promedio x: {prox}")
print(f"Promedio y: {proy}")
print(f"Pendiente: {a1}")
print(f"Intersección: {a0}")
print(f"Desviación estándar: {des}")
print(f"Error estándar: {err}")
print(f"Coeficiente de correlación: {r}")