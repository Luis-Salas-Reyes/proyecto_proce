from math import sqrt

lisx=[1, 2, 3, 4, 5, 6, 7]
lisy=[0.5, 2.5, 2.0, 4.0, 3.5, 6.0, 5.5]

def sumx():
    x=0
    for i in range(len(lisx)):
        x+=lisx[i]
    return x

def sumy():
    y=0
    for i in range(len(lisy)):
        y+=lisy[i]
    return y

def cuadrx():
    cuadx=0
    for i in range(len(lisx)):
        cuadx+=lisx[i]**2
    return cuadx

def sumxy():
    xy=0
    for i in range(len(lisx)):
        xy+=lisx[i]*lisy[i]
    return xy

def promedio(x,y):
    return (x/len(lisx),y/len(lisy))

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

x= sumx()
y= sumy()
cuadx= cuadrx()
sumxy= sumxy()

prox,proy=promedio(x,y)

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