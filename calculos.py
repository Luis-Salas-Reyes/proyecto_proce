from math import sqrt
import matplotlib.pyplot as plt # para instalar pip install matplotlib

lisx=[1, 2, 3, 4, 5, 6, 7]
lisy=[0.5, 2.5, 2.0, 4.0, 3.5, 6.0, 5.5]

def sumx():
    x=0
    for i in range(len(lisx)):
        x+=lisx[i]
    return round(x, 4)

def sumy():
    y=0
    for i in range(len(lisy)):
        y+=lisy[i]
    return round(y, 4)

def cuadrx():
    cuadx=0
    for i in range(len(lisx)):
        cuadx+=lisx[i]**2
    return round(cuadx, 4)

def sumxy():
    xy=0
    for i in range(len(lisx)):
        xy+=lisx[i]*lisy[i]
    return round(xy, 4)

def promedio(x,y):
    return (round(x/len(lisx), 4), round(y/len(lisy), 4))

def desviacion_estandar(proy): #(y-promedio y)^2
    desviacion_sum=0
    for i in range(len(lisx)):
        desviaciony= (lisy[i]-proy)**2
        desviacion_sum+=desviaciony
    return round(desviacion_sum, 4)

def error_estandar(a0, a1):#(y-a0-a1*x)^2
    error_sum=0
    for i in range(len(lisx)):
        error= (lisy[i]-a0-a1*lisx[i])**2
        error_sum+=error
    return round(error_sum, 4)



x= sumx()
y= sumy()
cuadx= cuadrx()
sumxy= sumxy()
prox,proy=promedio(x,y)

#pendiente y intersección
a1 = round((len(lisx)*sumxy-x*y)/(len(lisx)*cuadx-x**2), 4)
a0=round(proy-a1*prox, 4)

#desviacion estandar
desviacion_sum=desviacion_estandar(proy)
desviacion_est= round(sqrt(desviacion_sum/(len(lisy)-1)), 4)

# error estandar
error_sum=error_estandar(a0, a1)
error_est= round(sqrt(round(error_sum,4)/(len(lisy)-2)), 4)

print(f"desviacion estandar: {desviacion_sum}, error estandar: {error_sum}")

#coeficiente de correlación
r=round((desviacion_sum-error_sum)/desviacion_sum, 4)

print(f"Promedio x: {prox}")
print(f"Promedio y: {proy}")
print(f"Pendiente: {a1}")
print(f"Intersección: {a0}")
print(f"Desviación estándar: {desviacion_est}")
print(f"Error estándar: {error_est}")
print(f"Coeficiente de correlación: {r}")

#grafica de la regresion lineal
def recta(x):
    return [a1*xi + a0 for xi in x]

plt.scatter(lisx, lisy, label="Datos")

plt.plot(lisx, recta(lisx), label="Regresión lineal")

plt.xlabel("X")
plt.ylabel("Y")
plt.title("Regresión lineal")
plt.legend()
plt.grid()

plt.show()