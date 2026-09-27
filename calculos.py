from math import sqrt
import matplotlib.pyplot as plt # para instalar pip install matplotlib
import numpy as np

lisx=[0, 1, 2, 3, 4, 5]
lisy=[2.1, 7.7, 13.6, 27.2, 40.9, 61.1]

n=len(lisx)

def sumx():
    x=0
    for i in range(n):
        x+=lisx[i]
    return round(x, 4)

def sumy():
    y=0
    for i in range(n):
        y+=lisy[i]
    return round(y, 4)

def x_cuadrado():
    cuadx=0
    for i in range(n):
        cuadx+=lisx[i]**2
    return round(cuadx, 4)

def x_cubica():
    cubica=0
    for i in range(n):
        cubica+=lisx[i]**3
    return round(cubica, 4)

def x_cuarta():
    cuarta=0
    for i in range(n):
        cuarta+=lisx[i]**4
    return round(cuarta, 4)

def sumxy():
    xy=0
    for i in range(n):
        xy+=lisx[i]*lisy[i]
    return round(xy, 4)

def sumx_cuadradoy():
    x_cuadradoy=0
    for i in range(n):
        x_cuadradoy+=lisx[i]**2*lisy[i]
    return round(x_cuadradoy, 4)

def promedio(x,y):
    return (round(x/n, 4), round(y/n, 4))

def desviacion_estandar(proy): #(y-promedio y)^2
    desviacion_sum=0
    for i in range(n):
        desviaciony= (lisy[i]-proy)**2
        desviacion_sum+=desviaciony
    return round(desviacion_sum, 4)

def sum_error_lineal(a0, a1):#(y-a0-a1*x)^2
    error_sum=0
    for i in range(n):
        error= (lisy[i]-a0-a1*lisx[i])**2
        error_sum+=error
    return round(error_sum, 4)

def sum_error_cuadratica(a0, a1, a2):#(y-a0-a1*x-a2*x^2)^2
    error_sum=0
    for i in range(n):
        error= (lisy[i]-a0-a1*lisx[i]-a2*lisx[i]**2)**2
        error_sum+=error
    return round(error_sum, 4)



x= sumx()
y= sumy()
cuadx= x_cuadrado()
xcubica= x_cubica()
xcuarta= x_cuarta()
sumxy= sumxy()
x_cuadradoy= sumx_cuadradoy()
prox,proy=promedio(x,y)

#pendiente y intersección
a1 = round((n*sumxy-x*y)/(n*cuadx-x**2), 4)
a0=round(proy-a1*prox, 4)

#desviacion estandar
desviacion_sum=desviacion_estandar(proy)
desviacion_est= round(sqrt(desviacion_sum/(n-1)), 4)

# error estandar
error_sum_lineal=sum_error_lineal(a0, a1)
error_est= round(sqrt(round(error_sum_lineal,4)/(n-2)), 4)

print(f"sumatoria de desviacion estandar: {desviacion_sum}, sumatoria de error estandar: {error_sum_lineal}")

#coeficiente de determinación
r_cuadrado=round((desviacion_sum-error_sum_lineal)/desviacion_sum, 4)
r=round(sqrt(r_cuadrado), 4)

print(f"Promedio x: {prox}")
print(f"Promedio y: {proy}")
print(f"ecuacion lineal: y = {a1}x + {a0}")
print(f"Desviación estándar: {desviacion_est}")
print(f"Error estándar: {error_est}")
print(f"Coeficiente de determinación: {r_cuadrado}")
print(f"Coeficiente de correlación: {r}")

#regresion polinomial

A = [
    [n,    x,  cuadx],
    [x, cuadx, x_cubica()],
    [cuadx, x_cubica(), x_cuarta()]
]

b = [y, sumxy, x_cuadradoy]
a0_pol, a1_pol, a2_pol = np.linalg.solve(A, b)
a0_pol, a1_pol, a2_pol = round(a0_pol, 4), round(a1_pol, 4), round(a2_pol, 4)

print(f"ecuacion polinomial: y = {a0_pol} + {a1_pol}x + {a2_pol}x^2")

#error estandar polinomial
error_sum_polinomial = round(sum_error_cuadratica(a0_pol, a1_pol, a2_pol), 4)
error_est_polinomial = round(sqrt(error_sum_polinomial/(n-3)), 4)

print(f"sumatoria de error estandar polinomial: {error_sum_polinomial}")
print(f"Error estándar polinomial: {error_est_polinomial}")

# coeficiente de determinación polinomial
r_cuadrado_polinomial = round((desviacion_sum - error_sum_polinomial) / desviacion_sum, 4)
r_polinomial = round(sqrt(r_cuadrado_polinomial), 4)

print(f"Coeficiente de determinación polinomial: {r_cuadrado_polinomial}")
print(f"Coeficiente de correlación polinomial: {r_polinomial}")

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

#grafica de la regresion polinomial
def polinomio(x):
    return [a0_pol + a1_pol*xi + a2_pol*xi**2 for xi in x]

plt.scatter(lisx, lisy, label="Datos")
plt.plot(lisx, polinomio(lisx), label="Regresión polinomial")
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Regresión polinomial")
plt.legend()
plt.grid()
plt.show()