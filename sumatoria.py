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

x= sumx()
y= sumy()
cuadx= cuadrx()
sumxy= sumxy()