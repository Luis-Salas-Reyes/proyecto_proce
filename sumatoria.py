lisx=[1, 2, 3, 4, 5]
lisy=[2,2.5,3,3.5,4]

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

def cuadry():
    cuady=0
    for i in range(len(lisy)):
        cuady+=lisy[i]**2
    return cuady

def sumxy():
    xy=0
    for i in range(len(lisx)):
        xy+=lisx[i]*lisy[i]
    return xy

x= sumx()
y= sumy()
cuadx= cuadrx()
cuady= cuadry()
sumxy= sumxy()