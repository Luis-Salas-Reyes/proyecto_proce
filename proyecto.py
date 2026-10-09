# =============================================================================
# PROYECTO DE PROCESAMIENTO NUMERICO - CORTE 2
# Atenuacion acustica en pasillo interior (tono puro 1000 Hz)
#
# Mateo Enrique Pajaro Bossio - T00083783
# Luis Miguel Rodriguez Vega  - T00084082
# Luis Daniel Salas Reyes     - T00082453
# Jorge Isaac Hincapie Pautt  - T00082792
#
# Universidad Tecnologica de Bolivar
# Profesor: Oscar Porto Solano, M.Sc.
# Periodo: 2p2026
#
# Todos los resultados numericos estan redondeados a cuatro (4) cifras
# decimales y llevan unidades de medida segun corresponda.
# =============================================================================

from math import sqrt
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
import numpy as np

# =============================================================================
# PALETA DE COLORES
# =============================================================================
FONDO   = '#E8DEC8'
TEXTO   = '#3A3A3A'
PUNTO   = '#2A2A2A'
GRILLA  = '#7A7A7A'
LINEAL  = '#B85C38'
CUADRAT = '#4A6FA5'
COLORES_LAGRANGE = {1: '#1f77b4', 2: '#ff7f0e', 3: '#2ca02c', 4: '#d62728'}

# =============================================================================
# DATOS DEL EXPERIMENTO
# =============================================================================
proyecto = {
    "nombre": "Proyecto: Atenuacion acustica en pasillo interior",
    "titulo": "NPS (dB) vs. Distancia (m) - Tono puro 1000 Hz",
    "x_crudo": [0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50,
                2.75, 3.00, 3.25, 3.50, 3.75, 4.00, 4.25, 4.50, 4.75, 5.00],
    "y_crudo": [-22.6206, -28.9196, -26.4120, -28.7881, -36.9274,
                -28.9565, -39.2525, -40.5908, -33.8486, -33.8171,
                -32.6223, -30.1040, -37.5206, -37.3114, -39.6176,
                -33.8355, -28.7045, -37.7981, -33.6955, -37.4102],
    "xlabel": "Distancia d (m)",
    "ylabel": "NPS (dB)",
    "unidad_x": "m",
    "unidad_y": "dB",
    "unidad_a1": "dB/m",
    "unidad_a2": "dB/m^2",
}


# =============================================================================
# FUNCIONES AUXILIARES
# =============================================================================
def sumx(datos, n):
    s = 0
    for i in range(n):
        s += datos["x"][i]
    return s

def sumy(datos, n):
    s = 0
    for i in range(n):
        s += datos["y"][i]
    return s

def x_cuadrado(datos, n):
    s = 0
    for i in range(n):
        s += datos["x"][i] ** 2
    return s

def x_cubica(datos, n):
    s = 0
    for i in range(n):
        s += datos["x"][i] ** 3
    return s

def x_cuarta(datos, n):
    s = 0
    for i in range(n):
        s += datos["x"][i] ** 4
    return s

def sum_xy(datos, n):
    s = 0
    for i in range(n):
        s += datos["x"][i] * datos["y"][i]
    return s

def sumx_cuadradoy(datos, n):
    s = 0
    for i in range(n):
        s += datos["x"][i] ** 2 * datos["y"][i]
    return s

def promedio(x, y, n):
    return (x / n, y / n)

def desviacion_estandar(datos, n, proy):
    s = 0
    for i in range(n):
        s += (datos["y"][i] - proy) ** 2
    return s

def sum_error_lineal(datos, n, a0, a1):
    s = 0
    for i in range(n):
        s += (datos["y"][i] - a0 - a1 * datos["x"][i]) ** 2
    return s

def sum_error_cuadratica(datos, n, a0, a1, a2):
    s = 0
    for i in range(n):
        s += (datos["y"][i] - a0 - a1 * datos["x"][i] - a2 * datos["x"][i] ** 2) ** 2
    return s

def recta(x, a0, a1):
    return [a1 * xi + a0 for xi in x]

def polinomio(x, a0, a1, a2):
    return [a0 + a1 * xi + a2 * xi ** 2 for xi in x]


# =============================================================================
# FILTRO DE MEDIA MOVIL (ventana 5) - Chapra, capitulo 14
# =============================================================================
def media_movil(y, ventana=5):
    n = len(y)
    y_filtrado = []
    mitad = ventana // 2
    for i in range(n):
        ini = max(0, i - mitad)
        fin = min(n, i + mitad + 1)
        y_filtrado.append(sum(y[ini:fin]) / (fin - ini))
    return y_filtrado


# =============================================================================
# INTERPOLACION DE LAGRANGE
# =============================================================================
def lagrange(xp, yp, xe):
    n = len(xp)
    total = 0.0
    for i in range(n):
        Li = 1.0
        for j in range(n):
            if i != j:
                Li *= (xe - xp[j]) / (xp[i] - xp[j])
        total += yp[i] * Li
    return total


# =============================================================================
# VALIDACION CRUZADA LEAVE-ONE-OUT (LOOCV) CON LAGRANGE
# =============================================================================
def loocv_lagrange(x, y, grado):
    n = len(x)
    pred = [0.0] * n
    et = [0.0] * n
    er = [0.0] * n
    for i in range(n):
        vecinos = [j for j in range(n) if j != i]
        vecinos.sort(key=lambda j: abs(x[j] - x[i]))
        vecinos = vecinos[:grado + 1]
        vecinos.sort(key=lambda j: x[j])
        xp = [x[j] for j in vecinos]
        yp = [y[j] for j in vecinos]
        pred[i] = lagrange(xp, yp, x[i])
        et[i] = y[i] - pred[i]
        if abs(y[i]) > 1e-10:
            er[i] = 100.0 * et[i] / y[i]
    return pred, et, er


# =============================================================================
# PROGRAMA PRINCIPAL
# =============================================================================
print("\n" + "="*75)
print(f"  {proyecto['nombre']}")
print("="*75)

# Aplicar filtro de media movil
x = proyecto["x_crudo"]
y = media_movil(proyecto["y_crudo"], ventana=5)
n = len(x)

ux  = proyecto["unidad_x"]
uy  = proyecto["unidad_y"]
ua1 = proyecto["unidad_a1"]
ua2 = proyecto["unidad_a2"]

print(f"\nDatos de X (distancia): {x}")
print(f"Datos de Y crudos:      {proyecto['y_crudo']}")
print(f"Datos de Y filtrados:   {[round(v, 4) for v in y]}")

# Mostrar tabla de filtrado
print("\n--- FILTRADO POR MEDIA MOVIL (ventana 5) ---")
print(f"{'#':>3} {'d (m)':>7} {'NPS crudo':>12} {'NPS filtrado':>14}")
print("-"*42)
for i in range(n):
    print(f"{i+1:>3} {x[i]:>7.2f} {proyecto['y_crudo'][i]:>12.4f} {y[i]:>14.4f}")

# Actualizar datos del proyecto con y filtrado
datos = dict(proyecto)
datos["x"] = x
datos["y"] = y

# =============================================================================
# REGRESION LINEAL
# =============================================================================
sx        = sumx(datos, n)
sy        = sumy(datos, n)
cuadx     = x_cuadrado(datos, n)
xcubica   = x_cubica(datos, n)
xcuarta   = x_cuarta(datos, n)
sxy       = sum_xy(datos, n)
x_cuadradoy = sumx_cuadradoy(datos, n)
prox, proy = promedio(sx, sy, n)

# Coeficientes lineales
a1 = (n * sxy - sx * sy) / (n * cuadx - sx ** 2)
a0 = proy - a1 * prox

# Desviacion estandar y error estandar
desviacion_sum = desviacion_estandar(datos, n, proy)
desviacion_est = sqrt(desviacion_sum / (n - 1))
error_sum_lineal = sum_error_lineal(datos, n, a0, a1)
error_est = sqrt(error_sum_lineal / (n - 2))

# Coeficientes de determinacion y correlacion
r_cuadrado = (desviacion_sum - error_sum_lineal) / desviacion_sum
r = sqrt(r_cuadrado) if r_cuadrado >= 0 else float('nan')

print("\n" + "="*75)
print("  REGRESION LINEAL")
print("="*75)
print(f"Suma de cuadrados total (SST): {desviacion_sum:.4f} {uy}^2")
print(f"Suma de cuadrados del error (SSE) lineal: {error_sum_lineal:.4f} {uy}^2")
print(f"Promedio x: {prox:.4f} {ux}")
print(f"Promedio y: {proy:.4f} {uy}")
print(f"a1 lineal: {a1:.4f} {ua1}")
print(f"a0 lineal: {a0:.4f} {uy}")
print(f"Ecuacion lineal: y = {a0:.4f} {'+' if a1 >= 0 else '-'} {abs(a1):.4f}x")
print(f"Desviacion estandar (Sy):  {desviacion_est:.4f} {uy}")
print(f"Error estandar (Sy/x):     {error_est:.4f} {uy}")
print(f"Coeficiente de determinacion (R2): {r_cuadrado:.4f}")
print(f"Coeficiente de correlacion (r):    {r:.4f}")


# =============================================================================
# REGRESION POLINOMIAL GRADO 2
# =============================================================================
print("\n" + "="*75)
print("  REGRESION POLINOMIAL GRADO 2")
print("="*75)
print("Datos de la regresion polinomial:")
A = [
    [n, sx, cuadx],
    [sx, cuadx, xcubica],
    [cuadx, xcubica, xcuarta]
]
b = [sy, sxy, x_cuadradoy]

a0_pol, a1_pol, a2_pol = np.linalg.solve(A, b)

print(f"a2 polinomial: {a2_pol:.4f} {ua2}")
print(f"a1 polinomial: {a1_pol:.4f} {ua1}")
print(f"a0 polinomial: {a0_pol:.4f} {uy}")

if abs(a2_pol) >= 0.01:
    a2_texto = f"{abs(a2_pol):.4f}"
else:
    a2_texto = f"{abs(a2_pol):.4e}"

print(f"Ecuacion polinomial: y = {a0_pol:.4f} "
      f"{'+' if a1_pol >= 0 else '-'} {abs(a1_pol):.4f}x "
      f"{'+' if a2_pol >= 0 else '-'} {a2_texto}x^2")

# Error estandar polinomial
error_sum_polinomial = sum_error_cuadratica(datos, n, a0_pol, a1_pol, a2_pol)
error_est_polinomial = sqrt(error_sum_polinomial / (n - 3))

print(f"Sumatoria del error estandar polinomial: {error_sum_polinomial:.4f}")
print(f"Error estandar polinomial (Sy/x): {error_est_polinomial:.4f} {uy}")

# Coeficiente de determinacion polinomial
r_cuadrado_polinomial = (desviacion_sum - error_sum_polinomial) / desviacion_sum
r_polinomial = sqrt(r_cuadrado_polinomial) if r_cuadrado_polinomial >= 0 else float('nan')

print(f"Coeficiente de determinacion polinomial (R2): {r_cuadrado_polinomial:.4f}")
print(f"Coeficiente de correlacion polinomial (r):    {r_polinomial:.4f}")

if r_cuadrado_polinomial > r_cuadrado:
    print("Mejor modelo: CUADRATICO")
else:
    print("Mejor modelo: LINEAL")


# =============================================================================
# INTERPOLACION DE LAGRANGE + LOOCV (grados 1 a 4)
# =============================================================================
print("\n" + "="*75)
print("  INTERPOLACION DE LAGRANGE + LOOCV (grados 1 a 4)")
print("="*75)

resultados_lagrange = {}
for g in [1, 2, 3, 4]:
    pred, et, er = loocv_lagrange(x, y, g)

    er_abs = [abs(e) for e in er]
    er_mean = sum(er_abs) / n
    var = sum((e - er_mean) ** 2 for e in er_abs) / (n - 1)
    er_std = sqrt(var)

    rmse = sqrt(sum((pred[i] - y[i]) ** 2 for i in range(n)) / n)
    mae = sum(abs(pred[i] - y[i]) for i in range(n)) / n
    sr = sum((y[i] - pred[i]) ** 2 for i in range(n))
    r2 = 1 - sr / desviacion_sum

    resultados_lagrange[g] = {
        "pred": pred, "et": et, "er": er,
        "er_mean": er_mean, "er_std": er_std,
        "rmse": rmse, "mae": mae, "r2": r2
    }

    print(f"\n--- GRADO {g} ---")
    print(f"Er% promedio (Ēr%): {er_mean:.4f} %")
    print(f"Desviacion estandar (S_Er%): {er_std:.4f} %")
    print(f"RMSE: {rmse:.4f} {uy}")
    print(f"MAE:  {mae:.4f} {uy}")
    print(f"R2:   {r2:.4f}")

print("\n--- TABLA RESUMEN LOOCV ---")
print(f"{'Grado':>6} {'Ēr%':>10} {'S_Er%':>10} {'RMSE(dB)':>12} {'MAE(dB)':>10} {'R2':>10}")
print("-"*62)
for g in [1, 2, 3, 4]:
    r_ = resultados_lagrange[g]
    print(f"{g:>6} {r_['er_mean']:>10.4f} {r_['er_std']:>10.4f} "
          f"{r_['rmse']:>12.4f} {r_['mae']:>10.4f} {r_['r2']:>10.4f}")


# =============================================================================
# FIGURA 1: Dispersion + Lagrange grados 1 a 4
# =============================================================================
xs = np.linspace(min(x), max(x), 500)

fig1, ax1 = plt.subplots(figsize=(10, 6))
fig1.patch.set_facecolor(FONDO)
ax1.set_facecolor(FONDO)

ax1.scatter(x, y, s=110, color=PUNTO, edgecolors=FONDO,
            linewidths=1.5, zorder=6, label='Datos experimentales')

for g in [1, 2, 3, 4]:
    ys_curva = []
    for xi in xs:
        vecinos = list(range(n))
        vecinos.sort(key=lambda j: abs(x[j] - xi))
        vecinos = vecinos[:g + 1]
        vecinos.sort(key=lambda j: x[j])
        xp = [x[j] for j in vecinos]
        yp = [y[j] for j in vecinos]
        ys_curva.append(lagrange(xp, yp, xi))
    ax1.plot(xs, ys_curva, color=COLORES_LAGRANGE[g], lw=2.2,
             label=f'Lagrange grado {g}')

ax1.set_xlabel(proyecto["xlabel"], fontsize=13, fontweight='bold', color=TEXTO)
ax1.set_ylabel(proyecto["ylabel"], fontsize=13, fontweight='bold', color=TEXTO)
ax1.set_title('Figura 1: Interpolacion de Lagrange (grados 1 a 4)',
              fontsize=14, fontweight='bold', color=TEXTO, pad=15)
ax1.grid(alpha=0.3, linestyle='--', color=GRILLA)
ax1.set_axisbelow(True)
ax1.tick_params(axis='both', colors=TEXTO, labelsize=11)
for spine in ax1.spines.values():
    spine.set_color(TEXTO)
    spine.set_linewidth(1.2)
legend1 = ax1.legend(fontsize=11, loc='lower right', frameon=True,
                     facecolor=FONDO, edgecolor=TEXTO)
for text in legend1.get_texts():
    text.set_color(TEXTO)

plt.tight_layout()
plt.savefig('figura1_lagrange.png', dpi=220, bbox_inches='tight', facecolor=FONDO)
plt.show()


# =============================================================================
# FIGURA 2: Perfil espacial del error (LOOCV)
# =============================================================================
fig2, ax2 = plt.subplots(figsize=(10, 6))
fig2.patch.set_facecolor(FONDO)
ax2.set_facecolor(FONDO)

ax2.axhline(0, color=GRILLA, ls='--', alpha=0.5, lw=1)

for g in [1, 2, 3, 4]:
    er_abs_g = [abs(e) for e in resultados_lagrange[g]["er"]]
    ax2.plot(x, er_abs_g, 'o-', color=COLORES_LAGRANGE[g],
             lw=2.0, ms=8, label=f'Grado {g}',
             markeredgecolor=FONDO, markeredgewidth=1.5, zorder=3)

ax2.set_xlabel(proyecto["xlabel"], fontsize=13, fontweight='bold', color=TEXTO)
ax2.set_ylabel('Error relativo porcentual |Er%|',
               fontsize=13, fontweight='bold', color=TEXTO)
ax2.set_title('Figura 2: Perfil espacial del error (LOOCV)',
              fontsize=14, fontweight='bold', color=TEXTO, pad=15)
ax2.grid(alpha=0.3, linestyle='--', color=GRILLA)
ax2.set_axisbelow(True)
ax2.tick_params(axis='both', colors=TEXTO, labelsize=11)
for spine in ax2.spines.values():
    spine.set_color(TEXTO)
    spine.set_linewidth(1.2)
legend2 = ax2.legend(fontsize=11, loc='upper right', frameon=True,
                     facecolor=FONDO, edgecolor=TEXTO)
for text in legend2.get_texts():
    text.set_color(TEXTO)

plt.tight_layout()
plt.savefig('figura2_residuos.png', dpi=220, bbox_inches='tight', facecolor=FONDO)
plt.show()


# =============================================================================
# FIGURA 3: Dispersion (solo datos)
# =============================================================================
fig3, ax3 = plt.subplots(figsize=(10, 6))
fig3.patch.set_facecolor(FONDO)
ax3.set_facecolor(FONDO)

ax3.scatter(x, y, s=110, color=PUNTO, edgecolors=FONDO,
            linewidths=1.5, zorder=3)
ax3.set_xlabel(proyecto["xlabel"], fontsize=13, fontweight='bold', color=TEXTO)
ax3.set_ylabel(proyecto["ylabel"], fontsize=13, fontweight='bold', color=TEXTO)
ax3.set_title('Datos experimentales: NPS vs. Distancia',
              fontsize=14, fontweight='bold', color=TEXTO, pad=15)
ax3.grid(alpha=0.3, linestyle='--', color=GRILLA)
ax3.set_axisbelow(True)
ax3.tick_params(axis='both', colors=TEXTO, labelsize=11)
for spine in ax3.spines.values():
    spine.set_color(TEXTO)
    spine.set_linewidth(1.2)

plt.tight_layout()
plt.savefig('grafica1_dispersion.png', dpi=220, bbox_inches='tight', facecolor=FONDO)
plt.show()


# =============================================================================
# FIGURA 4: Regresion lineal vs cuadratica
# =============================================================================
fig4, ax4 = plt.subplots(figsize=(10, 6))
fig4.patch.set_facecolor(FONDO)
ax4.set_facecolor(FONDO)

ax4.scatter(x, y, s=110, color=PUNTO, edgecolors=FONDO,
            linewidths=1.5, zorder=5, label='Datos experimentales')

xs_graf = np.linspace(min(x), max(x), 300)
ax4.plot(xs_graf, recta(xs_graf, a0, a1), '--', color=LINEAL, lw=2.5,
         label=f'Lineal (R2 = {r_cuadrado:.3f})')
ax4.plot(xs_graf, polinomio(xs_graf, a0_pol, a1_pol, a2_pol), '-',
         color=CUADRAT, lw=2.8, label=f'Cuadratica (R2 = {r_cuadrado_polinomial:.3f})')

ax4.set_xlabel(proyecto["xlabel"], fontsize=13, fontweight='bold', color=TEXTO)
ax4.set_ylabel(proyecto["ylabel"], fontsize=13, fontweight='bold', color=TEXTO)
ax4.set_title('Regresion por Minimos Cuadrados: Lineal vs. Cuadratica',
              fontsize=14, fontweight='bold', color=TEXTO, pad=15)
ax4.grid(alpha=0.3, linestyle='--', color=GRILLA)
ax4.set_axisbelow(True)
ax4.tick_params(axis='both', colors=TEXTO, labelsize=11)
for spine in ax4.spines.values():
    spine.set_color(TEXTO)
    spine.set_linewidth(1.2)
legend4 = ax4.legend(fontsize=11, loc='lower right', frameon=True,
                     facecolor=FONDO, edgecolor=TEXTO)
for text in legend4.get_texts():
    text.set_color(TEXTO)

plt.tight_layout()
plt.savefig('grafica2_regresiones.png', dpi=220, bbox_inches='tight', facecolor=FONDO)
plt.show()


# =============================================================================
# ANALISIS CRITICO - Grado optimo
# =============================================================================
print("\n" + "="*75)
print("  ANALISIS CRITICO Y GRADO OPTIMO")
print("="*75)

rmse_min = min(resultados_lagrange[g]["rmse"] for g in [1, 2, 3, 4])
ser_min  = min(resultados_lagrange[g]["er_std"] for g in [1, 2, 3, 4])
r2_max   = max(resultados_lagrange[g]["r2"] for g in [1, 2, 3, 4])

print(f"\nGrado con menor RMSE:  ", end="")
for g in [1, 2, 3, 4]:
    if abs(resultados_lagrange[g]["rmse"] - rmse_min) < 1e-6:
        print(f"{g} ({rmse_min:.4f} dB)")

print(f"Grado con menor S_Er%: ", end="")
for g in [1, 2, 3, 4]:
    if abs(resultados_lagrange[g]["er_std"] - ser_min) < 1e-6:
        print(f"{g} ({ser_min:.4f} %)")

print(f"Grado con mayor R2:    ", end="")
for g in [1, 2, 3, 4]:
    if abs(resultados_lagrange[g]["r2"] - r2_max) < 1e-6:
        print(f"{g} ({r2_max:.4f})")

print("\nGrado optimo: 1 (interpolacion lineal local)")
print("Justificacion: principio de parsimonia (Chapra, ejemplos 18.6 y 18.7)")
print("Los grados altos sobreajustan el ruido experimental.")

print("\n" + "="*75)
print("  PROCESO COMPLETO")
print("="*75)