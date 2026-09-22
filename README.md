# Efectos de Superficie en Resonancia Magnética Nuclear a Nanoescala

**Trabajo Fin de Máster (TFM)** — Máster Oficial en Ingeniería Matemática y Computación.

> **Título:** *Efectos de superficie en Resonancia Magnética Nuclear a Nanoescala*  
> Modelización de la función de correlación Nano-NMR con difusión inhomogénea D(z) en geometría cilíndrica.

---

## Contexto físico

Las moléculas de líquido confinadas cerca de superficies sólidas experimentan una **reducción del coeficiente de difusión** D(z) debida a interacciones de Van der Waals, puentes de hidrógeno e hidrofobicidad. Este fenómeno modifica la señal de correlación en experimentos de Resonancia Magnética Nuclear a nanoescala (Nano-NMR), donde un defecto NV (nitrogen-vacancy) en diamante actúa como sensor a profundidad d = 15 nm.

El modelo resuelve la **ecuación de difusión en coordenadas cilíndricas** (ρ, z) con D = D(z) inhomogéneo:

$$\frac{\partial P}{\partial t} = \frac{1}{\rho}\frac{\partial}{\partial \rho}\!\left(\rho D(z)\frac{\partial P}{\partial \rho}\right) + \frac{\partial}{\partial z}\!\left(D(z)\frac{\partial P}{\partial z}\right)$$

con perfil de difusión de transición suave:

$$D(z) = D_s + (D_0 - D_s)\!\left[\left(1 - e^{-\lambda(z-d)}\right)(1 - S) + S\right], \quad S = \tfrac{1}{2}\left[1 + \tanh\!\left(\tfrac{z - (z_c + d)}{\delta}\right)\right]$$

y parámetros físicos: d = 15 nm, D₀ = 2×10³ nm²/μs, Dₛ = 2×10¹ nm²/μs, R = L = 50d = 750 nm, z_c = 2.5d = 37.5 nm, δ = 2d/3 ≈ 10 nm, λ = 0.75/d = 0.05 nm⁻¹.

La **función de correlación G(t)** — observable experimental de Nano-NMR — se calcula como:

$$G(t) = \frac{5}{16\pi} \sum_{n,m} F_{nm}^2\, e^{-\Lambda_{nm}\, t},\qquad F_{nm}=\int_0^R\!\!\int_d^{d+L} f(\rho,z)\,\mathcal{R}_n(\rho)\,Z_{nm}(z)\,\rho\,d\rho\,dz$$

donde Λ_nm son los autovalores del operador de difusión (unidades μs⁻¹) y F_nm la proyección del kernel dipolar sobre cada modo. El kernel es el armónico esférico Y₂⁰:

$$f(\rho,z)=\frac{3\cos^2\theta-1}{r^3}=\frac{2z^2-\rho^2}{(\rho^2+z^2)^{5/2}},\qquad \cos\theta=\frac{z}{r},\quad r^2=\rho^2+z^2$$

Por la identidad de Parseval, esta forma da **G(0) = B_rms² = (5/16π) Σ F_nm²**, y la diferencia respecto a (5/16π)∫f² ρ dρ dz mide directamente el error de truncamiento modal.

> **Revisión física y matemática:** las correcciones aplicadas a los cuadernos, con las comprobaciones numéricas que las respaldan, están documentadas en [`CORRECTIONS_LOG.md`](CORRECTIONS_LOG.md).

---

## Resultados clave

Cifras medidas sobre la última ejecución de los cuadernos (τ = D_eff·t/d², D_eff = 1792.6 nm²/μs):

- **Régimen de tiempos largos** (10 ≲ τ ≲ 10³): la pendiente ajustada de log G frente a log τ es **−1.43**, consistente con la ley de potencias τ^{−3/2} del confinamiento 3D.
- **Meseta a τ ≳ 10³**: G(τ) se aplana (pendiente ≈ 0). No es un fallo de la ley anterior: el cilindro es **cerrado**, de modo que el modo uniforme (n=0, m=0, Λ=0) da un valor estacionario G(∞) > 0. La ley τ^{−3/2} sólo es observable antes de que la molécula explore todo el volumen.
- **Régimen de tiempos cortos** (τ ≪ 1): el cuaderno compara con exp(−6τ). Véase la nota J2 de [`CORRECTIONS_LOG.md`](CORRECTIONS_LOG.md): esta referencia y la de `CurvaTeórica.ipynb` no coinciden y conviene fijar la convención.
- **Convergencia modal** (comprobación de Parseval, impresa por el propio cuaderno): con 100×100 modos la suma modal captura el **99.93 %** de (5/16π)∫f² ρ dρ dz; con 50×50, el 94.80 %. Este control es el que detecta el truncamiento; comparar curvas ya normalizadas G/G(0) es una prueba mucho más débil.
- **Contraste entre los dos modelos**: `TFM_V5.ipynb` (perfil D(z) continuo, volumen finito) y `TFM_TwoRegionModel_z'.ipynb` (dos regiones, raíces trascendentales y autofunciones analíticas) usan maquinaria completamente distinta y coinciden en el 99.93 % capturado y en (5/16π)∫f² ρ dρ dz a cinco cifras (3.6829×10⁻⁶ frente a 3.6831×10⁻⁶).
- **Conservación de masa** (comprobación del operador axial, impresa por el cuaderno): Λ₀,₀ = 1.0×10⁻⁷ μs⁻¹ ≈ 0 con κ = 0, y la norma axial cumple ∫Z² dz = 1 con desviación máxima 4×10⁻¹⁵.
- **Sensibilidad paramétrica**: los barridos de (λ, δ, z_c, κ) muestran que la transición entre regímenes es sensible al perfil de interfaz, lo que permite inferir parámetros microscópicos a partir de datos experimentales de G(t).
- **Modelo de dos regiones**: alternativa analítica con autovalores transcendentales, mucho más rápida que el modelo de perfil continuo. Su espectro incluye **modos superficiales evanescentes** (Λ² < α_n²D₂), que son los más lentos del problema: para n = 1…4 hay 1, 2, 2 y 3 de ellos respectivamente.
- **Paredes pegajosas**: la medida armónica G_surface(∞) cuantifica la contribución estacionaria de moléculas adsorbidas y depende sistemáticamente de la relación de aspecto R/L del cilindro.

---

## Estructura del repositorio

### Notebooks principales

| Notebook | Descripción | Método |
|----------|-------------|--------|
| [`TFM_V5.ipynb`](TFM_V5.ipynb) | Modelo principal: D(z) inhomogéneo, 100×100 modos, malla sinh, G(t) + análisis de truncamiento | Expansión modal + ARPACK |
| [`TFM_TwoRegionModel.ipynb`](TFM_TwoRegionModel.ipynb) | Modelo de dos regiones (50×50 modos, D₁/D₂ piecewise) | Ecuación trascendental |
| [`TFM_TwoRegionModel_z'.ipynb`](TFM_TwoRegionModel_z'.ipynb) | Dos regiones en coordenada desplazada z', heatmap de convergencia | Ecuación trascendental |
| [`DiagramaEfectosSuperficie.ipynb`](DiagramaEfectosSuperficie.ipynb) | Figura conceptual: capas moleculares + perfil D(z) | matplotlib |

### Barridos paramétricos — Modelo D(z) continuo

| Notebook | Parámetro | Valores |
|----------|-----------|---------|
| [`TFM_V5_LoopG_delta.ipynb`](TFM_V5_LoopG_delta.ipynb) | δ — ancho de transición tanh | {0.01, 0.1, 1, 10, 100} nm |
| [`TFM_V5_LoopG_kappa.ipynb`](TFM_V5_LoopG_kappa.ipynb) | κ — permeabilidad Robin BC | {0, 0.1, 1, 10, 100} nm/μs |
| [`TFM_V5_LoopG_lambda.ipynb`](TFM_V5_LoopG_lambda.ipynb) | λ — tasa de recuperación exponencial | {0.001, 0.01, 0.05, 0.1, 1, 10} nm⁻¹ |
| [`TFM_V5_LoopG_zc.ipynb`](TFM_V5_LoopG_zc.ipynb) | z_c — altura crítica de la zona superficial | {25, 75, 125, 250, 500, 1000} nm |

### Barridos paramétricos — Modelo de dos regiones (z')

| Notebook | Parámetro | Valores |
|----------|-----------|---------|
| [`TFM_TwoRegionModel_z'_LoopG_D1D2.ipynb`](TFM_TwoRegionModel_z'_LoopG_D1D2.ipynb) | D₁/D₂ | {0.001, 0.01, 0.1, 1} |
| [`TFM_TwoRegionModel_z'_LoopG_kappa.ipynb`](TFM_TwoRegionModel_z'_LoopG_kappa.ipynb) | κ | {0, 0.1, 1, 10, 100} nm/μs |
| [`TFM_TwoRegionModel_z'_LoopG_zc.ipynb`](TFM_TwoRegionModel_z'_LoopG_zc.ipynb) | z_c | {0, 25, 75, 125, 250, 500, 1000} nm |

### Versión adimensional

| Notebook | Descripción |
|----------|-------------|
| [`Adimensional/TFM_v5_adim.ipynb`](Adimensional/TFM_v5_adim.ipynb) | Modelo principal en variables adimensionales (L_ref=100 nm, D_ref=2×10⁹ nm²/s) |

### Solución Numérica

| Notebook | Descripción |
|----------|-------------|
| [`Solución Numérica/CurvaTeórica.ipynb`](Solución%20Numérica/CurvaTeórica.ipynb) | Curva analítica G(τ) ~ erfc + τ^{−3/2} como referencia de validación |
| [`Solución Numérica/SoluciónP.ipynb`](Solución%20Numérica/SoluciónP.ipynb) | Solución FEniCS/Crank–Nicolson (requiere FEniCS; elemento P₁) |
| [`Solución Numérica/EsquemaparaP.ipynb`](Solución%20Numérica/EsquemaparaP.ipynb) | Esquema FD+CN completo (Nr=5000, Nz=10000; descripción del método) |

### Sticky Walls Correlation

| Notebook | Descripción |
|----------|-------------|
| [`Sticky Walls Correlation/sticky-correlation-modal.ipynb`](Sticky%20Walls%20Correlation/sticky-correlation-modal.ipynb) | G_bulk + G_surface: modos de Bessel + Gauss–Legendre + Numba |
| [`Sticky Walls Correlation/G_Matrix_Loop.ipynb`](Sticky%20Walls%20Correlation/G_Matrix_Loop.ipynb) | Heatmap G_surface(∞) vs (R, L) en paralelo (joblib) |

---

## Métodos numéricos

| Componente | Técnica |
|------------|---------|
| Modos radiales | Raíces de la ecuación trascendental Robin vía `scipy.optimize.brentq` |
| Modos axiales (D(z) continuo) | Volumen finito sobre malla sinh-adaptativa (N_z = 90 000) → problema generalizado simétrico K Z = Λ M Z, simetrizado con A = M^(−1/2) K M^(−1/2) y diagonalizado con `scipy.linalg.eigh_tridiagonal` |
| Modos axiales (dos regiones) | Raíces de la ecuación trascendental analítica vía `brentq`, con rama **trigonométrica o hiperbólica** según el signo de (Λ² − α²D_i)/D_i (los modos superficiales tienen la región 2 evanescente) |
| Paralelización | `joblib.Parallel` sobre modos radiales n = 0..99 (4 núcleos) |
| Cuadratura modal | Gauss–Legendre (`numpy.polynomial.legendre.leggauss`) |
| Suma espectral G(t) | Numpy vectorizado; Numba `@njit(parallel=True)` en sticky-correlation |

---

## Instalación

```bash
pip install -r requirements.txt
```

Probado con Python 3.12, NumPy 2.5, SciPy 1.18 y pandas 3.

### Qué se ejecuta y qué no

Catorce de los dieciocho cuadernos corren de principio a fin sin errores. Los cuatro restantes:

| Cuaderno | Estado |
|---|---|
| `Solución Numérica/SoluciónP.ipynb` | requiere FEniCS, que no está en el entorno base de Anaconda |
| `Solución Numérica/EsquemaparaP.ipynb` | `MemoryError`: con `Nr = 5000`, `Nz = 10000` son 5×10⁷ incógnitas y la matriz de Crank–Nicolson pide ~2 GB sólo en las diagonales. Hay que bajar la resolución o pasar a un esquema ADI |
| `Sticky Walls Correlation/G_Matrix_Loop.ipynb` | barrido de 501×501 puntos (R, L) con 500×500 modos cada uno; no se ejecutó por coste |
| `Adimensional/TFM_TwoRegionModel_adim.ipynb` | está vacío (cero celdas) |

### Requisitos de memoria

`TFM_V5.ipynb` mantiene los 100×100 modos axiales sobre la malla de 90 000 nodos: **7.2 GB** de pico, y tarda unos 9 minutos con 4 núcleos. Los barridos `TFM_V5_LoopG_*` usan `N_z = 45 000` (≈3.6 GB) y tardan 15–18 minutos cada uno. Bajar `N_z` a 5 000 no cambia el espectro en 5 cifras y reduce el pico en más de un orden de magnitud.
