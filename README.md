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

y parámetros físicos: d = 15 nm, D₀ = 2×10³ nm²/μs, Dₛ = 2×10¹ nm²/μs, R = L = 50d = 750 nm.

La **función de correlación G(t)** — observable experimental de Nano-NMR — se calcula como:

$$G(t) = \frac{5}{16\pi} \sum_{n,m} C_{nm}^2\, e^{-\Lambda_{nm}^2\, t}$$

donde Λ²_nm son los autovalores del operador de difusión y C_nm los coeficientes de proyección del kernel dipolar.

---

## Resultados clave

- **Régimen de tiempos cortos** (τ ≪ 1): G(τ) ~ exp(−6τ), difusión libre antes de sentir la pared.
- **Régimen de tiempos largos** (τ ≫ 1): G(τ) ~ τ^{−3/2}, ley de potencias universal del confinamiento 3D.
- **Convergencia modal**: el error relativo entre truncamientos consecutivos es < 1% para ≥ 50 modos axiales y ≥ 50 modos radiales.
- **Sensibilidad paramétrica**: los barridos de (λ, δ, z_c, κ) muestran que la transición entre los dos regímenes asintóticos es sensible al perfil de interfaz, lo que permite inferir parámetros microscópicos a partir de datos experimentales de G(t).
- **Modelo de dos regiones**: alternativa analítica con autovalores transcendentales, 100× más rápida que el modelo de perfil continuo y cualitativamente equivalente.
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
| Modos axiales (D(z) continuo) | Diagonalización ARPACK (`eigsh`) de matriz FD tridiagonal sobre malla sinh-adaptativa (N_z = 90 000) |
| Modos axiales (dos regiones) | Raíces de ecuación trascendental analítica (β₁, β₂) vía `brentq` |
| Paralelización | `joblib.Parallel` sobre modos radiales n = 0..99 (4 núcleos) |
| Cuadratura modal | Gauss–Legendre (`numpy.polynomial.legendre.leggauss`) |
| Suma espectral G(t) | Numpy vectorizado; Numba `@njit(parallel=True)` en sticky-correlation |

---

## Instalación

```bash
pip install -r requirements.txt
```

> **Nota:** `Solución Numérica/SoluciónP.ipynb` requiere FEniCS, que no está disponible en el entorno base de Anaconda. Los demás notebooks funcionan con el entorno estándar.
