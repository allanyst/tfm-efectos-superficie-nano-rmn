# Efectos de Superficie en Resonancia Magnética Nuclear a Nanoescala

**Trabajo Fin de Máster (TFM)** — Máster Oficial en Ingeniería Matemática y Computación (PER 11087, Octubre 2024).

> **Título completo:** *Efectos de superficie en Resonancia Magnética Nuclear a Nanoescala*

---

## Contexto físico

Las moléculas de agua confinadas cerca de superficies sólidas experimentan una **reducción del coeficiente de difusión** D(z) debido a interacciones de van der Waals, puentes de hidrógeno e hidrofobicidad. Este fenómeno afecta directamente la señal de relajación en experimentos de NMR a nanoescala.

El modelo resuelve la **ecuación de difusión en coordenadas cilíndricas** con D = D(z) inhomogéneo:

$$\frac{\partial P}{\partial t} = \nabla \cdot [D(z)\, \nabla P]$$

con perfil de difusión de transición suave:
$$D(z) = D_0 \left[(1 - e^{-\lambda z})(1 - S(z)) + S(z)\right], \quad S(z) = \frac{1}{2}\left[1 + \tanh\!\left(\frac{z - z_c}{\delta}\right)\right]$$

La **función de correlación** G(t) — observable experimental de NMR — se obtiene integrando la densidad de probabilidad P(ρ, z, t) sobre el dominio cilíndrico.

---

## Estructura del repositorio

### Notebooks principales

| Archivo | Descripción |
|---------|-------------|
| [`TFM_V5.ipynb`](TFM_V5.ipynb) | Modelo principal en coordenadas físicas (nm); perfil D(z), modos de Bessel radiales, autovalores axiales, P(ρ,z,t) y G(t) |
| [`TFM_TwoRegionModel.ipynb`](TFM_TwoRegionModel.ipynb) | Extensión con dos regiones de difusión diferenciada |
| [`TFM_TwoRegionModel_z'.ipynb`](TFM_TwoRegionModel_z'.ipynb) | Modelo de dos regiones con cambio de variable z' |
| [`DiagramaEfectosSuperficie.ipynb`](DiagramaEfectosSuperficie.ipynb) | Figura conceptual: capas moleculares (adsorbida, transición, bulk) + perfil D(z) |

### Barridos paramétricos (Loop G)

| Archivo | Parámetro variado |
|---------|-----------------|
| [`TFM_V5_LoopG_delta.ipynb`](TFM_V5_LoopG_delta.ipynb) | δ — ancho de la zona de transición |
| [`TFM_V5_LoopG_kappa.ipynb`](TFM_V5_LoopG_kappa.ipynb) | κ — permeabilidad de la superficie (Robin BC) |
| [`TFM_V5_LoopG_lambda.ipynb`](TFM_V5_LoopG_lambda.ipynb) | λ — tasa de recuperación exponencial de D |
| [`TFM_V5_LoopG_zc.ipynb`](TFM_V5_LoopG_zc.ipynb) | z_c — altura crítica de la zona superficial |
| [`TFM_TwoRegionModel_z'_LoopG_*.ipynb`](.) | Equivalentes para el modelo de dos regiones |

### Versión adimensional

| Archivo | Descripción |
|---------|-------------|
| [`Adimensional/TFM_v5_adim.ipynb`](Adimensional/TFM_v5_adim.ipynb) | Modelo principal reescrito en variables adimensionales (L_ref=100 nm, D_ref=2×10⁹ nm²/s) |
| [`Adimensional/TFM_TwoRegionModel_adim.ipynb`](Adimensional/TFM_TwoRegionModel_adim.ipynb) | Modelo de dos regiones adimensional |

### Solución Numérica (exploratorio)

| Archivo | Descripción |
|---------|-------------|
| [`Solución Numérica/SoluciónP.ipynb`](Solución%20Numérica/SoluciónP.ipynb) | Solución numérica directa de la densidad P |
| [`Solución Numérica/EsquemaparaP.ipynb`](Solución%20Numérica/EsquemaparaP.ipynb) | Esquema de discretización para P |
| [`Solución Numérica/CurvaTeórica.ipynb`](Solución%20Numérica/CurvaTeórica.ipynb) | Curva teórica de referencia G(t) ~ t^{-3/2} |

### Sticky Walls Correlation

| Archivo | Descripción |
|---------|-------------|
| [`Sticky Walls Correlation/sticky-correlation-modal.ipynb`](Sticky%20Walls%20Correlation/sticky-correlation-modal.ipynb) | Modelo de correlación con paredes "pegajosas" (κ > 0) |

---

## Métodos numéricos

- **Modos radiales:** ceros de J₀ (funciones de Bessel de orden 0) via `scipy.special.jn_zeros`
- **Modos axiales:** problema de autovalores sobre malla sinh-adaptada (concentrada cerca de z=0) via `scipy.sparse.linalg.eigsh`
- **Paralelización:** `joblib.Parallel` sobre modos radiales n=0..15
- **Integración:** `scipy.integrate.quad` para coeficientes A_{nm}; regla trapezoidal no uniforme para G(t)

---

## Instalación

```bash
pip install -r requirements.txt
```
