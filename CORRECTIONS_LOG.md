# Log de correcciones — efectos-superficie-nano-RMN

Revisión física y matemática. Cada entrada indica el archivo, la celda, el error
y la corrección.

---

# Segunda revisión (2026-09-21)

Revisión centrada en el operador de difusión discreto, el kernel dipolar y la
ejecutabilidad de los cuadernos. Los hallazgos están ordenados por impacto.

## A. El solucionador axial de `TFM_V5.ipynb` no era fiable

### A1 — Celda 19: índice equivocado en la subdiagonal (`lower_diag[-2]`)

**Error.** El ensamblado hacía

```python
lo[jp] = -D_prev/(h_m*h_j)        # el bucle interior llena lo[0..N-3]
...
lower_diag[-2] = -D_3_2/hL/h_3_2  # frontera superior
```

`lower_diag` tiene longitud `N-1` (índices `0..N-2`) y, con `offsets=[0,1,-1]`,
`A[N-1, N-2]` lo da `lower_diag[N-2]`, es decir `lower_diag[-1]`. El bucle
interior sólo llega a `lower_diag[N-3]`, así que **`lower_diag[N-2]` nunca se
asignaba y quedaba en cero**: el último nodo quedaba desacoplado del resto, como
un bloque 1×1 con autovalor espurio `D/h_L²`.

**Consecuencia medible.** El operador dejaba de conservar la masa. Con
`alpha = 0` y `varkappa = 0`, `d/dz(D(z) d/dz ·)` debe anular la función
constante para cualquier `D(z)`; en cambio

```
||A @ 1||_inf = 5.69e+04   (el único nodo que falla es el N-1)
```

**Corrección.** Ensamblado en volumen finito con la subdiagonal completa; el
residuo pasa a `2.8e-09` (relativo `2e-17`).

### A2 — Celda 19: las filas de frontera usaban celdas de ancho `h_0` en vez de `h_0/2`

Las celdas de los extremos son **medias** celdas. El balance de flujos correcto es

```
(2 D_{1/2}/h_0²)(Z_0 - Z_1) + (2κ/h_0) Z_0 + α² D_0 Z_0 = Λ Z_0
```

y el cuaderno tenía el factor 2 ausente en ambos términos. Es una discretización
**inconsistente** de la condición de contorno (no converge al límite continuo:
para un modo `cos(kπ(z-d)/L)` la fila daba `Λ/2` en vez de `Λ`), y además reducía
la permeabilidad efectiva de Robin a `κ/2` en todos los barridos de κ.

### A3 — Celda 19: `varkappa` aplicado también en `z = d+L`

El código sumaba `varkappa/hL` en la fila superior. El propio markdown del
cuaderno dice que ahí va **Neumann**: `z = d+L` es la truncación artificial del
líquido semiinfinito y es reflejante. Sólo la superficie (`z = d`) lleva Robin.

### A4 — Celda 19: la matriz no era simétrica en malla no uniforme

`M⁻¹K` tiene `-D_{j+1/2}/(h_j w_j)` arriba y `-D_{j+1/2}/(h_j w_{j+1})` abajo:
en la malla sinh no son iguales. Se pasaba a `eigsh`, que asume simetría.

### A5 — Celda 19: `eigsh(..., sigma=0.0)` sobre una matriz singular

Con `varkappa = 0` y `alpha_0 = 0` el modo uniforme tiene autovalor exactamente
cero, de modo que la factorización *shift-invert* en `sigma = 0` es singular.
Combinado con A1–A4, **el solucionador no convergía**. Dos comprobaciones que no
necesitan ninguna solución analítica:

*El mismo operador, con distintos desplazamientos, daba espectros distintos:*

| `sigma` | 5 autovalores más bajos |
|---|---|
| `0` | 0.006627, 0.191691, 0.467544, 0.816077, 1.293817 |
| `-1e-3` | 0.008129, 0.173615, 0.204620, 0.453166, 0.726121 |
| `-1` | −0.009139, 0.042120, 0.212444, 0.538582, 0.729174 |
| `-10` | −0.359763, 0.000213, 0.158922, 0.359161, 0.403911 |
| volumen finito | **0, 0.035081, 0.140188, 0.314919, 0.558607** |

*Y no se estabilizaba al refinar la malla:*

| `N_z` | cuaderno | volumen finito |
|---|---|---|
| 5 000 | 0.009812, 0.050641, 0.100912, 0.407782 | 0, 0.0350807, 0.1401880, 0.3149181 |
| 10 000 | 0.008009, 0.130267, 0.348652, 0.442452 | 0, 0.0350807, 0.1401881, 0.3149187 |
| 20 000 | 0.007312, 0.200117, 0.222952, 0.622856 | 0, 0.0350807, 0.1401881, 0.3149188 |
| 40 000 | 0.008628, 0.091794, 0.192920, 0.416316 | 0, 0.0350808, 0.1401881, 0.3149189 |

Como comprobación adicional, en el caso degenerado `D(z) = D₀` constante y
`α = 0` —donde la respuesta se conoce, `Λ_m = D₀(mπ/L)²`— el cuaderno erraba
un 125–306 % y el volumen finito un 0.0004 %.

### Corrección A (celda 19 reescrita)

Discretización en **volumen finito** sobre la celda `[z_{j-1/2}, z_{j+1/2}]` de
ancho `w_j`, que da el problema generalizado **simétrico** `K Z = Λ M Z` con
`M = diag(w_j)`:

```
K[j,j+1] = K[j+1,j] = -D_{j+1/2}/h_j
K[j,j]   = D_{j-1/2}/h_{j-1} + D_{j+1/2}/h_j + α² D_j w_j
w_0 = h_0/2 (Robin: + κ)      w_{N-1} = h_{N-2}/2 (Neumann)
```

Se simetriza con `A = M^{-1/2} K M^{-1/2}` y se diagonaliza con
`scipy.linalg.eigh_tridiagonal`, sin ARPACK. Los pesos `w_j` son exactamente los
de la regla del trapecio, así que la ortonormalidad discreta coincide con
`∫ Z² dz = 1`.

La celda imprime ahora `Λ_(0,0)`, que debe ser 0: es la comprobación que detecta
cualquier error de ensamblado en las fronteras.

La misma corrección se aplicó a `build_system` en los cuatro cuadernos
`TFM_V5_LoopG_*.ipynb`.

---

## B. El kernel dipolar estaba mal en 11 cuadernos

**Error.** `f_mat = (3*Z2 - R2) / (R2 + Z2)**2.5`, es decir
`f = (3z² − ρ²)/r⁵`.

**Problema.** El factor dipolar es el armónico esférico `Y_2^0`:

```
f = (3cos²θ − 1)/r³ ,  cosθ = z/r ,  r² = ρ² + z²
  = (3z² − r²)/r⁵
  = (2z² − ρ²)/(ρ² + z²)^{5/2}
```

La resta va contra `r²`, no contra `ρ²`. Tres indicios coincidentes:

1. El prefactor `5/(16π)` de `G(t)` es `[(1/4)√(5/π)]²`, el coeficiente de
   `Y_2^0` al cuadrado, así que `f` tiene que ser ese mismo armónico.
2. El markdown escribía `Y_2^0 = (1/4)√(5/π)(cos²θ − 1)`, al que le falta el 3.
3. **El propio repositorio ya lo tenía bien** en
   `Sticky Walls Correlation/sticky-correlation-modal.ipynb`:
   `W = ρ[3z²/r⁵ − 1/r³]`, que es exactamente `ρ(2z² − ρ²)/r⁵`.

El kernel del cuaderno se anulaba en θ = 60° en vez de en el ángulo mágico
θ = 54.7°, y valía `3/d³` en vez de `2/d³` sobre el eje.

**Corrección.** `(2*Z2 - R2)` en `TFM_V5.ipynb`, `TFM_TwoRegionModel.ipynb`,
`TFM_TwoRegionModel_z'.ipynb`, `Adimensional/TFM_v5_adim.ipynb` y los 7
cuadernos `*_LoopG_*.ipynb`; markdown de `TFM_V5.ipynb` corregido.

---

## C. El modelo de dos regiones omitía sus modos más lentos

**Error.** En `find_m_roots` / `find_lambda_roots`:

```python
Lambda_min = alpha * np.sqrt(D2)
```

y en la ecuación trascendental

```python
beta_1_sq = np.abs(Lambda**2 - alpha_n**2 * D1) / D1
beta_2_sq = np.abs(Lambda**2 - alpha_n**2 * D2) / D2
```

**Problema.** El espectro de `-(D Z')' + α²D Z` está acotado por abajo por
`α² min(D) = α² D1`, no por `α² D2`. En la banda `α²D1 ≤ Λ² < α²D2` la región 2
es **evanescente** (`β₂` imaginario, soluciones hiperbólicas): son los modos
localizados en la capa superficial lenta. El `np.abs()` los convertía en
trigonométricos (rama equivocada) y el barrido ni siquiera llegaba a esa banda.

**Consecuencia medible** (contra volumen finito con `D` constante a trozos):

| n | Λ² más bajos, volumen finito | Λ² más bajos, cuaderno | omitidos |
|---|---|---|---|
| 1 | **0.03050**, 0.05494, 0.09357 | 0.05494, 0.09357, 0.20354 | 1 |
| 2 | **0.03471, 0.17401**, 0.20919 | 0.20919, 0.29567, 0.35468 | 2 |
| 3 | **0.03736, 0.29388**, 0.37303 | 0.36800, 0.41609, 0.52560 | 2 |
| 4 | **0.04034, 0.30933**, 0.62918 | 0.63118, 0.66414, 0.76849 | 3 |
| 5 | **0.04388, 0.31636, 0.83309** | 0.96456, 1.01936, 1.12817 | 3 |

Los modos que faltaban son los **más lentos**, es decir los que dominan la cola
de `G(t)` a tiempos largos — justo el régimen donde se verifica la ley `τ^{-3/2}`.
Para n ≥ 3 la primera «raíz» del cuaderno era además `α²D2` exacto, un artefacto
del punto singular donde arranca el barrido.

**Corrección.** Ecuación trascendental con las dos ramas, escrita como el empalme
de `Z` y del flujo `D Z'` en la interfaz:

```
F(Λ²) = D₁ Z₁'(z_c) g(z_c) − D₂ g'(z_c) Z₁(z_c)
```

con `Z₁` trigonométrica o hiperbólica según el signo de `(Λ² − α²D₁)/D₁`, y `g`
igual en la región 2, normalizada por `cosh(q(L−z_c))` para evitar desbordes.
El barrido arranca en `α²D₁` y añade `Λ² = 0` explícitamente cuando `α = 0` y
`κ = 0` (ahí `F` se anula sin cambiar de signo).

**Validación** contra volumen finito con `N_z = 200 000`:

- error relativo de los autovalores: `≤ 8.6e-05`
- salto de la autofunción en la interfaz: `≤ 2.9e-08`
- ortonormalidad, `|Gram − I|`: `≤ 3.3e-09`

Aplicado a `TFM_TwoRegionModel.ipynb`, `TFM_TwoRegionModel_z'.ipynb` y los tres
`TFM_TwoRegionModel_z'_LoopG_*.ipynb` (conservando en cada uno las firmas que ya
usaban sus bucles).

---

## D. Cuadernos que no llegaban a ejecutarse

### D1 — `np.trapz` no existe en NumPy ≥ 2.0

Se usaba en 8 cuadernos. El entorno tiene NumPy 2.5.1. Sustituido por
`np.trapezoid` (que el resto del repositorio ya usaba).

### D2 — `TFM_V5_LoopG_kappa.ipynb`: celda 10 rota

Faltaba el bloque `Parallel(...)` de resolución axial, dejando
`eigenvectors_ax[n] = evecs` con sangría suelta → `IndentationError`. Además
`C_nm`, `rho0` y `z0` no estaban definidos. El cuaderno no corría ni una celda.
Restaurado a partir de la estructura de `TFM_V5_LoopG_delta.ipynb`.

### D3 — `TFM_V5_LoopG_lambda.ipynb`: `C_nm`, `rho0`, `z0` sin definir

Quedaba el bloque de la condición inicial tipo delta (`A_nm = R0_pt ⊗ Z0`) de la
versión antigua, mientras `G_values` ya usaba `C_nm`. Sustituido por la
proyección dipolar `C_nm = ∫ f R_n Z_nm ρ dρ dz`.

### D4 — `TFM_V5_LoopG_zc.ipynb`: `del A_nm`, `del R0_pt`, `del Z_tens`

Nombres que ya no existían: `NameError` al final de la primera iteración del
barrido, después de todo el cálculo. Eliminados / corregidos a `Z_tens_int`.

### D5 — Memoria: 29 GB en splines

`TFM_V5.ipynb` construía `n_modes_radial × m_max = 10 000` objetos
`CubicSpline` sobre 90 000 nodos. Cada uno guarda `4 × N_z` coeficientes:
**28.8 GB**. Los barridos hacían lo mismo con `N_z = 45 000` (14 GB).

Sustituido por **un spline vectorial por modo radial** (sus `m_max` modos a la
vez), evaluado y descartado: 288 MB de pico y 3× más rápido. Resultado numérico
idéntico.

También se eliminaron dos duplicaciones de 7.2 GB:

- el apilado `np.array([[...] for ...])`, que mantiene viva a la vez la lista y
  el resultado → ahora se copia modo a modo liberando el diccionario;
- `np.trapezoid(eigenvectors_ax**2, z, axis=2)`, un temporal del tamaño del
  array completo → ahora la normalización sólo se **comprueba**, por bloques.

### D6 — `eigenvectors_ax.astype(np.float32)` en la norma

Perdía ~7 cifras significativas al normalizar los modos. Eliminado (era un
paliativo para el problema de memoria de D5). La desviación de `∫Z² dz` respecto
a 1 es ahora `4.0e-15`.

---

## E. El kernel dipolar usaba z′ en vez de la distancia al NV

En la familia de cuadernos en coordenada desplazada `z′ = z − d` (origen en la
superficie), el factor dipolar se construía con `z_nodes`, que corre sobre
`[0, L]`:

```python
z_nodes = np.linspace(0.0, L, N_quad)
Z2 = (z_nodes + 1e-30)**2   # "z_nodes empieza en d; no se suma d de nuevo"
```

El comentario es falso en estos cuadernos: la distancia NV→molécula es `z′ + d`,
no `z′`. El `+ 1e-30` estaba precisamente para tapar la división por cero en
`(ρ, z′) = (0, 0)`. Consecuencias:

- `TFM_TwoRegionModel_z'.ipynb` producía directamente `B_rms² = nan`
  (`0/0` en el origen);
- los tres barridos en `z′` daban un valor finito pero muy equivocado: cerca de
  la superficie `f ~ 2/z′³` se dispara y domina toda la integral, como si el
  NV estuviese pegado al líquido en vez de a 15 nm.

**Corrección:** `Z2 = (z_nodes + d)**2` en los cuatro cuadernos en `z′`.
(En `TFM_TwoRegionModel.ipynb` y en la familia `TFM_V5` la malla ya corre sobre
`[d, d+L]`, así que ahí `Z2 = z_nodes**2` es lo correcto.)

---

## F. La cuadratura de F_nm no resolvía el integrando

**Error.** Malla uniforme de 1000 nodos y **suma rectangular**:

```python
dr = rho_nodes[1] - rho_nodes[0]
dz = z_nodes[1]  - z_nodes[0]
F_nm[_n, :] = dr * dz * B_nj[_n, :] @ Z_tens[:, _n, :]
```

**Problema.** El integrando está concentrado cerca de la superficie: `f` decae
como `z⁻³` y los modos axiales tienen capa límite en `z = d`. Con `Δz = 0.75 nm`
y regla rectangular, el error de cuadratura es grande justo donde está toda la
señal.

**Consecuencia medible.** Con la identidad de Parseval ahora disponible, la suma
modal truncada **no puede superar** a la integral directa. Medido sobre la salida
de `TFM_V5.ipynb`:

```
(5/16π) ∫ f² ρ dρ dz  = 3.964e-06
(5/16π) Σ F_nm²       = 4.420e-06      ->  111.5 %   (imposible)
```

Con una cuadratura adecuada (malla en ρ refinada con pesos del trapecio y la
malla sinh en z) el mismo cálculo en una réplica reducida de 40×40 modos da un
**87 %**, que sí es un defecto de truncamiento.

**Resultado tras corregir.** El control de Parseval pasa a dar valores admisibles,
y además los dos modelos —que no comparten nada de maquinaria— coinciden:

| Cuaderno | modos | capturado | (5/16π)∫f² ρ dρ dz |
|---|---|---|---|
| `TFM_V5.ipynb` (volumen finito, D(z) continuo) | 100×100 | 99.93 % | 3.682878e-06 |
| `TFM_TwoRegionModel_z'.ipynb` (raíces trascendentales) | 100×100 | 99.93 % | 3.683072e-06 |
| `TFM_TwoRegionModel.ipynb` | 50×50 | 94.80 % | 3.683072e-06 |

**Corrección.** En los 11 cuadernos que calculan `F_nm`:

- `N_quad` de 1000 a 2000;
- `z` con la misma transformación sinh de la malla axial, agrupando nodos en la
  superficie;
- **pesos de la regla del trapecio** (`w_rho`, `w_z`) en lugar de `dr·dz`, lo que
  además corrige los medios pesos de los extremos;
- se añade el **control de Parseval** impreso en la propia celda:

```
(5/16pi) * int f^2 rho drho dz = ...
capturado por 100x100 modos = ... %  (debe ser <= 100 %)
```

---

## G. Incompatibilidades con el entorno

| Problema | Cuadernos | Corrección |
|---|---|---|
| `np.trapz` eliminado en NumPy ≥ 2.0 | 8 | `np.trapezoid` |
| `pd.read_csv(delim_whitespace=True)` eliminado en pandas 3 | 7 | `sep=r'\s+'` |
| `varkappa = 0*L/D_0` → `NameError` | `Adimensional/TFM_v5_adim.ipynb` | `D0` (el valor es 0 de todas formas) |

---

## H. Errores de dimensiones y de unidades

| Archivo | Antes | Ahora |
|---|---|---|
| `TFM_V5.ipynb`, celda 4 | `z_c = lambda_param*L` (nm⁻¹·nm = adimensional) | `z_c = 2.5*d` — mismo valor, 37.5 nm |
| `TFM_TwoRegionModel*.ipynb` | `z_c = 0.75*R/d` / `0.75/d*L` (adimensional) | `z_c = 2.5*d` — mismo valor |
| `TFM_V5.ipynb`, celda 15 | `α_n [nm⁻²]` | `α_n (nm⁻¹)` |
| `TFM_V5.ipynb`, celda 21 | `Autovalor (s⁻¹)` | `Λ_{n,m} (µs⁻¹)` |
| `TFM_TwoRegionModel*.ipynb` | `Λ²_{n,m} (s⁻¹)` | `Λ²_{n,m} (µs⁻¹)` |
| varios | `varkappa # nm/s` | `nm/µs` (coherente con `D` en nm²/µs) |
| `TFM_V5.ipynb`, celda 4 | `R = L = 50*d # 1e-6 m -> 1500 nm` | `50 d = 750 nm` |
| `TFM_V5.ipynb`, celda 4 | `delta = 2*d/3 # 1 nm` | `# 10 nm` |
| `TFM_V5.ipynb`, celda 12 | «Promedio armónico de 1/D(z) en [0, L]» | «Promedio armónico de D(z) sobre [d, d+L]» |
| `TFM_V5.ipynb`, celda 8 | docstring partido a la mitad (`D0 [ex\np(-λz)...`) | fórmula completa y correcta |

---

## I. Correcciones de markdown

### I1 — `TFM_V5.ipynb`, celda 33: `Y_2^0`

`Y_2^0(θ) = (1/4)√(5/π)(cos²θ − 1)` → `(3cos²θ − 1)`. Añadida la reducción
explícita a coordenadas cilíndricas y una nota sobre por qué el numerador es
`2z² − ρ²`. Se eliminó además el texto obsoleto sobre subdominios
(`ρ ∈ [0.1R, 0.9R]`, `z ∈ [7,13] nm`, malla no uniforme en z para la cuadratura),
que no correspondía a lo que hace el código, y el párrafo residual sobre la
condición inicial delta.

### I2 — `TFM_V5.ipynb`, celda 18: discretización axial

La versión anterior tenía `h_{j+1}²` y `h_j²` donde van `h·w`, escribía
`a_0 = 2D_0/h_0² − 2κ/h_0` (signo de κ invertido, y `D_0` en vez de `D_{1/2}`),
y afirmaba que la matriz era simétrica. Reescrita como derivación en volumen
finito, con las dos condiciones de contorno y la simetrización explícitas.

### I3 — Notación `Λ` vs `Λ²`

El markdown definía `A Z = Λ² Z` mientras el código usaba `exp(-Λ t)` con `Λ` el
autovalor. Unificado a `Λ_{n,m}` (autovalor, unidades µs⁻¹) en las celdas 18, 24,
28 y 33 de `TFM_V5.ipynb`.

### I4 — `TFM_V5.ipynb`, celda 26

`M_{n,m} = √(∫₀^L Z² dz)` → `√(∫_d^{d+L} Z² dz)`, que es el dominio real.

### I5 — `TFM_TwoRegionModel.ipynb`, celda 16: dominio de normalización

`z_nodes = np.linspace(0.0, L, N_z)` mientras el dominio físico es `[d, d+L]` y
`fz_raw` usa `z-d`. La norma axial se calculaba sobre el intervalo equivocado.
Corregido a `np.linspace(Z_OFF, Z_OFF+L, N_z)` con `Z_OFF = d`.
(La variante `_z'` sí estaba bien: ahí la coordenada es `z' = z − d`.)

---

## J. Segunda pasada: cuadernos que había revisado por encima

Los cuadernos auxiliares se habían mirado sólo en las partes que se tocaron.
Al leerlos enteros, y al abrir el único módulo `.py` del repositorio, apareció
lo siguiente.

### J1 — `sticky-correlation-modal.ipynb` tenía una copia obsoleta de su propio módulo

`Sticky Walls Correlation/time_dependent_accumulation_surface.py` —que es lo que
importa `G_Matrix_Loop.ipynb`— **ya trae la versión correcta**:
`precompute_surface_modes` devuelve `F_even` y `F_odd` por separado, cada uno con
su autovalor (`lambda2_even`, `lambda2_odd`), y `sum_surface_modes_numeric`
aplica a cada familia su propio factor temporal.

El cuaderno, en cambio, tenía inline una versión anterior que **suma las dos
familias en un único peso** `F[m,p]` y luego le aplica sólo `lambda2_o`:

```python
F_caps = (B_delta / (J1*J1)) * (C_sin*Is + C_cos*Ic)   # par e impar juntos
F_wall = (beta / J1) * (Is*Zs/lam2_e + Ic*Zc/lam2_o)
...
acc += coeff * (1.0 - math.exp(-D*lam2*t))             # un solo lambda2_o
```

Las dos familias tienen autovalores distintos
(`eta_e = 2mπ/L` frente a `eta_o = (2m−1)π/L`), así que no pueden agruparse.

**Corrección:** el cuaderno se alineó con el módulo (celdas 27, 29 y 31), lo que
además incorpora `-expm1(-x)` en lugar de `1 - exp(-x)` para evitar la
cancelación cuando `x ≪ 1`. El markdown de la celda 28 se reescribió con las dos
familias y con el signo del prefactor que usa el código (`−8π²/V`, no `+`).

*Esto resuelve lo que la revisión anterior había dejado como pendiente: no hacía
falta el manuscrito, la respuesta estaba en el propio repositorio.*

### J2 — `Adimensional/TFM_v5_adim.ipynb`: α mal adimensionalizado

```python
alpha_adim = alpha_zeros/(R*L_ref)  # α' = α*L_ref
```

El comentario dice lo correcto y el código hace otra cosa. `R` ya es
adimensional (`R = R_fisico/L_ref = 10`), de modo que dividir entre
`R*L_ref = 1000 nm` devuelve α **en nm⁻¹**, no adimensional. En la diagonal
conviven entonces `D_adim/h_adim²` (adimensional) y `alpha**2 * D_adim` (nm⁻²):
el acoplamiento radial quedaba **10⁴ veces por debajo** de lo debido.

**Corrección:** `alpha_adim = alpha_zeros / R`. También se corrigió el signo del
término de Robin (`− 2κ/h` → `+ 2κ/h`) en las dos celdas.

**No corregido, documentado en el propio cuaderno:** los nodos internos usan
`D/h²` en vez de `D/(h·w_j)` (la malla es bizona sinh, así que el operador no es
el de volumen finito y la matriz no es simétrica); el borde superior lleva
`D_L/h_L²` mientras el inferior lleva `2D_0/h_0²`; la norma axial usa peso `1/D`
cuando la forma de Sturm–Liouville de este problema tiene `r(x)=1`; y la
condición radial es de Dirichlet, no de Neumann, así que **no** resuelve el mismo
problema que `TFM_V5.ipynb` y sus autovalores no sirven para validarlo, pese a lo
que afirman sus conclusiones.

### J3 — `Solución Numérica/SoluciónP.ipynb`: falta el jacobiano cilíndrico

La forma débil se monta sobre una malla **rectangular** en (ρ, z) con `dx`, sin
multiplicar por `x[0]`:

```python
F_time  = (u - u0)/dt * v * dx
F_diff  = 0.5*D*dot(grad(u), grad(v))*dx + ...
G[n]    = assemble(K * u * dx)
```

En coordenadas cilíndricas axisimétricas `dV = 2π ρ dρ dz`, de modo que tal como
está resuelve la difusión **cartesiana** en 2-D, no la cilíndrica: el término
`(1/ρ)∂_ρ(ρ D ∂_ρ·)` no está representado. Lo mismo vale para la integral de
correlación, a la que le falta el peso `ρ`.

Además: la celda 6 usa `r_edges` y `z_edges`, que **nunca se definen** en este
cuaderno (vienen de `EsquemaparaP.ipynb`) → `NameError`; y la condición de Robin
se aplica a las tres fronteras (`ds(1)+ds(2)+ds(3)`), incluida la superior, que
en el modelo es Neumann.

**No corregido a propósito:** FEniCS no está en el entorno, así que cualquier
cambio aquí sería código que no puedo ejecutar ni verificar. Sólo se le corrigió
el kernel dipolar de la `Expression`, que es una sustitución literal.

### J4 — `Solución Numérica/EsquemaparaP.ipynb`: falta dividir entre ρ

En el stencil interior:

```python
a_r_p = (rho_i + dr/2)*D_r/(dr**2)
a_r_m = (rho_i - dr/2)*D_r/(dr**2)
```

La forma conservativa es `(1/ρ_i)[ρ_{i+1/2}D(P_{i+1}−P_i) − ρ_{i−1/2}D(P_i−P_{i−1})]/dr²`:
**falta la división entre `ρ_i`**. Tal como está, el operador radial es ρ veces el
correcto, y además queda inconsistente con los términos axiales, que no llevan
factor ρ. Los términos axiales usan `D/h²` en vez de `D/(h·w_j)` sobre la malla
sinh, el mismo defecto del punto A2.

**No corregido:** el cuaderno no llega a ejecutarse (`MemoryError`, ver J4), de
modo que es código muerto hasta que se baje la resolución.

### J5 — `z_c` tiene dos valores distintos en el repositorio

| Cuadernos | L | z_c | D_eff |
|---|---|---|---|
| `TFM_V5`, modelo de dos regiones y sus barridos | 50 d = 750 nm | 2.5 d = 37.5 nm | 1792.6 |
| `CurvaTeórica`, `SoluciónP`, `EsquemaparaP` | 100 d = 1500 nm | 0.75 L = 1125 nm | 1884.3 |

El impacto sobre `D_eff` es pequeño (5 %) porque con λ⁻¹ = 20 nm la exponencial
ya satura mucho antes de `z_c`, pero `CurvaTeórica.ipynb` se presenta en el README
como «referencia de validación» de `TFM_V5.ipynb` y está resolviendo otra
geometría y otro perfil. O se unifican los parámetros o se deja de llamarla
referencia.

### J6 — Comprobado y correcto

Cosas que se revisaron en esta pasada y **no** tienen problema:

- Las etiquetas de columna de los `.dat` (`G_del1`, `G_lam50`, `G_k10`…) son
  enteros escalados a propósito (`int(delta*100)`, `int(lam*1000)`, `int(kp*10)`)
  para que sean identificadores válidos de pandas; las leyendas de las gráficas
  mapean de vuelta a los valores reales correctamente.
- Los barridos sí son sensibles a sus parámetros. Dispersión máxima de
  `G/G(0)` entre curvas: z_c 63 %, δ 68 %, κ 272 %, λ 457 %. La afirmación del
  README sobre sensibilidad paramétrica se sostiene.
- El markdown de los siete cuadernos de barrido (dos celdas cada uno) usa ya las
  fórmulas corregidas de `D(z)` y de la condición de contorno.
- `time_dependent_accumulation_surface.py` está completo y es sintácticamente
  correcto; su `W_rz_vec` usa el kernel dipolar correcto.

---

## K. Pendientes — no corregidos, requieren decisión del autor

### K1 — Los dos asintóticos «teóricos» se contradicen entre cuadernos

`TFM_V5.ipynb` dibuja `exp(-6τ)` como teoría a tiempos cortos y lo repite en las
conclusiones. `Solución Numérica/CurvaTeórica.ipynb` concluye que a tiempos
cortos `G(τ)` **diverge** como `τ^{-3/2}`. Las dos afirmaciones no pueden ser
ambas la referencia del mismo `G(τ)`; conviene fijar cuál es la convención y de
dónde sale el 6.

### K2 — Truncamiento modal y `G(0)`

Con la fórmula `G(0) = (5/16π) Σ F_nm²`, la identidad de Parseval da el error de
truncamiento directamente: en una réplica reducida (40×40 modos) la suma captura
un **87 %** de `(5/16π)∫f² ρ dρ dz`. El análisis de convergencia del cuaderno
compara curvas **normalizadas** `G/G(0)`, que es una prueba mucho más débil.
Conviene reportar también la fracción de `B_rms²` capturada.

### K3 — Cuadernos que siguen sin ejecutarse

Tres cuadernos quedan sin ejecutar, por razones ajenas a esta revisión:

| Cuaderno | Motivo |
|---|---|
| `Solución Numérica/SoluciónP.ipynb` | requiere FEniCS, que no está en el entorno (ya documentado en el README). Se le corrigió el kernel dipolar de la `Expression`. |
| `Solución Numérica/EsquemaparaP.ipynb` | `MemoryError` en la celda 7: con `Nr = 5000`, `Nz = 10000` son 5.0×10⁷ incógnitas y la matriz de Crank–Nicolson pide ~2 GB sólo en las diagonales. Habría que bajar la resolución o pasar a un esquema ADI / dirección alternada. Se le corrigió el kernel dipolar (código y markdown). |
| `Sticky Walls Correlation/G_Matrix_Loop.ipynb` | el barrido es de 501×501 puntos (R, L), cada uno con 500×500 modos y 500×500 nodos de Gauss. No se modificó y no se intentó ejecutar por coste. |

Además, `Adimensional/TFM_TwoRegionModel_adim.ipynb` está **vacío** (cero celdas):
o se borra o se completa, pero tal cual sólo añade ruido al repositorio.

### K4 — `N_z = 90 000` es innecesario

El esquema de volumen finito ya da 5 cifras estables con `N_z = 5 000`
(ver la tabla de A5). Bajar `N_z` reduciría el pico de 7.2 GB en más de un orden
de magnitud.

---
---

# Primera revisión

## TFM_V5.ipynb (modelo principal)

### Bug 1 — Cell 7 [markdown]: D(z) sin D_s ni desplazamiento d
**Error:** La fórmula mostraba
```
D(z) = D₀ [(1−e^{−λz})(1−S) + S],   S = tanh((z−z_c)/δ)
```
**Problema:** (a) no hay D_s; el valor en la superficie sería D(d)=D₀·(1−e^{−λd})·(1−S(d))+S(d) ≠ D_s; (b) el exponente usa z en vez de z−d (el dominio del fluido empieza en z=d, no z=0); (c) S usa z−z_c en vez de z−(z_c+d).

**Corrección:**
```
D(z) = D_s + (D₀−D_s)[(1−e^{−λ(z−d)})(1−S) + S],   S = tanh((z−(z_c+d))/δ)
```
D(d)=D_s exactamente; D→D₀ para z≫z_c+d.

---

### Bug 2 — Cell 11 [markdown]: BC radial Dirichlet en vez de Neumann
**Error:** Se aplicaba R(R)=0 (pared absorbente), obteniendo α_n = j_{0,n}/R (ceros de J₀).

**Problema:** La pared lateral del cilindro es **reflejante** (no hay flujo de espines a través de ella); la condición física es de no-flujo: dR/dρ|_{ρ=R}=0.

**Corrección:** BC Neumann → J₁(αR)=0 → α_n = j_{1,n}/R (ceros de J₁). El factor
dipolar f(ρ,z) no se anula en ρ=R con BC Neumann, así que la diferencia es significativa.

---

### Bug 3 — Cell 18 [markdown]: D_{j+1/2} media aritmética → media armónica
**Error 3a:** La fórmula mostraba D_{j+1/2} = (D(z_{j+1})−D(z_j))/2 (una diferencia, no un promedio).

**Error 3b:** Typo en la fórmula de diferencias finitas: Z_{j+1}+Z_j en vez de Z_{j+1}−Z_j.

**Corrección:**
```
D_{j+1/2} = 2 D_j D_{j+1} / (D_j + D_{j+1})
```
La media armónica conserva el flujo J=−D dZ/dz en cada interfaz entre nodos,
que es la discretización consistente para d/dz[D(z) d/dz].

---

### Bug 4 — Cell 26 [markdown]: norma radial J₁²(j_{0,n})
**Error:** Se usaba N_n = sqrt(R²/2 · J₁²(j_{0,n})) para BC Dirichlet.

**Corrección (BC Neumann):**
```
N_n = sqrt(R²/2 · [J₀²(α_n R) + J₁²(α_n R)])
```
Con J₁(α_n R)=0 (BC Neumann), esto se simplifica a N_n = (R/√2)|J₀(α_n R)|.

---

### Bug 5 — Cell 34 [code]: factor dipolar f_mat·f_mat en vez de f₀·f_mat
**Error:** `rho_factor = rho_nodes * f_mat * f_mat * dr * dz` — se usaba f² en la integral.

**Problema:** La función de correlación es G(t)=⟨B(t)B(0)⟩; en la aproximación
de condición inicial delta en r₀=(0,d), la integral correcta es f(r₀)·∫f(r)P(r,t|r₀)d³r,
no ∫f²(r)P d³r.

**Corrección:** `f0 = (3z₀²−ρ₀²)/(ρ₀²+z₀²)^{5/2} = 3/d³`; `rho_factor = rho_nodes * f0 * f_mat * dr * dz`.

---

### Bug 6 — Cell 34 [code]: doble suma de d en f₀ (z₀+d)
**Error:** `f0 = (3*(z0+d)**2 − rho0**2) / ((z0+d)**2 + rho0**2)**2.5`

**Problema:** `z_nodes` ya empieza en d (el NV está en z=0 y el fluido en [d, d+L]);
así que z₀=d ya es la distancia NV→molécula. Sumar d de nuevo da z=2d=30 nm en vez de d=15 nm,
reduciendo f por factor (1/2)³=1/8 y f² por 1/64.

**Corrección:** `f0 = (3*z0**2 − rho0**2) / (rho0**2 + z0**2)**2.5`

---

### Bug 7 — Cell 33 + 35 [markdown + code]: G(t) no garantiza G(0)=B_rms²
**Error:** La fórmula G(t) = (5/16π) Σ_i Σ_j f(rᵢ,zⱼ) P(rᵢ,zⱼ,t) ρᵢ Δρ wⱼ
(integral simple sobre r con condición inicial delta en r₀) no satisface G(0)=B_rms² porque
en t→0 P→δ(r−r₀)/(2πρ₀) y la integral converge a f(r₀)² ≠ ∫f²d³r = B_rms².

**Corrección — fórmula modal exacta:**
```
G(t) = (5/16π) Σ_{n,m} F_{nm}² exp(−Λ_{nm} t)

F_{nm} = ∫₀^R ∫_d^{d+L} f(ρ,z) R_n(ρ) Z_{nm}(z) ρ dρ dz
```
Esta fórmula garantiza **G(0) = B_rms² = (5/16π) Σ_{n,m} F_{nm}²** exactamente,
porque es la expansión de Parseval de B_rms² = (5/16π)∫f²d³r en la base {R_n ⊗ Z_nm}.

---

## TFM_TwoRegionModel.ipynb

### Bug 8 — Cell 24 [code]: Z2 con doble suma de d
**Error:** `Z2 = (z_nodes + d)**2` cuando z_nodes ya empieza en d.
**Corrección:** `Z2 = z_nodes**2`

### Bug 9 — Cell 25 [code]: G(t) sin fórmula F_nm²
**Igual que Bug 7.** Reemplazado con la misma implementación modal.

---

## TFM_TwoRegionModel_z'.ipynb

### Bug 10 — Cell 24 [code]: Z2 con doble d
Mismo que Bug 8. Corregido: `Z2 = z_nodes**2`.

### Bug 11 — Cell 25 [code]: G(t) sin fórmula F_nm²
Mismo que Bug 9. Corregido con fórmula modal.

---

## Adimensional/TFM_v5_adim.ipynb (prototipo exploratorio)

### Bug 12 — [code]: D(z) sin D_s
**Error:** La línea activa era `return D0 * (...)`, la correcta con D_s estaba comentada.
**Corrección:** Se activó `return D_s+(D0-D_s) * (...)` (D_s=D0·0.05 ya estaba definido).
