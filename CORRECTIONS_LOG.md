# Log de correcciones — efectos-superficie-nano-RMN

Revisión física y matemática completa. Cada entrada indica el archivo, la celda,
el error y la corrección.

---

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

**Implementación (Cell 35):**
```python
R_rho = R0_mat * rho_nodes[None, :]           # (n_modes, Nr)
B_nj = np.dot(R_rho, f_mat)                   # (n_modes, Nz)
F_nm_correct = np.zeros((n_modes_radial, m_max))
for _n in range(n_modes_radial):
    F_nm_correct[_n, :] = dr * dz * B_nj[_n, :] @ Z_tens[:, _n, :]
B_rms_sq = (5 / (16 * np.pi)) * np.sum(F_nm_correct**2)
```

---

## TFM_TwoRegionModel.ipynb

### Bug 8 — Cell 24 [code]: Z2 con doble suma de d
**Error:** `Z2 = (z_nodes + d)**2` cuando z_nodes ya empieza en d.
**Corrección:** `Z2 = z_nodes**2`
**Impacto:** Factor dipolar f~(3z²−ρ²)/r⁵ evaluado con z mínimo 2d=30 nm en vez de d=15 nm.

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
**Nota:** Este cuaderno usa Dirichlet BC (ceros de J₀) porque es un prototipo con geometría
diferente (R=1000 nm, L=100 nm, z₀=10). La G(t) usa aproximación de un solo integral —
no se actualizó a F_nm² ya que las conclusiones lo señalan como "prototipo exploratorio".

---

## Cuadernos sin cambios (revisados y correctos)

| Cuaderno | Estado |
|---|---|
| `Solución Numérica/CurvaTeórica.ipynb` | Curva analítica de referencia — sin G(t) |
| `Solución Numérica/EsquemaparaP.ipynb` | Esquema de discretización — sin G(t) |
| `Solución Numérica/SoluciónP.ipynb` | Solucionador Crank-Nicolson — sin G(t) |
| `Sticky Walls Correlation/sticky-correlation-modal.ipynb` | Usa fórmula correcta I_sin²+I_cos² — Dirichlet intencional (paredes absorbentes) |
| `DiagramaEfectosSuperficie.ipynb` | Solo visualización |
| `*_LoopG_*.ipynb` | Barren parámetros sobre G(t) — sin bugs Z2/f² |
