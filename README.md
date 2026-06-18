# TFG-MIGUEL-HABA-VERDEJO
Codis de MATLAB emprats en el meu TFG ''Estudi de l'aproximació de Born-Oppenheimer. Aplicacions'' 

## Arxius inclosos

### `separacio_variables.m`
**Descripció:**  
Programa principal que implementa el mètode de separació de variables per resoldre l'equació de Schrödinger de l'H₂⁺. Compara els resultats amb el mètode LCAO i estudia els estats vibracionals i rotacionals.

**Funcionament:**
- Defineix les constants físiques fonamentals (α, c, h, hc, me, a₀, E₁s)
- Calcula les energies LCAO per a diferents distàncies internuclears
- Per a cada distància, aplica el mètode de Newton-Raphson per trobar la constant de separació C
- Calcula les energies exactes mitjançant la fórmula E = αhc/(2R) - hc²C/(2meR²)
- Troba el mínim de l'energia i mostra els resultats per pantalla
- Genera un gràfic de convergència per a N=1,2,3,4
- Calcula la funció d'ona exacta en 2D
- Resol l'equació vibracional amb el mètode de Numerov per a v=0,1,2
- Calcula els nivells rotacionals per a L=0,1,2,3
- Genera diagrames d'energia amb els nivells vibracionals i rotacionals

**Variables principals:**
- `E_exacta`: energies calculades per separació de variables
- `C_opt`: constant de separació òptima per a cada distància
- `info`: estructura amb informació del mínim (posició, energia, constants)
- `resultats_num`: estructura amb els resultats del mètode Numerov

---

### `calcul_separacio_variables.m`
**Descripció:**  
Funció que implementa el nucli del mètode de separació de variables.

**Funcionament:**
- Per a cada distància internuclear, obté una estimació inicial de C a partir de l'energia LCAO
- Aplica el mètode de Newton-Raphson per resoldre l'equació frac_cont(C, B, N) = 0
- Calcula l'energia exacta corresponent a la C trobada
- Retorna el vector d'energies, les constants C òptimes i informació del mínim

**Variables principals:**
- `frac`: funció de fracció contínua que defineix la condició de quantització
- `B`: paràmetre que depèn de la distància internuclear (B = 4me·α·R/hc)
- `tol`: tolerància per al mètode de Newton-Raphson
- `max_iter`: nombre màxim d'iteracions

---

### `frac_cont.m`
**Descripció:**  
Funció que calcula la fracció contínua necessària per al mètode de separació de variables.

**Funcionament:**
- Rep els paràmetres C (constant de separació), B (paràmetre de distància) i N (ordre de l'aproximació)
- Calcula la constant A (valor propi angular) segons l'ordre N:
  - N=1: A = -C/3
  - N=2: fórmula analítica
  - N=3: diagonalització d'una matriu 3x3
  - N=4: diagonalització d'una matriu 4x4
- Utilitza l'algorisme de Lentz modificat per avaluar la fracció contínua
- La fracció contínua representa la condició de quantització radial

**Variables principals:**
- `A`: valor propi angular (constant de separació)
- `s`: paràmetre definit com s = B/(2√C) - 1
- `a(n)`, `b(n)`: coeficients de la fracció contínua
- `tiny`: valor petit per evitar divisions per zero

---

### `calcul_LCAO.m`
**Descripció:**  
Funció auxiliar que calcula les integrals necessàries per al mètode LCAO.

**Funcionament:**
- Rep com a paràmetres la distància internuclear (dosRa0), les constants físiques i l'energia E₁s
- Calcula l'integral de solapament S = e⁻ᴰ(1 + D + D²/3)
- Calcula l'element de matriu H₁₁ = E₁s + e⁻²ᴰ(1/D + 1)·αhc/a₀
- Calcula l'element de matriu H₁₂ = 2E₁s·e⁻ᴰ(1/2 - 1/D + 7D/6 + D²/6)
- Retorna l'energia de l'estat enllaçant: E_LCAO = (H₁₁ + H₁₂)/(1 + S)

**Variables principals:**
- `S`: integral de solapament
- `H11`, `H12`: elements de matriu hamiltoniana
- `E_LCAO`: energia de l'estat enllaçant

---

### `LCAO.m`
**Descripció:**  
Programa que implementa el mètode de combinació lineal d'orbitals atòmics (LCAO) per a l'H₂⁺. Calcula les energies dels estats enllaçant i antienllaçant, i representa gràficament les funcions d'ona.

**Funcionament:**
- Defineix les constants físiques fonamentals
- Crida la funció `calcul_LCAO` per obtenir les energies, l'integral de solapament S i els elements de matriu H₁₁ i H₁₂
- Calcula l'estat antienllaçant mitjançant la fórmula E_anti = (H₁₁ - H₁₂)/(1 - S)
- Troba el mínim de l'estat enllaçant i mostra els resultats
- Genera tres gràfics:
  1. Corbes d'energia per als estats enllaçant i antienllaçant
  2. Funció d'ona 2D per a l'estat enllaçant
  3. Tall 1D de la funció d'ona al llarg de l'eix internuclear

**Variables principals:**
- `dosRa0`: vector de distàncies internuclears en unitats de a₀
- `E_LCAO`: energies de l'estat enllaçant
- `E_anti`: energies de l'estat antienllaçant
- `Psi_enllacant`, `Psi_antienllacant`: funcions d'ona normalitzades

---

### `numerov.m`
**Descripció:**  
Implementa el mètode de Numerov per resoldre l'equació de Schrödinger vibracional.

**Funcionament:**
- Retalla el domini de la distància internuclear a l'interval [0.5, 4.5] a₀
- Construeix les matrius de Numerov (B i A_mat) per a la discretització
- Resol l'equació de Schrödinger: -K·B⁻¹·A·ψ + V·ψ = E·ψ
- Diagonalitza l'hamiltonià per obtenir els nivells vibracionals
- Per als estats rotacionals, afegeix el terme de barrera centrífuga V_centrifug = K·L(L+1)/r²
- Calcula els nivells d'energia per a L = 0, 1, 2, 3

**Variables principals:**
- `H`: hamiltonià discretitzat
- `En`: energies vibracionals
- `psi`: funcions d'ona vibracionals
- `Energies_vL`: matriu d'energies per a diferents (v, L)

---

### `funcio_ona_2D.m`
**Descripció:**  
Funció que calcula i representa la funció d'ona exacta en dues dimensions.

**Funcionament:**
- Calcula els coeficients angulars resolent el problema de valors propis:
  - Construeix la matriu M amb els elements d'acoblament angulars
  - Diagonalitza -M per obtenir els valors propis A₀
  - Selecciona l'autovector associat al valor propi mínim
- Calcula els coeficients radials:
  - Construeix la matriu A per al problema radial
  - Diagonalitza per obtenir els coeficients d
- Utilitza coordenades el·líptiques (ξ, η) per avaluar la funció d'ona
- Retorna la funció d'ona normalitzada i la informació dels coeficients

**Variables principals:**
- `c_coefs`: coeficients angulars
- `d_coefs`: coeficients radials
- `A0_val`: valor propi angular mínim
- `Psi_exacta`: funció d'ona exacta


