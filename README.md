# TFG-MIGUEL-HABA-VERDEJO
Codis de MATLAB emprats en el meu TFG ''Estudi de l'aproximació de Born-Oppenheimer. Aplicacions'' 

## Arxius inclosos

### `separacio_variables.m`
**Descripció:**  
Programa principal que implementa el mètode de separació de variables per resoldre l'equació de Schrödinger de l'H₂⁺. Compara els resultats amb el mètode LCAO i estudia els estats vibracionals i rotacionals.

**Funcionament:**
- Calcula les energies LCAO per a diferents distàncies internuclears
- Per a cada distància, aplica el mètode de Newton-Raphson per trobar la constant de separació C
- Calcula les energies exactes mitjançant la fórmula E = αhc/(2R) - hc²C/(2meR²)
- Troba el mínim de l'energia i mostra els resultats per pantalla
- Genera un gràfic de convergència per a N=1,2,3,4
- Calcula la funció d'ona exacta en 2D
- Resol l'equació vibracional amb el mètode de Numerov per a v=0,1,2
- Calcula els nivells rotacionals per a L=0,1,2,3
- Genera diagrames d'energia amb els nivells vibracionals i rotacionals

---

### `calcul_separacio_variables.m`
**Descripció:**  
Funció que implementa el mètode de separació de variables.

**Funcionament:**
- Per a cada distància internuclear, obté una estimació inicial de C a partir de l'energia LCAO
- Aplica el mètode de Newton-Raphson per resoldre l'equació frac_cont(C, B, N) = 0
- Calcula l'energia exacta corresponent a la C trobada
- Retorna el vector d'energies, les constants C i informació del mínim

---

### `frac_cont.m`
**Descripció:**  
Funció que calcula la fracció contínua necessària per al mètode de separació de variables.

**Funcionament:**
- Rep els paràmetres C, B i N (ordre de l'aproximació)
- Calcula la constant A segons l'ordre N
- Utilitza l'algorisme de Lentz modificat per avaluar la fracció contínua

---

### `calcul_LCAO.m`
**Descripció:**  
Funció auxiliar que calcula les integrals necessàries per al mètode LCAO.

**Funcionament:**
- Rep com a paràmetres la distància internuclear (dosRa0), les constants físiques i l'energia E_1s

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
  2. Funció d'ona 2D per als estats enllaçant i antienllaçant
  3. Tall 1D de la funció d'ona al llarg de l'eix internuclear

---

### `numerov.m`
**Descripció:**  
Implementa el mètode de Numerov per resoldre l'equació de Schrödinger vibracional.

**Funcionament:**
- Retalla el domini de la distància internuclear a l'interval [0.5, 4.5] a₀
- Construeix les matrius de Numerov (B i A_mat) per a la discretització
- Resol l'equació de Schrödinger
- Diagonalitza l'hamiltonià per obtenir els nivells vibracionals
- Per als estats rotacionals, afegeix el terme de barrera centrífuga V_centrifug
- Calcula els nivells d'energia per a L = 0, 1, 2, 3

---

### `funcio_ona_2D.m`
**Descripció:**  
Funció que calcula i representa la funció d'ona exacta en dues dimensions.

**Funcionament:**
- Calcula els coeficients de l'equació polar
- Calcula els coeficients radials
- Utilitza coordenades el·líptiques (ξ, η) per avaluar la funció d'ona
- Retorna la funció d'ona normalitzada

---

### `estat_senar.m`
**Descripció:**  
Aquest codi implementa el mètode de separació de variables per resoldre l'equació de Schrödinger de l'H₂⁺ utilitzant exclusivament moments angulars senars (l = 1, 3, 5, ...). L'estimació inicial per a la constant de separació C s'obté a partir de l'estat antienllaçant del mètode LCAO.

**Funcionament:**
- Calcula les energies LCAO per a l'estat antienllaçant
- Per a cada distància internuclear, obté una estimació inicial de C a partir de l'energia LCAO antienllaçant
- Aplica el mètode de Newton-Raphson per resoldre l'equació frac_cont_imparell(C, B, N) = 0
- Calcula les energies exactes per a N=1,2,3,4
- Genera un gràfic de convergència que mostra les corbes d'energia per a cada N
- Compara els resultats amb la corba LCAO antienllaçant i l'energia E₁s

---

### `frac_cont_imparell.m` (funció interna)
**Descripció:**  Funció equivalent a frac_cont.m però en lloc de l=0,2,4,3 es té l=1,3,5,7.
