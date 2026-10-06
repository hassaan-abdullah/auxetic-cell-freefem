# Auxetic Unit Cell: Inverse Homogenisation and Large-Strain Response in FreeFEM

![FreeFEM](https://img.shields.io/badge/FreeFEM-4.13%2B-blue) ![License: MIT](https://img.shields.io/badge/License-MIT-green)

A FreeFEM/Python workflow for designing a 2D periodic unit cell using **topology optimisation** to have negative effective Poisson's ratio, i.e. as an *auxetic* material, and then checking whether the optimised small-strain behaviour survives thresholding, finite strain, finite specimen size, and instability. Everything is implemented from the weak form in [FreeFEM](https://freefem.org) [1].

## 1. Motivation

Auxetic lattices contract laterally when subjected to compression [2, 3], making them attractive for impact-absorbing structures like cushioning layers in footwear. The design of these lattices is typically done under two idealisations: **small strain**, where the cell is optimised for its linearly elastic homogenised stiffness [5, 6], and **infinite lattice**, where the homogenisation assumes a periodic medium with infinitely many cells [7].

A real component typically works at 10–30 % compression and is only a few cells thick, while being bonded to stiffer parts; at such strains, instabilities, and large rotations change the lattice's response [8,9]. 

This project quantifies what survives each step away from the idealisations, i.e.:
- **Printability:** the optimiser's grey density is thresholded to a black-and-white geometry.
- **Large strain:** the cell is loaded far beyond the linear range.
- **Finite size:** the lattice becomes a block of a few cells between two plates.

Each stage is verified against a known or independent answer before its results are used.

## 2. Periodic Homogenisation: The Physics

Assume cell $Y = [0,1]^2$ to be one period of an infinite lattice, where, for an average (macroscopic) strain of $\boldsymbol\varepsilon^0$, the displacement in the cell would be

$$\mathbf u(\mathbf x) = \boldsymbol\varepsilon^0\mathbf x + \boldsymbol\chi(\mathbf x),$$

where $\boldsymbol\chi$ is a **periodic fluctuation**. $\boldsymbol\chi$ is periodic because every cell in an infinite lattice deforms identically, so the left edge of one cell must move exactly like the right edge of its neighbour [10]. 

One can use an analogy of claccial mechanics to a better understand this concept. If a truss is under uniform macroscopic strain, the joints follow the macroscopic strain on average, but each bay would also have a local (repeating) sag. Therefore, the fluctuation is that sag.

The equilibrium of the fluctuation, written in weak form, is to find periodic $\boldsymbol\chi$ such that for every periodic test field $\mathbf v$,

$$\int_Y \mathbf C\,\boldsymbol\varepsilon(\boldsymbol\chi) : \boldsymbol\varepsilon(\mathbf v)\,dY = -\int_Y \mathbf C\,\boldsymbol\varepsilon^0 : \boldsymbol\varepsilon(\mathbf v)\,dY .$$


The effective stiffness is then the average strain energy of the total strain [5, 11]:

$$Q_{ij} = \int_Y \big(\boldsymbol\varepsilon^0_i + \boldsymbol\varepsilon(\boldsymbol\chi_i)\big) : \mathbf C : \big(\boldsymbol\varepsilon^0_j + \boldsymbol\varepsilon(\boldsymbol\chi_j)\big)\,dY \qquad (|Y| = 1).$$

The three unit load cases are
- $\boldsymbol\varepsilon^0_1 = (1,0,0)$ (stretch in $x$), 
- $\boldsymbol\varepsilon^0_2 = (0,1,0)$ (stretch in $y$), and
- $\boldsymbol\varepsilon^0_3 = (0,0,1)$ (engineering shear $\gamma_{xy} = 1$).

For an uniaxial stress, the effective Poisson's ratio is $\nu_{yx} = Q_{12}/Q_{11}$, where a negative value means the lattice is auxetic.

### 2.1 Element Choice
Fluctuations are discretised with periodic P2 elements, which represent bending of thin walls accurately.
### 2.2 Rigid Translation
A periodic field is defined only up to a rigid translation, so a $10^{-9}$ mass term removes it without affecting strains.
### 2.3 Base Material
**$E_0 = 1$, $\nu_0 = 0.5$:** the plane-stress limit of the incompressible rubber is used at large strain, so all stages share one material. Plane stress is free of volumetric locking at $\nu_0 = 0.5$, because the material can thin out of plane.

## 3. Code/Script Reference

Every script documents its purpose, inputs, outputs and parameters in its header. Parameters are given on the command line as `-name value`; without them the verified defaults are used.

### 3.1 `freefem/check_solid.edp` 

**Purpose:** verify the periodic homogenisation on a problem with a known answer. A completely solid cell is just the base material, so its homogenised plane-stress (or effective) stiffness must be 
- $Q_{11} = 4/3$, 
- $Q_{12} = 2/3$, and 
- $Q_{33} = 1/3$. 

For a homogeneous cell the load term vanishes and the fluctuation is zero, so the answer is exact on any mesh. The check therefore tests the bookkeeping, i.e. the sign of the load term, engineering versus tensor shear, the periodic pairing, the energy formula, rather than the discretisation.

**Output:**
```
Q11 = 1.33333   (expected 1.33333)
Q12 = 0.666667   (expected 0.666667)
Q33 = 0.333333   (expected 0.333333)
```

## 4. Verification and testing

Codes at each stage are checked against an independent or known answer before its results are used.

| Check | Expected | Obtained |
|---|---|---|
| `check_solid.edp`: solid cell returns the base material | $Q_{11}, Q_{12}, Q_{33} = 4/3, 2/3, 1/3$ | 1.33333, 0.666667, 0.333333 |

## 5. Limitations

- **2D Plane Stress:** real cells are 3D and plane stress is a model of a thin sheet.
- **Normalised Units:** stresses are relative to the base material's Young's modulus.

## 6. References
1. Hecht F. 2012 New development in FreeFem++. *J. Numer. Math.* **20**, 251–265. (doi:10.1515/jnum-2012-0013)
2. Lakes R. 1987 Foam structures with a negative Poisson's ratio. *Science* **235**, 1038–1040. (doi:10.1126/science.235.4792.1038)
3. Evans KE, Nkansah MA, Hutchinson IJ, Rogers SC. 1991 Molecular network design. *Nature* **353**, 124. (doi:10.1038/353124a0)
4. Bertoldi K, Vitelli V, Christensen J, van Hecke M. 2017 Flexible mechanical metamaterials. *Nat. Rev. Mater.* **2**, 17066. (doi:10.1038/natrevmats.2017.66)
5. Sigmund O. 1994 Materials with prescribed constitutive parameters: an inverse homogenization problem. *Int. J. Solids Struct.* **31**, 2313–2329. (doi:10.1016/0020-7683(94)90154-6)
6. Xia L, Breitkopf P. 2015 Design of materials using topology optimization and energy-based homogenization approach in Matlab. *Struct. Multidiscip. Optim.* **52**, 1229–1241. (doi:10.1007/s00158-015-1294-0)
7. Bensoussan A, Lions J-L, Papanicolaou G. 1978 *Asymptotic analysis for periodic structures*. Amsterdam, The Netherlands: North-Holland.
8. Clausen A, Wang F, Jensen JS, Sigmund O, Lewis JA. 2015 Topology optimized architectures with programmable Poisson's ratio over large deformations. *Adv. Mater.* **27**, 5523–5527. (doi:10.1002/adma.201502485)
9. Bertoldi K, Reis PM, Willshaw S, Mullin T. 2010 Negative Poisson's ratio behavior induced by an elastic instability. *Adv. Mater.* **22**, 361–366. (doi:10.1002/adma.200901956)
10. Andreassen E, Andreasen CS. 2014 How to determine composite material properties using numerical homogenization. *Comput. Mater. Sci.* **83**, 488–495. (doi:10.1016/j.commatsci.2013.09.006)
11. Hill R. 1963 Elastic properties of reinforced solids: some theoretical principles. *J. Mech. Phys. Solids* **11**, 357–372. (doi:10.1016/0022-5096(63)90036-X)

## 7. License, citation and acknowledgement

Released under the MIT License (see `LICENSE`). If you use this code, please cite the repository:

> Abdullah, S. H. (2026). *Auxetic Unit Cell: Inverse Homogenisation and Large-Strain Stability in FreeFEM* [Computer software]. GitHub.

Code development was assisted by an AI coding assistant. The problem formulation, the verification runs, and the interpretation of the results are the author's.