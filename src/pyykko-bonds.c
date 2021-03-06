#include <gtk/gtk.h>
#include <GL/freeglut.h>
#include <math.h>

#include "gmol.h"

typedef struct {
  double sbond;
  double dbond;
  double tbond;
} pyykko_s;
#define pyykko_single(b)	(pyykkos[b].sbond)
#define pyykko_double(b)	(pyykkos[b].dbond)
#define pyykko_triple(b)	(pyykkos[b].tbond)

pyykko_s pyykkos[] = {
  {-1.00, -1.00, -1.00},		//   0 empty
  { 0.32, -1.00, -1.00},		//   1 H
  { 0.46, -1.00, -1.00},		//   2 He
  { 1.33,  1.24, -1.00},		//   3 Li
  { 1.02,  0.90,  0.85},		//   4 Be
  { 0.85,  0.78,  0.73},		//   5 B
  { 0.75,  0.67,  0.60},		//   6 C
  { 0.71,  0.60,  0.54},		//   7 N
  { 0.63,  0.57,  0.53},		//   8 O
  { 0.64,  0.59,  0.53},		//   9 F
  { 0.67,  0.96, -1.00},		//  10 Ne
  { 1.55,  1.60, -1.00},		//  11 Na
  { 1.39,  1.32,  1.27},		//  12 Mg
  { 1.26,  1.13,  1.11},		//  13 Al
  { 1.16,  1.07,  1.02},		//  14 Si
  { 1.11,  1.02,  0.94},		//  15 P
  { 1.03,  0.94,  0.85},		//  16 S
  { 0.99,  0.95,  0.93},		//  17 Cl
  { 0.96,  1.07,  0.96},		//  18 Ar
  { 1.96,  1.93, -1.00},		//  19 K
  { 1.71,  1.47,  1.33},		//  20 Ca
  { 1.48,  1.16,  1.14},		//  21 Sc
  { 1.36,  1.17,  1.08},		//  22 Ti
  { 1.34,  1.12,  1.06},		//  23 V
  { 1.22,  1.11,  1.03},		//  24 Cr
  { 1.19,  1.05,  1.03},		//  25 Mn
  { 1.16,  1.09,  1.02},		//  26 Fe
  { 1.11,  1.03,  0.96},		//  27 Co
  { 1.10,  1.01,  1.01},		//  28 Ni
  { 1.12,  1.15,  1.20},		//  29 Cu
  { 1.18,  1.20, -1.00},		//  30 Zn
  { 1.24,  1.17,  1.21},		//  31 Ga
  { 1.21,  1.11,  1.14},		//  32 Ge
  { 1.21,  1.14,  1.06},		//  33 As
  { 1.16,  1.07,  1.07},		//  34 Se
  { 1.14,  1.09,  1.10},		//  35 Br
  { 1.17,  1.21,  1.08},		//  36 Kr
  { 2.10,  2.02, -1.00},		//  37 Rb
  { 1.85,  1.57,  1.39},		//  38 Sr
  { 1.63,  1.30,  1.24},		//  39 Y
  { 1.54,  1.27,  1.21},		//  40 Zr
  { 1.47,  1.25,  1.16},		//  41 Nb
  { 1.38,  1.21,  1.13},		//  42 Mo
  { 1.28,  1.20,  1.10},		//  43 Tc
  { 1.25,  1.14,  1.03},		//  44 Ru
  { 1.25,  1.10,  1.06},		//  45 Rh
  { 1.20,  1.17,  1.12},		//  46 Pd
  { 1.28,  1.39,  1.37},		//  47 Ag
  { 1.36,  1.44, -1.00},		//  48 Cd
  { 1.42,  1.36,  1.46},		//  49 In
  { 1.40,  1.30,  1.32},		//  50 Sn
  { 1.40,  1.33,  1.27},		//  51 Sb
  { 1.36,  1.28,  1.21},		//  52 Te
  { 1.33,  1.29,  1.25},		//  53 I
  { 1.31,  1.35,  1.22},		//  54 Xe
  { 2.32,  2.09, -1.00},		//  55 Cs
  { 1.96,  1.61,  1.49},		//  56 Ba
  { 1.80,  1.39,  1.39},		//  57 La
  { 1.63,  1.37,  1.31},		//  58 Ce
  { 1.76,  1.38,  1.28},		//  59 Pr
  { 1.74,  1.37, -1.00},		//  60 Nd
  { 1.73,  1.35, -1.00},		//  61 Pm
  { 1.72,  1.34, -1.00},		//  62 Sm
  { 1.68,  1.34, -1.00},		//  63 Eu
  { 1.69,  1.35,  1.32},		//  64 Gd
  { 1.68,  1.35, -1.00},		//  65 Tb
  { 1.67,  1.33, -1.00},		//  66 Dy
  { 1.66,  1.33, -1.00},		//  67 Ho
  { 1.65,  1.33, -1.00},		//  68 Er
  { 1.64,  1.31, -1.00},		//  69 Tm
  { 1.70,  1.29, -1.00},		//  70 Yb
  { 1.62,  1.31,  1.31},		//  71 Lu
  { 1.52,  1.28,  1.22},		//  72 Hf
  { 1.46,  1.26,  1.19},		//  73 Ta
  { 1.37,  1.20,  1.15},		//  74 W
  { 1.31,  1.19,  1.10},		//  75 Re
  { 1.29,  1.16,  1.09},		//  76 Os
  { 1.22,  1.15,  1.07},		//  77 Ir
  { 1.23,  1.12,  1.10},		//  78 Pt
  { 1.24,  1.21,  1.23},		//  79 Au
  { 1.33,  1.42, -1.00},		//  80 Hg
  { 1.44,  1.42,  1.50},		//  8i Tl
  { 1.44,  1.35,  1.47},		//  82 Pb
  { 1.51,  1.41,  1.35},		//  83 Bi
  { 1.45,  1.35,  1.29},		//  84 Po
  { 1.47,  1.38,  1.38},		//  85 At
  { 1.42,  1.45,  1.33},		//  86 Rn
  { 2.23,  2.18, -1.00},		//  87 Fr
  { 2.01,  1.73,  1.59},		//  88 Ra
  { 1.86,  1.53,  1.40},		//  89 Ac
  { 1.75,  1.43,  1.36},		//  90 Th
  { 1.69,  1.38,  1.29},		//  91 Pa
  { 1.70,  1.34,  1.18},		//  92 U
  { 1.71,  1.36,  1.16},		//  93 Np
  { 1.72,  1.35, -1.00},		//  94 Pu
  { 1.66,  1.35, -1.00},		//  95 Am
  { 1.66,  1.36, -1.00},		//  96 Cm
  { 1.68,  1.39, -1.00},		//  97 Bk
  { 1.68,  1.40, -1.00},		//  98 Cf
  { 1.65,  1.40, -1.00},		//  99 Es
  { 1.67, -1.00, -1.00},		// 100 Fm
  { 1.73,  1.39, -1.00},		// 101 Md
  { 1.76, -1.00, -1.00},		// 102 No
  { 1.61,  1.41, -1.00},		// 103 Lr
  { 1.57,  1.40,  1.31},		// 104 Rf
  { 1.49,  1.36,  1.29},		// 105 Db
  { 1.43,  1.28,  1.21},		// 106 Sg
  { 1.41,  1.28,  1.19},		// 107 Bh
  { 1.34,  1.25,  1.18},		// 108 Hs
  { 1.29,  1.25,  1.13},		// 109 Mt
  { 1.28,  1.16,  1.12},		// 110 Ds
  { 1.21,  1.16,  1.18},		// 111 Rg
  { 1.22,  1.37,  1.30},		// 112 Cn
  { 1.36, -1.00, -1.00},		// 113 Nh
  { 1.43, -1.00, -1.00},		// 114 Fl
  { 1.62, -1.00, -1.00},		// 115 Mc
  { 1.75, -1.00, -1.00},		// 116 Lv
  { 1.65, -1.00, -1.00},		// 117 Ts
  { 1.57, -1.00, -1.00},		// 118 Og
};

//static size_t nr_py_elements = sizeof(pyykkos) / sizeof(pyykko_s);

#define MAGIC_NUMBER 	0.95		// fixme

bond_e
get_bond_type (atom_s *atoma, atom_s *atomb)
{
  bond_e bond = BOND_NONE;
  
  radii_s *rada = atom_data (atoma);
  radii_s *radb = atom_data (atomb);
  if (rada && radb) {
    double dx = atom_x (atoma) - atom_x (atomb);
    double dy = atom_y (atoma) - atom_y (atomb);
    double dz = atom_z (atoma) - atom_z (atomb);
    double sp = MAGIC_NUMBER * sqrt ((dx * dx) + (dy * dy) + (dz * dz));
    double cv_single =
      (pyykko_single (tp_nr (rada)) <= 0.0 ||
       pyykko_single (tp_nr (radb)) <= 0.0) ? -1.0 : 
      pyykko_single (tp_nr (rada)) + pyykko_single (tp_nr (radb));
    double cv_double =
      (pyykko_double (tp_nr (rada)) <= 0.0 ||
       pyykko_double (tp_nr (radb)) <= 0.0) ? -1.0 : 
      pyykko_double (tp_nr (rada)) + pyykko_double (tp_nr (radb));
    double cv_triple =
      (pyykko_triple (tp_nr (rada)) <= 0.0 ||
       pyykko_triple (tp_nr (radb)) <= 0.0) ? -1.0 : 
      pyykko_triple (tp_nr (rada)) + pyykko_triple (tp_nr (radb));
    if (sp < cv_single || sp < cv_double || sp < cv_triple) {
      bond = BOND_SINGLE;
      if (sp < cv_triple) bond = BOND_TRIPLE; 
      if (sp < cv_double) bond = BOND_DOUBLE;
    }
  }
  return bond;
}
