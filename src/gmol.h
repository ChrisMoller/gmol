#ifndef GMOL_H
#define GMOL_H

typedef struct {
  int atomic_nr;
  char *element_symbol;
  double atomic_radius;
  double ionic_radius;
  double covalent_radius;
  double van_der_waals_radius;
  double crystal_radius;
} radii_s;
#define tp_nr(p)  ((p)->atomic_nr)
#define tp_sym(p) ((p)->element_symbol)
#define tp_ar(p)  ((p)->atomic_radius)
#define tp_cvr(p) ((p)->covalent_radius)

typedef struct {
  double   x;
  double   y;
  double   z;
  char    *name;
  radii_s *data;
} atom_s;
#define atom_x(a) ((a)->x)
#define atom_y(a) ((a)->y)
#define atom_z(a) ((a)->z)
#define atom_name(a) ((a)->name)
#define atom_data(a) ((a)->data)

typedef enum {
  BOND_NONE,
  BOND_SINGLE,
  BOND_DOUBLE,
  BOND_TRIPLE,
} bond_e;

typedef struct {
  atom_s *a;
  atom_s *b;
  bond_e  type;
  double  strength;
  double  separation;
} bond_s;
#define bond_a(p)		((p)->a)
#define bond_b(p) 		((p)->b)
#define bond_type(p)		((p)->type)
#define bond_strength(p)	((p)->strength)
#define bond_separation(p)	((p)->separation)

typedef struct {
  GSList *atoms;
  GSList *bonds;
  char   *comment;
  char   *file;
  char   *shortname;
  GtkGLArea *area;
  GtkWidget *window;
  gdouble zoom;
  gdouble latitude;
  gdouble longitude;
  gdouble voff;
  gdouble hoff;
} molecule_s;
#define molecule_atoms(m)	((m)->atoms)
#define molecule_bonds(m)	((m)->bonds)
#define molecule_comment(m)	((m)->comment)
#define molecule_file(m)	((m)->file)
#define molecule_shortname(m)	((m)->shortname)
#define molecule_area(m)	((m)->area)
#define molecule_window(m)	((m)->window)
#define molecule_zoom(m)	((m)->zoom)
#define molecule_latitude(m)	((m)->latitude)
#define molecule_longitude(m)	((m)->longitude)
#define molecule_voff(m)	((m)->voff)
#define molecule_hoff(m)	((m)->hoff)

#endif /* GMOL_H */
