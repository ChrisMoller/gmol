#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#define _GNU_SOURCE
#include <gtk/gtk.h>
#include <glib/gi18n-lib.h>
#include <GL/freeglut.h>
#include <FTGL/ftgl.h>
#include <wand/MagickWand.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <values.h>


#include "gmol.h"
#include "radii.h"
#include "pyykko-bonds.h"
#include "jmol-colours.h"
#include "cvbond.h"

typedef struct {
  double theta;
  double x;
  double y;
  double z;
} rotation_s;

static void open_molecule_window (gpointer data, gpointer user_data);

static GtkWidget	*window		= NULL;
static GSList 		*molecules	= NULL;

#define AR_SCALE	0.5
#define INITIAL_ZOOM	1.0
#define INITIAL_LATITUDE	0.0
#define INITIAL_LONGITUDE	0.0

#define deg_to_rad(d) (((d) * M_PI) / 180.0)
#define rad_to_deg(r) (((r) * 180.0) / M_PI)

// for ftgl 
static FTGLfont *font[3] = {NULL, NULL, NULL};
static int fontindex = 1;
#define FONT_FILE	"/usr/share/fonts/dejavu/DejaVuSerif.ttf"

GtkListStore *store = NULL;
enum {
  COLUMN_NAME,
  COLUMN_COMMENT,
  COLUMN_MOLECULE,
  NR_COLUMNS
};

static void inline
cross_product (rotation_s *r, rotation_s *a, rotation_s *b)
{
  g_assert (r && a && b);
  /***b
      | i  j  k  |
      | ax ay az |
      | bx by bz |

      (ay * bz - az * by)i
      (az * bx - ax * bz)j
      (ax * by - ay * bx)k
  ***/
  r->x = (a->y * b->z   -  a->z * b->y);
  r->y = (a->z * b->x   -  a->x * b->z);
  r->z = (a->x * b->y   -  a->y * b->x);
}

static double inline
dot_product (rotation_s *a, rotation_s *b)
{
  /***
      (a . b) / |a| * |b| = cos (theta)
      (dx * 0 + dy * 0 + dz * s) / (s * s) = cos (theta)
  ***/

  double dotp = a->x * b->x   +   a->y * b->y  +  a->z * b->z;
  double maga = sqrt (a->x * a->x  +  a->y * a->y  +  a->z * a->z);
  double magb = sqrt (b->x * b->x  +  b->y * b->y  +  b->z * b->z);
  double theta = (maga > 0.0 && magb > 0.0) ?
    acos (dotp / (maga * magb)) : 0.0;
  return theta;
}

static void
gl_realise (GtkGLArea *area, gpointer user_data)
{
  GLfloat white[]       = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat black[]       = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat direction[]   = { 1.0, 1.0, 1.0, 0.0 };

  gtk_gl_area_make_current (area);

  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, white);
  glMaterialf(GL_FRONT,  GL_SHININESS, 5);
  glColorMaterial(GL_FRONT, GL_DIFFUSE); //    for FTGL
  

  glLightfv(GL_LIGHT0, GL_AMBIENT, black);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
  glLightfv(GL_LIGHT0, GL_POSITION, direction);

  glEnable(GL_LIGHTING);                // so the renderer considers light
  glEnable(GL_LIGHT0);                  // turn LIGHT0 on
  glEnable(GL_DEPTH_TEST);              // so the renderer considers depth
  glEnable(GL_COLOR_MATERIAL);
}

static void
gl_resize (GtkGLArea *area,
	   gint width,
	   gint height,
	   gpointer user_data)
{
  glViewport (0, 0, width, height);
  glMatrixMode (GL_PROJECTION);

  GLdouble aspect = (GLdouble)width / (GLdouble)height;
  glLoadIdentity ();

#define ORTHO_SCALE 2.5
  if (width <= height) {
    // width is smaller, so stretch out the height
    glOrtho (-ORTHO_SCALE, ORTHO_SCALE,
	     -ORTHO_SCALE/aspect, ORTHO_SCALE/aspect,
	     -10.0, 10.0);
  } else {
    // height is smaller, so stretch out the width
    glOrtho (-ORTHO_SCALE*aspect, ORTHO_SCALE*aspect,
	     -ORTHO_SCALE, ORTHO_SCALE,
	     -10.0, 10.0);
  }

}

static void
draw_bond (gpointer data, gpointer user_data)
{
  bond_s *bond = data;

  atom_s *atoma = bond_a (bond);
  atom_s *atomb = bond_b (bond);
  glPushMatrix();
  glTranslated(atom_x (atoma), 
	       atom_y (atoma), 
	       atom_z (atoma));
  
  double dx = atom_x (atomb) - atom_x (atoma);	// calculate
  double dy = atom_y (atomb) - atom_y (atoma);	// the distance
  double dz = atom_z (atomb) - atom_z (atoma);

  /***
      rotational axis is the cross product between the axis 0f
      the default cylinder ([0, 0, s]) and the desired cylinder.
      rotation angle is from the dot product.
  ***/

  rotation_s dist = {0.0, dx, dy, dz};
  rotation_s cyl  = {0.0, 0.0, 0.0, bond_separation (bond)};
  rotation_s cp;

  double theta = rad_to_deg (dot_product (&dist, &cyl));
  cross_product (&cp, &dist, &cyl);
  cp.theta = theta;
  glRotated (-cp.theta, cp.x, cp.y, cp.z);

  double br = 0.0, bg = 0.0, bb = 0.0;
  
  switch (bond_type (bond)) {
  case BOND_NONE:	// should nrver happen
    br = 0.0;
    bg = 0.0;
    bb = 0.0;
    break;
  case BOND_SINGLE:
    br = bond_strength (bond);
    bg = 0.0;
    bb = 0.0;
    break;
  case BOND_DOUBLE:
    br = 0.0;
    bg = bond_strength (bond);
    bb = 0.0;
    break;
  case BOND_TRIPLE:
    br = 0.0;
    bg = 0.0;
    bb = 1.0 + bond_strength (bond);
    break;
  }
  glColor3d (br, bg, bb);

#if 1		// glut version doesn't work
  glutSolidCylinder (0.1,		// radius
		      bond_separation (bond),	// height
		     30,			// slices
		     30);			// stacks
#else
  GLUquadric *quad = gluNewQuadric ();
  gluCylinder (quad,
	       0.1,		// base radius
	       0.1,		// top radius
	       s,			// height
	       30,			// slices
	       30);			// stacks
  gluDeleteQuadric (quad);
#endif  
  glPopMatrix();
}

static void
draw_atom (gpointer data, gpointer user_data)
{
  atom_s *atom = data;
  rotation_s *rotation = user_data;
  double red = 0.0;
  double green = 1.0;
  double blue = 0.0;
  radii_s *rad = atom_data (atom);

  get_jmol_colour (atom, &red, &green, &blue);
  
  glColor3d (red, green, blue);
  
  glPushMatrix();
  glTranslated (atom_x (atom), 
		atom_y (atom), 
		atom_z (atom));
  glutSolidSphere (rad ? (AR_SCALE * tp_ar (rad)) : 0.0, 30, 30);

#define FONT_SCALE 0.004
  
  glColor3d (0.0, 1.0, 1.0);
  glRotated (rotation->theta, rotation->x, rotation->y, rotation->z);
  glTranslated (0.0, 0.0, AR_SCALE * tp_ar (rad));
  glScaled (FONT_SCALE,  FONT_SCALE,  FONT_SCALE);

  float bounds[6] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
  ftglGetFontBBox (font[fontindex], atom_name (atom), -1, bounds);

  glTranslatef ((bounds[0] - bounds[3]) / 2.0f,
		(bounds[1] - bounds[4]) / 2.0f,
		0.0f);
  ftglRenderFont (font[fontindex], atom_name (atom), FTGL_RENDER_FRONT);
  glPopMatrix();
}

static void
draw_molecule (molecule_s *molecule, rotation_s *rotation)
{
  if (molecule_atoms (molecule))
    g_slist_foreach (molecule_atoms (molecule), draw_atom, rotation);
  if (molecule_bonds (molecule))
    g_slist_foreach (molecule_bonds (molecule), draw_bond, NULL);
}

static gboolean
gl_render (GtkGLArea *area,
	   GdkGLContext *context,
	   gpointer      user_data)
{
  molecule_s *molecule = user_data;
  glClearColor (0, 0, 0, 0);
  glClearDepth (1.0);

  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode (GL_MODELVIEW);
  glPushMatrix ();


  
  
  double zoom		= molecule_zoom (molecule);
  double latitude	= deg_to_rad (molecule_latitude (molecule));
  double longitude	= deg_to_rad (molecule_longitude (molecule));
  double voff		= molecule_voff (molecule);
  double hoff		= molecule_hoff (molecule);
  double dx = sin (longitude) * cos (latitude);
  double dy = sin (latitude);
  double dz = cos (longitude) * cos (latitude);

  if (molecule_comment (molecule) && *molecule_comment (molecule)) {
    glColor3d (1.0, 0.0, 0.0);
    glRasterPos2d (-2.2, 2.2);
    glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10,
		     (const unsigned char *)molecule_comment (molecule));
  }

  if (molecule_file (molecule) && *molecule_file (molecule)) {
    glColor3d (1.0, 0.0, 0.0);
    glRasterPos2d (-2.2, 2.0);
    if (molecule_shortname (molecule))
      glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10,
		       (const unsigned char *)molecule_shortname (molecule));
  }

  glScaled (zoom, zoom, zoom);

  /**
     compute rotational axis and rotation angle from [0 0 1] to
     latitude and longitude.
  **/
  rotation_s dist = {0.0, dx, dy, dz};
  rotation_s cyl  = {0.0, 0.0, 0.0, 1.0};
  rotation_s cp;

  double theta = rad_to_deg (dot_product (&dist, &cyl));
  cross_product (&cp, &dist, &cyl);
  cp.theta = theta;
  glRotated (-theta, cp.x, cp.y, cp.z);
  glTranslated (-hoff, voff, 0.0);
  
  draw_molecule (molecule, &cp);

  glPopMatrix();

  return TRUE;
}

static void inline
average_atoms (gpointer data, gpointer user_data)
{
  atom_s *atom = data;
  rotation_s *rotation = user_data;

  atom_x (atom) -= rotation->x;
  atom_y (atom) -= rotation->y;
  atom_z (atom) -= rotation->z;
}

static molecule_s *
read_file (gchar *fn)
{
  molecule_s *molecule = NULL;
  FILE *xyz = fopen (fn, "r");
  if (xyz) {
    char  *buffer = NULL;
    size_t buffer_len = 0;
    size_t char_ct;

    double min_x = MAXDOUBLE;
    double max_x = MINDOUBLE;
    double min_y = MAXDOUBLE;
    double max_y = MINDOUBLE;
    double min_z = MAXDOUBLE;
    double max_z = MINDOUBLE;

    if (buffer) *buffer = 0;
    char_ct = getline (&buffer, &buffer_len, xyz);	// read atom count
    if (char_ct > 0) {
      int count = atoi (buffer);
      if (count > 0) {
	char  *comment = NULL;
	size_t comment_len = 0;
	char_ct = getline (&comment, &comment_len, xyz);	// read comment
	if (comment) comment[char_ct - 1] = 0;			// kill the nl
	GSList *atoms = NULL;
	GSList *bonds = NULL;
	int i;
	for (i = 0; i < count; i++) {
#define SYM_LEN	8
	  char sym[SYM_LEN];		// place to hold the element symbol
	  atom_s *atom = malloc (sizeof(atom_s));	// create an atom record
	  atom_seq (atom) = i;
	  sym[0] = 0;				// initialise the symbol
	  int ct = fscanf (xyz, " %8s %lf %lf %lf%*[^\n]\n",
			   sym,		// read and parse the atom symbol
			   &atom_x (atom),	// and xyz data
			   &atom_y (atom),
			   &atom_z (atom));
	  if (min_x > atom_x (atom)) min_x = atom_x (atom); // extend 
	  if (max_x < atom_x (atom)) max_x = atom_x (atom); // bounding
	  if (min_y > atom_y (atom)) min_y = atom_y (atom); // box
	  if (max_y < atom_y (atom)) max_y = atom_y (atom);
	  if (min_z > atom_z (atom)) min_z = atom_z (atom);
	  if (max_z < atom_z (atom)) max_z = atom_z (atom);
	  if (ct >= 4) {
	     int k;
	     for (k = 0; k < strlen (sym) && g_ascii_isalpha (sym[k]); k++) {}
	     sym[k] = 0;
	     atom_name (atom) = strdup (sym);	// copy in the raw symbol
	     atom_data (atom) = get_rad_data (sym);
	     atoms = g_slist_append (atoms, atom);
	  }
	}
	if (i == count) {		// got the right number of atoms
	  if (count >= 2) {
	    rotation_s rotation;
	    rotation.x = (max_x + min_x) / 2.0;
	    rotation.y = (max_y + min_y) / 2.0;
	    rotation.z = (max_z + min_z) / 2.0;
	    g_slist_foreach (atoms, average_atoms, &rotation);
	    
	    for (int a = 0; a < count; a++) {
	      atom_s *atoma = g_slist_nth_data (atoms, a);
	      for (int b = a+1; b < count; b++) {
		atom_s *atomb = g_slist_nth_data (atoms, b);
		if (atoma && atomb) {
		  bond_e bond_type = get_bond_type (atoma, atomb);
		  
		  if (bond_type != BOND_NONE) {
		    double strength =  1.0;
		    radii_s *rada = atom_data (atoma);
		    radii_s *radb = atom_data (atomb);
		    bond_s *bond = malloc (sizeof(bond_s));
		    if (rada && radb)
		      strength = get_strength (tp_nr (rada), tp_nr (radb));
		    bond_a (bond) = atoma;
		    bond_b (bond) = atomb;
		    bond_type (bond) = bond_type;
		    bond_strength (bond) = strength;
		    double dx = atom_x (atomb) - atom_x (atoma);
		    double dy = atom_y (atomb) - atom_y (atoma);
		    double dz = atom_z (atomb) - atom_z (atoma);
		    bond_separation (bond)
		      = sqrt ((dx * dx) + (dy * dy) + (dz * dz));
		    bonds = g_slist_append (bonds, bond);
		  }	
		}
	      }
	    }
	  }
	  char *slash = strrchr (fn, '/');
	  slash = slash ? slash+1 : fn;
	  
	  molecule = malloc (sizeof(molecule_s));
	  molecule_atoms (molecule)	= atoms;
	  molecule_bonds (molecule)	= bonds;
	  molecule_comment (molecule)	= comment;
	  molecule_file (molecule)	= strdup (fn);
	  molecule_shortname (molecule)	= strdup (slash);
	  char *dot = strrchr (molecule_shortname (molecule), '.');
	  if (dot) *dot = 0;
	  molecule_area (molecule)	= NULL;
	  molecule_window (molecule)	= NULL;
	  molecule_zoom (molecule)	= INITIAL_ZOOM;
	  molecule_latitude (molecule)	= INITIAL_LATITUDE;
	  molecule_longitude (molecule)	= INITIAL_LONGITUDE;
	  molecule_voff (molecule)	= 0.0;
	  molecule_hoff (molecule)	= 0.0;
	  molecules = g_slist_append (molecules, molecule);
	  GtkTreeIter   iter;
	  gtk_list_store_append (store, &iter);
	  gtk_list_store_set (store, &iter,
			      COLUMN_NAME, molecule_shortname (molecule),
			      COLUMN_COMMENT,
			      molecule_comment (molecule) ?: "",
			      COLUMN_MOLECULE, molecule,
			      -1);
	}
	else {
	  // fixme -- error
	}
      }
    }
    fclose (xyz);
  }
  else {
    // fixme -- handle errors
  }
  return molecule;
}

static void
print_atom (gpointer data, gpointer user_data)
{
  atom_s *atom = data;
  FILE *file = user_data;
  fprintf (file, "    Atom: [%2d] %s\n", atom_seq (atom), atom_name (atom));
}

static void
print_bond (gpointer data, gpointer user_data)
{
  bond_s *bond = data;
  FILE *file = user_data;
  char   *btype = "";
  atom_s *atoma = bond_a (bond);
  atom_s *atomb = bond_b (bond);
  switch (bond_type (bond)) {
  case BOND_NONE:   btype = " x "; break;
  case BOND_SINGLE: btype = "<->"; break;
  case BOND_DOUBLE: btype = "<=>"; break;
  case BOND_TRIPLE: btype = "<#>"; break;
  }
  fprintf (file, "    Bond: [%2d] %s %s [%2d] %s\n",
	  atom_seq (atoma), atom_name (atoma),
	  btype,
	  atom_seq (atomb), atom_name (atomb));
}

static void
print_molecule (molecule_s *molecule, gchar *fn)
{
  FILE *file = fopen (fn, "w");
  if (file) {
    fprintf (file, "  Comment: %s\n", molecule_comment (molecule) ?: "<>");
    if (molecule_atoms (molecule))
      g_slist_foreach (molecule_atoms (molecule), print_atom, file);
    if (molecule_bonds (molecule))
      g_slist_foreach (molecule_bonds (molecule), print_bond, file);
    fclose (file);
  }
  else {
    // fixme -- error
  }
}

static void
print_view (GtkWidget *widget, gpointer data)
{
  molecule_s *molecule = data;
   GtkWidget *dialog =
     gtk_file_chooser_dialog_new (_("Print molecule"),
				 GTK_WINDOW (molecule_window (molecule)),
				 GTK_FILE_CHOOSER_ACTION_SAVE,
				  _ ("Accept"), GTK_RESPONSE_ACCEPT,
				  _ ("Cancel"), GTK_RESPONSE_CANCEL,
				 NULL);
  gtk_window_set_position (GTK_WINDOW (dialog), GTK_WIN_POS_MOUSE);
  gtk_dialog_set_default_response (GTK_DIALOG (dialog),
                                   GTK_RESPONSE_ACCEPT);
  gtk_widget_show_all (dialog);
  gint response = gtk_dialog_run (GTK_DIALOG (dialog));
  if (response == GTK_RESPONSE_ACCEPT) {
    gchar *file = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (dialog));
    gboolean doit = TRUE;
    if (g_file_test (file, G_FILE_TEST_EXISTS)) {
      GtkWidget *dialogq =
	gtk_message_dialog_new (GTK_WINDOW (molecule_window (molecule)),
				GTK_DIALOG_MODAL,
				GTK_MESSAGE_QUESTION,
				GTK_BUTTONS_YES_NO,
				_ ("File %s exists.  Overwrite it?"),
				file);
      gtk_widget_show_all (dialogq);
      response = gtk_dialog_run (GTK_DIALOG (dialogq));
      if (response != GTK_RESPONSE_YES) doit = FALSE;
      gtk_widget_destroy (dialogq);
    }
    if (doit) print_molecule (molecule, file);
    g_free (file);
  }
  gtk_widget_destroy (dialog);
}

static gboolean
zoom_changed (GtkRange *range,
	      GtkScrollType scroll,
	      gdouble       value,
	      gpointer user_data)
{
  molecule_s *molecule = user_data;
  if (molecule) {
    molecule_zoom (molecule) = value;
    gtk_gl_area_queue_render (molecule_area (molecule));
  }
  return FALSE;		// do other stuff if needed
}

static gboolean
latitude_changed (GtkRange *range,
		  GtkScrollType scroll,
		  gdouble       value,
		  gpointer user_data)
{
  molecule_s *molecule = user_data;
  if (molecule) {
    if (value < -180.0) value = -180.0;
    if (value >  180.0) value =  180.0;
    molecule_latitude (molecule) = -value;
    gtk_gl_area_queue_render (molecule_area (molecule));
  }
  return FALSE;		// do other stuff if needed
}

static gboolean
voff_changed (GtkRange *range,
	      GtkScrollType scroll,
	      gdouble       value,
	      gpointer user_data)
{
  molecule_s *molecule = user_data;
  if (molecule) {
    if (value < -5.0) value = -5.0;
    if (value >  5.0) value =  5.0;
    molecule_voff (molecule) = value;
    gtk_gl_area_queue_render (molecule_area (molecule));
  }
  return FALSE;		// do other stuff if needed
}

static gboolean
hoff_changed (GtkRange *range,
	      GtkScrollType scroll,
	      gdouble       value,
	      gpointer user_data)
{
  molecule_s *molecule = user_data;
  if (molecule) {
    if (value < -5.0) value = -5.0;
    if (value >  5.0) value =  5.0;
    molecule_hoff (molecule) = value;
    gtk_gl_area_queue_render (molecule_area (molecule));
  }
  return FALSE;		// do other stuff if needed
}

static gboolean
longitude_changed (GtkRange *range,
		   GtkScrollType scroll,
		   gdouble       value,
		   gpointer user_data)
{
  molecule_s *molecule = user_data;
  if (molecule) {
    if (value < -180.0) value = -180.0;
    if (value >  180.0) value =  180.0;
    molecule_longitude (molecule) = value;
    gtk_gl_area_queue_render (molecule_area (molecule));
  }
  return FALSE;		// do other stuff if needed
}

#if 0
static gboolean
button_press (GtkWidget *widget,
	      GdkEvent  *event,
	      gpointer   user_data)
{
  g_print ("button press\n");
}
#endif

static void
free_bond (gpointer data)
{
  bond_s *bond = data;
  if (bond) free (bond);
}

static void
free_atom (gpointer data)
{
  atom_s *atom = data;
  if (atom) {
    if (atom_name (atom)) free (atom_name (atom));
    free (atom);
  }
}
  
static void
free_molecule (gpointer data)
{
  molecule_s *molecule = data;
  if (molecule) {
    if (molecule_file (molecule))    free (molecule_file (molecule));
    if (molecule_shortname (molecule)) free (molecule_shortname (molecule));
    if (molecule_comment (molecule)) free (molecule_comment (molecule));
    if (molecule_bonds (molecule))
      g_slist_free_full (molecule_bonds (molecule), free_bond);
    if (molecule_atoms (molecule))
      g_slist_free_full (molecule_atoms (molecule), free_atom);
    free (molecule);
  }
}

static void
gmol_quit (GtkWidget *widget, gpointer data)
{
  g_slist_free_full (molecules, free_molecule);

  for (int i = 0; i < 3; i++)
    if (font[i]) ftglDestroyFont (font[i]);
  
  gtk_main_quit ();
}

static void
dump_view (molecule_s *molecule, const char *fname)
{
  GtkGLArea *area = molecule_area (molecule);
  gtk_gl_area_make_current (area);
  int ww = gtk_widget_get_allocated_width (GTK_WIDGET (area));
  int wh = gtk_widget_get_allocated_height (GTK_WIDGET (area));

  unsigned char *buffer = g_malloc0 (4 * ww * wh);
  glReadBuffer (GL_BACK_LEFT);
  glReadPixels (0,                     // GLint x,
	        0,                     // GLint y,
	        ww,                    // GLsizei width,
	        wh,                    // GLsizei height,
	        GL_BGRA,               // enum format,
	        GL_UNSIGNED_BYTE,      // GLenum type,
	        buffer);


  //  MagickBooleanType status;
  MagickWand *magick_wand;
  DrawingWand *d_wand;
  PixelWand *p_wand;


  magick_wand = NewMagickWand ();
  d_wand = NewDrawingWand ();
  p_wand = NewPixelWand ();

  MagickWandGenesis ();

  /*status = */MagickConstituteImage (magick_wand,
				  ww,  // const size_t columns,
				  wh,  // const size_t rows,
				  "BGRA",      // const char *map,
				  CharPixel,   // const StorageType storage,
				  buffer);     // void *pixels)
  /***
      Don't know why, but the image gets dumped upside down, so this
      flips it.
  ***/
  MagickFlipImage (magick_wand);

  MagickDrawImage (magick_wand, d_wand);
  /*status = */MagickWriteImages (magick_wand, fname, MagickTrue);

  if(magick_wand) magick_wand = DestroyMagickWand (magick_wand);
  if(d_wand) d_wand = DestroyDrawingWand (d_wand);
  if(p_wand) p_wand = DestroyPixelWand (p_wand);

  MagickWandTerminus(); 
  
  g_free (buffer);
}

static void
open_files (gpointer data, gpointer user_data)
{
  gchar *file = data;
  molecule_s *molecule = read_file (file);
  if (molecule) open_molecule_window (molecule, NULL);
}

static void
import_view (GtkWidget *widget, gpointer data)
{
  GtkWidget *window = data;

  GtkFileFilter *filter = gtk_file_filter_new ();
  gtk_file_filter_add_pattern (filter, "*.xyz");
  
  GtkWidget *dialog =
    gtk_file_chooser_dialog_new (_ ("Import XYZ file"),
				 GTK_WINDOW (window),
				 GTK_FILE_CHOOSER_ACTION_OPEN,
				 _ ("Accept"), GTK_RESPONSE_ACCEPT,
				 _ ("Cancel"), GTK_RESPONSE_CANCEL,
				 NULL);
  gtk_file_chooser_set_select_multiple (GTK_FILE_CHOOSER (dialog), TRUE);
  gtk_window_set_position (GTK_WINDOW (dialog), GTK_WIN_POS_MOUSE);
  gtk_dialog_set_default_response (GTK_DIALOG (dialog),
                                   GTK_RESPONSE_ACCEPT);
  gtk_file_chooser_add_filter (GTK_FILE_CHOOSER (dialog), filter);
  gtk_widget_show_all (dialog);
  gint response = gtk_dialog_run (GTK_DIALOG (dialog));
  if (response == GTK_RESPONSE_ACCEPT) {
#if 1
    GSList *filenames =
      gtk_file_chooser_get_filenames (GTK_FILE_CHOOSER (dialog));
    g_slist_foreach (filenames, open_files, NULL);
    g_slist_free_full (filenames, g_free);
#else
    gchar *file = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (dialog));
    molecule_s *molecule = read_file (file);
    if (molecule) open_molecule_window (molecule, NULL);
    g_free (file);
#endif
  }
  gtk_widget_destroy (dialog);
}

static void
export_view (GtkWidget *widget, gpointer data)
{
  molecule_s *molecule = data;

  GtkWidget *dialog =
    gtk_file_chooser_dialog_new (_ ("Export image"),
				 GTK_WINDOW (molecule_window (molecule)),
				 GTK_FILE_CHOOSER_ACTION_SAVE,
				 _ ("Accept"), GTK_RESPONSE_ACCEPT,
				 _ ("Cancel"), GTK_RESPONSE_CANCEL,
				 NULL);
  gtk_window_set_position (GTK_WINDOW (dialog), GTK_WIN_POS_MOUSE);
  gtk_dialog_set_default_response (GTK_DIALOG (dialog),
                                   GTK_RESPONSE_ACCEPT);
  gtk_widget_show_all (dialog);
  gint response = gtk_dialog_run (GTK_DIALOG (dialog));
  if (response == GTK_RESPONSE_ACCEPT) {
    gchar *file = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (dialog));
    gboolean doit = TRUE;
    if (g_file_test (file, G_FILE_TEST_EXISTS)) {
      GtkWidget *dialogq =
	gtk_message_dialog_new (GTK_WINDOW (molecule_window (molecule)),
				GTK_DIALOG_MODAL,
				GTK_MESSAGE_QUESTION,
				GTK_BUTTONS_YES_NO,
				_ ("File %s exists.  Overwrite it?"),
				file);
      gtk_widget_show_all (dialogq);
      response = gtk_dialog_run (GTK_DIALOG (dialogq));
      if (response != GTK_RESPONSE_YES) doit = FALSE;
      gtk_widget_destroy (dialogq);
    }
    if (doit) dump_view (molecule, file);
    g_free (file);
  }
  gtk_widget_destroy (dialog);
}

static void
export_sgi (gpointer widget, gpointer data)
{
  molecule_s *molecule = data;

  /***
      being lazy and sneaky here--if export_sgi is called from a menu,
      widget holds a widget we don't need and the molecule is in the data.
      if called from foreach, data is null and widget is repurposed to hold
      the molecule.
  ***/
  
  if (!molecule) molecule = (molecule_s *)widget;

  if (molecule_file (molecule)) {
    time_t t = time (NULL);
#define TBUF_LEN        128
    char *tbuf = alloca (TBUF_LEN);
    strftime (tbuf, TBUF_LEN, "%Y-%m-%d-%H-%m-%S", localtime (&t));
    char *fname = g_strdup_printf ("%s-%s.sgi",
				   molecule_shortname (molecule),
				   tbuf);
    dump_view (molecule, fname);
    g_free (fname);
  }
}

 
static void
export_all (GtkWidget *widget, gpointer data)
{
  if (molecules) {
    g_slist_foreach (molecules, export_sgi, NULL);
  }
}

static void
close_view (GtkWidget *widget, gpointer data)
{
  molecule_s *molecule = data;
  GtkWidget *window = molecule_window (molecule);
  molecule_window (molecule) = NULL;
  molecule_area (molecule) = NULL;
  if (window) gtk_widget_destroy (window);
}

static void
open_molecule_window (gpointer data, gpointer user_data)
{
  GtkWidget *thingy;
  GtkAdjustment *adj;
  molecule_s *molecule = data;
  
  window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
  molecule_window (molecule) = window;
  g_signal_connect (window, "destroy", G_CALLBACK (close_view), molecule);
  gtk_window_set_title (GTK_WINDOW (window), molecule_file (molecule));

  GtkWidget *vbox = gtk_box_new (GTK_ORIENTATION_VERTICAL, 8);
  gtk_container_add (GTK_CONTAINER (window), vbox);

  {
    GtkWidget *menubar;
    GtkWidget *menu;
    GtkWidget *item;

    menubar = gtk_menu_bar_new();

    /********* file menu ********/

    menu = gtk_menu_new();
    item = gtk_menu_item_new_with_label (_ ("File"));
    gtk_menu_item_set_submenu (GTK_MENU_ITEM (item), menu);
    gtk_menu_shell_append (GTK_MENU_SHELL (menubar), item);

    item = gtk_menu_item_new_with_label (_ ("Import XYZ..."));
    g_signal_connect (G_OBJECT (item), "activate",
		      G_CALLBACK (import_view), molecule);
    gtk_menu_shell_append (GTK_MENU_SHELL (menu), item);

    item = gtk_separator_menu_item_new();
    gtk_menu_shell_append (GTK_MENU_SHELL (menu), item);

    item = gtk_menu_item_new_with_label (_ ("Export view..."));
    g_signal_connect (G_OBJECT (item), "activate",
		      G_CALLBACK (export_view), molecule);
    gtk_menu_shell_append (GTK_MENU_SHELL (menu), item);

    item = gtk_menu_item_new_with_label (_ ("Snapshot SGI"));
    g_signal_connect (G_OBJECT (item), "activate",
		      G_CALLBACK (export_sgi), molecule);
    gtk_menu_shell_append (GTK_MENU_SHELL (menu), item);

    item = gtk_menu_item_new_with_label (_ ("Print view..."));
    g_signal_connect (G_OBJECT (item), "activate",
		      G_CALLBACK (print_view), molecule);
    gtk_menu_shell_append (GTK_MENU_SHELL (menu), item);

    item = gtk_separator_menu_item_new();
    gtk_menu_shell_append (GTK_MENU_SHELL (menu), item);

    item = gtk_menu_item_new_with_label (_ ("Close view"));
    g_signal_connect (G_OBJECT (item), "activate",
		      G_CALLBACK (close_view), molecule);
    gtk_menu_shell_append (GTK_MENU_SHELL (menu), item);

    gtk_box_pack_start (GTK_BOX (vbox), GTK_WIDGET (menubar),
			FALSE, FALSE, 2);

  }

  GtkWidget *grid = gtk_grid_new ();
  gtk_box_pack_start (GTK_BOX (vbox), grid, TRUE, TRUE, 2);

  /***** zoom *****/
  adj = gtk_adjustment_new (INITIAL_ZOOM,	// gdouble value,
			    0.1,	// gdouble lower,
			    10.0,	// gdouble upper,
			    0.1,	// gdouble step_increment,
			    0.1,	// gdouble page_increment,
			    0.1);	// gdouble page_size);
  thingy = gtk_scrollbar_new (GTK_ORIENTATION_HORIZONTAL, adj);
  gtk_widget_set_tooltip_text (thingy, "Zoom");
  g_signal_connect (thingy, "change-value",
		    G_CALLBACK (zoom_changed), molecule);
  gtk_grid_attach (GTK_GRID(grid), thingy, 0, 0, 1, 1);


  /***** open gl area *****/
  GtkWidget *gl_area = gtk_gl_area_new ();

  gtk_gl_area_set_has_depth_buffer (GTK_GL_AREA (gl_area), TRUE);
  molecule_area (molecule) = GTK_GL_AREA (gl_area);
  gtk_widget_set_size_request (gl_area, 400, 400);
  gtk_widget_set_hexpand (gl_area, TRUE);
  gtk_widget_set_vexpand (gl_area, TRUE);
#if 0
  gtk_widget_add_events (gl_area, GDK_BUTTON_PRESS_MASK);
  g_signal_connect (gl_area, "button-press-event",
		    G_CALLBACK (button_press), molecule);
#endif
  g_signal_connect (gl_area, "render",
		    G_CALLBACK (gl_render), molecule);
  g_signal_connect (gl_area, "realize",
		    G_CALLBACK (gl_realise), molecule);
  g_signal_connect (gl_area, "resize",
		    G_CALLBACK (gl_resize), molecule);
  gtk_grid_attach (GTK_GRID(grid), gl_area, 0, 1, 1, 1);

  /***** latitude *****/
  adj = gtk_adjustment_new (molecule_latitude (molecule),
			    -180.0,	// gdouble lower,
			     180.0,	// gdouble upper,
			    1.0,	// gdouble step_increment,
			    10.0,	// gdouble page_increment,
			    10.0);	// gdouble page_size);
  thingy = gtk_scrollbar_new (GTK_ORIENTATION_VERTICAL, adj);
  gtk_widget_set_tooltip_text (thingy, "Latitude");
  g_signal_connect (thingy, "change-value",
		    G_CALLBACK (latitude_changed), molecule);
  gtk_grid_attach (GTK_GRID(grid), thingy, 1, 1, 1, 1);

  /***** voff *****/
  adj = gtk_adjustment_new (molecule_voff (molecule),
			    -5.0,	// gdouble lower,
			     5.0,	// gdouble upper,
			    0.1,	// gdouble step_increment,
			    1.0,	// gdouble page_increment,
			    1.0);	// gdouble page_size);
  thingy = gtk_scrollbar_new (GTK_ORIENTATION_VERTICAL, adj);
  gtk_widget_set_tooltip_text (thingy, "Vertical offset");
  g_signal_connect (thingy, "change-value",
		    G_CALLBACK (voff_changed), molecule);
  gtk_grid_attach (GTK_GRID(grid), thingy, 2, 1, 1, 1);

  /***** longitude *****/
  adj = gtk_adjustment_new (molecule_longitude (molecule),
			    -180.0,	// gdouble lower,
			     180.0,	// gdouble upper,
			    1.0,	// gdouble step_increment,
			    10.0,	// gdouble page_increment,
			    10.0);	// gdouble page_size);
  thingy = gtk_scrollbar_new (GTK_ORIENTATION_HORIZONTAL, adj);
  gtk_widget_set_tooltip_text (thingy, "Longitude");
  g_signal_connect (thingy, "change-value",
		    G_CALLBACK (longitude_changed), molecule);
  gtk_grid_attach (GTK_GRID(grid), thingy, 0, 2, 1, 1);

  /***** hoff *****/
  adj = gtk_adjustment_new (molecule_hoff (molecule),
			    -5.0,	// gdouble lower,
			     5.0,	// gdouble upper,
			    0.1,	// gdouble step_increment,
			    1.0,	// gdouble page_increment,
			    1.0);	// gdouble page_size);
  thingy = gtk_scrollbar_new (GTK_ORIENTATION_HORIZONTAL, adj);
  gtk_widget_set_tooltip_text (thingy, "Horizontalal offset");
  g_signal_connect (thingy, "change-value",
		    G_CALLBACK (hoff_changed), molecule);
  gtk_grid_attach (GTK_GRID(grid), thingy, 0, 3, 1, 1);

  gtk_widget_show_all (window);
}

static void
row_activated (GtkTreeView *tree_view,
	       GtkTreePath *path,
	       GtkTreeViewColumn *column,
	       gpointer user_data)
{
  GtkTreeIter iter;

  if (gtk_tree_model_get_iter (GTK_TREE_MODEL (store), &iter, path)) {
    molecule_s *molecule;
    gtk_tree_model_get (GTK_TREE_MODEL (store), &iter,
			COLUMN_MOLECULE, &molecule,
			-1);
    if (molecule_window (molecule)) close_view (NULL, molecule);
    else open_molecule_window (molecule, NULL);
  }
}

int
main (int ac, char *av[])
{
  GError *error = NULL;
  GOptionEntry entries[] = {
    { NULL }
  };
  GOptionContext *context =
    g_option_context_new (_ ("string string string..."));
  g_option_context_add_main_entries (context, entries, NULL);
  g_option_context_add_group (context, gtk_get_option_group (TRUE));
  
  gtk_init (&ac, &av);
  
  store = gtk_list_store_new (NR_COLUMNS,
			      G_TYPE_STRING,
			      G_TYPE_STRING,
			      G_TYPE_POINTER);

  if (!g_option_context_parse (context, &ac, &av, &error)) {
    g_warning (_ ("option parsing failed: %s\n"), error->message);
    g_clear_error (&error);
  }
  
  glutInit (&ac, av);

  font[0] = ftglCreateExtrudeFont (FONT_FILE);
  ftglSetFontFaceSize (font[0], 80.0, 72.0);
  ftglSetFontDepth (font[0], 10.0);
  //  ftglSetFontOutset (font[0], 0, 3);
  ftglSetFontCharMap (font[0], ft_encoding_unicode);

#if 0
  font[1] = ftglCreateBufferFont (FONT_FILE);
  ftglSetFontFaceSize (font[1], 80.0, 72.0);
  ftglSetFontCharMap (font[1], ft_encoding_unicode);

  font[2] = ftglCreateOutlineFont (FONT_FILE);
  ftglSetFontFaceSize (font[2], 80.0, 72.0);
  ftglSetFontCharMap (font[2], ft_encoding_unicode);
#endif

  fontindex = 0;

  for (int i = 1; i < ac; i++) read_file (av[i]);
#if 0
  if (molecules) {
    g_print ("Molecules:\n\n");
    g_slist_foreach (molecules, dump_molecule, NULL);
  }
#endif

  
  if (molecules) {
    g_slist_foreach (molecules, open_molecule_window, NULL);
  }

  window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
  gtk_window_set_default_size (GTK_WINDOW (window), 256, 256);
  gtk_window_set_title (GTK_WINDOW (window), "gmol");

  g_signal_connect (window, "destroy", G_CALLBACK (gmol_quit), NULL);

  GtkWidget *vbox = gtk_box_new (GTK_ORIENTATION_VERTICAL, 8);
  gtk_container_add (GTK_CONTAINER (window), vbox);

   {
    GtkWidget *menubar;
    GtkWidget *menu;
    GtkWidget *item;

    menubar = gtk_menu_bar_new();

    /********* file menu ********/

    menu = gtk_menu_new();
    item = gtk_menu_item_new_with_label (_ ("File"));
    gtk_menu_item_set_submenu (GTK_MENU_ITEM (item), menu);
    gtk_menu_shell_append (GTK_MENU_SHELL (menubar), item);

    item = gtk_menu_item_new_with_label (_ ("Import XYZ..."));
    g_signal_connect (G_OBJECT (item), "activate",
		      G_CALLBACK (import_view), window);
    gtk_menu_shell_append (GTK_MENU_SHELL (menu), item);

    item = gtk_separator_menu_item_new();
    gtk_menu_shell_append (GTK_MENU_SHELL (menu), item);

    item = gtk_menu_item_new_with_label (_ ("Export all"));
    g_signal_connect (G_OBJECT (item), "activate",
		      G_CALLBACK (export_all), NULL);
    gtk_menu_shell_append (GTK_MENU_SHELL (menu), item);

    item = gtk_separator_menu_item_new();
    gtk_menu_shell_append (GTK_MENU_SHELL (menu), item);

    item = gtk_menu_item_new_with_label (_ ("Quit"));
    g_signal_connect (G_OBJECT (item), "activate",
		      G_CALLBACK (gmol_quit), window);
    gtk_menu_shell_append (GTK_MENU_SHELL (menu), item);

    gtk_box_pack_start (GTK_BOX (vbox), GTK_WIDGET (menubar),
			FALSE, FALSE, 2);
   }

   {
     GtkWidget *sw = gtk_scrolled_window_new (NULL, NULL);
     gtk_box_pack_start (GTK_BOX (vbox), GTK_WIDGET (sw), TRUE, TRUE, 2);

     GtkWidget *tree = gtk_tree_view_new_with_model (GTK_TREE_MODEL (store));
     //     gtk_tree_view_set_headers_visible (GTK_TREE_VIEW (tree), FALSE);
     g_signal_connect (tree, "row-activated",
		       G_CALLBACK (row_activated), NULL);
     
     gtk_tree_view_set_activate_on_single_click  (GTK_TREE_VIEW (tree), TRUE);

     GtkCellRenderer *renderer;
     GtkTreeViewColumn *column;
     
     renderer = gtk_cell_renderer_text_new ();
     column = gtk_tree_view_column_new_with_attributes (_ ("Name"),
							renderer,
							"text",
							COLUMN_NAME,
							NULL);
     gtk_tree_view_append_column (GTK_TREE_VIEW (tree), column);
     
     renderer = gtk_cell_renderer_text_new ();
     column = gtk_tree_view_column_new_with_attributes (_ ("Description"),
							renderer,
							"text",
							COLUMN_COMMENT,
							NULL);
     gtk_tree_view_append_column (GTK_TREE_VIEW (tree), column);

     
     gtk_container_add (GTK_CONTAINER (sw), tree);

   }

  gtk_widget_show_all (window);
  gtk_main ();
 
}
