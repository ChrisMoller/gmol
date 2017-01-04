 bin_PROGRAMS = gmol
gmol_SOURCES = cvbond.c cvbond.h \
               gmol.c gmol.h \
               jmol-colours.c jmol-colours.h \
               pyykko-bonds.c pyykko-bonds.h \
               radii.c radii.h
DISTCLEANFILES = .deps/*.P


gmol_CFLAGS = -std=c99 -Wall -Werror \
              $(GTK_CFLAGS) $(GL_CFLAGS) $(GLU_CFLAGS) \
              $(FTGL_CFLAGS) $(MAGICKWAND_CFLAGS)

gmol_LDFLAGS = -lm -lglut $(GTK_LIBS) $(GL_LIBS) $(GLU_LIBS) \
                $(FTGL_LIBS) $(MAGICKWAND_LIBS)

