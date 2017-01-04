# gmol

An GTK+/OpenGL-based XYZ-format molecule visualiser/format-converter.

gmol is a simple molecule viewer for molecules described in XYZ-format
files.  Using OpenGL, it provides 3-space visualisation of the molecules,
allowing the user to rotate them on all axes and to zoom in and out of them.

XYZ-format molecular representation contains no information concerning
atomic bonds, bonds being inferred by the proximity of the atoms--atoms
physically within the sum of the covalent radii of the respective atoms
are considered bonded.  To identify bonds, this application uses the table
of covalent radii devised by Pekka Pyyko; this table provides radii for
single, double, and triple bonds, allowing the application to visually
distinguish those bonds.  Bond strength is also visually displayed both
by the apparent diameter of the bond representation and by its colour, where
the strength is mapped to the colour spectrum, red corresponding to weaker
bonds through violet corresponding to stronger bonds.  At present, double
bonds are represented by bright yellow, triple bonds by bright cyan, though
this may change in the future.

The application distinguishes atoms using the jmol colour scheme.

Comments, questions, and requests are welcome.

Chris Moller
moller@mollerware.com
