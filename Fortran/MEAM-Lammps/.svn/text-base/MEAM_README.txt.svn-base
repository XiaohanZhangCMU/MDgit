
Instructions for building and using the MEAM potential in LAMMPS
----------------------------------------------------------------


1) Download and install the latest LAMMPS installation from
http://lammps.sandia.gov.  Make sure you can build LAMMPS on its own
-- if you can't get this to work, adding the MEAM library is only
going to make things harder!  You may need to create a new Makefile
under the MAKE director that works on your system.


2) Compile the MEAM library:

   a) In the lib directory under LAMMPS, type "mkdir meam; cd meam"

   b) Copy the MEAM.tar file into this directory

   c) Untar it: "tar -xvf MEAM.tar"

   d) Type "make".

The Makefile currently uses the gfortran compiler.  You may need to
edit it to use a different compiler on your platform.


3) Compile lammps with the MEAM pair potential:

   a) In the src director under LAMMPS, type "mkdir MEAM; cd MEAM"

   b) Copy the MEAM-lammps.tar file into this directory

   c) Untar it: "tar -xvf MEAM-lammps.tar"

   d) Copy the Makefile into the src director: "cp Makefile .."
      (Alternatively, if you have a new-ish version of LAMMPS, you may
      be able to simply uncomment some MEAM-related lines in the
      existing src/Makefile.

   e) Copy the style_manybody.h file into src: "cp style_manybody.h .."

   f) Edit the machine-specific Makefile under src/MAKE to include the
      MEAM library.  Add "-L../../lib/meam" to the LINKFLAGS, and
      "-lmeam" to the USRLIB.

   f) Back in the src directory, type "make yes-meam"

   g) Try to compile LAMMPS.  If everything works, great -- but you'll
      probably get link errors complaining of undefined references.
      On to step 4...


4) Resolve link errors

There are a few things that might be going wrong if you get link
errors.  If you see undefined references to meam-related subroutines
(like "meam_setup_param_"), it may be that you're not properly linking
to the MEAM library file meam.a -- check your link flags.  Also
possible is that there's a problem with leading or trailing
underscores in symbol names.  This can usually be handled by changing
the compile flags in lib/meam/Makefile.  It may be helpful to use the
"nm" tool on the libmeam.a file to see what format the symbol names in
that library take.  In the worst case, you may need to edit subroutine
calls in pair_meam.cpp to get the underscores right.

More difficult are errors brought on by trying to link fortran and C
objects.  These are undefined references to symbols that live in
libraries that get linked in automatically if you use a fortran
compiler to link, but are missed because LAMMPS uses a C compiler to
link.  Which library you need depends on which fortran compiler you're
using.  On my system (a Redhat Linux box), when I use the gfortran
compiler I need to link to libgfortran.a (and include the path to it
in the link flags).  Your results may vary.


5) Test

CD into the examples directory, copy the file MEAM-examples.tar file
here, and untar it.  The resulting meam directory includes two example
files, but of which run NVE systems.  One is a nickel fcc lattice, the
other a Si/C system read in from an input file.  The main command
files are of the form "in.*".  The other files provide other input
data.  

