c dyn.inc
c$Author: baskes $
c$Date: 2003/12/06 04:19:37 $
c$Source: /home/baskes/programs/MEAM/v88/RCS/dyn.inc,v $
c$Revision: 1.12 $
c$Log: dyn.inc,v $
cRevision 1.12  2003/12/06 04:19:37  baskes
cupdating
c.'
c
cRevision 1.11  2003/07/01 03:39:34  baskes
cgetting ready to change dyn.inc
c
cRevision 1.10  2002/04/26 20:56:00  baskes
cadd vcm averaging
c
cRevision 1.9  2002/01/04 23:15:33  baskes
cmodified averaging
c
cRevision 1.8  2002/01/03 23:25:21  baskes
cadded sheardot
c
cRevision 1.7  2001/06/12 15:34:25  baskes
cadded acf
c
cRevision 1.6  2001/05/15 21:34:18  baskes
cadded table
c
cRevision 1.5  2001/01/17 01:03:33  baskes
caverage slocal added
c
cRevision 1.4  2001/01/08 23:31:47  baskes
cfixed header
c
cRevision 1.3  2001/01/08 23:02:36  baskes
cfixed header
c
	include 'utils.inc'
      logical nn
      common /twonn/ a2nn(nelmax),b2nn(nelmax),nn
      common /grids/ phitab(ntable,nelmax,nelmax),itable
      common /mea/ meamf,meafile,forcef,printf,rstrtf,conff,timestr
      character*80 printf,rstrtf,conff,meamf,forcef,meafile
      character*24 timestr
