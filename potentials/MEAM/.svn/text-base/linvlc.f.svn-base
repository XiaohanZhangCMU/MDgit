        subroutine linvlc
c     *****************
c
c sets up the link cell map
c
      include 'moldy.h'

        hmeps = 0.5 - 1.e-9
	nn1 = nm

c******** Position assignment
c******** x0, y0, z0 non-dimensional
c
      fnlcx = float(nlcx)
      fnlcy = float(nlcy)
      fnlcz = float(nlcz)

      do 100 i=1,nlc
      ltop(i)=0
100   continue

c determine in which link cell atom i is situated.
 
      do 110 i=1,nn1
      ix = int( (x0(i)+hmeps)*fnlcx     )
      iy = int( (y0(i)+hmeps)*fnlcy     )
      iz = int( (z0(i)+hmeps)*fnlcz     )
      ip = 1 + ix + nlcx*iy + nlcx*nlcy*iz
c
c assign atom i to link cell ip
c
      j = ltop(ip)
      ltop(ip) = i
      linkmp(i) = j
110   continue

        do 310 ic = 1,nlc
        j = ltop(ic)

         if(j.eq.0) go to 310
311      continue

        j = linkmp(j)
         if(j.gt.0) go to 311
310     continue
c
        return
        end

