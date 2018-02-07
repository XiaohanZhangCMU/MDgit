
 $prntcard
  printf='test4.p'
  rstrtf='test4.r'
  ipatoms=3
  ipitera=-1
  iconst=0
  meamf='../meamf'
 $
 $headcard
  header='Ni '
 $
 $meacard 
        ntypes=1,
        enames(1)='Ni4',
        kodes(1)='library',
        rcut=4.,cmin=0.8
 $
 $initcard
  initf='ni.atm'
  genvel=.t.,
  scale=1.,1,1.0
 $
 $latcard $
 $velcard temp=600,icm=.f. $
 $bndcard ibdtyp=2,idynper=+0,+0,+1,bndmas=1.,1.,1.,dpress=1e6 $
 $neicard nmeth=-2$
 $defcard $
 $fixcard $
 $tmpcard ifxtmp=2 follow=.t. $
 $regcard destmp=600.,tmptim=0.1 $
 $regcard $
 $avecard  eqtim=0.5 $
 $intcard inte=+1, nfmax=1,dt=0.002,outtim=0.01,tottim=0.002  $
 $continue contin=.f. $
