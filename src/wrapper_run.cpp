//Only works for ni-cs
void MDFrame::wrapper_run(){

  double n6 = 450;

  double sigma_xy  = 0.0;
  double sigma_yz  = 300;
  double sigma_zz  = -2015.254326;
  double sigma_yy  = 2015.254326;
  double sigma_xz  = 848.528137;
  double sigma_xx  = 0.0;


  double newT = 300;
  int maxiter = totalsteps/200;
  double delta_T = (n6-300)*2 / maxiter;
  for (int iter = 0; iter < maxiter; iter++) { 

    double sig_avg_xx = 0;
    double sig_avg_yy = 0;
    double sig_avg_zz = 0;
    double sig_avg_xy = 0;
    double sig_avg_yz = 0;
    double sig_avg_xz = 0;

    double factor = 0;
    
    if ( iter < 50 ) {
       factor = 400.0e3;
    } else if (iter < 100 ) {
       factor = 800.0e3;
    } else if (iter < 150 ) {
       factor = 1200.0e3;
    } else if (iter <= 200)  {
       factor = 1600.0e3;
    } else {
       factor = 2000.0e3;
    }
    
    if ( newT < n6 ) {
    	newT += delta_T;
       _TDES = newT;
    }
    

    for (int subiter = 0; subiter < 20; subiter ++) {
       continue_curstep = 1;
       totalsteps = 10 ;
       run();
       double sig_xx = _TOTSTRESSinMPa[0][0];
       double sig_yy = _TOTSTRESSinMPa[1][1];
       double sig_zz = _TOTSTRESSinMPa[2][2];
       double sig_xy = _TOTSTRESSinMPa[0][1];
       double sig_yz = _TOTSTRESSinMPa[1][2];
       double sig_xz = _TOTSTRESSinMPa[0][2];
       printf("sig_xx = %e sig_yy = %e sig_zz = %e\n", sig_xx, sig_yy, sig_zz);
       printf("sig_yz = %e sig_xz = %e sig_xy = %e\n", sig_yz, sig_xz, sig_xy);

       if ( subiter > 10 ) {
          sig_avg_xx += sig_xx;
          sig_avg_yy += sig_yy;
          sig_avg_zz += sig_zz;
          sig_avg_xy += sig_xy;
          sig_avg_yz += sig_yz;
          sig_avg_xz += sig_xz;
       }
    }
    sig_avg_xx /= 10.0;
    sig_avg_yy /= 10.0;
    sig_avg_zz /= 10.0;
    sig_avg_xy /= 10.0;
    sig_avg_yz /= 10.0;
    sig_avg_xz /= 10.0;

       printf("sig_avg_xx = %e sig_avg_yy = %e sig_avg_zz = %e\n", sig_avg_xx, sig_avg_yy, sig_avg_zz);
       printf("sig_avg_yz = %e sig_avg_xz = %e sig_avg_xy = %e\n", sig_avg_yz, sig_avg_xz, sig_avg_xy);

    if ( iter < maxiter ) {
       double xx_err = sig_avg_xx - sigma_xx;
       double yy_err = sig_avg_yy - sigma_yy;
       double zz_err = sig_avg_zz - sigma_zz;
       double xz_err = sig_avg_xz - sigma_xz;
       double yz_err = sig_avg_yz - sigma_yz;
       double delta_xx = xx_err / factor; 
       double delta_yy = yy_err / factor; 
       double delta_zz = zz_err / factor; 
       double delta_xz = xz_err / factor; 
       input[0] =1; input[1] =1; input[2] = delta_xx; shiftbox(); // ={ 1, 1, delta_xx }; shiftbox();
       input[0] =1; input[1] =1; input[2] = delta_xx; shiftbox(); //input ={ 2, 2, delta_yy }; shiftbox();
       input[0] =1; input[1] =1; input[2] = delta_xx; shiftbox(); //input ={ 3, 3, delta_zz }; shiftbox();
       input[0] =1; input[1] =1; input[2] = delta_xx; shiftbox(); //input ={ 3, 1, delta_xz }; shiftbox();
       swap_velocities(20);
    }
  }
}

void MDFrame::swap_velocities ( double v_ratio ) { 
  int maxAtomsPerturbed = _NP*1.0/v_ratio;
  Vector3 v = _VSR[0];
  for (int i = 0; i < maxAtomsPerturbed ; i++)  {
     int atom =  int(rand()*_NP);
     int ind = atom*3;
     Vector3 v0 = _VSR[ind];
     _VSR[ind] = v;
    // if (i == 0 ) {
    //   puts "vx = $vx, vx0 = $vx0"
    //   puts "vy = $vy, vy0 = $vy0"
    //   puts "vz = $vz, vz0 = $vz0"
    //   puts "..."
    // }
     v = v0;
  }
}

