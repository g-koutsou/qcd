/* qcd_particle_data.c
 *  
 * Collection of some particle data needed for threep_all.c
 * 
 * June 2012
 * Christos Kallidonis
 **************************************/

char *particles_parts[53][4] = {
     {"u" ,   "d",  "z", "z" }, // DUMMY		
     {"u" ,   "d",  "z", "z" }, // PROTON
     {"u" ,   "d",  "z", "z" }, // NEUTRON
     {"u" ,   "d",  "s", "z" }, // LAMBDA
     {"u" ,   "z",  "s", "z" }, // SIGMA_PLUS
     {"u" ,   "d",  "s", "z" }, // SIGMA_ZERO
     {"z" ,   "d",  "s", "z" }, // SIGMA_MINUS
     {"u" ,   "z",  "s", "z" }, // KSI_ZERO
     {"z" ,   "d",  "s", "z" }, // KSI_MINUS
  
     {"u" ,   "z",  "z", "z" }, // DELTA_PLUS_PLUS
     {"u" ,   "d",  "z", "z" }, // DELTA_PLUS
     {"u" ,   "d",  "z", "z" }, // DELTA_ZERO
     {"z" ,   "d",  "z", "z" }, // DELTA_MINUS
     {"u" ,   "z",  "s", "z" }, // SIGMA_STAR_PLUS
     {"u" ,   "d",  "s", "z" }, // SIGMA_STAR_ZERO
     {"z" ,   "d",  "s", "z" }, // SIGMA_STAR_MINUS
     {"u" ,   "z",  "s", "z" }, // KSI_STAR_ZERO
     {"z" ,   "d",  "s", "z" }, // KSI_STAR_MINUS
     {"z" ,   "z",  "s", "z" }, // OMEGA

     {"u" ,   "d",  "z", "c" }, // LAMBDA_C
     {"u" ,   "z",  "z", "c" }, // SIGMA_PLUS_PLUS_C
     {"u" ,   "d",  "z", "c" }, // SIGMA_PLUS_C
     {"z" ,   "d",  "z", "c" }, // SIGMA_ZERO_C
     {"u" ,   "z",  "s", "c" }, // KSI_PLUS_C
     {"u" ,   "z",  "s", "c" }, // KSI_PRIME_PLUS_C
     {"z" ,   "d",  "s", "c" }, // KSI_ZERO_C
     {"z" ,   "d",  "s", "c" }, // KSI_PRIME_ZERO_C
     {"z" ,   "z",  "s", "c" }, // OMEGA_ZERO_C

     {"u" ,   "z",  "z", "c" }, // SIGMA_STAR_PLUS_PLUS_C
     {"u" ,   "d",  "z", "c" }, // SIGMA_STAR_PLUS_C
     {"z" ,   "d",  "z", "c" }, // SIGMA_STAR_ZERO_C
     {"u" ,   "z",  "s", "c" }, // KSI_STAR_PLUS_C
     {"z" ,   "d",  "s", "c" }, // KSI_STAR_ZERO_C
     {"z" ,   "z",  "s", "c" }, // OMEGA_STAR_ZERO_C

     {"u" ,   "z",  "z", "c" }, // KSI_PLUS_PLUS_C_C
     {"z" ,   "d",  "z", "c" }, // KSI_PLUS_C_C
     {"z" ,   "z",  "s", "c" }, // OMEGA_PLUS_C_C

     {"u" ,   "z",  "z", "c" }, // KSI_STAR_PLUS_PLUS_C_C
     {"z" ,   "d",  "z", "c" }, // KSI_STAR_PLUS_C_C
     {"z" ,   "z",  "s", "c" }, // OMEGA_STAR_PLUS_C_C

     {"z" ,   "z",  "z", "c" }, // OMEGA_PLUS_PLUS_C_C_C

//-EXTRA PARTICLES     
     {"u" ,   "d",  "z", "z" }, // PROTON EXTRA
     {"u" ,   "d",  "z", "z" }, // NEUTRON EXTRA
     
     {"u" ,   "z",  "s", "z" }, // KSI_STAR_ZERO  EXTRA
     {"z" ,   "d",  "s", "z" }, // KSI_STAR_MINUS EXTRA
     {"u" ,   "z",  "s", "c" }, // KSI_PLUS_C EXTRA
     {"z" ,   "d",  "s", "c" }, // KSI_ZERO_C EXTRA
     {"u" ,   "z",  "s", "c" }, // KSI_STAR_PLUS_C EXTRA
     {"z" ,   "d",  "s", "c" }, // KSI_STAR_ZERO_C EXTRA
     {"z" ,   "z",  "s", "c" }, // OMEGA_STAR_ZERO_C EXTRA          
     {"u" ,   "z",  "z", "c" }, // KSI_STAR_PLUS_PLUS_C_C EXTRA
     {"z" ,   "d",  "z", "c" }, // KSI_STAR_PLUS_C_C EXTRA
     {"z" ,   "z",  "s", "c" }  // OMEGA_STAR_PLUS_C_C EXTRA
};

//------------------------------------------------------------------
 int particles_pnum[53][4] = {
     {1 ,   1,   0 , 0 }, // DUMMY	 
     {1 ,   1,   0 , 0 }, // PROTON
     {1 ,   1,   0 , 0 }, // NEUTRON
     {1 ,   1,   1 , 0 }, // LAMBDA
     {1 ,   0,   1 , 0 }, // SIGMA_PLUS
     {1 ,   1,   1 , 0 }, // SIGMA_ZERO
     {0 ,   1,   1 , 0 }, // SIGMA_MINUS
     {1 ,   0,   1 , 0 }, // KSI_ZERO
     {0 ,   1,   1 , 0 }, // KSI_MINUS
  
     {1 ,   0,   0 , 0 }, // DELTA_PLUS_PLUS
     {1 ,   1,   0 , 0 }, // DELTA_PLUS
     {1 ,   1,   0 , 0 }, // DELTA_ZERO
     {0 ,   1,   0 , 0 }, // DELTA_MINUS
     {1 ,   0,   1 , 0 }, // SIGMA_STAR_PLUS
     {1 ,   1,   1 , 0 }, // SIGMA_STAR_ZERO
     {0 ,   1,   1 , 0 }, // SIGMA_STAR_MINUS
     {1 ,   0,   1 , 0 }, // KSI_STAR_ZERO
     {0 ,   1,   1 , 0 }, // KSI_STAR_MINUS
     {0 ,   0,   1 , 0 }, // OMEGA

     {1 ,   1,   0 , 1 }, // LAMBDA_C
     {1 ,   0,   0 , 1 }, // SIGMA_PLUS_PLUS_C
     {1 ,   1,   0 , 1 }, // SIGMA_PLUS_C
     {0 ,   1,   0 , 1 }, // SIGMA_ZERO_C
     {1 ,   0,   1 , 1 }, // KSI_PLUS_C
     {1 ,   0,   1 , 1 }, // KSI_PRIME_PLUS_C
     {0 ,   1,   1 , 1 }, // KSI_ZERO_C
     {0 ,   1,   1 , 1 }, // KSI_PRIME_ZERO_C
     {0 ,   0,   1 , 1 }, // OMEGA_ZERO_C

     {1 ,   0,   0 , 1 }, // SIGMA_STAR_PLUS_PLUS_C
     {1 ,   1,   0 , 1 }, // SIGMA_STAR_PLUS_C
     {0 ,   1,   0 , 1 }, // SIGMA_STAR_ZERO_C
     {1 ,   0,   1 , 1 }, // KSI_STAR_PLUS_C
     {0 ,   1,   1 , 1 }, // KSI_STAR_ZERO_C
     {0 ,   0,   1 , 1 }, // OMEGA_STAR_ZERO_C

     {1 ,   0,   0 , 1 }, // KSI_PLUS_PLUS_C_C
     {0 ,   1,   0 , 1 }, // KSI_PLUS_C_C
     {0 ,   0,   1 , 1 }, // OMEGA_PLUS_C_C

     {1 ,   0,   0 , 1 }, // KSI_STAR_PLUS_PLUS_C_C
     {0 ,   1,   0 , 1 }, // KSI_STAR_PLUS_C_C
     {0 ,   0,   1 , 1 }, // OMEGA_STAR_PLUS_C_C

     {0 ,   0,   0 , 1 }, // OMEGA_PLUS_PLUS_C_C_C 

//-EXTRA PARTICLES     
     {1 ,   1,   0 , 0 }, // PROTON EXTRA
     {1 ,   1,   0 , 0 }, // NEUTRON EXTRA
     
     {1 ,   0,   1 , 0 }, // KSI_STAR_ZERO  EXTRA
     {0 ,   1,   1 , 0 }, // KSI_STAR_MINUS EXTRA
     {1 ,   0,   1 , 1 }, // KSI_PLUS_C EXTRA
     {0 ,   1,   1 , 1 }, // KSI_ZERO_C EXTRA
     {1 ,   0,   1 , 1 }, // KSI_STAR_PLUS_C EXTRA
     {0 ,   1,   1 , 1 }, // KSI_STAR_ZERO_C EXTRA
     {0 ,   0,   1 , 1 }, // OMEGA_STAR_ZERO_C EXTRA          
     {1 ,   0,   0 , 1 }, // KSI_STAR_PLUS_PLUS_C_C EXTRA
     {0 ,   1,   0 , 1 }, // KSI_STAR_PLUS_C_C EXTRA
     {0 ,   0,   1 , 1 }  // OMEGA_STAR_PLUS_C_C EXTRA        
	  
 };

//------------------------------------------------------------------

char *particle_names[53] = {
	"dummy",
    "proton",
    "neutron",
    "lambda",
    "sigmaplus",
    "sigmazero",
    "sigmaminus",
    "xizero",
    "ximinus",
//------------------------------
     "deltaplusplus",
     "deltaplus",
    "deltazero",
    "deltaminus",
    "sigmastarplus",
    "sigmastarzero",
    "sigmastarminus",
    "xistarzero",
    "xistarminus",
    "omega",
//------------------------------
    "lambda_cplus",
    "sigma_cplusplus",
    "sigma_cplus",
    "sigma_czero",
    "xi_cplus",
    "xi_prime_cplus",
    "xi_czero",
    "xi_prime_czero",
    "omega_czero",
//------------------------------
    "sigmastar_cplusplus",
    "sigmastar_cplus",
    "sigmastar_czero",
    "xistar_cplus",
    "xistar_czero",
    "omegastar_czero",
//------------------------------
    "xi_ccplusplus",
    "xi_ccplus",
    "omega_ccplus",
//------------------------------
    "xistar_ccplusplus",
    "xistar_ccplus",
    "omegastar_ccplus",
//------------------------------
    "omega_cccplusplus",
//------------------------------

    "proton_extra",
    "neutron_extra",    
    "xistarzero_extra",
    "xistarminus_extra",
    "xi_cplus_extra",
    "xi_czero_extra",    
    "xistar_cplus_extra",
    "xistar_czero_extra",
    "omegastar_czero_extra",
    "xistar_ccplusplus_extra",
    "xistar_ccplus_extra",
    "omegastar_ccplus_extra"
};

//------------------------------------------------------------------
int particles32[53] = {0,
						0,0,0,0,0,0,0,0,
						1,1,1,1,1,1,1,1,1,1,
						0,0,0,0,0,0,0,0,0,
						1,1,1,1,1,1,
						0,0,0,
						1,1,1,
						1,
						0,0,
						1,1,0,0,1,1,1,1,1,1
					  };
