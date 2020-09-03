/** includes */
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <chrono>

/*********************************************************
 *
 * Program starts here, Main function:
 *
 *********************************************************/

int main(int argc, char *argv[] )
{
    //Initial conditions
    double nfkb_n = 0.25;
    double nlrp3_i = 0;
    double nlrp3_a = 0; 
    double nlrp3_b = 0;
    double asc_b = 0;
    double caspase1_b = 0; 
    double gsdmd_c = 0; 
    double il_1b_p = 0; 
    double il_1b_c = 0; 
    double il_1b_e = 0; 
    double il_18_c = 0; 
    double il_18_e =0; 
    double volume_c = 1; 

    //Model constants
    double k_nfkb_ctn = 0.3;
    double k_nfkb_ntc =  0.03;
    double k_nlrp3_ita =  0.07;
    double k_nlrp3_atb = 0.07;
    double k_asc_ftb = 0.02;
    double k_c1_ftb = 0.04;
    double k_il1b_cte = 0.8;
    double k_il18_cte = 0.8;
    double k_vol_c = 0.1;
    // Decay constants
    double d_nlrp3 = 0.002;
    double d_il =0.004;
    // Hill function rates
    double a_nlrp3 = 0.025;
    double a_il1b_p = 0.007;
    double a_gsdmd =0.08;
    double a_il1b_c = 0.8;
    double a_il18 = 0.8;
    double hm_nfkb = 0.3;
    double hm_c1 = 0.3;
    double hex_nfkb = 2.0;
    double hex_c1 = 2.0 ;
   
    double dt=6.0;

    int t=0;

     FILE *result_file;
     result_file=fopen("ODE_data_cpp","a");

    while(volume_c<1.5)
    {

        double F_ib = 1;
        if( nlrp3_b >= 1)
        {F_ib=0;}

        //Update nuclear NFkB (updated backward)
        //double dNFkB_n = (nfkb_n+k_nfkb_ctn*F_ib*dt)/(1+dt*k_nfkb_ntc+k_nfkb_ctn*F_ib*dt);
        // nfkb_n  = dNFkB_n; 
         nfkb_n  = ( nfkb_n +k_nfkb_ctn*F_ib*dt)/(1+dt*k_nfkb_ntc+k_nfkb_ctn*F_ib*dt);
        //printf("NFkBn= %lf .\n", nfkb_n );

        //Set Hill function 1 (same)
        double hill_nfkb = (pow( nfkb_n ,hex_nfkb))/(pow(hm_nfkb,hex_nfkb)+pow( nfkb_n ,hex_nfkb));
        //hill_nfkb = 0.5;
        printf("hill= %lf .\n",hill_nfkb);
        //Update NLRP3 (inactive, active and bound) (updated backward)
        //double dNLRP3_i = (nlrp3_i+dt*a_nlrp3*hill_nfkb)/(1+dt*k_nlrp3_ita+dt*d_nlrp3);
        // nlrp3_i  = dNLRP3_i;
         nlrp3_i  = ( nlrp3_i +dt*a_nlrp3*hill_nfkb)/(1+dt*k_nlrp3_ita+dt*d_nlrp3);
        //printf("NLRP3I= %lf .\n", nlrp3_i );

        //double dNLRP3_a = (nlrp3_a+k_nlrp3_ita*dt*( nlrp3_i ))/(1+dt*k_nlrp3_atb+dt*d_nlrp3);
        // nlrp3_a  = dNLRP3_a;
         nlrp3_a  = ( nlrp3_a +k_nlrp3_ita*dt*( nlrp3_i ))/(1+dt*k_nlrp3_atb+dt*d_nlrp3);

        //double dNLRP3_b = dt * k_nlrp3_atb * F_ib *  nlrp3_a ;
        // nlrp3_b  += dNLRP3_b; //(remains same with backward method)
         nlrp3_b  =  nlrp3_b  + dt * k_nlrp3_atb * F_ib *  nlrp3_a ;
        //printf("NLRP3b= %lf .\n", nlrp3_b );

        //Update bound ASC (updated backward)
        //double dASC_b = (asc_b+dt*k_asc_ftb*(1-F_ib)*( nlrp3_b ))/(1+dt*k_asc_ftb*(1-F_ib)*( nlrp3_b ));
        // asc_b  = dASC_b; 
         asc_b  = ( asc_b  + dt*k_asc_ftb*(1-F_ib)*( nlrp3_b ))/(1+dt*k_asc_ftb*(1-F_ib)*( nlrp3_b ));


        //Update bound caspase1 (updated backward)
        //double dCasp1 = (caspase1_b +dt*k_c1_ftb*( asc_b ))/(1+dt*k_c1_ftb*( asc_b ));
        // caspase1_b  = dCasp1; 
         caspase1_b  = ( caspase1_b  + dt*k_c1_ftb*( asc_b ))/(1+dt*k_c1_ftb*( asc_b )); 

        //Set Hill function 2 (same)
        double hill_caspase1 = (pow( caspase1_b ,hex_c1))/(pow(hm_c1,hex_c1)+pow( caspase1_b ,hex_c1));
        //hill_caspase1 = 0.5; 
        //Update cleaved GSDMD (updated backward)
        //double dGSDMD_c = (gsdmd_c+dt*a_gsdmd*hill_caspase1)/(1+dt*a_gsdmd*hill_caspase1);
        // gsdmd_c  = dGSDMD_c;
         gsdmd_c  = ( gsdmd_c +dt*a_gsdmd*hill_caspase1)/(1+dt*a_gsdmd*hill_caspase1);

        //Set G function (same now that total GSDMD concentration is 1 au of concentration)
        double g_gsdmd =  gsdmd_c /1;

        //Update IL1b (pro, cytoplasmic, external)  We want to relate this to secreted cytokine IL1b (updated backward)
        //double dIL1b_p = (il_1b_p+dt*a_il1b_p*hill_nfkb)/(1+dt*a_il1b_c*hill_caspase1+dt*d_il);
        // il_1b_p  = dIL1b_p;
         il_1b_p  = ( il_1b_p +dt*a_il1b_p*hill_nfkb)/(1+dt*a_il1b_c*hill_caspase1+dt*d_il);

        //double dIL1b_c = (il_1b_c+dt*a_il1b_c*hill_caspase1*( il_1b_p ))/(1+dt*d_il+dt*k_il1b_cte*g_gsdmd);             
        // il_1b_c  = dIL1b_c;
         il_1b_c  = ( il_1b_c +dt*a_il1b_c*hill_caspase1*( il_1b_p ))/(1+dt*d_il+dt*k_il1b_cte*g_gsdmd);      

        //double dIL1b_e = dt * (k_il1b_cte*g_gsdmd* il_1b_c );
        // il_1b_e  += dIL1b_e; //We want to relate this to secreted cytokine IL1b (remains same backward)
         il_1b_e  =  il_1b_e  + dt * (k_il1b_cte*g_gsdmd* il_1b_c );

        //Update IL18 (cytoplasmic, external)(updated backward)
        //double dIL18_c = (il_18_c+dt*a_il18*hill_caspase1*(1-il_18_e))/((1+dt*a_il18*hill_caspase1)*(1+dt*k_il18_cte*g_gsdmd));
        // il_18_c  = dIL18_c;
         il_18_c  = ( il_18_c +dt*a_il18*hill_caspase1*(1- il_18_e ))/((1+dt*a_il18*hill_caspase1)*(1+dt*k_il18_cte*g_gsdmd));

        //double dIL18_e = dt * k_il18_cte*g_gsdmd* il_18_c ;
        // il_18_e  += dIL18_e;//We want to relate this to secreted cytokine IL18 (stays same backward)
         il_18_e  =  il_18_e  +  dt * k_il18_cte*g_gsdmd* il_18_c ;

        //Update cytoplasmic volume (updated backward)
        //double dVol_c = 1/(1-dt * k_vol_c * g_gsdmd);
        // volume_c  *= dVol_c;//We want to relate this to cell cytoplasm volume
         volume_c  =  volume_c /(1-dt * k_vol_c * g_gsdmd);
        //printf("volume %lf \n",  volume_c );

        //if(  volume_c >1.5)
        //{
            // (Adrianne) cell lyses
          //  std::cout<<"Pyroptosis cell burst!"<<std::endl;
        //}
         fprintf(result_file,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",t, nfkb_n, nlrp3_i, nlrp3_a, nlrp3_b, asc_b ,caspase1_b, gsdmd_c, il_1b_p, il_1b_c, il_1b_e, il_18_c, il_18_e, volume_c  );
    
         t=t+dt;
     }

    printf("Cell rupture at t=%d.\n",t );

 
    fclose(result_file);
   return 0;
}


