#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>

#define IPOPN       2000
#define DT          0.01923
#define TS          0
#define DISCARD     52*10
#define DISCARD2    52*100
#define TE          52*200
#define minf(a,b)     (a)>(b)? (a):(b)
#define YEAR        52


// Several random number generator functions listed below are not included in this code. Please specify by yourself.
// Uniform should be random number followed uniform distibution ranged from 0 to 1, and, double Uniform(void)
// poisn should be random number followed poisson distribution, and, int poisn(double par) where par = mean
// rand_gamma should be gamma distributed random number, and, double rand_gamma(double par1, double par2) where par1 = variance/mean; par2 = pow(mean,2.0)/variance;

int main(){
  int popn,popn2, t, agen, cnt, popcnt, popcnt2, popcnt5, popt[20], cnt2, temp, sum, sum4, expn, expn2, expn3, expn4, expn5, expn6, *tempm;
  double mu[100], muHIV[4], tpHIV[4], transHIV[4], tpHSV[4], transHSV[4], transgono[4], transchla[4], transsyp[9], tpgonof, tpchlaf, tpsypf, tempd, marriage, gam1, gam2, gam3, gam4, mean1, var1, mean2, var2;
  double removelink, probcoit_age[51], acqage[51], probcoit_HIV[5], probcoit_HSV[5], marriage_acq, marriage_div;
  char  *stHIV, *stHSV, *stsyp, *stchla, *stgono, *stHIV2, *stHSV2, *stsyp2, *stchla2, *stgono2, *survivem, *mm;
  double *agem, *probacq, *nposclu;
  char **neighm, **neighm2;
  int    *cas_degree, activ_p, classn, div1, remainedge, totalnewedge, bab_cnt, lastm;
  double mr[20], beki, beki2, cur_newedge, denom, random1, ndcdnm;
  int par_cnt, par_cnt2, par_cnt3, par_cnt4;
  int *newedge, newedge1, newedge2;
  double phi_gono_m, phi_gono_f, phi_chla_m, phi_chla_f;
  double tmp20180115, tmp201801152;
  
  neighm = (char **)malloc(IPOPN * sizeof(char *));
  for (cnt = 0; cnt<IPOPN; cnt++){
    neighm[cnt] = (char *)malloc(IPOPN * sizeof(char));
  }
  neighm2 = (char **)malloc(IPOPN * sizeof(char *));
  for (cnt = 0; cnt<IPOPN; cnt++){
    neighm2[cnt] = (char *)malloc(IPOPN * sizeof(char));
  }
  
  stHIV	= (char *)malloc((sizeof(char)) * IPOPN);
  stHSV	= (char *)malloc((sizeof(char)) * IPOPN);
  stgono = (char *)malloc((sizeof(char)) * IPOPN);
  stchla = (char *)malloc((sizeof(char)) * IPOPN);
  stsyp  = (char *)malloc((sizeof(char)) * IPOPN);
  stHIV2 = (char *)malloc((sizeof(char)) * IPOPN);
  stHSV2 = (char *)malloc((sizeof(char)) * IPOPN);
  stgono2 = (char *)malloc((sizeof(char)) * IPOPN);
  stchla2 = (char *)malloc((sizeof(char)) * IPOPN);
  stsyp2 = (char *)malloc((sizeof(char)) * IPOPN);
  survivem	= (char *)malloc((sizeof(char)) * IPOPN);
  mm	= (char *)malloc((sizeof(char)) * IPOPN);
  tempm = (int *)malloc((sizeof(int)) * IPOPN);
  agem = (double *)malloc((sizeof(double)) * IPOPN);
  probacq	= (double *)malloc((sizeof(double)) * IPOPN);
  cas_degree	= (int *)malloc((sizeof(int)) * IPOPN);
  newedge	= (int *)malloc((sizeof(int)) * IPOPN);
  nposclu	= (double *)malloc((sizeof(double)) * IPOPN);
  agen = 100;
  marriage_acq = 0.07*DT;
  marriage_div = 1.0/5.0*DT;
  classn = 20;
  
  marriage 	= 0.53;
  beki    = -10.0+20.0*Uniform();
  beki2   = 20.0*Uniform();
  mean1    = Uniform()*5.0;
  var1     = Uniform()*5.0;
  mean2    = Uniform()*2.5;
  var2     = Uniform()*2.5;
  gam1 = var1/mean1;
  gam2 = pow(mean1,2.0)/var1;
  gam3 = var2/mean2;
  gam4 = pow(mean2,2.0)/var2;
  
  for (cnt = 0; cnt < 5; cnt++){
    mu[cnt] = 0.04 *DT;
  }
  for (cnt = 5; cnt < 69; cnt++){
    mu[cnt] = 0.0026 *DT;
  }
  for (cnt = 69; cnt < 100; cnt++){
    mu[cnt] = 0.0998 *DT;
  }
  muHIV[0] = 0.0;
  muHIV[1] = 0.0;
  muHIV[2] = 0.5 * DT;
  tpHIV[0] = 0.0;
  tpHIV[1] = 0.0360;
  tpHIV[2] = 0.0008;
  tpHIV[3] = 0.0042;
  tpHSV[0] = 0.0;
  tpHSV[1] = 0.004;
  tpHSV[2] = 0.0;
  tpHSV[3] = 0.004;
  tpgonof = 0.46;
  phi_gono_m = 0.64;
  tpchlaf = 0.17;
  phi_chla_m = 0.3;
  tpsypf = 0.2;
  transHIV[0] = 365.0 / 49.0*DT;
  transHIV[1] = 1.0 / 9.0*DT;
  transHIV[2] = 0.0;
  transHSV[0] = 365.0 / 20.0*DT;
  transHSV[1] = 365.0 / 78.5*DT;
  transHSV[2] = 365.0 / 12.8*DT;
  transgono[0] = 6.08*0.5*DT;
  transgono[1] = (365.0/(20.0*7.0)+6.08*0.5)*DT;
  transgono[2] = 365.0 / (20.0*7.0)*DT;
  transgono[3] = 365.0 / (52.0*7.0)*DT;
  transchla[0] = 7.6*0.5*DT;
  transchla[1] = (365.0/(16.0*7.0)+7.6*0.5)*DT;
  transchla[2] = 365.0 / (90.0*7.0)*DT;
  transchla[3] = 365.0 / (520.0*7.0)*DT;
  transsyp[0] = 365.0 / (4.4*7.0) * DT;
  transsyp[1] = 2.86 * DT;
  transsyp[2] = 4.29 * DT;
  transsyp[3] = 365.0 / (6.6*7.0) * DT;
  transsyp[4] = 2.50 * DT;
  transsyp[5] = 365.0 / (15.6*7.0) * DT;
  transsyp[6] = 365.0 / (520.0*7.0) * DT;
  transsyp[7] = 365.0 / (26.0*7.0) * DT;
  transsyp[8] = 365.0 / (52.0*7.0) * DT;
  removelink = 1.0/(14.0/365.0)*DT;
   
  for (cnt = 0; cnt < 5; cnt++){
    probcoit_age[cnt]   = 1.6338;
    probcoit_age[cnt + 5 * 1]   = 1.49296;
    probcoit_age[cnt + 5 * 2]   = 1.35211;
    probcoit_age[cnt + 5 * 3]   = 1.21127;
    probcoit_age[cnt + 5 * 4]   = 1.07042;
    probcoit_age[cnt + 5 * 5]   = 0.929577;
    probcoit_age[cnt + 5 * 6]   = 0.788734;
    probcoit_age[cnt + 5 * 7]   = 0.647886;
    probcoit_age[cnt + 5 * 8]   = 0.507042;
    probcoit_age[cnt + 5 * 9]   = 0.366197;
  }
  
  for (cnt = 0; cnt < 5; cnt++){
    acqage[cnt] = 2.25;
    acqage[cnt + 5 * 1] = 3.18;
    acqage[cnt + 5 * 2] = 1.62;
    acqage[cnt + 5 * 3] = 0.87;
    acqage[cnt + 5 * 4] = 0.66;
    acqage[cnt + 5 * 5] = 0.50;
    acqage[cnt + 5 * 6] = 0.37;
    acqage[cnt + 5 * 7] = 0.25;
    acqage[cnt + 5 * 8] = 0.17;
    acqage[cnt + 5 * 9] = 0.12;
  }
  
  probcoit_HIV[0] = 11.0 * 365.0 / 30.0*DT;
  probcoit_HIV[1] = 11.0 * 365.0 / 30.0*DT;
  probcoit_HIV[2] = 11.0 * 365.0 / 30.0*DT;
  probcoit_HIV[3] = 7.0  * 365.0 / 30.0*DT;
  
  probcoit_HSV[0] = 11.0 * 365.0 / 30.0*DT;
  probcoit_HSV[1] = 11.0 * 365.0 / 30.0*DT;
  probcoit_HSV[2] = 11.0 * 365.0 / 30.0*DT;
  probcoit_HSV[3] = 11.0 * 365.0 / 30.0*DT;
  
  for (popcnt = 0; popcnt < IPOPN; popcnt++){
    survivem[popcnt] = 0;
    mm[popcnt] = 0;
    agem[popcnt] = 0.0;
    cas_degree[popcnt] = 0;
    for (popcnt2 = 0; popcnt2 < IPOPN; popcnt2++){
      neighm[popcnt][popcnt2] = 0;
      neighm[popcnt2][popcnt] = 0;
      neighm2[popcnt][popcnt2] = 0;
      neighm2[popcnt2][popcnt] = 0;
    }
  }
  
  for (cnt = 0; cnt < IPOPN; cnt++){
    survivem[cnt] = 1;
    stHIV[cnt] = 0;
    stHSV[cnt] = 0;
    stgono[cnt] = 0;
    stchla[cnt] = 0;
    stsyp[cnt] = 0;
    stHIV2[cnt] = 0;
    stHSV2[cnt] = 0;
    stgono2[cnt] = 0;
    stchla2[cnt] = 0;
    stsyp2[cnt] = 0;
    probacq[cnt] = 2.0*rand_gamma(gam1, gam2)/(1.0/removelink+1.0/DT);
    
  }
  
  popt[0] = 138;
  popt[1] = 272;
  popt[2] = 406;
  popt[3] = 538;
  popt[4] = 668;
  popt[5] = 796;
  popt[6] = 924;
  popt[7] = 1050;
  popt[8] = 1174;
  popt[9] = 1296;
  popt[10] = 1416;
  popt[11] = 1536;
  popt[12] = 1652;
  popt[13] = 1768;
  popt[14] = 1844;
  popt[15] = 1896;
  popt[16] = 1930;
  popt[17] = 1954;
  popt[18] = 1970;
  popt[19] = 2000;
  
  mr[0] = 0.0;
  mr[1] = 0.0;
  mr[2] = 0.0;
  mr[3] = 0.23;
  mr[4] = 0.50;
  mr[5] = 0.62;
  mr[6] = 0.67;
  mr[7] = 0.69;
  mr[8] = 0.70;
  mr[9] = 0.70;
  mr[10] = 0.70;
  mr[11] = 0.70;
  mr[12] = 0.70;
  mr[13] = 0.70;
  mr[14] = 0.70;
  mr[15] = 0.70;
  mr[16] = 0.70;
  mr[17] = 0.70;
  mr[18] = 0.70;
  mr[19] = 0.70;
  
  for (cnt = 0; cnt < popt[1]; cnt++){
    agem[cnt] = Uniform() * 5.0;
  }
  
  for (cnt2 = 1; cnt2<20; cnt2++){
    for (cnt = popt[cnt2 - 1]; cnt < popt[cnt2]; cnt++){
      agem[cnt] = Uniform() * 5.0 + (double)(cnt2)* 5.0;
    }
  }
    
  for (t = TS; t < TE; t++){
    
    for (popcnt = 0; popcnt < IPOPN; popcnt++){
      if (survivem[popcnt] == 1){
        ndcdnm = 0.0;
        for (popcnt2 = 0; popcnt2 < IPOPN; popcnt2++){
          if (survivem[popcnt2] == 1){
            ndcdnm += (double)(neighm[popcnt][popcnt2]);
          }
        }
        cas_degree[popcnt] = (int)(ndcdnm);
      }
    }
    
    reselectnewedge:;
    sum = 0;
    popn = 0;
    popn2 = 0;
    bab_cnt = 0;
    totalnewedge=0;
    for (popcnt = 0; popcnt < IPOPN; popcnt++){
      newedge[popcnt] = 0;
      if (survivem[popcnt] == 1){
        popn++;
        if(agem[popcnt]>15.0 &&agem[popcnt] < 65.0){
          popn2++;
          newedge[popcnt] = poisn(probacq[popcnt]*acqage[(int)(floor(agem[popcnt]-15.0))]);
          totalnewedge += newedge[popcnt];
          if (newedge[popcnt] > 0) {
            bab_cnt++;
          }
        }
      }
    }
    if((totalnewedge%2)==1) goto reselectnewedge;
    
    remainedge = totalnewedge;
    while (remainedge > 1){
      denom = (double)(totalnewedge);
      cur_newedge = 0.0;
      for (popcnt = 0; popcnt < IPOPN; popcnt++){
        if(survivem[popcnt] == 1 &&agem[popcnt] > 15.0 &&agem[popcnt] < 65.0 && newedge[popcnt]>0){
          cur_newedge += (double)(newedge[popcnt])/denom;
          if(Uniform() <= cur_newedge){
            newedge1 = popcnt;
            goto firstend;
          }
          lastm = popcnt;
        }
      }
      newedge1 = lastm;
      firstend:;
      denom = 0.0;
      for (popcnt = 0; popcnt < IPOPN; popcnt++){
        if(survivem[popcnt] == 1 &&agem[popcnt] > 15.0 &&agem[popcnt] < 65.0 && popcnt != newedge1 && newedge[popcnt]>0){
          nposclu[popcnt] = 0.0;
          for (popcnt2 = 0; popcnt2 < IPOPN; popcnt2++){
            if(survivem[popcnt2] == 1 &&agem[popcnt2] > 15.0 &&agem[popcnt2] < 65.0 && popcnt2 != newedge1 && popcnt != popcnt2 && newedge[popcnt2]>0 && (neighm[newedge1][popcnt2]+neighm2[newedge1][popcnt2])>0&& (neighm[popcnt][popcnt2]+neighm2[popcnt][popcnt2])>0){
              nposclu[popcnt] += 1.0;
            }
          }
          denom += (double)(newedge[popcnt])*pow((fabs((double)(cas_degree[popcnt]-cas_degree[newedge1]))+1.0),beki)*pow(((double)(nposclu[popcnt])+1.0),beki2);
        }
      }
      cur_newedge = 0.0;
      for (popcnt = 0; popcnt < IPOPN; popcnt++){
        if(survivem[popcnt] == 1 &&agem[popcnt] > 15.0 &&agem[popcnt] < 65.0 && popcnt != newedge1 && newedge[popcnt]>0){
          cur_newedge += (double)(newedge[popcnt])*pow((fabs((double)(cas_degree[popcnt]-cas_degree[newedge1]))+1.0),beki)*pow(((double)(nposclu[popcnt])+1.0),beki2)/denom;
          if(Uniform() <= cur_newedge){
            newedge2 = popcnt;
            goto secondend;
          }
          lastm = popcnt;
        }
      }
      newedge2 = lastm;
      secondend:;
      if(newedge1==newedge2) goto reselectnewedge;
      
      neighm[newedge1][newedge2] = 1;
      neighm[newedge2][newedge1] = 1;
      newedge[newedge1] -= 1;
      newedge[newedge2] -= 1;
      remainedge -= 2;
    }
    
    for (popcnt = 0; popcnt < IPOPN - 1; popcnt++){
      if (survivem[popcnt] == 1){
        for (popcnt2 = popcnt + 1; popcnt2 < IPOPN; popcnt2++){
          if (survivem[popcnt2] == 1){
            if (neighm[popcnt][popcnt2] == 1 || neighm[popcnt2][popcnt] == 1){
              if (Uniform() < removelink){
                neighm[popcnt][popcnt2] = 0;
                neighm[popcnt2][popcnt] = 0;
              }
            }
          }
        }
      }
    }
    
    for (popcnt = 0; popcnt < IPOPN; popcnt++){
      if(survivem[popcnt] == 1 && mm[popcnt] == 0 && agem[popcnt] >= 15.0){
        for (popcnt5 = 0; popcnt5 < IPOPN; popcnt5++){
          if (popcnt != popcnt5 && survivem[popcnt5] == 1 && mm[popcnt5] == 0 && agem[popcnt5] >= 15.0){
            if(Uniform()<(marriage_acq)){
              neighm2[popcnt][popcnt5] = 1;
              neighm2[popcnt5][popcnt] = 1;
              mm[popcnt] = 1;
              mm[popcnt5] = 1;
              tempd = rand_gamma(gam3, gam4);
              probacq[popcnt] = 2.0*tempd/(1.0/removelink+1.0/DT);;
              tempd = rand_gamma(gam3, gam4);
              probacq[popcnt5] = 2.0*tempd/(1.0/removelink+1.0/DT);;
              
            }
          }
        }
      }
    }
    
    for (popcnt = 0; popcnt < IPOPN; popcnt++){
      if(survivem[popcnt] == 1 && mm[popcnt] == 1){
        if(Uniform()<(marriage_div/2.0)){
          mm[popcnt]  = 0;
          tempd = rand_gamma(gam1, gam2);
          probacq[popcnt] = 2.0*tempd/(1.0/removelink+1.0/DT);
          for (popcnt2 = 0; popcnt2 < IPOPN; popcnt2++){
            if(neighm2[popcnt][popcnt2] == 1){
              mm[popcnt2]  = 0;
              tempd = rand_gamma(gam1, gam2);
              probacq[popcnt2] = 2.0*tempd/(1.0/removelink+1.0/DT);
              neighm2[popcnt][popcnt2] = 0;
              neighm2[popcnt2][popcnt] = 0;
            }
          }
        }
      }
    }
    
    for (popcnt = 0; popcnt < IPOPN; popcnt++){
      if (survivem[popcnt] == 1){
        if (floor(agem[popcnt]) < 100.0){
          if (Uniform() < mu[(int)(floor(agem[popcnt]))]){
            if(mm[popcnt]==1){
              mm[popcnt] = 0;
              for (popcnt2 = 0; popcnt2 < IPOPN; popcnt2++){
                if(neighm2[popcnt][popcnt2] == 1){
                  mm[popcnt2]  = 0;
                  tempd = rand_gamma(gam1, gam2);
                  probacq[popcnt2] = 2.0*tempd/(1.0/removelink+1.0/DT);
                }
              }
            }
            survivem[popcnt] = 1;
            stHSV[popcnt] = 0;
            stHIV[popcnt] = 0;
            stgono[popcnt] = 0;
            stchla[popcnt] = 0;
            stsyp[popcnt] = 0;
            stHIV2[popcnt] = 0;
            stHSV2[popcnt] = 0;
            stgono2[popcnt] = 0;
            stchla2[popcnt] = 0;
            stsyp2[popcnt] = 0;
            
            tempd = rand_gamma(gam1, gam2);
            probacq[popcnt] = 2.0*tempd/(1.0/removelink+1.0/DT);
            agem[popcnt] = 0.0;
            cas_degree[popcnt] =0;
            for (popcnt2 = 0; popcnt2 < IPOPN; popcnt2++){
              if(neighm[popcnt][popcnt2]==1){
                neighm[popcnt][popcnt2] = 0;
                neighm[popcnt2][popcnt] = 0;
              }
              neighm2[popcnt][popcnt2] = 0;
              neighm2[popcnt2][popcnt] = 0;
            }
          }
          else{
            agem[popcnt] += DT;
          }
        }
        else{
          if (Uniform() < mu[99]){
            survivem[popcnt] = 0;
            if(mm[popcnt]==1){
              mm[popcnt] = 0;
              for (popcnt2 = 0; popcnt2 < IPOPN; popcnt2++){
                if(neighm2[popcnt][popcnt2] == 1){
                  mm[popcnt2]  = 0;
                  tempd = rand_gamma(gam1, gam2);
                  probacq[popcnt2] = 2.0*tempd/(1.0/removelink+1.0/DT);
                }
              }
            }
            survivem[popcnt] = 1;
            stHSV[popcnt] = 0;
            stHIV[popcnt] = 0;
            stgono[popcnt] = 0;
            stchla[popcnt] = 0;
            stsyp[popcnt] = 0;
            stHIV2[popcnt] = 0;
            stHSV2[popcnt] = 0;
            stgono2[popcnt] = 0;
            stchla2[popcnt] = 0;
            stsyp2[popcnt] = 0;
            
            tempd = rand_gamma(gam1, gam2);
            probacq[popcnt] = 2.0*tempd/(1.0/removelink+1.0/DT);
            agem[popcnt] = 0.0;
            cas_degree[popcnt] =0;
            for (popcnt2 = 0; popcnt2 < IPOPN; popcnt2++){
              if(neighm[popcnt][popcnt2]==1){
                neighm[popcnt][popcnt2] = 0;
                neighm[popcnt2][popcnt] = 0;
              }
              neighm2[popcnt][popcnt2] = 0;
              neighm2[popcnt2][popcnt] = 0;
            }
          }
          else{
            agem[popcnt] += DT;
          }
        }
      }
    }
    
    if(t == DISCARD){
      
      activ_p = 0;
      while(activ_p < 10){
        popcnt2 = (int)(Uniform()*(double)(IPOPN));
        if (survivem[popcnt2] == 1 && agem[popcnt2]>15.0 && agem[popcnt2]<65.0){
          stHIV[popcnt2] = 1;
          activ_p++;
        }
      }
      activ_p = 0;
      while(activ_p < 10){
        popcnt2 = (int)(Uniform()*(double)(IPOPN));
        if (survivem[popcnt2] == 1 && agem[popcnt2]>15.0 && agem[popcnt2]<65.0){
          stHSV[popcnt2] = 1;
          activ_p++;
        }
      }
      activ_p = 0;
      while(activ_p < 10){
        popcnt2 = (int)(Uniform()*(double)(IPOPN));
        if (survivem[popcnt2] == 1 && agem[popcnt2]>15.0 && agem[popcnt2]<65.0){
          stgono[popcnt2] = 2;
          activ_p++;
        }
      }
      activ_p = 0;
      while(activ_p < 10){
        popcnt2 = (int)(Uniform()*(double)(IPOPN));
        if (survivem[popcnt2] == 1 && agem[popcnt2]>15.0 && agem[popcnt2]<65.0){
          stchla[popcnt2] = 2;
          activ_p++;
        }
      }
      activ_p = 0;
      while(activ_p < 10){
        popcnt2 = (int)(Uniform()*(double)(IPOPN));
        if (survivem[popcnt2] == 1 && agem[popcnt2]>15.0 && agem[popcnt2]<65.0){
          stsyp[popcnt2] = 1;
          activ_p++;
        }
      }
    }
    
    if(t > DISCARD){
      expn2 = poisn(1.0*DT);
      expn3 = poisn(1.0*DT);
      expn4 = poisn(1.0*DT);
      expn5 = poisn(1.0*DT);
      expn6 = poisn(1.0*DT);
      for (popcnt = 0; popcnt<expn2; popcnt++){
        popcnt2 = (int)(Uniform()*(double)(IPOPN));
        if (survivem[popcnt2] == 1 && agem[popcnt2]>15.0 && agem[popcnt2] < 65.0){
          stHIV[popcnt2] = 1;
          stHIV2[popcnt2] = 1;
        }
      }
      for (popcnt = 0; popcnt<expn3; popcnt++){
        popcnt2 = (int)(Uniform()*(double)(IPOPN));
        if (survivem[popcnt2] == 1 && agem[popcnt2]>15.0 && agem[popcnt2] < 65.0){
          stHSV[popcnt2] = 1;
          stHSV2[popcnt2] = 1;
        }
      }
      for (popcnt = 0; popcnt<expn4; popcnt++){
        popcnt2 = (int)(Uniform()*(double)(IPOPN));
        if (survivem[popcnt2] == 1 && agem[popcnt2]>15.0 && agem[popcnt2] < 65.0){
          stgono[popcnt2] = 2;
          stgono2[popcnt2] = 1;
        }
      }
      for (popcnt = 0; popcnt<expn5; popcnt++){
        popcnt2 = (int)(Uniform()*(double)(IPOPN));
        if (survivem[popcnt2] == 1 && agem[popcnt2]>15.0 && agem[popcnt2] < 65.0){
          stchla[popcnt2] = 2;
          stchla2[popcnt2] = 1;
        }
      }
      for (popcnt = 0; popcnt<expn6; popcnt++){
        popcnt2 = (int)(Uniform()*(double)(IPOPN));
        if (survivem[popcnt2] == 1 && agem[popcnt2]>15.0 && agem[popcnt2] < 65.0){
          stsyp[popcnt2] = 1;
          stsyp2[popcnt2] = 1;
        }
      }
    }
    
    if (t > DISCARD){
      for (popcnt = 0; popcnt < IPOPN; popcnt++){
        if (survivem[popcnt] == 1 && agem[popcnt]>15.0 && agem[popcnt]<65.0){
          for (popcnt2 = 0; popcnt2 < IPOPN; popcnt2++){
            if (survivem[popcnt2] == 1 && agem[popcnt2]>15.0 && agem[popcnt2] < 65.0){
              if (neighm[popcnt][popcnt2] == 1 || neighm2[popcnt][popcnt2] == 1){
                tempd = minf(minf(minf(probcoit_HIV[stHIV[popcnt]], probcoit_HIV[stHIV[popcnt2]]), probcoit_HSV[stHSV[popcnt]]), probcoit_HSV[stHSV[popcnt2]])*(probcoit_age[(int)(floor(agem[popcnt]-15.0))]+probcoit_age[(int)(floor(agem[popcnt2]-15.0))])/2.0;
                if (stHIV[popcnt] != 0 && stHIV[popcnt2] == 0){
                  if (Uniform() < (tempd * tpHIV[stHIV[popcnt]])){
                    stHIV[popcnt2] = 1;
                    stHIV2[popcnt2] = 2;
                  }
                }
                if (stHSV[popcnt] != 0 && stHSV[popcnt2] == 0){
                  if (Uniform() < (tempd * tpHSV[stHSV[popcnt]])){
                    stHSV[popcnt2] = 1;
                    stHSV2[popcnt2] = 2;
                  }
                }
                if (((stgono[popcnt] == 1)||(stgono[popcnt] == 2)) && stgono[popcnt2] == 0){
                  if (Uniform() < (tempd * tpgonof)){
                    if(Uniform() < phi_gono_m){
                      stgono[popcnt2] = 1;
                      stgono2[popcnt2] = 2;
                    }
                    else{
                      stgono[popcnt2] = 2;
                      stgono2[popcnt2] = 2;
                    }
                  }
                }
                if (((stchla[popcnt] == 1)||(stchla[popcnt] == 2)) && stchla[popcnt2] == 0){
                  if (Uniform() < (tempd * tpchlaf)){
                    if(Uniform() < phi_chla_m){
                      stchla[popcnt2] = 1;
                      stchla2[popcnt2] = 2;
                    }
                    else{
                      stchla[popcnt2] = 2;
                      stchla2[popcnt2] = 2;
                    }
                  }
                }
                if (((stsyp[popcnt] == 2)||(stsyp[popcnt] == 3)) && stsyp[popcnt2] == 0){
                  if (Uniform() < (tempd * tpsypf)){
                    stsyp[popcnt2] = 1;
                    stsyp2[popcnt2] = 2;
                  }
                }
              }
            }
          }
        }
      }
      
      for (popcnt = 0; popcnt < IPOPN; popcnt++){
        if (survivem[popcnt] == 1){
          tmp201801152 = stHIV[popcnt];
          if (stHIV[popcnt] == 1){
            if (Uniform() < transHIV[0]) tmp201801152 = 2;
          }
          if (stHIV[popcnt] == 2){
            if (Uniform() < transHIV[1]) tmp201801152 = 3;
          }
          if (stHIV[popcnt] == 3){
            if (Uniform() < muHIV[2]){
              survivem[popcnt] = 1;
              tmp201801152 = 0;
              stHIV[popcnt] = 0;
              stHSV[popcnt] = 0;
              stgono[popcnt] = 0;
              stchla[popcnt] = 0;
              stsyp[popcnt] = 0;
              stHIV2[popcnt] = 0;
              stHSV2[popcnt] = 0;
              stgono2[popcnt] = 0;
              stchla2[popcnt] = 0;
              stsyp2[popcnt] = 0;
              agem[popcnt] = 0.0;
              tempd = rand_gamma(gam1, gam2);
              probacq[popcnt] = 2.0*tempd/(1.0/removelink+1.0/DT);
              cas_degree[popcnt] =0;
              for (popcnt2 = 0; popcnt2 < IPOPN; popcnt2++){
                neighm[popcnt][popcnt2] = 0;
                neighm[popcnt2][popcnt] = 0;
              }
              if(mm[popcnt]==1){
                mm[popcnt] = 0;
                for (popcnt2 = 0; popcnt2 < IPOPN; popcnt2++){
                  if(neighm2[popcnt][popcnt2] == 1){
                    mm[popcnt2]  = 0;
                    tempd = rand_gamma(gam1, gam2);
                    probacq[popcnt2] = 2.0*tempd/(1.0/removelink+1.0/DT);
                    neighm2[popcnt][popcnt2] = 0;
                    neighm2[popcnt2][popcnt] = 0;
                  }
                }
              }
            }
          }
          stHIV[popcnt] = tmp201801152;
          tmp201801152 = stHSV[popcnt];
          if (stHSV[popcnt] == 1){
            if (Uniform() < transHSV[0]) tmp201801152 = 2;
          }
          if (stHSV[popcnt] == 2){
            if (Uniform() < transHSV[1]) tmp201801152 = 3;
          }
          if (stHSV[popcnt] == 3){
            if (Uniform() < transHSV[2]) tmp201801152 = 2;
          }
          stHSV[popcnt] = tmp201801152;
          
          tmp201801152 = stgono[popcnt];
          if (stgono[popcnt] == 1){
            if (Uniform() < (transgono[0] + transgono[1])){
              if(Uniform() < (transgono[0]/(transgono[0] + transgono[1]))) tmp201801152 = 0;
              else tmp201801152 = 3;
            }
          }
          if (stgono[popcnt] == 2){
            if (Uniform() < transgono[2]) tmp201801152 = 3;
          }
          if (stgono[popcnt] == 3){
            if (Uniform() < transgono[3]) tmp201801152 = 0;
          }
          stgono[popcnt] = tmp201801152;
          
          tmp201801152 = stchla[popcnt];
          if (stchla[popcnt] == 1){
            if (Uniform() < (transchla[0] + transchla[1])){
              if(Uniform() < (transchla[0]/(transchla[0] + transchla[1]))) tmp201801152 = 0;
              else tmp201801152 = 3;
            }
          }
          if (stchla[popcnt] == 2){
            if (Uniform() < transchla[2]) tmp201801152 = 3;
          }
          if (stchla[popcnt] == 3){
            if (Uniform() < transchla[3]) tmp201801152 = 0;
          }
          stchla[popcnt] = tmp201801152;
          
          tmp201801152 = stsyp[popcnt];
          if (stsyp[popcnt] == 1){
            if (Uniform() < transsyp[0]) tmp201801152 = 2;
          }
          if (stsyp[popcnt] == 2){
            if (Uniform() < (transsyp[1] + transsyp[2] + transsyp[3])){
              tmp20180115 = Uniform();
              if(tmp20180115 < (transsyp[1]/(transsyp[1] + transsyp[2] + transsyp[3]))) tmp201801152 = 0;
              else if(((transsyp[1]/(transsyp[1] + transsyp[2] + transsyp[3])) <= tmp20180115) && (tmp20180115 < ((transsyp[1] + transsyp[2])/(transsyp[1] + transsyp[2] + transsyp[3])))) tmp201801152 = 5;
              else tmp201801152 = 3;
            }
          }
          if (stsyp[popcnt] == 3){
            if (Uniform() < (transsyp[4] + transsyp[5])){
              tmp20180115 = Uniform();
              if(tmp20180115 < (transsyp[4]/(transsyp[4] + transsyp[5]))) tmp201801152 = 5;
              else tmp201801152 = 4;
            }
          }
          if (stsyp[popcnt] == 4){
            if (Uniform() < transsyp[6]) tmp201801152 = 6;
          }
          if (stsyp[popcnt] == 5){
            if (Uniform() < transsyp[7]) tmp201801152 = 0;
          }
          if (stsyp[popcnt] == 6){
            if (Uniform() < transsyp[8]) tmp201801152 = 0;
          }
          stsyp[popcnt] = tmp201801152;
        }
      }
    }
  }
}

