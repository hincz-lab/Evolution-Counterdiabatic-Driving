// Implementing the non-vectorized ABM to Fokker Planck mapping in Mike's notes.
//Code for varying s(t) and implementing sCD(t)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

float Th(float r) //Defining Heaviside Theta function
{
  float theta=(r>=0.) ? (1.) : (0.);
  return theta;
}

float sVarCD(double x) //Defining tanh based gen varying s function
{
  float s,ds,scd;
  s=(float)(0.00075+0.00075*tanh((x-500.)/270.));
  ds=0.00075/(270.*pow(cosh((x-500.)/270.),2.))/0.05;
  scd=s+(ds/pow((pow((0.0008-s),2.)+4.*0.0004*s),0.5));
  return scd;
}

float sVar(double x) //Defining tanh based gen varying s function
{
  float s;
  s=(float)(0.00075+0.00075*tanh((x-500.)/270.));
  return s;
}

float X0() //Using Box Muller transform for inverse sampling of initial eq. distribution p_eq(t=0)
{
  float U1,U2,Z,X;
  float mu=0.5,sigma=0.1816;


  U1= ((float)rand()/(RAND_MAX));
  U2= ((float)rand()/(RAND_MAX));

  Z= pow(-(2.*log(U1)),0.5)*cos(2.*3.14159*U2);
  X=Z*sigma + mu;
  return X;
}

int main()
{

	/* Defining simulation constants and variables
	   K ~ carrying capacity
       d ~ death rate (constant)
       b0~ maximum birth rate
       mu, nu ~mutation rates A->B, B->A
       s: 1+s is relative fitness of A over B
       Ai ~ pop.no of mutant A in i-th generation
       Bi ~ pop.no of mutant B in i-th gen
       Ni ~ total population no Ai+Bi in i-th gen (constrained by K)
bi ~ birth rate of mutant A in i-th gen (will determine Ai+1)

	*/




	//FILE *fp1;
	//FILE *fp2;
  FILE *fp3;
	//FILE *fp4;
	//fp1 = fopen("VxVar.csv","w+");
	//fp2 = fopen("DxVar.csv","w+");
  fp3 = fopen("popAsVar.csv","w+");
	//fp4 = fopen("popBVarCD.csv","w+");


	float tAi,tBi,Ai,Bi,K=10000.;
	float b0=2.,d=0.05,s,mu=0.0004,nu=0.0004,x;


  float bi,sma,sna,smb,snb,rm[3],rn[3];
	int i,j,m,n,pos,k,l;//i,m,n,l counters. pos variable stores dx in the correct bin
  //no.of bins to divide span of x 0->1





  int imax; //no.of generations to run for
  printf("Enter imax \n");
  scanf("%d", &imax);
  //printf("Enter s \n");
  //scanf("%f",&s);
        //printf("\nyour choice of imax= %d \n",imax);

  float **trajxA = (float **)malloc((imax+1) * sizeof(float *));
  for (i=0; i<(imax+1); i++)
      trajxA[i] = (float *)malloc(1000 * sizeof(float));

  //float *trajA=(float*)malloc((imax+1)*sizeof(float));// arrays storing population nos. of each mutant  in each generation */
  //float *trajB=(float*)malloc((imax+1)*sizeof(float));
  //double *trajA=(double*)malloc((imax+1)*sizeof(double));//xi=Ai/(Ai+Bi)


	srand(time(NULL));

  clock_t begin=clock();

	//trajA[0]=Ai; trajB[0]=Bi;
  //trajxA[0][0]=(double)(Ai/(Ai+Bi));
	k=0;

while(k<1000){
  i=0;
  x=X0();
  Ai=(float)floor(x*K);
  Bi=(float)floor((1.-x)*K);
	while(i<=imax){
        trajxA[i][k]=(double)(Ai/(Ai+Bi));
        s=sVar((double)i);
        tAi=Ai;tBi=Bi;
        bi=( (Ai+Bi)<=K) ? (b0*(1.-(Ai+Bi)/K)) : (0.);//calculating birth rate of mut. A in i-th gen


        sma=0.;sna=0.;smb=0.;snb=0.;
        for ( m=1; m<=tAi; m++ ) {
                for(j=0;j<=2;j++){rm[j]= ((float)rand()/(RAND_MAX));};
                sma=(Th(rm[0]-d)*(1.+Th(bi-rm[1])*Th(rm[2]-mu)))+sma;
                smb=(Th(rm[0]-d)*Th(bi-rm[1])*Th(mu-rm[2]))+smb;
                };

        for ( n=1; n<=tBi; n++ ) {
                for(j=0;j<=2;j++){rn[j]= ((float)rand()/(RAND_MAX));};
                sna=(Th(rn[0]-d)*Th(bi/(1.+s)-rn[1])*Th(nu-rn[2]))+sna;
                snb=(Th(rn[0]-d)*(1.+Th(bi/(1.+s)-rn[1])*Th(rn[2]-nu)))+snb;
                };

        Ai=sna+sma;
        Bi=snb+smb;

        //printf("bi= %lf  Ai+1= %lf Bi+1= %lf \n",bi,Ai,Bi);

            //rajA[i][k] =Ai;
            //trajB[i] =Bi;
            //trajxA[i][k]=(double)(Ai/(Ai+Bi));
            //pos=(int)(floor(trajx[i-1]/intv));
            //dx[pos]+=trajx[i]-trajx[i-1];
            //dx2[pos]+=pow((trajx[i]-trajx[i-1]),2.);
            //count[pos]++;
            i++;
		};k++;
  };



		/* ending traj evolution*/
		 //for( l=0; l<=imax; l++){
		//fprintf(fp1,"%lf \n",trajA[l]);
                //fprintf(fp2,"%lf \n",trajB[l]);
                 //};

    //for( j=0; j<50; j++){
      //dx[j]=(count[j]>0) ?(dx[j]/(double)count[j]) : (0.);//calc v(x)= <dx>
      //dx2[j]=(count[j]>0) ? (((dx2[j]/(double)count[j])-pow(dx[j],2.0))/2.) : (0.);//calc D(x)=(<dx^2>-<dx>^2 )/2
      //printf("dx=%.5f , dx^2=%.8f, count=%.1f \n",dx[j],dx2[j],count[j]);
      //fprintf(fp1,"%.2f,%.5f\n",intv*((float)j+1.),dx[j]);
      //fprintf(fp2,"%.2f,%.8f\n",intv*((float)j+1),dx2[j]);

      //printf(fp2,"%lf \n",trajB[l]);
       //};


    for( j=0; j<=imax; j++){
      for(l=0;l<1000;l++){
      fprintf(fp3,"%.5f,",trajxA[j][l]);
      //fprintf(fp4,"%.5f\n",trajB[j]);
       };
       fprintf(fp3,"\n");
     };

     clock_t end=clock();
     double t=((double)(end-begin)/CLOCKS_PER_SEC)/60.;
     printf("run time=%lf mins \n",t);

     //fclose(fp1);
     //fclose(fp2);
     fclose(fp3);
     //fclose(fp4);
     free(trajxA);
     //free(trajA);
     //free(trajB);


    }
