#ifndef GUARD_LocationClass_h
#define GUARD_LocationClass_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

double Uniform( void );
double Exponential( void );
double Exponential_g( double lambda );
double Normal( void );
double Gauss( double mu, double sigma );
double ChiSquared( void );
double ChiSquared_g( int k );
double Cauchy_g( double mu, double gamma);
double Gamma( double kappa );
double Gamma_g( double theta, double kappa );
double LogNormal_g( double mu, double sigma );
double LogNormal( double sigma );
double InvGauss_g( double mu, double lambda );
double InvGauss( double lambda );

void init_genrand(unsigned long s);
unsigned long genrand_int32(void);
double genrand_real3(void);

double Uniform( void ){
    
    static int check=0;
    if( check==0 ){ init_genrand(10000); check=1; }
    
    return genrand_real3();
}

double rand_exp( double lambda ){
   return -log(Uniform())/lambda;
}

double rand_normal( double mu, double sigma ){
    double z=sqrt( -2.0*log(Uniform()) ) * sin( 2.0*M_PI*Uniform() );
	return mu + sigma*z;
}

double rand_chi( int k ){
    int i;
    double z,w=0;
     
    for(i=0;i<k;i++){
        z=sqrt( -2.0*log(Uniform()) ) * sin( 2.0*M_PI*Uniform() );
        w+=z*z;
    }
     
    return w;
}

double rand_cauchy( double mu, double gamma){
    return mu + gamma*tan(M_PI*( Uniform()-0.5 ));
}

double rand_gamma( double theta, double kappa ){
    
   int int_kappa;
   double frac_kappa;
    
   int_kappa  = (int)kappa;
   frac_kappa = kappa - (double)int_kappa;
    
   double u,uu;
   double b,p,x_frac,x_int;
   int i;
    
   x_int=0;
   for(i=0;i<int_kappa;i++){
       x_int+=-log(Uniform()); 
   }
    
   if( fabs(frac_kappa) < 0.01 ) x_frac=0;
 
   else{
       b=(exp(1.0)+frac_kappa)/exp(1.0);
       while(1){
        
           u=Uniform();
           p=b*u;
            
           uu=Uniform();
            
           if(p<=1.0){
               x_frac=pow(p,1.0/frac_kappa);
               if(uu<=exp(-x_frac)) break;
           }
            
           else{
               x_frac=-log((b-p)/frac_kappa);
               if(uu<=pow(x_frac,frac_kappa-1.0)) break;
           }
        
       }
   }
    
   return (x_int+x_frac)*theta;
}

double rand_Lnormal( double mu, double sigma ){
   double z= mu + sigma*sqrt(-2.0*log(Uniform()))*sin(2.0*M_PI*Uniform());
   return exp(z);
}

double rand_Igauss( double mu, double lambda ){
   double x,y,w,z;
   x=sqrt(-2.0*log(Uniform()))*sin(2.0*M_PI*Uniform());
   y=x*x;
   w= mu + 0.5*y*mu*mu/lambda -(0.5*mu/lambda)*sqrt(4.0*mu*lambda*y+mu*mu*y*y); 
   z=Uniform();
    
   if( z< mu/(mu+w) )   return w;
   else                return mu*mu/w;
}

#define MT_N 624
#define MT_M 397
#define MATRIX_A 0x9908b0dfUL   
#define UPPER_MASK 0x80000000UL 
#define LOWER_MASK 0x7fffffffUL 

static unsigned long mt[MT_N]; 
static int mti=MT_N+1; 

void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<MT_N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        mt[mti] &= 0xffffffffUL;
    }
}

void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (MT_N>key_length ? MT_N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j;
        mt[i] &= 0xffffffffUL;
        i++; j++;
        if (i>=MT_N) { mt[0] = mt[MT_N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=MT_N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i;
        mt[i] &= 0xffffffffUL;
        i++;
        if (i>=MT_N) { mt[0] = mt[MT_N-1]; i=1; }
    }

    mt[0] = 0x80000000UL;
}

unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};

    if (mti >= MT_N) {
        int kk;

        if (mti == MT_N+1)   
            init_genrand(5489UL); 

        for (kk=0;kk<MT_N-MT_M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<MT_N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(MT_M-MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[MT_N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[MT_N-1] = mt[MT_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
}

double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
}

double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
}

double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 

#endif