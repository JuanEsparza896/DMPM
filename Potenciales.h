#ifndef HEADER_POT
#define HEADER_POT
void CalculoParamLJ(double sig,double eps,double *param)
{
    //u(r)=4eps*((sig/r)12-(sig/r)6)
    double d2=sig*sig;
    double d6=d2*d2*d2;
    param[0]=4.0*eps*d6;
    param[1]=param[0]*d6;

}

__device__ __forceinline__ double2 InteraccionLJ(int i, int j,double dis, double *param,bool nconf,bool eshift,double r2c=0.0)
{
    double apot=0;
    double2 fuepotshift,val;
    bool gidj=i-j;

    double d2=gidj?(1/dis):0;
    double d6=d2*d2*d2;
    double d12=d6*d6;
    double fue=6.0*(param[1]*2.0*d12-param[0]*d6)*d2;
    if(nconf){
        apot=param[1]*d12-param[0]*d6;
    }
    
    val.x=fue;
    val.y=apot;
    if(eshift){
        d2=1/r2c;
        d6=d2*d2*d2;
        d12=d6*d6;
        fuepotshift.x=6.0*(param[1]*2.0*d12-param[0]*d6)*d2;
        fuepotshift.y=param[1]*d12-param[0]*d6;
        val.x+=fuepotshift.x;
        val.y+=fuepotshift.y;
    }
    return val;
}
#endif
