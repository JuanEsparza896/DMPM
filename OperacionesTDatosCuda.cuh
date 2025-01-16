#ifndef DATOSCUDA
#define DATOSCUDA

/*
No Usamos clases o estructuras por lo que no podemos recurrir a sobrecarga de operadores
así que definimos las operaciones para los tipos de datos vectorizados que se usan que
proporciona cuda aquí
*/
template<typename T,typename K>
__device__ __host__ T InitDataType2(K a,K b)
{
    T temp;
    temp.x=a;
    temp.y=b;
    return temp;
}

template<typename T>
__device__ __host__ T InvDataType2(T var)
{
    T temp;
    temp.x=1.0/var.x;
    temp.y=1.0/var.y;
    return temp;
}

template<typename T>
__device__ __host__ T InitDataType3(auto a,auto b,auto c)
{
    T temp;
    temp.x=a;
    temp.y=b;
    temp.z=c;
    return temp;
}

template<typename T>
T InvDataType3(T var)
{
    T temp;
    temp.x=1.0/var.x;
    temp.y=1.0/var.y;
    temp.z=1.0/var.z;
    return temp;
}
#endif