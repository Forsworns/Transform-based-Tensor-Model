//
// Created by jcfei on 4/30/18.
//

#ifndef TENSOR_TENSOR_H
#define TENSOR_TENSOR_H
#include <fftw3.h>

#include "string.h"
#include <iostream>
#include <armadillo>
using namespace std;
using namespace arma;

template <class T>
class Tensor{
public:

    Tensor<T>(int,int,int); //创建三维张量
    Tensor<T>();
    ~Tensor<T>();

    Tensor<T> & operator=(const Tensor&);
    Tensor<T>& operator+=(const Tensor&);
    Tensor<T>& operator-=(const Tensor&);
    Tensor<T>& operator*=(const Tensor&);
    template <class T1> Tensor<T>& operator*=(T1);
    template <class T1> Tensor<T>& operator/=(T1);

    inline T& operator()(int x, int y, int z) {
        if (x<n1 && y<n2 && z<n3) {
            return p[x][y][z];
        }
        else {
            cout<<"Warning: Index beyond range."<<endl;
            exit(0);
        }
    }

    static Tensor<T> rand(int, int, int); //零张量
    static Tensor<T> Identity(int, int, int); //零张量
//    Tensor<T>(int,int,int,int); //创建三维张量

    //    private:
    int n1,n2,n3;
    T ***p;
    T*** allocSpace();
};

//get the size of tensor
template <class T>
int * getsize(const Tensor<T> &a);

//fiber
template <class T>
vec fiber(const Tensor<T> & t, int m, int n ,int order);

//slice
template <class T>
fmat slice(const Tensor<T> & t, int m, int order);

//tensor 转 matrix
template<class T>
mat ten2mat(Tensor<T> & a , int b);

//tensor 转 vector
template<class T>
vec ten2vec(Tensor<T> & a);

//norm
template<class T>
double norm(Tensor<T> &a);       //写在类外时，不能在前面加friend.   //为什么要写两个Tensor,第一个tensor是表示类型

//转置运算；不支持复数暂时
template<class T>
Tensor<T> Transpose(Tensor<T> &a);

//Inner product
template<class T1, class T2>
double dotProduct(Tensor <T1>a, Tensor<T2> b);

//t-prod of tensor
template<class T, class T1, class T2>
Tensor<T> tprod(Tensor<T1> &a, Tensor<T2> &b);

//返回结构
template<class T>
struct tuckercore{
    Tensor<T> & core;
    fmat u1,u2,u3;
};

//HOSVD
template<class T>
tuckercore<T> HOSVD(Tensor <T> & a, int  r1, int  r2, int r3);

//n-mode product
template<class T>
mat ttm(Tensor<T> & a, mat & b, int c);

/**************************************
*******Realization of functions********
**************************************/
//定义三维张量，初始化为0
template <class T>
Tensor<T>::Tensor(int n1, int n2, int n3) : n1(n1), n2(n2),n3(n3)
{
    allocSpace();
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            for(int k=0; k<n3; ++k) {
                p[i][j][k] = 0;  //randi...
            }
        }
    }
}

//template<class T>
//Tensor<T>::Tensor(int n1, int n2, int n3, int n4): n1(n1), n2(n2),n3(n3) {
//    if(n4==0) {
//        allocSpace();
//        for (int i = 0; i < n1; ++i) {
//            for (int j = 0; j < n2; ++j) {
//                for(int k=0; k<n3; ++k) {
//                    p[i][j][k] = 0;  //randi...
//                }
//            }
//        }
//    }
//    if(n4==1){
//        if(n1==n2){
//            allocSpace();
//            for (int i = 0; i < n1; ++i) {
//                for (int j = 0; j < n2; ++j) {
//                    for(int k=0; k<n3; ++k) {
//                        if(k==0 && i==j) {
//                            p[i][j][k] = 1;
//                        }
//                        else {
//                            p[i][j][k]=0;
//                        }
//                    }
//                }
//            }
//        }
//        else{
//            cout<<"Warning: no identity tensor";
//            exit(0);
//        }
//    }
//}

//析构函数
template <class T>
Tensor<T>::~Tensor() {
    for (int i=0; i<n1;i++) {
        for (int j=0;j<n2;j++){
            delete[] p[i][j];
        }
        delete[] p[i];
    }
    delete[] p;
    p=NULL;
//    cout<<"destructing T"<<endl;
}

template<class T>
Tensor<T>::Tensor() {
    allocSpace();
    p[0][0][0]=1;
}

//分配内存空间（三维动态数组）
template <class T>
T*** Tensor<T>::allocSpace()
{
    p =new T**[n1];
    for (int i=0; i<n1;i++) {
        p[i] = new T *[n2];
        for (int j=0;j<n2;j++){
            p[i][j]=new T [n3];
        }
    }
    return p;
}

template<class T>
Tensor<T> &Tensor<T>::operator=(const Tensor & a) {
    if(this == &a){
        return *this;
    }
    if(n1 != a.n1 || n2 != a.n2 || n3 != a.n3){
        for (int i=0; i<n1;i++) {
            for (int j=0;j<n2;j++){
                delete[] p[i][j];
            }
            delete[] p[i];
        }
        delete[] p;
        n1=a.n1;n2=a.n2;n3=a.n3;
        allocSpace();
    }
    for (int i=0; i<n1; ++i){
        for (int j=0; j<n2;++j){
            for(int k=0;k<n3;++k){
                p[i][j][k]=a.p[i][j][k];
            }
        }
    }
    return *this;
}

template<class T>
Tensor<T> &Tensor<T>::operator+=(const Tensor &a) {
    if (n1 == a.n1 && n2 == a.n2 && n3 == a.n3){
        for (int i=0; i<a.n1; ++i){
            for (int j=0; j<a.n2;++j){
                for(int k=0;k<a.n3;++k){
                    p[i][j][k]+=a.p[i][j][k];
                }
            }
        }
        return *this;}
    else{
        cout<<"check size match"<<endl;
        exit(0);
    }
}

template<class T>
Tensor<T> &Tensor<T>::operator-=(const Tensor &a) {
    if (n1 == a.n1 && n2 == a.n2 && n3 == a.n3){
        for (int i=0; i<a.n1; ++i){
            for (int j=0; j<a.n2;++j){
                for(int k=0;k<a.n3;++k){
                    p[i][j][k]-=a.p[i][j][k];
                }
            }
        }
        return *this;}
    else{
        cout<<"check size match"<<endl;
        exit(0);
    }
}

template<class T>
Tensor<T> &Tensor<T>::operator*=(const Tensor &a) {
    if (n1 == a.n1 && n2 == a.n2 && n3 == a.n3){
        for (int i=0; i<a.n1; ++i){
            for (int j=0; j<a.n2;++j){
                for(int k=0;k<a.n3;++k){
                    p[i][j][k]*=a.p[i][j][k];
                }
            }
        }
        return *this;}
    else{
        cout<<"check size match"<<endl;
        exit(0);
    }
}

template<class T>
template<class T1>
Tensor<T> &Tensor<T>::operator/=(T1 num) {
    if (num != 0){
        for (int i=0; i<n1; ++i){
            for (int j=0; j<n2;++j){
                for(int k=0;k<n3;++k){
                    p[i][j][k]/=num;
                }
            }
        }
    }
    else {
        cout << "divide number is zero."<< endl;
        exit(0);
    }
}

template<class T>
template<class T1>
Tensor<T> &Tensor<T>::operator*=(T1 num) {
    for (int i=0; i<n1; ++i){
        for (int j=0; j<n2;++j){
            for(int k=0;k<n3;++k){
                p[i][j][k]/=num;
            }
        }
    }
}


template<class T>
Tensor<T> Tensor<T>::rand(int n1, int n2, int n3) {
    Tensor<T> tem(n1,n2,n3);
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            for(int k=0; k<n3; ++k) {
                tem.p[i][j][k] = randn();
            }
        }
    }
    return tem;
}

template<class T>
Tensor<T> Tensor<T>::Identity(int n1, int n2, int n3) {
    if(n1==n2){
        Tensor tem(n1,n2,n3);
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
                for(int k=0; k<n3; ++k) {
                    if(k==0 && i==j) {
                        tem.p[i][j][k] = 1;
                    }
                    else {
                        tem.p[i][j][k]=0;
                    }
                }
            }
        }
        return tem;
    }
    else{
        cout<<"Warning: no identity tensor";
        exit(0);
    }
}

template<class T>
int *getsize(const Tensor<T> &a) {
    static int p[3]; //为什么加了static就行了。。。
    p[0]=a.n1;
    p[1]=a.n2;
    p[2]=a.n3;
    return p;
}

template<class T>
vec fiber(const Tensor<T> &t, int m, int n, int order) {
    if(order==1){
        vec c(t.n1);
        for (int i=0; i<t.n1;++i){
            c(i)=t.p[i][m][n];
        }
        return c;
    }
    if(order==2){
        vec c(t.n2);
        for (int i=0; i<t.n1;++i){
            c(i)=t.p[i][m][n];
        }
        return c;
    }
    if(order==3){
        vec c(t.n3);
        for (int i=0; i<t.n1;++i){
            c(i)=t.p[i][m][n];
        }
        return c;
    }
}

template<class T>
fmat slice(const Tensor<T> &t, int m, int order) {
    if(order==1){
        fmat c(t.n2,t.n3);
        for (int i=0;i<t.n2;++i){
            for (int j=0;j<t.n3;++j){
                c(i,j)=t.p[m][i][j];
            }
        }
        return c;
    }
    if(order==2){
        fmat c(t.n1,t.n3);
        for (int i=0;i<t.n1;++i){
            for (int j=0;j<t.n3;++j){
                c(i,j)=t.p[i][m][j];
            }
        }
        return c;
    }
    if(order==3){
        fmat c(t.n1,t.n2);
        for (int i=0;i<t.n1;++i){
            for (int j=0;j<t.n2;++j){
                c(i,j)=t.p[i][j][m];
            }
        }
        return c;
    }
}

template<class T>
vec ten2vec(Tensor<T> &a) {
    int len = a.n1*a.n2*a.n3;
    vec result(len);
    for (int i = 0;  i < a.n1;i++){
        for (int j = 0; j < a.n2; j++){
            for (int k = 0; k < a.n3; k++ ){
                int I = a.n1;
                int J = a.n1*a.n2;
                result(i+j*I+k*J) = a.p[i][j][k];
            }
        }
    }
    return result;
}

template<class T>
mat ten2mat(Tensor<T> &a, int b) {
    int m,len;
    if (b==1){ m=a.n1;len = a.n2*a.n3;
        mat result = zeros(m,len);
        for (int k = 0;  k < a.n3;k++){
            for (int i = 0; i < a.n1; i++){
                for (int j = 0; j < a.n2; j++ ){
                    int J = a.n2;
                    result(i,j+k*J) = a.p[i][j][k];
                }
            }
        }
        return result;
    }
    else if (b == 2){ m=a.n2;len = a.n1*a.n3;
        mat result = zeros(m,len);
        for (int i = 0; i < a.n1; i++){
            for (int k = 0;  k < a.n3;k++){
                for (int j = 0; j < a.n2; j++ ){
                    int J = a.n1;
                    result(j,i+k*J) = a.p[i][j][k];
                }
            }
        }
        return result;
    }
    else { m=a.n3;len = a.n1*a.n2;
        mat result = zeros(m,len);
        for (int j = 0; j < a.n2; j++ ){
            for (int k = 0;  k < a.n3;k++){
                for (int i = 0; i < a.n1; i++){
                    int J = a.n1;
                    result(k,i+j*J) = a.p[i][j][k];
                }
            }
        }
        return result;
    }
}

template<class T>
double norm(Tensor<T> &a) {
    double norm=0;
    for (int i = 0; i < a.n1; ++i) {
        for (int j = 0; j < a.n2; ++j) {
            for(int k=0; k< a.n3; ++k) {
                norm=norm+a.p[i][j][k]*a.p[i][j][k];
            }
        }
    }
    return norm;
}

template<class T>
Tensor<T> Transpose(Tensor<T> &a) {
    Tensor<T> tem=a;
    for (int i = 0; i < a.n1; ++i) {
        for (int j = 0; j < a.n2; ++j) {
            for(int k=0; k< a.n3; ++k) {
                if (k==0) tem.p[i][j][k]=a.p[j][i][k];
                else tem.p[i][j][k]=a.p[j][i][a.n3-k];
            }
        }
    }
    return tem;
}

template<class T1, class T2>
double dotProduct(Tensor<T1> a, Tensor<T2> b) {
    double sum = 0;
    if(a.n1==b.n1 && a.n2==b.n2 && a.n3==b.n3) {
        for (int i = 0; i < a.n1; ++i) {
            for (int j = 0; j < a.n2; ++j) {
                for (int k = 0; k < a.n3; ++k) {
                    sum = sum + a.p[i][j][k] * b.p[i][j][k];
                }
            }
        }
    }
    else{
        cout<<"Warning: The size is not match."<<endl;
        exit(0);
    }
    return sum;
}

template<class T1, class T2>
Tensor<T1> tprod(Tensor<T1> &a, Tensor<T2> &b) {
    Tensor<T1> tmp(a.n1,b.n2,a.n3);
    tmp=tmp.zeros(a.n1,b.n2,a.n3);
    if (a.n2==b.n1 && a.n3==b.n3){
        for(int k=0;k<a.n3;++k){
            for (int i=0;i<a.n1;i++){
                for (int j=0;j<b.n2;++j){
                    double s=0;
                    for (int l=0;l<a.n2;++l){
                        s=s+a.p[i][l][k]*b.p[l][j][k];
                    }
                    tmp.p[i][j][k]=s;
                }
            }
        }
        return tmp;
    }

}

template<class T>
tuckercore <T> HOSVD(Tensor<T> &a, int r1, int r2, int r3) {
    int n1=a.n1; int n2=a.n2; int n3=a.n3;
    fmat tmp(n3,n3);
    fmat cal(n1,n3);
    tmp.zeros();

    for (int i=0; i< n2;i++){
        cal = slice(a,i,2);
        tmp = tmp + cal.t()*cal; //slice computing
    }

    fmat U;//U,V均为正交矩阵
    fvec S;//S为奇异值构成的列向量
    fmat U3(n3,r3);
    eig_sym(S,U,tmp);
    U3 = U.cols(n3-r3,n3-1);

    tmp.resize(n1,n1);
    cal.resize(n1,n2);
    tmp.zeros();
//    m1.reset(); //直接改变矩阵大小 不需要删除堆栈？？内存管理好神奇。。。

    for (int i=0; i< n3;i++){
        cal = slice(a,i,3);
        tmp = tmp + cal*cal.t();
    }

    eig_sym(S,U,tmp);
    fmat U1(n1,r1);
    U1 = U.cols(n1-r1,n1-1);

    tmp.resize(n2,n2);
    tmp.zeros();
    for (int i=0; i< n3;i++){
        cal = slice(a,i,3);
        tmp = tmp + cal.t()*cal;
    }

    eig_sym(S,U,tmp);
    fmat U2(n2,r2);
    U2 = U.cols(n2-r2,n2-1);

    U.reset();S.reset();
//    S.reset();U.reset();
    Tensor<T>  g_tmp(r1,r2,n3);
    tmp.resize(r1,r2);
    for (int i=0; i< n3;i++){
        cal = slice(a,i,3);
        tmp = U1.t()*cal*U2;
        for(int j=0;j<r1;j++){
            for(int k=0;k<r2;k++){
                g_tmp(j,k,i) = tmp(j,k);
            }
        }
    }
//    tmp.reset();

    Tensor<T> g(r1,r2,n3);
    cal.resize(r1,r2);
    tmp.resize(r1,r3);
    for (int i=0; i< r2;i++){
        cal = slice(g_tmp,i,2);
        tmp = cal.t()*U3.t();
        for(int j=0;j<r1;j++){
            for(int k=0;k<r3;k++){
                g(j,k,i) = tmp(j,k);
            }
        }
    }
//    cout << U1 << endl;
//    cout << U2 << endl;
//    cout << U3 << endl;

    tuckercore<T> A{g,U1,U2,U3};
    return A;
}

template <class T>
fmat cp_als(Tensor<T> &a, int r){
    int n1=a.n1; int n2=a.n2; int n3=a.n3;
    fmat A=randu<fmat>(n1,r); fmat B=randu<fmat>(n2,r); fmat C=randu<fmat>(n3,r); //randu(fmat)

    fmat cal(n1,n2); fmat tmp(n1,r); fmat krtmp(n2,r);

    for (int turn=0; turn<1;turn++) {
        tmp.zeros();
        for (int j = 0; j < n3; j++) {
            for (int i = 0; i < r; i++) {
                krtmp.col(i) = C(j, i) * B.col(i);
            }
            cal = slice(a, j, 3);
            tmp = tmp + cal * krtmp; //slice computing
        }
        A = tmp * inv((C.t() * C) % (B.t() * B));
        A = normalise(A);

        tmp.set_size(n2, r);
        tmp.zeros();
        krtmp.set_size(n1, r);
        for (int j = 0; j < n3; j++) {
            for (int i = 0; i < r; i++) {
                krtmp.col(i) = C(j, i) * A.col(i);
            }
            cal = slice(a, j, 3);
            tmp = tmp + cal.t() * krtmp; //slice computing
        }
        B = tmp * inv((C.t() * C) % (A.t() * A));
        B = normalise(B);

        tmp.set_size(n3, r);
        tmp.zeros();
        krtmp.set_size(n1, r);
        for (int j = 0; j < n2; j++) {
            for (int i = 0; i < r; i++) {
                krtmp.col(i) = B(j, i) * A.col(i);
            }
            cal = slice(a, j, 2);
            tmp = tmp + cal.t() * krtmp; //slice computing
        }
        C = tmp * inv((B.t() * B) % (A.t() * A));
        if (turn != 1){
            C = normalise(C);
        }
    }

//    Tensor<T> t_con(n1,n2,n3);
//
//    for (int i=0;i<n1;i++){
//        for (int j=0;j<n2;j++){
//            for(int k=0;k<n3;k++){
//                T sum = 0;
//                for (int l=0;l<r;l++){
//                    sum = sum + A(i,l)*B(j,l)*C(k,l);
//                }
//                t_con(i,j,k) = sum;
//            }
//        }
//    }
//    mat m1=ten2mat(t_con,1);
//    mat m2=ten2mat(a,1);
//    cout << m1-m2 << endl;

//    cout << A << endl;
//    cout << B << endl;
//    cout << C << endl;

    return B;
}


//double *** al(int &n1,int &n2,int &n3){
//    double ***p;
//    p=new double**[n1];
//    for (int i=0;i<n1;i++){
//        p[i] = new double *[n2];
//        for (int j=0;j<n2;j++) {
//            p[i][j] = new double[n3];
//        }
//    }
//    return p;
//}

////template <typename T>
//void del(double *** p, int &n1,int &n2){
//    for(int i=0; i<n1; i++) {
//        for (int j = 0; j < n2; j++) {
//            delete [] p[i][j];
//        }
//        delete [] p[i];
//    }
//    delete  [] p;
//    p=NULL;
//}

//double *** loadfile(int &n1, int &n2, int &n3, char *path){
//
//    int count = 0;
//    FILE* fp;
//    char str[100];
//
//    double ***v = al(n1,n2,n3);
//
//    fp = fopen(path,"r");
//
//    string tmp;
//    while (fscanf(fp, "%s", str) != EOF)
//    {
//        int NUM=count;
//        int k=NUM%(n1*n2);
//        k=(NUM-k)/(n1*n2);
//        int j=(NUM-k*n1*n2)%n2;
//        int i=(NUM-k*n1*n2)/n2;
//        tmp=str;
//        v[i][j][k]=(double)atof(tmp.c_str());
//        count++;
//        if(count==n1*n2*n3){break;}
//    }
//
//    fclose(fp);
//
////释放空间 内存泄漏。。。还会影响运行时间
//    return v;
//}

//void display(double *** v, int &n1, int &n2, int &n3){
//
//    cout << "Tensor: " << endl;
//
//    for (int k = 0; k < n3; k++) {
//        for (int i = 0; i < n1; i++) {
//            for (int j = 0; j < n2; j++) {
//                cout << v[i][j][k] << " ";
//            }
//            cout << endl;
//        }
//        cout << endl;
//    }
//}

//template <class T>
//int tsvd(Tensor<T> &a){
//    int n1=a.n1; int n2=a.n2; int n3=a.n3;int N0 = floor(n3/2.0) + 1;
//
//    Tensor<T> v_t(n1,n2,N0);
//    Tensor<T> v_t1(n1,n2,N0);
//
//    fftw_complex out[N0];
//    T *in = fftw_alloc_real(n3);
//
//    fftw_plan p_fft;
//    p_fft = fftw_plan_dft_r2c_1d(n3, in, out, FFTW_ESTIMATE);
////    p=fftw_plan_dft_1d(n3,in,out,FFTW_FORWARD,FFTW_MEASURE);
//
//    for (int i = 0; i < n1; i++) {
//        for (int j = 0; j < n2; j++) {
//            in = fiber(a,i,j,3);   //第三维赋值
//            fftw_execute_dft_r2c(p_fft, in, out);
//            for (int k = 0; k < N0; k++) {
//                v_t[i][j][k] = out[k][0];
//                v_t1[i][j][k] = out[k][1];
//            }               //局部变量归零了？？
//        }
//    }
////    del(M,n1,n2);
//
//    Tensor<T> uf(n1,n1,N0);
//    Tensor<T> uf1(n1,n1,N0);
//    Tensor<T> theta(n1,n2,N0);
//    Tensor<T> vf(n2,n2,N0);
//    Tensor<T> vf1(n2,n2,N0);
//
//    cx_mat TMP = zeros<cx_mat>(n1, n2);
//    cx_mat TMPU = zeros<cx_mat>(n1, n1);
//    cx_mat TMPV = zeros<cx_mat>(n2, n2);
//    colvec TMPT;
//
////    t2=gettime();
////    cout<<"alloc space: "<<t2-t1<<endl;
//
////#pragma omp parallel for num_threads(8) //如何写多线程的svd 会不会覆盖同一个svd，的结果？？
//    for (int k = 0; k < N0; k++) {
//        for (int i = 0; i < n1; i++) {
//            for (int j = 0; j < n2; j++) {
//                TMP(i, j).real(v_t[i][j][k]);
//                TMP(i, j).imag(v_t1[i][j][k]);
//            }
//        }
//
//        svd(TMPU, TMPT, TMPV, TMP, "dc");
////        svd(TMPU,TMPT,TMPV,TMP,"std");
//
//        for (int i = 0; i < n1; i++) {
//            for (int j = 0; j < n1; j++) {
//                uf[i][j][k] = TMPU(i, j).real();
//                uf1[i][j][k] = TMPU(i, j).imag();
//            }
//        }
//        for (int i = 0; i < n2; i++) {
//            for (int j = 0; j < n2; j++) {
//                vf[i][j][k] = TMPV(i, j).real();
//                vf1[i][j][k] = TMPV(i, j).imag();
//            }
//        }
//        if (n1 <= n2) {
//            for (int i = 0; i < n1; i++) {
//                theta[i][i][k] = TMPT(i);
//            }
//        } else {
//            for (int i = 0; i < n2; i++) {
//                theta[i][i][k] = TMPT(i);
//            }
//        }
//    }
//
//    fftw_complex out1[N0]; //fftw_alloc_real()
//    double *in1 = fftw_alloc_real(n3);
//    p_fft = fftw_plan_dft_c2r_1d(n3, out1, in1, FFTW_ESTIMATE);
//
//    Tensor<T> U(n1,n1,n3);
//
//    for (int i = 0; i < n1; i++) {
//        for (int j = 0; j < n1; j++) {
//            for (int k = 0; k < N0; k++) {
//                out1[k][0] = uf[i][j][k];
//                out1[k][1] = uf1[i][j][k];
//            }
//            fftw_execute_dft_c2r(p_fft, out1, in1); //为什么执行不了？？和申请Theta有什么关系？？
//
//            for (int k = 0; k < n3; k++) {
//                U[i][j][k] = 1.0 / n3 * in1[k];
//            }
//        }
//    }
//    Tensor<T> V(n2,n2,n3);
//    for (int i = 0; i < n2; i++) {
//        for (int j = 0; j < n2; j++) {
//            for (int k = 0; k < N0; k++) {
//                out1[k][0] = vf[i][j][k];
//                out1[k][1] = vf1[i][j][k];
//            }
//            fftw_execute_dft_c2r(p_fft, out1, in1);
//            for (int k = 0; k < n3; k++) {
//                V[i][j][k] = 1.0 / n3 * in1[k];
//            }
//        }
//    }
//
//    Tensor<T> Theta(n1,n2,3);
//    for (int i = 0; i < n1; i++) {
//        for (int j = 0; j < n2; j++) {
//            for (int k = 0; k < N0; k++) {
//                out1[k][0] = theta[i][j][k];
//                out1[k][1] = 0;
//            }
//            fftw_execute_dft_c2r(p_fft, out1, in1);
//            for (int k = 0; k < n3; k++) {
//                Theta[i][j][k] = 1.0 / n3 * in1[k];
//            }
//        }
//    }
//    del(theta, n1, n2);
//
//    fftw_destroy_plan(p_fft);
//
////    t4 = gettime();
////    cout<<"ifft time: "<<t4-t3<<endl;
////    cout << "Total time: " << t4 - t0 << endl;
//    return 0;
//}

template<class T>
mat ttm(Tensor<T> &a, mat &b, int c) {
    mat result;
    mat a_c = ten2mat(a,c);
    result = b*a_c;
    return result;
}



#endif //TENSOR_TENSOR_H
