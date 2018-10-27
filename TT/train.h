#ifndef TENSOR_TRAIN
#define TENSOR_TRAIN

#include "tensor.h"

namespace yph{
    template<class T>
    class tensorTrain{
    private:
        Tensor<T> tensor;
        Tensor<T> *cores;
        int modeNum;
        int *modeSize;
        int *ttRank;
        double delta;
        void getTTform();
        int numel(){ 
            int temp = 1;
            for(int i = 0;i<modeNum;++i) temp*=modeSize[i];
            return temp;
        }
    public:
        tensorTrain(){}
        tensorTrain(Tensor<T> t,double epsilon);
        ~tensorTrain();
        Tensor<T> getTensor() const { return tensor; }
        int getSize(){
            // 返回维度数目
            return modeNum;
        }
        int getSize(int mode){
            // 返回选中维度大小
            return modeSize[mode-1];
        }
    };

    template<class T>
    tensorTrain<T>::tensorTrain(Tensor<T> t,double epsilon):tensor(t),modeNum(3){
        modeSize = new int(modeNum);
        modeSize[0] = t.n1;
        modeSize[1] = t.n2;
        modeSize[2] = t.n3;
        delta = epsilon * norm(t) / sqrt(modeNum-1);
        cores = new Tensor<T>(modeNum);
        ttRank = new int(modeNum + 1);
        getTTform();
    }
    template<class T>
    tensorTrain<T>::~tensorTrain(){
        delete []modeSize;
        delete []cores;
        delete []ttRank;
    }

    template<class T>
    void tensorTrain<T>::getTTform(){
        ttRank[0] = 1;
        ttRank[modeNum] = 1;
        Tensor<T> tempTensor = tensor;
        Mat<T> tempMat;
        for(int i = 1;i<modeNum;++i){
            int rowSize = ttRank[i-1]*modeSize[i-1];
            tempMat = ten2mat(tempTensor,rowSize,numel(tensor)); // need to implement
            Mat<T> U,S,V;
            int r;
            deltaSVD(tempMat,U,S,V,r,delta);
            ttRank[i] = r;
            cores[i-1] = reshape(U*S*V,ttRank[i-1],modeSize[i-1],ttRank[i]);// need to implement
        }
        cores[modeNum-1] = tempTensor;
    }

    // -------------------------------------------------------------
    template<class T>
    void deltaSVD(Mat<T> mat, Mat<T> &U, Mat<T> &S, Mat<T> &V, int r,double delta){
       
    }
    
    template<class T>
    void product(){

    }

    template<class T>
    void rounding(){

    }
}

#endif
