#ifndef TENSOR_TRAIN
#define TENSOR_TRAIN

#include "tensor.h"

namespace yph
{
template <class T>
class TensorTrain
{
  private:
    Cube<T> tensor;
    Cube<T> *cores;
    int modeNum;
    int *modeSize;
    int *ttRank;
    double delta;
    void getTTform();
    int numel()
    {
        int temp = 1;
        for (int i = 0; i < modeNum; ++i)
            temp *= modeSize[i];
        return temp;
    }
    inline void rounding();
    inline void copy();
    inline Mat<T> ten2mat(Cube<T> &a, int row, int col);

  public:
    TensorTrain() {}
    TensorTrain(const Cube<T> &t, const double &epsilon);
    TensorTrain(const TensorTrain<T> &t, bool onlyCopy);
    TensorTrain(const cp_mats<T> &t);
    ~TensorTrain();
    Cube<T> &getTensor() const { return tensor; }
    int getSize()
    {
        // 返回维度数目
        return modeNum;
    }
    int getSize(int mode)
    {
        // 返回选中维度大小
        return modeSize[mode - 1];
    }
    friend TensorTrain<T> operator+(const TensorTrain<T> &t1, const TensorTrain<T> &t2);
};

template <class T>
TensorTrain<T>::TensorTrain(const Cube<T> &t, const double &epsilon) : tensor(t), modeNum(3)
{
    modeSize = new int(modeNum);
    modeSize[0] = t.n_rows;
    modeSize[1] = t.n_cols;
    modeSize[2] = t.n_slices;
    delta = epsilon * norm(t) / sqrt(modeNum - 1); //norm is in tensor.cpp
    cores = new Cube<T>(modeNum);
    ttRank = new int(modeNum + 1);
    getTTform();
}

template <class T>
TensorTrain<T>::TensorTrain(const TensorTrain<T> &t, bool onlyCopy)
{
    if (onlyCopy)
    {
        copy();
    }
    else
    {
        // TT-rounding
        rounding();
    }
}

template <class T>
TensorTrain<T>::~TensorTrain()
{
    delete[] modeSize;
    delete[] cores;
    delete[] ttRank;
}

template <class T>
void TensorTrain<T>::getTTform()
{
    ttRank[0] = 1;
    ttRank[modeNum] = 1;
    Cube<T> tempTensor = tensor;
    Mat<T> tempMat;
    for (int i = 1; i < modeNum; ++i)
    {
        int rowSize = ttRank[i - 1] * modeSize[i - 1];
        tempMat = ten2mat(tempTensor, rowSize, numel(tensor)); // 这里仅适用于三维
        Mat<T> U, S, V;
        int r;
        deltaSVD(tempMat, U, S, V, r, delta);
        ttRank[i] = r;
        cores[i - 1] = reshape(U * S * V.t(), ttRank[i - 1], modeSize[i - 1], ttRank[i]); // need to implement
    }
    cores[modeNum - 1] = tempTensor;
}

template <class T>
inline void TensorTrain<T>::rounding()
{
}
template <class T>
inline void TensorTrain<T>::copy()
{
}
template <class T>
inline Mat<T> TensorTrain<T>::ten2mat(Cube<T> &t, int row, int col)
{
}
template <class T>
inline Mat<T> TensorTrain<T>::reshape(Mat<T> &m, Cube<T> &t, int row, int col, int slice)
{
}
// -------------------------------------------------------------
template <class T>
void deltaSVD(Mat<T> mat, Mat<T> &U, Mat<T> &S, Mat<T> &V, int &r, const double &delta)
{
    Col<T> s = zeros<Col<T>>(min(mat.n_rows, mat.n_cols));
    svd(U, s, V, mat); // use V.t() to get trans
    truncation(s, r, delta);
    S = diagmat(s);
}

template <class T>
void truncation(Col<T> &s, int &r, const double &delta)
{
    T temp = 0;
    for (uword i = s.n_rows - 1; i >= 0; --i)
    {
        temp += s(i); // 利用的是norm(A) == norm(S)
        if (temp > delta)
        {
            for (uword j = s.n_rows - 1; j > i; --j)
            {
                s(j) = 0;
            }
            r = i + 1;
            break;
        }
    }
}
template <class T>
TensorTrain<T> operator+(const TensorTrain<T> &t1, const TensorTrain<T> &t2)
{
    if (t1.modeNum != t2.modeNum)
        return new TensorTrain<T>();
    TensorTrain<T> answer(t1, true);
    answer.cores[0] = new Cube<T>(1, t1.modeSize[t1.modeNum], t1.cores[0].n3 + t2.cores[0].n3);
    answer.cores[t1.modeNum] = new Cube<T>(t1.cores[0].n1 + t2.cores[0].n1, t1.modeSize[t1.modeNum], 1);
    for (int i = 1; i < t1.modeNum; ++i)
    {
        answer.cores[i]; // need 2 implement
    }
    return answer;
}

template <class T>
void dotProduct()
{
}

template <class T>
void multiContraction()
{
}

} // namespace yph

#endif
