/*!
* \file matrix.h
* \brief 概述 
* 
*详细概述 
* 
* \author 作者名字
* \version 版本号(maj.min，主版本.分版本格式) 
* \date 日期 
*/


#ifndef MATRIX_H
#define MATRIX_H


#include <vector>
#include <complex>

using namespace std;

/// \brief 命名空间的简单概述 
/// 
///命名空间的详细概述
namespace matrixlab
{

    template <typename Type>
	/// \brief Ctext的doxygen测试
	///
	/// 作doxygen测试用
    class Matrix
    {

    public:
        Matrix();
        Matrix( const Matrix<Type> &A );
        Matrix( int rows, int columns, const Type &x = Type(0) );
        Matrix( int rows, int columns, const Type *v );
        ~Matrix();

        Matrix<Type>& operator=( const Matrix<Type> &A );
        Matrix<Type>& operator=( const Type &x );

        Type* operator[]( int i );
        const Type* operator[]( int i ) const;
        Type& operator()( int row, int column );
        const Type& operator()( int  row, int column ) const;

        operator Type*();
        operator const Type*() const;

        long size() const;
        int dim( int dimension ) const;
        int rows() const;
        int cols() const;
        Matrix<Type>& resize( int rows, int columns );
        vector<Type> getRow( int row ) const;
        vector<Type> getColumn( int column ) const;
		void setRow( const vector<Type> &v, int row );
        void setColumn( const vector<Type> &v, int column );

        Matrix<Type>& operator+=( const Type& );
        Matrix<Type>& operator+=( const Matrix<Type>& );
        Matrix<Type>& operator-=( const Type& );
        Matrix<Type>& operator-=( const Matrix<Type>& );
        Matrix<Type>& operator*=( const Type& );
        Matrix<Type>& operator*=( const Matrix<Type>& );
        Matrix<Type>& operator/=( const Type& );
        Matrix<Type>& operator/=( const Matrix<Type>& );

    private:

        Type *pv0;///< 矩阵'0'基指针
        Type **prow0;///<矩阵'0'基行指针
		Type **prow1;///<矩阵'1'基行指针

        int	 nRow;///< 矩阵行数
        int	 nColumn;///< 矩阵列数
        long nTotal;///< 矩阵总数

        void init( int rows, int columns );
        void copyFromArray( const Type *v );
        void setByScalar( const Type &x );
        void destroy();

    };

    template<typename Type>
    ostream& operator<<( ostream&, const Matrix<Type>& );
	template<typename Type>
	ostream& operator<<( ostream&, const vector<Type>& );
    template<typename Type>
    istream& operator>>( istream&, Matrix<Type>& );

    template<typename Type>
    Matrix<Type> operator-( const Matrix<Type>& );
    template<typename Type>
    Matrix<Type> operator+( const Matrix<Type>&, const Type& );
    template<typename Type>
    Matrix<Type> operator+( const Type&, const Matrix<Type>& );
   template<typename Type>
    Matrix<Type> operator+( const Matrix<Type>&, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type> operator-( const Matrix<Type>&, const Type& );
    template<typename Type>
    Matrix<Type> operator-( const Type&, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type> operator-( const Matrix<Type>&, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type> operator*( const Matrix<Type>&, const Type& );
    template<typename Type>
    Matrix<Type> operator*( const Type&, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type> operator*( const Matrix<Type>&, const Matrix<Type>& );
    template<typename Type>
    vector<Type> operator*( const Matrix<Type>&, const vector<Type>& );
    template<typename Type>
    Matrix<Type> operator/( const Matrix<Type>&, const Type& );
	template<typename Type>
	vector<Type> operator/( const vector<Type>&, const Type& );
    template<typename Type>
    Matrix<Type> operator/( const Type&, const Matrix<Type>& );

    template<typename Type>
    Matrix<Type>& optMult( const Matrix<Type>&, const Matrix<Type>&, Matrix<Type>& );
    template<typename Type>
    vector<Type>& optMult( const Matrix<Type>&, const vector<Type>&, vector<Type>& );
	
    template<typename Type>
    Matrix<Type> elemMult( const Matrix<Type>&, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type> elemDivd( const Matrix<Type>&, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type>& elemMultEq( Matrix<Type>&, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type>& elemDivdEq( Matrix<Type>&, const Matrix<Type>& );

    template<typename Type> Matrix<Type> trT( const Matrix<Type>& );
    template<typename Type> Matrix<Type> trH( const Matrix<Type>& );

    template<typename Type>
    Matrix<Type> trMult( const Matrix<Type>&, const Matrix<Type>& );
    template<typename Type>
    vector<Type> trMult( const Matrix<Type>&, const vector<Type>& );
    template<typename Type>
    Matrix<Type> multTr( const Matrix<Type>&, const Matrix<Type>& );
 /* template<typename Type>
    Matrix<Type> multTr( const vector<Type>&, const vector<Type>& );*/
    template<typename Type>
    Matrix<complex<Type> > trMult( const Matrix<complex<Type> >&,
                                   const Matrix<complex<Type> >& );
  /*  template<typename Type>
    vector<complex<Type> > trMult( const Matrix<complex<Type> >&,
                                   const vector<complex<Type> >& );*/
    template<typename Type>
    Matrix<complex<Type> > multTr( const Matrix<complex<Type> >&,
                                   const Matrix<complex<Type> >& );
 /*   template<typename Type>
    Matrix<complex<Type> > multTr( const vector<complex<Type> >&,
                                   const vector<complex<Type> >& );

    template<typename Type> Matrix<Type> eye( int, const Type& );*/
    template<typename Type> vector<Type> diag( const Matrix<Type>& );
     template<typename Type> Matrix<Type> diag( const vector<Type>& );

   /*template<typename Type> Type norm( const Matrix<Type>& );
    template<typename Type> Type norm( const Matrix<complex<Type> >& );
    template<typename Type> void swap( Matrix<Type>&, Matrix<Type>& );*/
    template<typename Type> vector<Type> sum( const Matrix<Type>& );
    template<typename Type> vector<Type> min( const Matrix<Type>& );
    template<typename Type> vector<Type> max( const Matrix<Type>& );
    template<typename Type> vector<Type> mean( const Matrix<Type>& );
    template<typename Type> Matrix<Type> abs( const Matrix<complex<Type> >& );
    template<typename Type> Matrix<Type> arg( const Matrix<complex<Type> >& );
    template<typename Type> Matrix<Type> real( const Matrix<complex<Type> >& );
    template<typename Type> Matrix<Type> imag( const Matrix<complex<Type> >& );
    template<typename Type>
    Matrix<complex<Type> > complexMatrix( const Matrix<Type>& );
/*  template<typename Type>
    Matrix<complex<Type> > complexMatrix( const Matrix<Type>&,
                                          const Matrix<Type>& );*/


    #include "matrix_impl.h"

}

#endif
// MATRIX_H