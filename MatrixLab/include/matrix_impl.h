/// \brief 矩阵初始化(pv0, prow0, prow1)
/// \param rows 矩阵行数
/// \param columns 矩阵列数
template <typename Type>
void Matrix<Type>::init( int rows, int columns )
{
	nRow = rows;
	nColumn = columns;
	nTotal = nRow * nColumn;

	pv0 = new Type[nTotal];
	prow0 = new Type*[nRow];
	prow1 = new Type*[nRow];
	Type *p = pv0;
	for( int i=0; i<nRow; ++i )
	{
		prow0[i] = p;
		prow1[i] = p-1;
		p += nColumn;
	}

	prow1--;
}

/// \brief 利用数组数据初始化矩阵数据pv0
/// \param v 数据指针
template <typename Type>
inline void Matrix<Type>::copyFromArray( const Type *v )
{
	for( long i=0; i<nTotal; ++i )
		pv0[i] = v[i];
}
/// \brief 设置矩阵数据为同一数x
/// \param x 矩阵系数
template <typename Type>
inline void Matrix<Type>::setByScalar( const Type &x )
{
	for( long i=0; i<nTotal; ++i )
		pv0[i] = x;
}
/// \brief 矩阵销毁(pv0, prow0, prow1)
template <typename Type>
void Matrix<Type>::destroy()
{
	if( pv0 == NULL )
		return ;
	else
		delete []pv0;

	if( prow0 != NULL )
		delete []prow0;

	prow1++;
	if( prow1 != NULL )
		delete []prow1;
}

/// \brief 默认构造函数
template <typename Type>
Matrix<Type>::Matrix()
: pv0(0), prow0(0), prow1(0), nRow(0), nColumn(0), nTotal(0)
{
}
/// \brief 构造函数,利用已知矩阵初始化
/// \param A 已知矩阵
/// \see init(),copyFromArray()
template <typename Type>
Matrix<Type>::Matrix( const Matrix<Type> &A )
{
	init( A.nRow, A.nColumn );
	copyFromArray( A.pv0 );
}

/// \brief 构造函数,已知矩阵行数、列数、系数,初始化矩阵
/// \param rows 矩阵行
/// \param columns 矩阵列
/// \param x 矩阵系数,默认为0
/// \see init(),setByScalar()
template <typename Type>
Matrix<Type>::Matrix( int rows, int columns, const Type &x = Type(0) )
{
	init( rows,columns );
	setByScalar(x);
}
/// \brief 构造函数,已知矩阵行数、列数、数组,初始化矩阵
/// \param rows 矩阵行
/// \param columns 矩阵列
/// \param array 矩阵数据指针
/// \see init(),copyFromArray()
template <typename Type>
Matrix<Type>::Matrix( int rows, int columns, const Type *arrays )
{
	init( rows,columns );
	copyFromArray( arrays );
}
/// \brief 类析构函数
template <typename Type>
Matrix<Type>::~Matrix()
{
	destroy();
}
/// \brief 重载"=",矩阵赋值,such as: matrixB = matrixA
/// \param A 赋值右参数
/// \return 赋值左参数
/// \remarks 如果maxtrixB 与 maxtrix A的维数不等,则将原先matrixB销毁,赋值为matrixA
template <typename Type>
Matrix<Type>& Matrix<Type>::operator=( const Matrix<Type> &A )
{
	if( pv0 == A.pv0 )
		return *this;

	if( nRow == A.nRow && nColumn == A.nColumn )
		copyFromArray( A.pv0 );
	else
	{
		destroy();
		init( A.nRow, A.nColumn );
		copyFromArray( A.pv0 );
	}

	return *this;
}
/// \brief 重载"=",系数赋值,such as: matrixA = matrix(x)
/// \param x 矩阵系数
/// \return 赋值左参数
template <typename Type>
inline Matrix<Type>& Matrix<Type>::operator=( const Type &x )
{
	setByScalar( x );
	return *this;
}
/// \brief 重载"[]",矩阵'0'基行数据访问. 
///	such as:matrixA =
///							 [1,  1, 1
///		*matrixA[1]=>首地址  (2), 2, 2
///							  3,  3, 3]
/// \param x 行数
/// \return 第i行的首地址
template <typename Type>
inline Type* Matrix<Type>::operator[]( int i )
{
	return prow0[i];
}
template <typename Type>
inline const Type* Matrix<Type>::operator[]( int i ) const
{
	return prow0[i];
}
/// \brief 重载"()",矩阵'1'基数据访问. 
///	such as:matrixA =
///							 [1,  1, 1
///		matrixA[2][1]=>      (2), 2, 2
///							  3,  3, 3]
/// \param row 行数
/// \param column 列数
/// \return 第i行第j列数据
template <typename Type>
inline Type& Matrix<Type>::operator()( int row, int column )
{
	return  prow1[row][column];
}

template <typename Type>
inline const Type& Matrix<Type>::operator()( int row, int column ) const
{
	return  prow1[row][column];
}
/// \brief 重载"*",矩阵数据首地址访问. 
///	such as: matrixA =
///		    *matrixA=>	[1,  1, 1
///						(2), 2, 2
///						 3,  3, 3]
/// \return 矩阵数据首地址pv0
template <typename Type>
inline Matrix<Type>::operator Type*()
{
	return pv0;
}

template <typename Type>
inline Matrix<Type>::operator const Type*() const
{
	return pv0;
}
/// \brief 返回矩阵大小
template <typename Type>
inline long Matrix<Type>::size() const
{
	return nTotal;
}
/// \brief 返回矩阵维数
/// \param dimension '1'代表行数,'2'代表列数
template <typename Type>
int Matrix<Type>::dim( int dimension ) const
{

	if( dimension == 1 )
		return nRow;
	else if( dimension == 2 )
		return nColumn;
	else
		return 0;
}
/// \brief 返回矩阵行数
template <typename Type>
inline int Matrix<Type>::rows() const
{
    return nRow;
}
/// \brief 返回矩阵行数
template <typename Type>
inline int Matrix<Type>::cols() const
{
    return nColumn;
}
/// \brief 重新设定矩阵大小
/// \param rows 行数
/// \param columns 列数
template <typename Type>
Matrix<Type>& Matrix<Type>::resize( int rows, int columns )
{
	if(  rows == nRow && columns == nColumn )
		return *this;

	destroy();
	init( rows, columns );

	return *this;
}
/// \brief 获得矩阵行数据
/// \param row 行数(0基)
/// \return 矩阵行向量
template <typename Type>
vector<Type> Matrix<Type>::getRow( int row ) const
{

	vector<Type> tmp( nColumn );
	for( int j=0; j<nColumn; ++j )
		tmp[j] = prow0[row][j];

	return tmp;
}
/// \brief 获得矩阵列数据
/// \param column 列数(0基)
/// \return 矩阵列向量
template <typename Type>
vector<Type> Matrix<Type>::getColumn( int column ) const
{
	vector<Type> tmp( nRow );
	for( int i=0; i<nRow; ++i )
		tmp[i] = prow0[i][column];

	return tmp;
}
/// \brief 设定矩阵行数据
/// \param row 行数(0基)
/// \param v 行向量
template <typename Type>
void Matrix<Type>::setRow( const vector<Type> &v, int row )
{
	for( int j=0; j<nColumn; ++j )
		prow0[row][j] = v[j];
}
/// \brief 设定矩阵列数据
/// \param column 列数(0基)
/// \param v 列向量
template <typename Type>
void Matrix<Type>::setColumn( const vector<Type> &v, int column )
{
	for( int i=0; i<nRow; ++i )
		prow0[i][column] = v[i];
}
/// \brief 重载'+=', 矩阵各元素加系数
/// such as :	 A = [ 1, 2, 3
///					   2, 3, 4
///					   3, 4, 5]
///		A += 2:	 A =  [3, 4, 5
///					   4, 5, 6
///					   5, 6, 7]
/// \param x 系数
template <typename Type>
Matrix<Type>& Matrix<Type>::operator+=( const Type &x )
{
    Type **rowPtr = prow0;
    Type *colPtr = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtr = *rowPtr++;
        for( int j=0; j<nColumn; ++j )
            *colPtr++ += x;
    }
	return *this;
}
/// \brief 重载'+=', 矩阵A各元素加对应矩阵B各元素系数
/// such as :	 A = [ 1, 2, 3		B = [1,	2, 3
///					   2, 3, 4			 2, 3, 4
///					   3, 4, 5]			 3, 4, 5]
///		A += B:	 A =  [2, 4, 6
///					   4, 6, 8
///					   6, 8, 10]
/// \param rhs 矩阵B
template <typename Type>
Matrix<Type>& Matrix<Type>::operator+=( const Matrix<Type> &rhs )
{

    Type **rowPtrL = prow0;
    Type *colPtrL = 0;
    Type **rowPtrR = rhs.prow0;
    const Type *colPtrR = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtrL = *rowPtrL++;
        colPtrR = *rowPtrR++;
        for( int j=0; j<nColumn; ++j )
            *colPtrL++ += *colPtrR++;
    }

	return *this;
}
/// \brief 重载'-=', 矩阵各元素减系数
/// such as :	 A = [ 1, 2, 3
///					   2, 3, 4
///					   3, 4, 5]
///		A -= 1:	 A =  [0, 1, 2
///					   1, 2, 3
///					   2, 3, 4]
/// \param x 系数
template <typename Type>
Matrix<Type>& Matrix<Type>::operator-=( const Type &x )
{
    Type **rowPtr = prow0;
    Type *colPtr = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtr = *rowPtr++;
        for( int j=0; j<nColumn; ++j )
            *colPtr++ -= x;
    }

	return *this;
}
/// \brief 重载'-=', 矩阵A各元素减对应矩阵B各元素系数
/// such as :	 A = [ 1, 2, 3		B = [1,	2, 3
///					   2, 3, 4			 2, 3, 4
///					   3, 4, 5]			 3, 4, 5]
///		A -= B:	 A =  [0, 0, 0
///					   0, 0, 0
///					   0, 0, 0]
/// \param rhs 矩阵B
template <typename Type>
Matrix<Type>& Matrix<Type>::operator-=( const Matrix<Type> &rhs )
{
    Type **rowPtrL = prow0;
    Type *colPtrL = 0;
    Type **rowPtrR = rhs.prow0;
    const Type *colPtrR = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtrL = *rowPtrL++;
        colPtrR = *rowPtrR++;
        for( int j=0; j<nColumn; ++j )
            *colPtrL++ -= *colPtrR++;
    }

	return *this;
}
/// \brief 重载'*=', 矩阵各元素乘系数
/// such as :	 A = [ 1, 2, 3
///					   2, 3, 4
///					   3, 4, 5]
///		A *= 2:	 A =  [2, 4, 6
///					   4, 6, 8
///					   6, 8, 10]
/// \param x 系数
template <typename Type>
Matrix<Type>& Matrix<Type>::operator*=( const Type &x )
{
    Type **rowPtr = prow0;
    Type *colPtr = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtr = *rowPtr++;
        for( int j=0; j<nColumn; ++j )
            *colPtr++ *= x;
    }

	return *this;
}
/// \brief 重载'*=', 矩阵A各元素乘对应矩阵B各元素系数
/// such as :	 A = [ 1, 2, 3		B = [1,	2, 3
///					   2, 3, 4			 2, 3, 4
///					   3, 4, 5]			 3, 4, 5]
///		A *= B:	 A =  [1, 4, 9
///					   4, 9, 16
///					   9, 16, 25]
/// \param rhs 矩阵B
template <typename Type>
Matrix<Type>& Matrix<Type>::operator*=( const Matrix<Type> &rhs )
{
    Type **rowPtrL = prow0;
    Type *colPtrL = 0;
    Type **rowPtrR = rhs.prow0;
    const Type *colPtrR = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtrL = *rowPtrL++;
        colPtrR = *rowPtrR++;
        for( int j=0; j<nColumn; ++j )
            *colPtrL++ *= *colPtrR++;
    }

	return *this;
}
/// \brief 重载'/=', 矩阵各元素除系数
/// such as :	 A = [ 1, 2, 3
///					   2, 3, 4
///					   3, 4, 5]
///		A /= 2:	 A =  [0.5, 1, 1.5
///					   1, 1.5, 2
///					   1.5, 2, 2.5]
/// \param x 系数
template <typename Type>
Matrix<Type>& Matrix<Type>::operator/=( const Type &x )
{
    Type **rowPtr = prow0;
    Type *colPtr = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtr = *rowPtr++;
        for( int j=0; j<nColumn; ++j )
            *colPtr++ /= x;
    }

	return *this;
}
/// \brief 重载'/=', 矩阵A各元素除对应矩阵B各元素系数
/// such as :	 A = [ 1, 2, 3		B = [1,	2, 3
///					   2, 3, 4			 2, 3, 4
///					   3, 4, 5]			 3, 4, 5]
///		A /= B:	 A =  [1, 1, 1
///					   1, 1, 1
///					   1, 1, 1]
/// \param rhs 矩阵B
template <typename Type>
Matrix<Type>& Matrix<Type>::operator/=( const Matrix<Type> &rhs )
{

    Type **rowPtrL = prow0;
    Type *colPtrL = 0;
    Type **rowPtrR = rhs.prow0;
    const Type *colPtrR = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtrL = *rowPtrL++;
        colPtrR = *rowPtrR++;
        for( int j=0; j<nColumn; ++j )
            *colPtrL++ /= *colPtrR++;
    }

	return *this;
}
/// \brief 重载'<<', 输出矩阵A各元素
/// such as :	size: 3 by 3 
///				1, 2, 3		
///				2, 3, 4		
///				3, 4, 5			 
/// \param A 矩阵A
template <typename Type>
ostream& operator<<( ostream &out, const Matrix<Type> &A )
{
	int rows = A.rows();
	int columns = A.cols();

	out << "size: " << rows << " by " << columns << "\n";
	for( int i=0; i<rows; ++i )
	{
		for( int j=0; j<columns; ++j )
			out << A[i][j] << "\t";
		out << "\n";
	}

	return out;
}
/// \brief 重载'<<', 输出向量v各元素
/// such as :	size: 3 by 1 
///				1	
///				2		
///				3		 
/// \param v 向量v
template<typename Type>
ostream& operator<<( ostream &out, const vector<Type> &v )
{
	int rows = v.size();

	out << "size: " << rows << " by " << 1 << "\n";
	for( int i=0; i<rows; ++i )
	{
			out << v.at(i) << "\n";
	}
	return out;
}
/// \brief 重载'>>', 输入矩阵A各元素
/// such as :	3 3 
///				1 1 1	
///				2 2 2		
///				3 3 3		 
/// \param A 输入矩阵
template <typename Type>
istream& operator>>( istream &in, Matrix<Type> &A )
{
	int rows, columns;
	in >> rows >> columns;

	if( !( rows == A.rows() && columns == A.cols() ) )
		A.resize( rows, columns );

	for( int i=0; i<rows; ++i )
	{
		for( int j=0; j<columns; ++j )
			in >> A[i][j];
		cout<<endl;
	}
	return in;
}
/// \brief 重载'-', 矩阵A各元素取反
/// such as :	 A = [ 1, 2, 3		
///					   2, 3, 4			
///					   3, 4, 5]			
///				 -A = [ -1, -2, -3		
///					    -2, -3, -4			
///					    -3, -4, -5]		
/// \param A 矩阵A
template<typename Type>
Matrix<Type> operator-( const Matrix<Type> &A )
{
	int rows = A.rows();
	int columns = A.cols();

	Matrix<Type> tmp( rows, columns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			tmp[i][j] = -A[i][j];

	return tmp;
}
/// \brief 重载'+', 矩阵A各元素加系数
/// such as :	A = [ 1, 2, 3		
///					   2, 3, 4		
///					   3, 4, 5]		
///		A+2:	      [3, 4, 5
///					   4, 5, 6
///					   5, 6, 7]
/// \param A 矩阵A
/// \param x 系数
template<typename Type>
inline Matrix<Type> operator+( const Matrix<Type> &A, const Type &x )
{
	Matrix<Type> tmp( A );
	return tmp += x;
}
/// \brief 重载'+', 系数加矩阵A各元素
/// such as :	A = [ 1, 2, 3		
///					   2, 3, 4		
///					   3, 4, 5]		
///		2+A:	      [3, 4, 5
///					   4, 5, 6
///					   5, 6, 7]
/// \param A 矩阵A
/// \param x 系数
template<typename Type>
inline Matrix<Type> operator+( const Type &x, const Matrix<Type> &A )
{
	return A + x;
}
/// \brief 重载'+', 矩阵A1各元素加对应矩阵A2各元素系数
/// such as :	 A1 = [ 1, 2, 3		A2 = [1,2, 3
///					    2, 3, 4			  2, 3, 4
///					    3, 4, 5]		  3, 4, 5]
///		A1 + A2:	  [ 2, 4, 6
///					    4, 6, 8
///					    6, 8, 10]
/// \param A1 矩阵A1
/// \param A2 矩阵A2
template<typename Type>
inline Matrix<Type> operator+( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	Matrix<Type> tmp( A1 );
	return tmp += A2;
}
/// \brief 重载'-', 矩阵A各元素减系数
/// such as :	A = [ 1, 2, 3		
///					   2, 3, 4		
///					   3, 4, 5]		
///		A-1:	      [0, 1, 2
///					   1, 2, 3
///					   2, 3, 4]
/// \param A 矩阵A
/// \param x 系数
template<typename Type>
inline Matrix<Type> operator-( const Matrix<Type> &A, const Type &x )
{
	Matrix<Type> tmp( A );
	return tmp -= x;
}
/// \brief 重载'-', 系数减矩阵A各元素
/// such as :	A = [ 1, 2, 3		
///					   2, 3, 4		
///					   3, 4, 5]		
///		1-A:	      [ 0, -1, -2
///					   -1, -2, -3
///					   -2, -3, -4]
/// \param A 矩阵A
/// \param x 系数
template<typename Type>
inline Matrix<Type> operator-( const Type &x, const Matrix<Type> &A )
{
	Matrix<Type> tmp( A );
	return -tmp += x;
}
/// \brief 重载'-', 矩阵A1各元素减对应矩阵A2各元素系数
/// such as :	 A1 = [ 1, 2, 3		A2 = [1,2, 3
///					    2, 3, 4			  2, 3, 4
///					    3, 4, 5]		  3, 4, 5]
///		A1 - A2:	  [ 0, 0, 0
///					    0, 0, 0
///					    0, 0, 0]
/// \param A1 矩阵A1
/// \param A2 矩阵A2
template<typename Type>
inline Matrix<Type> operator-( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	Matrix<Type> tmp( A1 );
	return tmp -= A2;
}
/// \brief 重载'*', 矩阵各元素乘系数
/// such as :	 A =  [1, 2, 3
///					   2, 3, 4
///					   3, 4, 5]
///		A * 2:	      [2, 4, 6
///					   4, 6, 8
///					   6, 8, 10]
/// \param A 矩阵
/// \param x 系数
template <typename Type>
inline Matrix<Type> operator*( const Matrix<Type> &A, const Type &x )
{
	Matrix<Type> tmp( A );
	return tmp *= x;
}
/// \brief 重载'*', 系数乘矩阵各元素
/// such as :	 A =  [1, 2, 3
///					   2, 3, 4
///					   3, 4, 5]
///		2 * A:	      [2, 4, 6
///					   4, 6, 8
///					   6, 8, 10]
/// \param A 矩阵
/// \param x 系数
template <typename Type>
inline Matrix<Type> operator*( const Type &x, const Matrix<Type> &A )
{
	return A * x;
}
/// \brief 重载'*', 矩阵A1各元素乘矩阵A2
/// such as :	 A1 = [1, 2, 3		A2 = [1, 2, 3
///					   2, 3, 4			  2, 3, 4
///					   3, 4, 5]			  3, 4, 5]
///		A1 .* A2:	  [14, 20, 26
///					   20, 26, 38
///					   26, 38, 50]
/// \param A1 矩阵A1
/// \param A1 矩阵A2
template <typename Type>
Matrix<Type> operator*( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	int rows = A1.rows();
	int columns = A2.cols();

	Matrix<Type> tmp( rows, columns );
    optMult( A1, A2, tmp );

	return tmp;
}
/// \brief 重载'*', 矩阵A各元素乘向量v
/// such as :	 A = [1, 2, 3		v = [1, 2, 3]T
///					   2, 3, 4			  
///					   3, 4, 5]			 
///		A1 .* v:	  [14, 20, 26]T
/// \param A 矩阵 
/// \param v 向量 
template <typename Type>
vector<Type> operator*( const Matrix<Type> &A, const vector<Type> &b )
{

	int rows = A.rows();

	vector<Type> tmp(rows);

    optMult( A, b, tmp );

	return tmp;
}
/// \brief 重载'/', 矩阵各元素除系数
/// such as :	 A =  [1, 2, 3
///					   2, 3, 4
///					   3, 4, 5]
///		A / 2:	      [0.5, 1, 1.5
///					   1, 1.5, 2
///					   1.5, 2, 2.5]
/// \param A 矩阵
/// \param x 系数
template <typename Type>
inline Matrix<Type> operator/( const Matrix<Type> &A, const Type &x )
{
	Matrix<Type> tmp( A );
	return tmp /= x;
}
/// \brief 重载'/', 向量各元素除系数
/// such as :	 v =  [1,2,3]
///		v / 2:	      [0.5, 1, 1.5]
/// \param v 向量
/// \param x 系数
template<typename Type>
vector<Type> operator/( const vector<Type> &v, const Type &x )
{
	int msize = v.size();
	vector<Type> tmp(msize);
	for (int i = 0; i < msize; i++)
	{
		tmp.at(i) = v.at(i)/x;
	}
	return tmp;
}
/// \brief 重载'/', 系数除矩阵各元素
/// such as :	 A =  [1, 2, 3
///					   2, 3, 4
///					   3, 4, 5]
///		15 / A:	      [15, 7.5, 5
///					   7.5, 5, 3.75
///					   5, 3.75, 3]
/// \param A 矩阵
/// \param x 系数
template <typename Type>
Matrix<Type> operator/( const Type &x, const Matrix<Type> &A )
{
	int rows = A.rows();
	int clumns = A.cols();

	Matrix<Type> tmp( rows,clumns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<clumns; ++j )
			tmp[i][j] = x / A[i][j];

	return tmp;
}
/// \brief 重载'*', 矩阵A各元素乘矩阵B
/// such as :	 A  = [1, 2, 3		B  = [1, 2, 3
///					   2, 3, 4			  2, 3, 4
///					   3, 4, 5]			  3, 4, 5]
///		C = A.* B:	  [14, 20, 26
///					   20, 26, 38
///					   26, 38, 50]
/// \param A 矩阵
/// \param B 矩阵
/// \param C 结果矩阵
template <typename Type>
Matrix<Type>& optMult( const Matrix<Type> &A, const Matrix<Type> &B,
                    Matrix<Type> &C )
{
    int M = A.rows();
    int N = B.cols();
    int K = A.cols();

    C.resize( M, N );
    Type        sum;
    const Type  *pRow,
                *pCol;

    for( int i=0; i<M; i++ )
        for( int j=0; j<N; ++j )
        {
            pRow  = &A[i][0];
            pCol  = &B[0][j];
            sum = 0;

            for( int k=0; k<K; ++k )
            {
                sum += (*pRow) * (*pCol);
                pRow++;
                pCol += N;
            }
            C[i][j] = sum;
        }
    return C;
}
/// \brief 重载'*', 矩阵A各元素乘向量b
/// such as :	 A = [1, 2, 3		b = [1, 2, 3]T
///					   2, 3, 4			  
///					   3, 4, 5]			 
///		A .* b:	  [14, 20, 26]T
/// \param A 矩阵 
/// \param b 向量 
/// \param c 
template <typename Type>
vector<Type>& optMult( const Matrix<Type> &A, const vector<Type> &b,
                    vector<Type> &c )
{
    int M = A.rows();
    int N = A.cols();

    c.resize( M );
    Type        sum;
    const Type  *pRow,
                *pCol;

    for( int i=0; i<M; i++ )
    {
        pRow  = &A[i][0];
        pCol  = &b[0];
        sum = 0;

        for( int j=0; j<N; ++j )
        {
            sum += (*pRow) * (*pCol);
            pRow++;
            pCol++;
        }
        c[i] = sum;
    }
    return c;
}


template<typename Type>
inline Matrix<Type> elemMult( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	Matrix<Type> tmp( A1 );
	return tmp *= A2;
}

template <typename Type>
inline Matrix<Type>& elemMultEq( Matrix<Type> &A1, const Matrix<Type> &A2 )
{
    return A1 *= A2;
}


template <typename Type>
inline Matrix<Type> elemDivd( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	Matrix<Type> tmp( A1 );
	return tmp /= A2;
}
template <typename Type>
inline Matrix<Type>& elemDivdEq( Matrix<Type> &A1, const Matrix<Type> &A2 )
{
    return A1 /= A2;
}


template <typename Type>
Matrix<Type> trT( const Matrix<Type> &A )
{
	int rows = A.cols();
	int clumns = A.rows();

	Matrix<Type> tmp( rows, clumns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<clumns; ++j )
			tmp[i][j] = A[j][i];

	return tmp;
}


template <typename Type>
Matrix<Type> trH( const Matrix<Type> &A )
{
	int rows = A.cols();
	int clumns = A.rows();

	Matrix<Type> tmp( rows, clumns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<clumns; ++j )
			tmp[i][j] = conj(A[j][i]);

	return tmp;
}


template <typename Type>
Matrix<Type> trMult( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	int rows = A1.cols();
	int columns = A2.cols();
	int K = A1.rows();

	Matrix<Type> tmp( rows, columns );

    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			for( int k=0; k<K; ++k )
			   tmp[i][j] += A1[k][i] * A2[k][j];

	return tmp;
}

template <typename Type>
vector<Type> trMult( const Matrix<Type> &A, const vector<Type> &v )
{
	int rows = A.rows();
	int columns = A.cols();

	vector<Type> tmp( columns );
    for( int i=0; i<columns; ++i )
		for( int j=0; j<rows; ++j )
			tmp[i] += A[j][i] * v[j];

	return tmp;
}

template <typename Type>
Matrix<Type> multTr( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{

	int rows = A1.rows();
	int columns = A2.rows();
	int K = A1.cols();

	Matrix<Type> tmp( rows, columns );
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			for( int k=0; k<K; ++k )
			   tmp[i][j] += A1[i][k] * A2[j][k];

	return tmp;
}

/*
template <typename Type>
Matrix<Type> multTr( const vector<Type> &a, const vector<Type> &b )
{
	int rows = a.dim();
	int columns = b.dim();

	Matrix<Type> tmp( rows, columns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			tmp[i][j] = a[i]*b[j];

	return tmp;
}

*/
template <typename Type>
Matrix<complex<Type> > trMult( const Matrix<complex<Type> > &A1,
                               const Matrix<complex<Type> > &A2 )
{

	int rows = A1.cols();
	int columns = A2.cols();
	int K = A1.rows();

	Matrix<complex<Type> > tmp( rows, columns );
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			for( int k=0; k<K; ++k )
			   tmp[i][j] += conj(A1[k][i]) * A2[k][j];

	return tmp;
}
/*
template <typename Type>
vector<complex<Type> > trMult( const Matrix<complex<Type> > &A, const vector<complex<Type> > &v )
{
	int rows = A.rows();
	int columns = A.cols();

	vector<complex<Type> > tmp( columns );

    for( int i=0; i<columns; ++i )
		for( int j=0; j<rows; ++j )
			tmp[i] += conj(A[j][i]) * v[j];

	return tmp;
}
*/

template <typename Type>
Matrix<complex<Type> > multTr( const Matrix<complex<Type> > &A1,
                               const Matrix<complex<Type> > &A2 )
{

	int rows = A1.rows();
	int columns = A2.rows();
	int K = A1.cols();

	Matrix<complex<Type> > tmp( rows, columns );
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			for( int k=0; k<K; ++k )
			   tmp[i][j] += A1[i][k] * conj(A2[j][k]);

	return tmp;
}


/*
template <typename Type>
Matrix<complex<Type> > multTr( const vector<complex<Type> > &a,
                               const vector<complex<Type> > &b )
{
	int rows = a.dim();
	int columns = b.dim();

	Matrix<complex<Type> > tmp( rows, columns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			tmp[i][j] = a[i]*conj(b[j]);

	return tmp;
}


template <typename Type>
Matrix<Type> eye( int N, const Type &x )
{
    Matrix<Type> tmp( N, N );
	for( int i=0; i<N; ++i )
		tmp[i][i] = x;

	return tmp;
}
*/


template <typename Type>
vector<Type> diag( const Matrix<Type> &A )
{
	int nColumn = A.rows();
	if( nColumn > A.cols() )
		nColumn = A.cols();

	vector<Type> tmp( nColumn );
	for( int i=0; i<nColumn; ++i )
		tmp[i] = A[i][i];

	return tmp;
}


template <typename Type>
Matrix<Type> diag( const vector<Type> &d )
{
	int N = d.size();

	Matrix<Type> tmp( N, N );
	for( int i=0; i<N; ++i )
		tmp[i][i] = d[i];

	return tmp;
}

/*
template <typename Type>
Type norm( const Matrix<Type> &A )
{
	int m = A.rows();
	int n = A.cols();

	Type sum = 0;
	for( int i=1; i<=m; ++i )
		for( int j=1; j<=n; ++j )
            sum += A(i,j) * A(i,j);

	return sqrt(sum);
}

template <typename Type>
Type norm( const Matrix<complex<Type> > &A )
{
	int m = A.rows();
	int n = A.cols();

	Type sum = 0;
	for( int i=1; i<=m; ++i )
		for( int j=1; j<=n; ++j )
            sum += norm(A(i,j));

	return sqrt(sum);
}


template <typename Type> void swap( Matrix<Type> &lhs, Matrix<Type> &rhs )
{
    int m = lhs.rows();
	int n = lhs.cols();

	assert( m == rhs.rows() );
	assert( n == rhs.cols() );

	for( int i=1; i<=m; ++i )
		for( int j=1; j<=n; ++j )
            swap( lhs(i,j), rhs(i,j) );
}

*/
template <typename Type>
vector<Type> sum( const Matrix<Type> &A )
{
	int m = A.rows();
	int n = A.cols();
	vector<Type> sum(n);

	for( int j=1; j<=n; ++j )
		for( int i=1; i<=m; ++i )
            sum.at(j-1) += A(i,j);

	return sum;
}


template <typename Type>
vector<Type> min( const Matrix<Type> &A )
{
	int m = A.rows();
	int n = A.cols();
	vector<Type> sum(n);

	for( int j=1; j<=n; ++j )
	{
	    Type tmp = A(1,j);
        for( int i=2; i<=m; ++i )
            if( tmp > A(i,j) )
                tmp = A(i,j);
        sum.at(j-1) = tmp;
	}

	return sum;
}


template <typename Type>
vector<Type> max( const Matrix<Type> &A )
{
	int m = A.rows();
	int n = A.cols();
	vector<Type> sum(n);

	for( int j=1; j<=n; ++j )
	{
	    Type tmp = A(1,j);
        for( int i=2; i<=m; ++i )
            if( tmp < A(i,j) )
                tmp = A(i,j);
        sum.at(j-1) = tmp;
	}

	return sum;
}


template <typename Type>
inline vector<Type> mean( const Matrix<Type> &A )
{
	return sum(A) / Type(A.rows());
}


template <typename Type>
Matrix<complex<Type> > complexMatrix( const Matrix<Type> &rA )
{
	int rows = rA.rows();
	int columns = rA.cols();

    Matrix<complex<Type> > cA( rows, columns );
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			cA[i][j] = rA[i][j];

    return cA;
}
/*
template <typename Type>
Matrix<complex<Type> > complexMatrix( const Matrix<Type> &mR,
                                      const Matrix<Type> &mI )
{
	int rows = mR.rows();
	int columns = mR.cols();

	assert( rows == mI.rows() );
	assert( columns == mI.cols() );

    Matrix<complex<Type> > cA( rows, columns );
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			cA[i][j] = complex<Type>( mR[i][j], mI[i][j] );

    return cA;
}
*/
template <typename Type>
Matrix<Type> abs( const Matrix<complex<Type> > &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = abs( A[i][j] );

    return tmp;
}


template <typename Type>
Matrix<Type> arg( const Matrix<complex<Type> > &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = arg( A[i][j] );

    return tmp;
}


template <typename Type>
Matrix<Type> real( const Matrix<complex<Type> > &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = A[i][j].real();

    return tmp;
}


template <typename Type>
Matrix<Type> imag( const Matrix<complex<Type> > &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = A[i][j].imag();

    return tmp;
}