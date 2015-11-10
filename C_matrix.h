/*
 #
 #  File            : C_matrix.h
 #                    ( C++ header file )
 #
 #  Description     : pas besoin
 #
 #  Project manager : moi
 #
 #  Licenses        : This file is 'dual-licensed', you have to choose one
 #                    of the two licenses below to apply.
 #
 #                    CeCILL-C
 #                    The CeCILL-C license is close to the GNU LGPL.
 #                    ( http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html )
 #
 #                or  CeCILL v2.0
 #                    The CeCILL license is compatible with the GNU GPL.
 #                    ( http://www.cecill.info/licences/Licence_CeCILL_V2-en.html )
 #
 #  This software is governed either by the CeCILL or the CeCILL-C license
 #  under French law and abiding by the rules of distribution of free software.
 #  You can  use, modify and or redistribute the software under the terms of
 #  the CeCILL or CeCILL-C licenses as circulated by CEA, CNRS and INRIA
 #  at the following URL: "http://www.cecill.info".
 #
 #  As a counterpart to the access to the source code and  rights to copy,
 #  modify and redistribute granted by the license, users are provided only
 #  with a limited warranty  and the software's author,  the holder of the
 #  economic rights,  and the successive licensors  have only  limited
 #  liability.
 #
 #  In this respect, the user's attention is drawn to the risks associated
 #  with loading,  using,  modifying and/or developing or reproducing the
 #  software by the user in light of its specific status of free software,
 #  that may mean  that it is complicated to manipulate,  and  that  also
 #  therefore means  that it is reserved for developers  and  experienced
 #  professionals having in-depth computer knowledge. Users are therefore
 #  encouraged to load and test the software's suitability as regards their
 #  requirements in conditions enabling the security of their systems and/or
 #  data to be ensured and,  more generally, to use and operate it in the
 #  same conditions as regards security.
 #
 #  The fact that you are presently reading this means that you have had
 #  knowledge of the CeCILL and CeCILL-C licenses and that you accept its terms.
 #
*/

#ifndef C_MATRIX_H
#define C_MATRIX_H

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <C_vector.h>
#include <cmath>
#include <float.h>
#include <cstring>
#include <vector>


const double epsilon=0.0000000001;

#ifndef ABS
#define ABS(a) ((a) >= 0 ? (a) : -(a))
#endif
#ifndef SQR
#define SQR(a) ((a)*(a))
#endif
#ifndef SIGN
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#endif
#ifndef SIGN_
#define SIGN_(a) ((a) >= 0.0 ? 1 : -1)
#endif
#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a) >= (b) ? (b) : (a))
#endif
#ifndef MINMOD
#define MINMOD(a,b) (0.5*(SIGN_(a)+SIGN_(b))*MIN(ABS(a),ABS(b)))
#endif
#ifndef Pi
#define Pi 3.141592653589793238462643383279502884197169399375105820974944592L
#endif

#define SWAP(x,y) do \
   { unsigned char swap_temp[sizeof(x) == sizeof(y) ? (signed)sizeof(x) : -1]; \
     memcpy(swap_temp,&y,sizeof(x)); \
     memcpy(&y,&x,       sizeof(x)); \
     memcpy(&x,swap_temp,sizeof(x)); \
    } while(0)


#define NEAREST 0
#define LINEAR  1
#define SMALL_NUM_F 1e-37
#define SMALL_NUM_D 1e-307

//this is a template, meaning that you need to specify the type when you instanciate a object.
template<class dataType> class C_matrix
{
public:
    //constructors and destructors
    //default ctor: create empty matrix (use resize to set a new size and allocate the image container)
    C_matrix();
    //ctor: create mtrix of size _L rows, _C columns
    C_matrix(int _L, int _C);
    //copy ctor
    C_matrix(const C_matrix &M);
    //ctor: create a matrix copying _M
    C_matrix(dataType **_M, int _L, int _C);
    //dtor: free memory
    virtual ~C_matrix();


    //IO methods
    dataType& operator()(const int l, const int c);
    const dataType& operator()(const int l, const int c) const;
    dataType& operator()(const long l, const long c);
    const dataType& operator()(const long l, const long c) const;
    dataType& operator()(const short l, const short c);
    const dataType& operator()(const short l, const short c) const;
    dataType& operator()(const unsigned int l, const unsigned int c);
    const dataType& operator()(const unsigned int l, const unsigned int c) const;
    dataType& operator()(const unsigned short l, const unsigned short c);
    const dataType& operator()(const unsigned short l, const unsigned short c) const;
    dataType& operator()(const unsigned long l, const unsigned long c);
    const dataType& operator()(const unsigned long l, const unsigned long c) const;
    dataType& operator()(const float l, const float c);
    const dataType& operator()(const float l, const float c, unsigned short INTERPOLATION=NEAREST) const;
    dataType& operator()(const double l, const double c);
    const dataType& operator()(const double l, const double c, unsigned short INTERPOLATION=NEAREST) const;

    //overload operators
    C_matrix<dataType> operator= (C_matrix const& c);
    C_matrix<dataType> operator= (dataType const& x);
    C_matrix<dataType> operator+ (C_matrix const& c);
    C_matrix<dataType> operator+ (dataType const& x);
    C_matrix<dataType> operator- (C_matrix const& c);
    C_matrix<dataType> operator- (dataType const& x);
    C_matrix<dataType> operator* (C_matrix const& c);
    //C_vector<dataType>& operator* (C_vector<dataType> const& c); //obselete
    C_matrix operator* (const dataType& x);

    C_matrix<dataType> operator== (C_matrix const& c);
    C_matrix<dataType> operator== (dataType const& x);

    C_matrix<dataType> operator> (C_matrix const& c);
    C_matrix<dataType> operator> (dataType const& x);
    C_matrix<dataType> operator>= (C_matrix const& c);
    C_matrix<dataType> operator>= (dataType const& x);

    C_matrix<dataType> operator< (C_matrix const& c);
    C_matrix<dataType> operator< (dataType const& x);
    C_matrix<dataType> operator<= (C_matrix const& c);
    C_matrix<dataType> operator<= (dataType const& x);

    //offset
    int offsetL, offsetC;

    //cast operator (later)
    //operator int();

    //usual operations
    C_matrix<dataType> m_abs(void); //must be multithreaded
    dataType maxVal(void);
    dataType minVal(void);
    double sum(void);
    double mean(void);
    double var(void);
    C_matrix<dataType> SQRT(void); //must be multithreaded
    C_matrix<dataType> dotPower(double alpha);//octave .^ //must be multithreaded
    C_matrix<dataType> dotExp(double alpha=1.0); //must be multithreaded
    C_matrix<dataType> dotLog(double alpha=1.0); //must be multithreaded
    C_matrix<dataType> dotProduct(C_matrix const& B);//octave .* //must be multithreaded
    C_matrix<dataType> dotDiv(C_matrix const& B);//octave ./ //must be multithreaded
    C_matrix<dataType> dotProduct(C_matrix& B); //must be multithreaded
    C_matrix<dataType> dotDiv(C_matrix& B); //must be multithreaded
    C_matrix<dataType> subset(int lbegin, int lend, int cbegin, int cend);
    void subset(C_matrix<dataType> M, int lbegin, int lend, int cbegin, int cend);

    //distance map
    C_matrix<double> bwdistEuclidean(void); //Algorithme de Danielson

    //convolution //must be multithreaded
    C_matrix<dataType> conv2(C_matrix<dataType> const& h);
    C_matrix<dataType> conv2(C_matrix<dataType>& h);
    C_matrix<dataType> gradX(void);
    C_matrix<dataType> gradY(void);

    //random generator
    void random(void);
    void randomf(void);
    void randomGauss(std::vector<double> mu, std::vector<double> sigma);
    void randomGauss(double mu, double sigma);


    //mesh
    void meshRow(dataType minValue, dataType maxValue);
    void meshColumn(dataType minValue, dataType maxValue);
    //linespace (later)

    //matrix inversion
    C_matrix<dataType> inv(void);
    C_matrix<dataType> pseudoInv(void);

    //fft (later)

    //tools for linear system solving
    //LU
    C_matrix<dataType> Transpose(void);
    C_matrix<dataType> LU(void);
    C_matrix<dataType> LUP(C_vector<int> &Indice);
    C_matrix<dataType> LMU(void);
    C_vector<dataType> LineAlgEq_LU(C_vector<dataType> &B);
    //svd
    std::vector<C_matrix<dataType> > svd(void);
    //diagonalisation: sym sys
    std::vector<C_matrix<dataType> > eigSym(bool yesvec);
protected:
    void tred2(C_matrix<dataType>* z, C_matrix<dataType>* d, C_matrix<dataType>* e, bool yesvecs);
    void tqli(C_matrix<dataType> *z, C_matrix<dataType> *d, C_matrix<dataType> *e, bool yesvecs);
    void eigsrt(C_matrix<dataType> *d, C_matrix<dataType> *z);
    void sortEig(C_matrix<dataType> *d, C_matrix<dataType> *z, bool yesvecs);
public:


    //other tools
    int getNbRow(void)const {return m_L;}
    int getNbColumn(void)const {return m_C;}
    int numel(void)const{return m_L*m_C;}
    void show(void);
    void save(std::string fileName);

    int getIndex(int l, int c){return l*m_C+c;}
    int getRow(int idx){return floor(idx/m_C);}//floor();
    int getColumn(int idx){return idx - (floor(idx/m_C)*m_C);}// return idx%m_C;

    //to use carefully
    bool resize(int newL, int newC);
    int endL;
    int endC;

protected:
    int m_L;
    int m_C;
    dataType* m_A;

    //privates tools to allocate and free memory
    dataType *allocate(int M, int N);
    void deallocation(dataType* B, int M, int N);

    double pythag(const double a, const double b);

};
template<class dataType> C_matrix<dataType>::C_matrix()
{
    m_L=0;
    m_C=0;
    m_A=NULL;
    endL = m_L-1;
    endC = m_C-1;
}
template<class dataType> C_matrix<dataType>::C_matrix(int _L, int _C) : m_L(_L),m_C(_C)
{
    m_A = allocate(m_L,m_C);
    if(m_A==NULL)
    {
        m_L=0;
        m_C=0;
    }
    endL = m_L-1;
    endC = m_C-1;
    *this = 0.0;
}

template<class dataType> C_matrix<dataType>::C_matrix(dataType **_M, int _L, int _C) : m_L(_L),m_C(_C)
{
    m_A = allocate(m_L, m_C);
    if(m_A==NULL)
    {
        m_L=0;
        m_C=0;
    }
    else
    {
        for(int i=0 ; i<m_L ; i++)
        {
            for(int j=0 ; j<m_C ; j++)
                m_A[getIndex(i,j)] = _M[i][j];
        }
    }
    endL = m_L-1;
    endC = m_C-1;
}

template<class dataType> C_matrix<dataType>::C_matrix(const C_matrix &X)
{
    m_L = X.getNbRow();
    m_C = X.getNbColumn();
    m_A = allocate(m_L, m_C);
    if(m_A==NULL)
    {
        m_L=0;
        m_C=0;
    }
    else
    {
        for(int i=0 ; i<m_L ; i++)
        {
            for(int j=0 ; j<m_C ; j++)
            {
                m_A[getIndex(i,j)] = X(i,j);
            }
        }
    }
    endL = m_L-1;
    endC = m_C-1;

    offsetL = X.offsetL;
    offsetC = X.offsetC;

}

template<class dataType> C_matrix<dataType>::~C_matrix()
{
    deallocation(m_A,m_L,m_C);
}

//to use carefully
template<class dataType> bool C_matrix<dataType>::resize(int newL, int newC)
{
    deallocation(m_A,m_L,m_C);
    m_A = NULL;
    m_L = 0;
    m_C = 0;
    m_A = allocate(newL, newC);
    if(m_A==NULL) return false;
    m_L = newL;
    m_C = newC;
    endL = m_L-1;
    endC = m_C-1;
    *this = 0.0;
    return true;
}

template<class dataType> dataType* C_matrix<dataType>::allocate(int M, int N)
{
    return new dataType[M*N];
//    dataType* A = new dataType[M*N];
//    return A
//    if(A==NULL) return A;

//    for(int i=0 ; i<M ; i++)
//        A[i] = new dataType[N];

//    //should check if all allocation is OK
//    bool cond=true;
//    for(int i=0 ; i<M ; i++)
//    {
//        if(A[i]==NULL)
//        {
//            cond=false;
//            i=M;
//        }
//    }
//    if(!cond)
//    {
//        deallocation(A,M,N);
//    }
//    return A;
}

template<class dataType> void C_matrix<dataType>::deallocation(dataType *B, int M, int N)
{
    if(B!=NULL)
    {
//        for(int i=0 ; i<M ; i++)
//        {
//            if(B[i]!=NULL)
//            {
//                delete [] B[i];
//            }
//        }
        delete [] B;
        B=NULL;
    }
    return;
}


template<class dataType> void C_matrix<dataType>::show(void)
{
    std::cout << "--------------------------------------------------------" << std::endl;
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            std::cout << m_A[getIndex(i,j)] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << std::endl;
}

template<class dataType> void C_matrix<dataType>::save(std::string fileName)
{
    std::ofstream myfileX;
    myfileX.open (fileName.data(), std::ios::out );

    if(myfileX.is_open())
    {
        for(unsigned short i=0 ; i<m_L ; i++)
        {
            for(unsigned short j=0 ; j<m_C ; j++)
            {
                myfileX << m_A[getIndex(i,j)] << "\t";
            }
            myfileX << std::endl;
        }
        myfileX.close();
    }
}



template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator= (C_matrix const& c)
{
    if(this==&c)return *this;

    if(c.getNbColumn()!=this->getNbColumn() || c.getNbRow()!=this->getNbRow())
    {
        this->resize(c.getNbRow(),c.getNbColumn());
        //throw "dimension matrix must agree";
    }

    for(unsigned short i=0 ; i<this->m_L ; i++)
    {
        for(unsigned short j=0 ; j<this->m_C ; j++)
        {
            this->m_A[getIndex(i,j)] = c(i,j);
        }
    }
    return *this;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator= (dataType const& x)
{
    for(unsigned short i=0 ; i<this->m_L ; i++)
    {
        for(unsigned short j=0 ; j<this->m_C ; j++)
        {
            m_A[getIndex(i,j)] = x;
        }
    }

    return *this;
}




template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator+ (C_matrix const& c)
{
    if(c.getNbColumn()!=this->getNbColumn() || c.getNbRow()!=this->getNbRow()) throw "dimension matrix must agree";
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            B(i,j) = m_A[getIndex(i,j)] + c(i,j);
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator+ (dataType const& x)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<this->m_L ; i++)
    {
        for(unsigned short j=0 ; j<this->m_C ; j++)
        {
            B(i,j) = m_A[getIndex(i,j)] + x;
        }
    }

    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator- (C_matrix const& c)
{
    if(c.getNbColumn()!=this->getNbColumn() || c.getNbRow()!=this->getNbRow()) throw "dimension matrix must agree";
    C_matrix<dataType> B(m_L,m_C);// = new C_matrix(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            B(i,j) = m_A[getIndex(i,j)] - c(i,j);
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator- (dataType const& x)
{
    C_matrix<dataType> B(m_L,m_C);// = new C_matrix
    for(unsigned short i=0 ; i<this->m_L ; i++)
    {
        for(unsigned short j=0 ; j<this->m_C ; j++)
        {
            B(i,j) = m_A[getIndex(i,j)] - x;
        }
    }

    return B;
}


template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator* (C_matrix const& c)
{
    if(m_C!=c.getNbRow()) throw "mismatch dimension matrix";

    C_matrix B(m_L,c.getNbColumn());
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<c.getNbColumn() ; j++)
        {
            dataType S = (dataType) 0;
            for(unsigned short k=0 ; k<m_C ; k++)
            {
                S += m_A[getIndex(i,j)]*c(k,j);
            }
            B(i,j) = S;
        }
    }
    return B;
}


template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator* (const dataType& x)
{
    C_matrix B(m_L,m_C);
    for(unsigned short i=0 ; i<this->m_L ; i++)
    {
        for(unsigned short j=0 ; j<this->m_C ; j++)
        {
            B(i,j) = m_A[getIndex(i,j)]*x;
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator== (C_matrix const& c)
{
    if(c.getNbColumn()!=this->getNbColumn() || c.getNbRow()!=this->getNbRow()) throw "dimension matrix must agree";
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[getIndex(i,j)]==c(i,j)) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator== (dataType const& x)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[getIndex(i,j)]==x) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator> (C_matrix const& c)
{
    if(c.getNbColumn()!=this->getNbColumn() || c.getNbRow()!=this->getNbRow()) throw "dimension matrix must agree";
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[getIndex(i,j)]>c(i,j)) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator> (dataType const& x)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[getIndex(i,j)]>x) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator>= (C_matrix const& c)
{
    if(c.getNbColumn()!=this->getNbColumn() || c.getNbRow()!=this->getNbRow()) throw "dimension matrix must agree";
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[getIndex(i,j)]>=c(i,j)) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator>= (dataType const& x)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[getIndex(i,j)]>=x) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}




template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator< (C_matrix const& c)
{
    if(c.getNbColumn()!=this->getNbColumn() || c.getNbRow()!=this->getNbRow()) throw "dimension matrix must agree";
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[getIndex(i,j)]<c(i,j)) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator< (dataType const& x)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[getIndex(i,j)]<x) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator<= (C_matrix const& c)
{
    if(c.getNbColumn()!=this->getNbColumn() || c.getNbRow()!=this->getNbRow()) throw "dimension matrix must agree";
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[getIndex(i,j)]<=c(i,j)) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator<= (dataType const& x)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[getIndex(i,j)]<=x) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}






template<class dataType> dataType& C_matrix<dataType>::operator()(const int l, const int c)
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[getIndex(l,c)];
}
template<class dataType> const dataType& C_matrix<dataType>::operator()(const int l, const int c) const
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[l*m_C+l];//TODO explain why getIndex(l,c) generates an error
}

template<class dataType> dataType& C_matrix<dataType>::operator()(const long l, const long c)
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[getIndex(l,c)];
}
template<class dataType> const dataType& C_matrix<dataType>::operator()(const long l, const long c) const
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[getIndex(l,c)];
}

template<class dataType> dataType& C_matrix<dataType>::operator()(const short l, const short c)
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[getIndex(l,c)];
}
template<class dataType> const dataType& C_matrix<dataType>::operator()(const short l, const short c) const
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[getIndex(l,c)];
}

template<class dataType> dataType& C_matrix<dataType>::operator()(const unsigned int l, const unsigned int c)
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[getIndex(l,c)];
}
template<class dataType> const dataType& C_matrix<dataType>::operator()(const unsigned int l, const unsigned int c) const
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[getIndex(l,c)];
}
template<class dataType> dataType& C_matrix<dataType>::operator()(const unsigned short l, const unsigned short c)
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[getIndex(l,c)];
}
template<class dataType> const dataType& C_matrix<dataType>::operator()(const unsigned short l, const unsigned short c) const
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[l*m_C+c];//TODO explain why getIndex(l,c) generates an error
}

template<class dataType> dataType& C_matrix<dataType>::operator()(const unsigned long l, const unsigned long c)
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[getIndex(l,c)];
}
template<class dataType> const dataType& C_matrix<dataType>::operator()(const unsigned long l, const unsigned long c) const
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[getIndex(l,c)];
}

template<class dataType> dataType& C_matrix<dataType>::operator()(const float l, const float c)
{
    if( (l<0.0) || (l>=(float) (this->m_L)) || (c<0.0) || (c>= (float) (this->m_C)))
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    unsigned int ll = (unsigned int) floor((double) l);
    unsigned int cc = (unsigned int) floor((double) c);
    return m_A[getIndex(ll,cc)];

}
template<class dataType> const dataType& C_matrix<dataType>::operator()(const float l, const float c, unsigned short INTERPOLATION) const
{
    if( (l<0.0) || (l>=(float) (this->m_L)) || (c<0.0) || (c>= (float) (this->m_C)))
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    if(INTERPOLATION==NEAREST)
    {
        unsigned int ll = (unsigned int) floor((double) l);
        unsigned int cc = (unsigned int) floor((double) c);
        return m_A[getIndex(ll,cc)];
    }
    else
    {
        //not yet available
        unsigned int ll = (unsigned int) floor((double) l);
        unsigned int cc = (unsigned int) floor((double) c);
        //return m_A[ll][cc];
        if(l<SMALL_NUM_F && c<SMALL_NUM_F)
        {
            return m_A[getIndex(ll,cc)];
        }
        if(l<SMALL_NUM_F)//column interpolation
        {
            unsigned int cc1 = cc+1;
            return (dataType) (((double)m_A[getIndex(ll,cc)])*( ((double)cc1) - ((double)c) ) + ((double)m_A[getIndex(ll,cc1)])*( ((double)c) - ((double)cc) ));
        }
        if(c<SMALL_NUM_F)//row interpolation
        {
            unsigned int ll1 = ll+1;
            return (dataType) (((double)m_A[getIndex(ll,cc)])*( ((double)ll1) - ((double)l) ) + ((double)m_A[getIndex(ll1,cc)])*( ((double)l) - ((double)ll) ));
        }
        //main case: bilinear interpolation
        unsigned int ll1 = ll+1;
        unsigned int cc1 = cc+1;
        double valuell = (((double)m_A[getIndex(ll,cc)])*( ((double)cc1) - ((double)c) ) + ((double)m_A[getIndex(ll,cc1)])*( ((double)c) - ((double)cc) ));
        double valuell1 = (((double)m_A[getIndex(ll1,cc)])*( ((double)cc1) - ((double)c) ) + ((double)m_A[getIndex(ll1,cc1)])*( ((double)c) - ((double)cc) ));
        return (dataType) (valuell*( ((double)ll1) - ((double)l) ) + valuell1*( ((double)l) - ((double)ll) ));
    }
}

template<class dataType> dataType& C_matrix<dataType>::operator()(const double l, const double c)
{
    if( (l<0.0) || (l>=(double) (this->m_L)) || (c<0.0) || (c>= (double) (this->m_C)))
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    unsigned int ll = (unsigned int) floor(l);
    unsigned int cc = (unsigned int) floor(c);
    return m_A[getIndex(ll,cc)];

}
template<class dataType> const dataType& C_matrix<dataType>::operator()(const double l, const double c, unsigned short INTERPOLATION) const
{
    if( (l<0.0) || (l>=(double) (this->m_L)) || (c<0.0) || (c>= (double) (this->m_C)))
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    if(INTERPOLATION==NEAREST)
    {
        unsigned int ll = (unsigned int) floor(l);
        unsigned int cc = (unsigned int) floor(c);
        return m_A[getIndex(ll,cc)];
    }
    else
    {
        //not yet available
        unsigned int ll = (unsigned int) floor(l);
        unsigned int cc = (unsigned int) floor(c);
        //return m_A[ll][cc];
        if(l<SMALL_NUM_F && c<SMALL_NUM_F)
        {
            return m_A[getIndex(ll,cc)];
        }
        if(l<SMALL_NUM_F)//column interpolation
        {
            unsigned int cc1 = cc+1;
            return (dataType) (((double)m_A[getIndex(ll,cc)])*( ((double)cc1) - c ) + ((double)m_A[getIndex(ll,cc1)])*( c - ((double)cc) ));
        }
        if(c<SMALL_NUM_F)//row interpolation
        {
            unsigned int ll1 = ll+1;
            return (dataType) (((double)m_A[getIndex(ll,cc)])*( ((double)ll1) - l ) + ((double)m_A[getIndex(ll1,cc)])*( l - ((double)ll) ));
        }
        //main case: bilinear interpolation
        unsigned int ll1 = ll+1;
        unsigned int cc1 = cc+1;
        double valuell = (((double)m_A[getIndex(ll,cc)])*( ((double)cc1) - c ) + ((double)m_A[getIndex(ll,cc1)])*( c - ((double)cc) ));
        double valuell1 = (((double)m_A[getIndex(ll1,cc)])*( ((double)cc1) - c ) + ((double)m_A[getIndex(ll1,cc1)])*( c - ((double)cc) ));
        return (dataType) (valuell*( ((double)ll1) - l ) + valuell1*( l - ((double)ll) ));
    }
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::subset(int lbegin, int lend, int cbegin, int cend)
{
    if(lend<lbegin || cend<cbegin) throw "end must be larger than beginning";
    if(lbegin<0 || cbegin<0) throw "begining must be larger than 0";
    if(lend>this->endL || cend>this->endC) throw "end must be smaller size-1";

    C_matrix<dataType> B(lend-lbegin+1,cend-cbegin+1);
    for(int l=0 ; l<lend-lbegin+1 ; l++)
    {
        for(int c=0 ; c<cend-cbegin+1 ; c++)
        {
            B(l,c) = m_A[getIndex(lbegin+l,cbegin+c)];
        }
    }
    return B;
}

template<class dataType> void C_matrix<dataType>::subset(C_matrix<dataType> M, int lbegin, int lend, int cbegin, int cend)
{
    if(lend<lbegin || cend<cbegin) throw "end must be larger than beginning";
    if(lbegin<0 || cbegin<0) throw "begining must be larger than 0";
    if(lend>this->endL || cend>this->endC) throw "end must be smaller size-1 k";
    if(M.getNbRow()!=(lend-lbegin+1) || M.getNbColumn()!=(cend-cbegin+1)) throw "dimension matrix must agree";
    for(int l=0 ; l<M.getNbRow() ; l++)
    {
        for(int c=0 ; c<M.getNbColumn() ; c++)
        {
            m_A[getIndex(lbegin+l,cbegin+c)] = M(l,c);
        }
    }
    return;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::m_abs(void)
{
    C_matrix<dataType> B(this->m_L,this->m_C);// = new C_matrix
    for(unsigned short i=0 ; i<this->m_L ; i++)
    {
        for(unsigned short j=0 ; j<this->m_C ; j++)
        {
            B(i,j) = ABS(m_A[i][j]);
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::SQRT(void)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            B(i,j) = sqrt(m_A[getIndex(i,j)]);
        }
    }
    return B;
}

template<class dataType> C_matrix<double> C_matrix<dataType>::bwdistEuclidean(void)
{
    C_matrix<double> B(m_L,m_C), B2(m_L,m_C), B3(m_L,m_C), B4(m_L,m_C);
    C_matrix<double> dL(m_L,m_C), dC(m_L,m_C), dL2(m_L,m_C), dC2(m_L,m_C), dL3(m_L,m_C), dC3(m_L,m_C), dL4(m_L,m_C), dC4(m_L,m_C);
    double BIGNUM = ((double)m_L)*((double)m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(ABS(m_A[i][j])<SMALL_NUM_F)//if not in the ocean
            {
                B(i,j) = 0.0; B2(i,j) = 0.0; B3(i,j) = 0.0; B4(i,j) = 0.0;
                dL(i,j) = 0.0; dL2(i,j) = 0.0; dL3(i,j) = 0.0; dL4(i,j) = 0.0;
                dC(i,j) = 0.0; dC2(i,j) = 0.0; dC3(i,j) = 0.0; dC4(i,j) = 0.0;
            }
            else
            {
                B(i,j) = BIGNUM; B2(i,j) = BIGNUM; B3(i,j) = BIGNUM; B4(i,j) = BIGNUM;
                dL(i,j) = BIGNUM; dL2(i,j) = BIGNUM; dL3(i,j) = BIGNUM; dL4(i,j) = BIGNUM;
                dC(i,j) = BIGNUM; dC2(i,j) = BIGNUM; dC3(i,j) = BIGNUM; dC4(i,j) = BIGNUM;
            }
        }
    }

    double v, vi, vj;
    bool firstIslandMet=false;
    //first sweep
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(!firstIslandMet && B(i,j)<BIGNUM)
            {
                firstIslandMet = true;
            }
            if(B(i,j)>SMALL_NUM_F && firstIslandMet)
            {
                v = SQR(dL(i,j))+SQR(dC(i,j));
                if(i>0) vi = SQR(dL((unsigned short)(i-1),j)-1.0)+SQR(dC((unsigned short)(i-1),j));
                else vi = BIGNUM;

                if(j>0) vj = SQR(dL(i,(unsigned short)(j-1)))+SQR(dC(i,(unsigned short)(j-1))-1.0);
                else vj = BIGNUM;

                if(v<vj && v<vi)
                {
                    //no change
                }
                else if(vi<v && vi<vj)
                {
                    if(i>0)
                    {
                        dL(i,j) = dL((unsigned short)(i-1),j)-1.0;
                        dC(i,j) = dC((unsigned short)(i-1),j);
                        if(B(i,j)>sqrt(vi)) B(i,j) = sqrt(vi);
                    }
                }else
                {
                    if(j>0)
                    {
                        dL(i,j) = dL(i,(unsigned short)(j-1));
                        dC(i,j) = dC(i,(unsigned short)(j-1))-1.0;
                        if(B(i,j)>sqrt(vj)) B(i,j) = sqrt(vj);
                    }
                }
            }
        }
    }

    //second sweep
    firstIslandMet=false;
    for(long i=(long)(m_L-1) ; i>=0 ; i--)
    {
        for(long j=0 ; j<m_C ; j++)
        {
            if(!firstIslandMet && ABS(m_A[i][j])<SMALL_NUM_F)
            {
                firstIslandMet = true;
            }
            if(B2(i,j)>SMALL_NUM_F && firstIslandMet)
            {
                try
                {
                    v = SQR(dL2(i,j))+SQR(dC2(i,j));
                    if(i<(long)(m_L-1)) vi = SQR(dL2((long)(i+1),j)+1.0)+SQR(dC2((long)(i+1),j));
                    else vi = BIGNUM;

                    if(j>0) vj = SQR(dL2(i,(long)(j-1)))+SQR(dC2(i,(long)(j-1))-1.0);
                    else vj = BIGNUM;

                    if(v<vj && v<vi)
                    {
                        //no change
                    }
                    else if(vi<v && vi<vj)
                    {
                        if(i<(long)(m_L-1))
                        {
                            dL2(i,j) = dL2((long)(i+1),j)+1.0;
                            dC2(i,j) = dC2((long)(i+1),j);
                            B2(i,j) = sqrt(vi);
                        }
                    }
                    else
                    {
                        if(j>0)
                        {
                            dL2(i,j) = dL2(i,(long)(j-1));
                            dC2(i,j) = dC2(i,(long)(j-1))-1.0;
                            B2(i,j) = sqrt(vj);
                        }
                    }
                }
                catch(const char* a)
                {
                    std::cout << a << std::endl;
                }
            }
            if(B(i,j)>B2(i,j)) B(i,j)=B2(i,j);
        }
    }

    //third sweep
    firstIslandMet=false;
    for(long i=(long)(m_L-1) ; i>=0 ; i--)
    {
        for(long j=(m_C-1) ; j>=0 ; j-- )
        {
            if(!firstIslandMet && ABS(m_A[i][j])<SMALL_NUM_F)
            {
                firstIslandMet = true;
            }
            if(B3(i,j)>SMALL_NUM_F && firstIslandMet)
            {
                try
                {
                    v = SQR(dL3(i,j))+SQR(dC3(i,j));
                    if(i<(long)(m_L-1)) vi = SQR(dL3((long)(i+1),j)+1.0)+SQR(dC3((long)(i+1),j));
                    else vi = BIGNUM;

                    if(j<(long)(m_C-1)) vj = SQR(dL3(i,(long)(j+1)))+SQR(dC3(i,(long)(j+1))+1.0);
                    else vj = BIGNUM;

                    if(v<vj && v<vi)
                    {
                        //no change
                    }
                    else if(vi<v && vi<vj)
                    {
                        if(i<(long)(m_L-1))
                        {
                            dL3(i,j) = dL3((long)(i+1),j)+1.0;
                            dC3(i,j) = dC3((long)(i+1),j);
                            B3(i,j) = sqrt(vi);
                        }
                    }
                    else
                    {
                        if(j<(long)(m_C-1))
                        {
                            dL3(i,j) = dL3(i,(long)(j+1));
                            dC3(i,j) = dC3(i,(long)(j+1))+1.0;
                            B3(i,j) = sqrt(vj);
                        }
                    }
                }
                catch(const char* a)
                {
                    std::cout << a << std::endl;
                }
            }
            if(B(i,j)>B3(i,j)) B(i,j)=B3(i,j);
        }
    }


    //fourth sweep
    firstIslandMet=false;
    std::cout << (unsigned short)(m_L-1) << " " << (unsigned short)(m_C-1) << std::endl;
    for(long i=0 ; i<m_L ; i++)
    {
        for(long j=(long)(m_C-1) ; j>=0 ; j--)
        {
            if(!firstIslandMet && ABS(m_A[i][j])<SMALL_NUM_F)
            {
                firstIslandMet = true;
            }
            if(B4(i,j)>SMALL_NUM_F && firstIslandMet)
            {
                try
                {
                    v = SQR(dL4(i,j))+SQR(dC4(i,j));
                    if(i>0) vi = SQR(dL4((long)(i-1),j)-1.0)+SQR(dC4((long)(i-1),j));
                    else vi = BIGNUM;

                    if(j<(long)(m_C-1)) vj = SQR(dL4(i,(long)(j+1)))+SQR(dC4(i,(long)(j+1))+1.0);
                    else vj = BIGNUM;

                    if(v<vj && v<vi)
                    {
                        //no change
                    }
                    else if(vi<v && vi<vj)
                    {
                        if(i>0)
                        {
                            dL4(i,j) = dL4((long)(i-1),j)-1.0;
                            dC4(i,j) = dC4((long)(i-1),j);
                            B4(i,j) = sqrt(vi);
                        }
                    }
                    else
                    {
                        if(j<(long)(m_C-1))
                        {
                            dL4(i,j) = dL4(i,(long)(j+1));
                            dC4(i,j) = dC4(i,(long)(j+1))+1.0;
                            B4(i,j) = sqrt(vj);
                        }
                    }
                }
                catch(const char* a)
                {
                    std::cout << a << std::endl;
                }
            }
            if(B(i,j)>B4(i,j)) B(i,j)=B4(i,j);
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::dotPower(double alpha)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            B(i,j) = (dataType) pow((double) m_A[getIndex(i,j)], alpha);
        }
    }
    return B;
}
template<class dataType> C_matrix<dataType> C_matrix<dataType>::dotExp(double alpha)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            B(i,j) = (dataType) exp(alpha*(double) m_A[getIndex(i,j)]);
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::dotLog(double alpha)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            B(i,j) = (dataType) log(alpha*(double) m_A[getIndex(i,j)]);
        }
    }
    return B;
}

template<class dataType> dataType C_matrix<dataType>::maxVal(void)
{
    dataType B = m_A[0];
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(B<m_A[getIndex(i,j)]) B = m_A[getIndex(i,j)];
        }
    }
    return B;
}

template<class dataType> dataType C_matrix<dataType>::minVal(void)
{
    dataType B = m_A[0];
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(B>m_A[getIndex(i,j)]) B = m_A[getIndex(i,j)];
        }
    }
    return B;
}

template<class dataType>double C_matrix<dataType>::sum(void)
{
    double B = 0.0;
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            B += (double) m_A[getIndex(i,j)];
        }
    }
    return B;
}

template<class dataType>double C_matrix<dataType>::mean(void)
{
    return sum()/(((double)m_L)*((double)m_C));
}

template<class dataType>double C_matrix<dataType>::var(void)
{
    double B = mean(), V=0.0;
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            V += SQR(((double) m_A[getIndex(i,j)]) - B);
        }
    }
    return V/(((double)m_L)*((double)m_C) - 1.0);
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::Transpose(void)
{
    C_matrix<dataType> B(this->m_C,this->m_L);
    for(unsigned short i=0 ; i<this->m_L ; i++)
    {
        for(unsigned short j=0 ; j<this->m_C ; j++)
        {
            B(j,i) = m_A[getIndex(i,j)];
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::dotProduct(C_matrix const& B)
{
    if(this->getNbRow()!=B.getNbRow() || this->getNbColumn()!=B.getNbColumn()) throw "mismatch dimension matrix";
    C_matrix C(m_L,m_C);
    for(unsigned short i=0 ; i<C.getNbRow() ; i++)
    {
        for(unsigned short j=0 ; j<C.getNbColumn() ; j++)
        {
            C(i,j) = m_A[getIndex(i,j)]*B(i,j);
        }
    }
    return C;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::dotProduct(C_matrix& B)
{
    if(this->getNbRow()!=B.getNbRow() || this->getNbColumn()!=B.getNbColumn()) throw "mismatch dimension matrix";
    C_matrix C(m_L,m_C);
    for(unsigned short i=0 ; i<C.getNbRow() ; i++)
    {
        for(unsigned short j=0 ; j<C.getNbColumn() ; j++)
        {
            C(i,j) = (m_A[getIndex(i,j)])*(B(i,j));
        }
    }
    return C;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::dotDiv(C_matrix const& B)
{
    if(this->getNbRow()!=B.getNbRow() || this->getNbColumn()!=B.getNbColumn()) throw "mismatch dimension matrix";
    C_matrix<dataType> C(m_L,m_C);
    for(unsigned short i=0 ; i<C.getNbRow() ; i++)
    {
        for(unsigned short j=0 ; j<C.getNbColumn() ; j++)
        {
            C(i,j) = m_A[getIndex(i,j)]/(B(i,j));
        }
    }
    return C;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::dotDiv(C_matrix& B)
{
    if(this->getNbRow()!=B.getNbRow() || this->getNbColumn()!=B.getNbColumn()) throw "mismatch dimension matrix";
    C_matrix C(m_L,m_C);
    for(unsigned short i=0 ; i<C.getNbRow() ; i++)
    {
        for(unsigned short j=0 ; j<C.getNbColumn() ; j++)
        {
            C(i,j) = (m_A[getIndex(i,j)])/(B(i,j));
        }
    }
    return C;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::conv2(C_matrix<dataType> const& h)
{
    if((h.getNbRow() % 2 == 0) || (h.getNbColumn() % 2 == 0))
    {
        std::cout << "even kernel size are not handled. Try to change the size of the convolution kernel to an odd size" << std::endl;
    }

    //returned matrix
    C_matrix<dataType> tmp(m_L,m_C);
    tmp = 0.0;


    //compute convolution
    unsigned long K = h.getNbRow()/2;
    unsigned long L = h.getNbColumn()/2;
    for(unsigned long i=0 ; i<m_L ; i++)
    {
        for(unsigned long j=0 ; j<m_C ; j++)
        {
            for(unsigned long ih=0 ; ih<2*K+1 ; ih++)
            {
                for(unsigned long jh=0 ; jh<2*L+1 ; jh++)
                {
                    if((i+K)>=ih && (j+L)>=jh && (i-ih+K)<m_L && (j-jh+L)<m_C) //check image boundaries
                    {
                        tmp(i,j) += m_A[getIndex(i-ih+K,j-jh+L)]*h(ih,jh);
                    }
                }
            }
        }
    }
    return tmp;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::conv2(C_matrix<dataType>& h)
{
    if((h.getNbRow() % 2 == 0) || (h.getNbColumn() % 2 == 0))
    {
        std::cout << "even kernel size are not handled. Try to change the size of the convolution kernel to an odd size" << std::endl;
    }
    //returned matrix
    C_matrix<dataType> tmp(m_L,m_C);
    //std::cout << m_L << " " << m_C << std::endl;
    tmp = 0.0;


    //compute convolution
    unsigned long K = h.getNbRow()/2;
    unsigned long L = h.getNbColumn()/2;
    for(unsigned long i=0 ; i<m_L ; i++)
    {
        for(unsigned long j=0 ; j<m_C ; j++)
        {
            for(unsigned long ih=0 ; ih<2*K+1 ; ih++)
            {
                for(unsigned long jh=0 ; jh<2*L+1 ; jh++)
                {
                    if((i+K)>=ih && (j+L)>=jh && (i-ih+K)<m_L && (j-jh+L)<m_C) //check image boundaries
                    {
                        //std::cout << "doing the job " << i << " " << j << " " << ih << " " << jh << std::endl;
                        tmp(i,j) += m_A[getIndex(i-ih+K,j-jh+L)]*h(ih,jh);
                    }
                }
            }
        }
    }
    return tmp;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::gradX(void)
{
    //returned matrix
    C_matrix<dataType> tmp(m_L,m_C);
    tmp = 0.0;


    //compute convolution
    for(unsigned long i=0 ; i<m_L ; i++)
    {
        for(unsigned long j=1 ; j<m_C-1 ; j++)
        {
            tmp(i,j) = 0.5*(m_A[getIndex(i,j+1)]-m_A[getIndex(i,j-1)]);
        }
    }
    for(unsigned long i=0 ; i<m_L ; i++)
    {
        tmp(i,(unsigned long)0) = m_A[getIndex(i,1)]-m_A[getIndex(i,0)];
        tmp(i,(unsigned long)m_C-1) = m_A[getIndex(i,m_C-1)]-m_A[getIndex(i,m_C-2)];
    }
    return tmp;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::gradY(void)
{
    //returned matrix
    C_matrix<dataType> tmp(m_L,m_C);
    tmp = 0.0;


    //compute convolution
    for(unsigned long i=1 ; i<m_L-1 ; i++)
    {
        for(unsigned long j=0 ; j<m_C ; j++)
        {
            tmp(i,j) = 0.5*(m_A[getIndex(i+1,j)]-m_A[getIndex(i-1,j)]);
        }
    }
    for(unsigned long i=0 ; i<m_C ; i++)
    {
        tmp((unsigned long) 0,i) = m_A[getIndex(1,i)]-m_A[getIndex(0,1)];
        tmp((unsigned long) m_L-1,i) = m_A[getIndex(m_L-1,i)]-m_A[getIndex(m_L-2,i)];
    }
    return tmp;
}


template<class dataType> void C_matrix<dataType>::random(void)
{
    srand(time(NULL));
    for(unsigned long i=0 ; i<m_L ; i++)
    {
        for(unsigned long j=0 ; j<m_C ; j++)
        {
            m_A[getIndex(i,j)] = (dataType) rand();
        }
    }
    return;
}

template<class dataType> void C_matrix<dataType>::randomf(void)
{
    srand(time(NULL));
    for(unsigned long i=0 ; i<m_L ; i++)
    {
        for(unsigned long j=0 ; j<m_C ; j++)
        {
            m_A[getIndex(i,j)] = (dataType) (((double) rand())/((double) RAND_MAX));
        }
    }
    return;
}

template<class dataType> void C_matrix<dataType>::randomGauss(double mu, double sigma)
{
    double u,v;
    srand(time(NULL));
    for(unsigned long i=0 ; i<m_L ; i++)
    {
        for(unsigned long j=0 ; j<m_C ; j++)
        {
            u = (((double) rand())/((double) RAND_MAX));
            v = (((double) rand())/((double) RAND_MAX));
            m_A[getIndex(i,j)] = (dataType) mu + ( sigma * sqrt(-2*log(u)) ) * cos(2*Pi*v);
        }
    }
    return;
}

template<class dataType> void C_matrix<dataType>::randomGauss(std::vector<double> mu, std::vector<double> sigma)
{
    double u,v;
    srand(time(NULL));
    for(unsigned long i=0 ; i<m_L ; i++)
    {
        for(unsigned long j=0 ; j<m_C ; j++)
        {
            u = (((double) rand())/((double) RAND_MAX));
            v = (((double) rand())/((double) RAND_MAX));
            m_A[getIndex(i,j)] = (dataType) mu.at(j) + ( sigma.at(j) * sqrt(-2*log(u)) ) * cos(2*Pi*v);
        }
    }
    return;
}



template<class dataType> void C_matrix<dataType>::meshRow(dataType minValue, dataType maxValue)
{
    if(minValue>=maxValue)
    {
        std::cout << "min value is larger than max value" << std::endl;
        return;
    }
    double delta = (((double) maxValue) - ((double) minValue))/((double) (m_L-1));
    for(unsigned i=0 ; i<m_C ; i++)
    {
        for(unsigned j=0 ; j<m_L ; j++)
        {
            m_A[getIndex(j,i)] = (dataType) (((double) minValue) + ((double) j)*delta);
        }
    }
    return;
}

template<class dataType> void C_matrix<dataType>::meshColumn(dataType minValue, dataType maxValue)
{
    if(minValue>=maxValue)
    {
        std::cout << "min value is larger than max value" << std::endl;
        return;
    }
    double delta = (((double) maxValue) - ((double) minValue))/((double) (m_C-1));
    for(unsigned i=0 ; i<m_C ; i++)
    {
        for(unsigned j=0 ; j<m_L ; j++)
        {
            m_A[getIndex(j,i)] = (dataType) (((double) minValue) + ((double) i)*delta);
        }
    }
    return;
}



















//**********************
//Décomposition LU
//**********************
template<class dataType> C_matrix<dataType> C_matrix<dataType>::LU(void)
{
    C_matrix<dataType> MLU(*this);// = new C_matrix<dataType>(*this);
    C_matrix<dataType> tmp(m_L,m_C);


    std::cout<<"\nLU Decomposition Steps\n";
    std::cout<<"------------------------\n";

    if(m_L != m_C)
    {
        std::cout<<"\nColumn and Rows must be equal"<<std::endl;
    }
    else
    {
        for(int j=0; j<m_C; j++)
        {
            if(ABS(MLU(j,j))<epsilon)
            {
                std::cout<<"Pivot is equal to zero..."<<std::endl;
                break; //pivot nul
            }
            else
            {

                MLU.show();
                //Division colonne j par pivot
                for(int i=j+1; i<m_L; i++)
                {
                    MLU(i,j) = MLU(i,j)/MLU(j,j);
                }

                //Complément de Shur
                for(int jj=j+1; jj<m_C; jj++)
                {
                    for(int ii=j+1; ii<m_L; ii++)
                    {
                        tmp(ii,jj) = MLU(j,jj)*MLU(ii,j);
                    }
                }
                for(int jj=j+1; jj<m_C; jj++)
                {
                    for(int ii=j+1; ii<m_L; ii++)
                    {
                        MLU(ii,jj) = MLU(ii,jj)-tmp(ii,jj);
                    }
                }

            }
        }
    }
    std::cout<<"\t\tEnd LU Decomposition\n";
    //MLU.Affiche();
    return MLU;
}












//******************************
//Décomposition LU + Permutation
//******************************
template<class dataType> C_matrix<dataType> C_matrix<dataType>::LUP(C_vector<int> &Indice)
{
    C_matrix<dataType> MLUP(*this);
    C_matrix<dataType> tmp(m_L,m_C);
    dataType Tmp;

    std::cout<<"\nLUP Decomposition Steps\n";
    std::cout<<"------------------------\n";

    if(m_L != m_C){
        std::cout<<"\nColumn and Rows must be equal"<<std::endl;
    }
    else
    {

        //Indice contiendra les permutations
        for(int j=0; j<m_C; j++)
            Indice.set(j,(dataType) j);

        //Pour tous les pivots
        for(int j=0; j<m_C; j++)
        {
            if(ABS(MLUP(j,j))<epsilon)
            {
                std::cout<<"Pivot is equal to zero..."<<std::endl;
                break; //pivot nul
            }
            else
            {
                MLUP.show();

                //Recherche du pivot le plus grand
                dataType max=MLUP(j,j);
                int indmax=j;
                for(int i=j+1; i<m_L; i++)
                {
                    if(MLUP(i,j)>max)
                    {
                        max=MLUP(i,j);
                        indmax=i;
                    }
                }

                //Echange des lignes
                for(int jj=0;jj<m_C;jj++)
                {
                    Tmp = MLUP(j,jj);
                    MLUP(j,jj) = MLUP(indmax,jj);
                    MLUP(indmax,jj) = Tmp;
                }
                //Mise à jour du tableau Indice (=les perumtations)
                Indice.set(j,indmax);
                Indice.set(indmax,j);

                //Division colonne j par pivot
                for(int i=j+1; i<m_L; i++)
                {
                    MLUP(i,j) = MLUP(i,j)/MLUP(j,j);
                }

                //Complément de Shur
                for(int jj=j+1; jj<m_C; jj++)
                {
                    for(int ii=j+1; ii<m_L; ii++)
                    {
                        tmp(ii,jj) = MLUP(j,jj)*MLUP(ii,j);
                    }
                }
                for(int jj=j+1; jj<m_C; jj++)
                {
                    for(int ii=j+1; ii<m_L; ii++)
                    {
                        MLUP(ii,jj) = MLUP(ii,jj)-tmp(ii,jj);
                    }
                }

            }
        }
    }

    std::cout<<"\t\tEnd LUP Decomposition\n";
    //MLUP.Affiche();
    return MLUP;
}







//**********************
//Produit L*U
//**********************
template<class dataType> C_matrix<dataType> C_matrix<dataType>::LMU(void)
{
    C_matrix<dataType> L(m_L,m_C);
    C_matrix<dataType> U(m_L,m_C);
    C_matrix<dataType> LMU(m_L,m_C);
    L = 0.0;
    U = 0.0;
    LMU = 0.0;


    //Matrice L
    for(int i=0;i<m_L;i++)
    {
        L(i,i)=1.0;
        for (int j=0;j<i; j++)
        {
            L(i,j) = (*this)(i,j);
        }
    }

    //Matrice L
    for(int i=0;i<m_L;i++)
    {
        for (int j=i;j<m_C; j++)
        {
            U(i,j) = (*this)(i,j);
        }
    }

    L.show();
    U.show();

    LMU = L * U;
    return LMU;

}


//**********************************
//Résolution du système MX=B avec LU
//**********************************
template<class dataType> C_vector<dataType> C_matrix<dataType>::LineAlgEq_LU(C_vector<dataType> &B)
{

    C_matrix<dataType> MLU = LU();
    C_vector<dataType> X(m_C);
    C_vector<dataType> Y(m_C);


    //Solve for LY=B
    for(int i=0;i<m_L;i++)
    {
        Y.set(i,B[i]);
        for(unsigned short j=0;j<i;j++)
        {
            Y.set(i,Y[i] - MLU((unsigned short)i,j)*Y[j]);
        }
    }

    //Solve for UX=Y
    for(int i=m_L-1;i>=0;i--)
    {
        X.set(i, Y[i]);
        for(int j=m_L-1;j>i;j--)
        {
            X.set(i,X[i] - MLU(i,j)*X[j]);
        }
        X.set(i,X[i]/MLU(i,i));
    }
    return X;
}
/**/

//http://fr.wikipedia.org/wiki/%C3%89limination_de_Gauss-Jordan
template<class dataType> C_matrix<dataType> C_matrix<dataType>::inv(void)
{
    if(m_L!=m_C)
    {
        throw "matrix must be square";
    }
    //std::cout << "use me carefully, I will return a matrix even though the matrix is not invertible" << std::endl;

    //copy the current matrix push
    C_matrix<dataType> tmp(*this);

    //
    long r = -1; //row index

    C_matrix<dataType> I(m_L,m_C);
    I= (dataType) 0.0;
    for(long i=0 ; i<m_L ; i++)
        I(i,i) = (dataType) 1.0;

    for(long j=0 ; j<m_C ; j++) //go through all column
    {
        //search max abs value in column j starting at index r+1
        dataType piv = -1.0;
        long k = 0;
        for(long i=r+1 ; i<m_L ; i++)
        {
            if(piv<ABS(tmp(i,j)))
            {
                piv = tmp(i,j);
                k = i;
            }
        }

        //if piv is not null
        if(piv>SMALL_NUM_F)
        {
            //divide row k by piv
            for(long i=0 ; i<m_C ; i++)
            {
                tmp(k,i) = tmp(k,i)/piv;
                I(k,i) = I(k,i)/piv;
            }

            //
            r++;

            //swapp rows k and r
            for(long i=0 ; i<m_C ; i++)
            {
                SWAP(tmp(k,i),tmp(r,i));
                SWAP(I(k,i),I(r,i));
            }

            for(long i=0 ; i<m_C ; i++)
            {
                if(i!=r)
                {
                    //compute Li <- Li - m_A[i][j]*Lr
                    dataType tmpD = tmp(i,j);
                    for(long m=0 ; m<m_C ; m++)
                    {
                        tmp(i,m) = tmp(i,m) - tmpD*tmp(r,m);
                        I(i,m) = I(i,m) - tmpD*I(r,m);
                    }
                }
            }
        }
    }
    return I;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::pseudoInv(void)
{
    C_matrix<dataType> tmp = ((*this).Transpose())*(*this);
    return (tmp.inv())*((*this).Transpose());
}



template<class dataType> std::vector<C_matrix<dataType> > C_matrix<dataType>::svd(void)
{
    long int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;

    long int m = getNbRow();
    long int n = getNbColumn();

    C_matrix<dataType> a(*this);
    C_matrix<dataType> v(n,n);
    C_matrix<dataType> w(n,n);
    v=0.0;
    w=0.0;
    std::vector<C_matrix<dataType> > svdDecomp;

    if (m < n)
    {
        throw "#rows must be > #cols";
        return svdDecomp;
    }

    rv1 = new double[n];

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += ABS((double)a(k,i));
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    a(k,i) = (double)((double)a(k,i)/scale);
                    s += ((double)a(k,i) * (double)a(k,i));
                }
                f = (double)a(i,i);
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a(i,i) = (double)(f - g);
                if (i != n - 1)
                {
                    for (j = l; j < n; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += ((double)a(k,i) * (double)a(k,j));
                        f = s / h;
                        for (k = i; k < m; k++)
                            a(k,j) += (double)(f * (double)a(k,i));
                    }
                }
                for (k = i; k < m; k++)
                    a(k,i) = (double)((double)a(k,i)*scale);
            }
        }
        w(i,i) = (double)(scale * g);

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1)
        {
            for (k = l; k < n; k++)
                scale += ABS((double)a(i,k));
            if (scale)
            {
                for (k = l; k < n; k++)
                {
                    a(i,k) = (double)((double)a(i,k)/scale);
                    s += ((double)a(i,k) * (double)a(i,k));
                }
                f = (double)a(i,l);
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a(i,l) = (double)(f - g);
                for (k = l; k < n; k++)
                    rv1[k] = (double)a(i,k) / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < n; k++)
                            s += ((double)a(j,k) * (double)a(i,k));
                        for (k = l; k < n; k++)
                            a(j,k) += (double)(s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++)
                    a(i,k) = (double)((double)a(i,k)*scale);
            }
        }
        anorm = MAX(anorm, (ABS((double)w(i,i)) + ABS(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--)
    {
        if (i < n - 1)
        {
            if (g)
            {
                for (j = l; j < n; j++)
                    v(j,i) = (double)(((double)a(i,j) / (double)a(i,l)) / g);
                    /* double division to avoid underflow */
                for (j = l; j < n; j++)
                {
                    for (s = 0.0, k = l; k < n; k++)
                        s += ((double)a(i,k) * (double)v(k,j));
                    for (k = l; k < n; k++)
                        v(k,j) += (double)(s * (double)v(k,i));
                }
            }
            for (j = l; j < n; j++)
                v(i,j) = v(j,i) = 0.0;
        }
        v(i,i) = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--)
    {
        l = i + 1;
        g = (double)w(i,i);
        if (i < n - 1)
            for (j = l; j < n; j++)
                a(i,j) = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != n - 1)
            {
                for (j = l; j < n; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += ((double)a(k,i) * (double)a(k,j));
                    f = (s / (double)a(i,i)) * g;
                    for (k = i; k < m; k++)
                        a(k,j) += (double)(f * (double)a(k,i));
                }
            }
            for (j = i; j < m; j++)
                a(j,i) = (double)((double)a(j,i)*g);
        }
        else
        {
            for (j = i; j < m; j++)
                a(j,i) = 0.0;
        }
        ++(a(i,i));
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--)
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++)
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (ABS(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (ABS((double)w(nm,nm)) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (ABS(f) + anorm != anorm)
                    {
                        g = (double)w(i,i);
                        h = pythag(f, g);
                        w(i,i) = (double)h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = (double)a(j,nm);
                            z = (double)a(j,i);
                            a(j,nm) = (double)(y * c + z * s);
                            a(j,i) = (double)(z * c - y * s);
                        }
                    }
                }
            }
            z = (double)w(k,k);
            if (l == k)
            {                  /* convergence */
                if (z < 0.0)
                {              /* make singular value nonnegative */
                    w(k,k) = (double)(-z);
                    for (j = 0; j < n; j++)
                        v(j,k) = (-v(j,k));
                }
                break;
            }
            if (its >= 50) {
                delete rv1;
                throw "No convergence after 30,000! iterations";
                return svdDecomp;
            }

            /* shift from bottom 2 x 2 minor */
            x = (double)w(l,l);
            nm = k - 1;
            y = (double)w(nm,nm);
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = (double)w(i,i);
                h = s * g;
                g = c * g;
                z = pythag(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++)
                {
                    x = (double)v(jj,j);
                    z = (double)v(jj,i);
                    v(jj,j) = (double)(x * c + z * s);
                    v(jj,i) = (double)(z * c - x * s);
                }
                z = pythag(f, h);
                w(j,j) = (double)z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = (double)a(jj,j);
                    z = (double)a(jj,i);
                    a(jj,j) = (double)(y * c + z * s);
                    a(jj,i) = (double)(z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w(k,k) = (double)x;
        }
    }
    //free((void*) rv1);
    delete rv1;

    svdDecomp.push_back(a);
    svdDecomp.push_back(w);
    svdDecomp.push_back(v);
    return svdDecomp;
}


template<class dataType> double C_matrix<dataType>::pythag(const double a, const double b)
{
        double absa=ABS(a), absb=ABS(b);
        return (absa > absb ? absa*sqrt(1.0+SQR(absb/absa)) :
                (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb))));
}













template<class dataType> std::vector<C_matrix<dataType> > C_matrix<dataType>::eigSym(/**const double** a, const int nb_rows,*/ bool yesvec)
{
    if(getNbColumn()!=getNbRow()) throw "matrix dimensions : square matrix expected";

    std::vector<C_matrix<dataType> > eigDecomp;
    int n = getNbRow();
    C_matrix<dataType> z(*this);//each column correspond to an eigenvector

    C_matrix<dataType> d(n,n);// = new double[n]; //contain the eigenvalues sorted descendingly
    C_matrix<dataType> e(n,n);//e = new double[n];
    tred2(&z,&d,&e,yesvec);
    tqli(&z,&d,&e,yesvec);
    sortEig(&d,&z,yesvec);

    eigDecomp.push_back(z);
    eigDecomp.push_back(d);
    eigDecomp.push_back(e);
    return eigDecomp;
}

template<class dataType> void C_matrix<dataType>::tred2(C_matrix<dataType> *z, C_matrix<dataType> *d, C_matrix<dataType> *e, bool yesvecs)
{
    int l,k,j,i;
    double scale,hh,h,g,f;
    int n = getNbRow();
    for (i=n-1;i>0;i--) {
        l=i-1;
        h=scale=0.0;
        if (l > 0) {
            for (k=0;k<i;k++)
                                scale += ABS((*z)(i,k));

            if (scale == 0.0)
                (*e)(i,i)=(*z)(i,l);
            else {
                for (k=0;k<i;k++) {
                    (*z)(i,k) /= scale;
                    h += (*z)(i,k)*(*z)(i,k);
                }
                f=(*z)(i,l);
                g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
                (*e)(i,i)=scale*g;
                h -= f*g;
                (*z)(i,l)=f-g;
                f=0.0;
                for (j=0;j<i;j++) {
                    if (yesvecs)
                        (*z)(j,i)=(*z)(i,j)/h;
                    g=0.0;
                    for (k=0;k<j+1;k++)
                        g += (*z)(j,k)*(*z)(i,k);
                    for (k=j+1;k<i;k++)
                        g += (*z)(k,j)*(*z)(i,k);
                    (*e)(j,j)=g/h;
                    f += (*e)(j,j)*(*z)(i,j);
                }
                hh=f/(h+h);
                for (j=0;j<i;j++) {
                    f=(*z)(i,j);
                    (*e)(j,j)=g=(*e)(j,j)-hh*f;
                    for (k=0;k<j+1;k++)
                        (*z)(j,k) -= (f*(*e)(k,k)+g*(*z)(i,k));
                }
            }
        } else
        {
            (*e)(i,i)=(*z)(i,l);
        }
        (*d)(i,i)=h;
    }
    if (yesvecs) (*d)(0,0)=0.0;
    (*e)(0,0)=0.0;
    for (i=0;i<n;i++) {
        if (yesvecs) {
            if ((*d)(i,i) != 0.0) {
                for (j=0;j<i;j++) {
                    g=0.0;
                    for (k=0;k<i;k++)
                        g += (*z)(i,k)*(*z)(k,j);
                    for (k=0;k<i;k++)
                        (*z)(k,j) -= g*(*z)(k,i);
                }
            }
            (*d)(i,i)=(*z)(i,i);
            (*z)(i,i)=1.0;
            for (j=0;j<i;j++) (*z)(j,i)=(*z)(i,j)=0.0;
        } else {
            (*d)(i,i)=(*z)(i,i);
        }
    }
}



template<class dataType> void C_matrix<dataType>::tqli(C_matrix<dataType> *z, C_matrix<dataType> *d, C_matrix<dataType> *e, bool yesvecs)
{
    int m,l,iter,i,k;
    double s,r,p,g,f,dd,c,b;
    const double EPS=DBL_EPSILON;
    int n = getNbRow();

    for (i=1;i<n;i++) (*e)(i-1,i-1)=(*e)(i,i);
    (*e)(n-1,n-1)=0.0;
    for (l=0;l<n;l++) {
        iter=0;
        do {
            for (m=l;m<n-1;m++) {
                                dd=ABS((*d)(m,m))+ABS((*d)(m+1,m+1));
                                if (ABS((*e)(m,m)) <= EPS*dd) break;
            }
            if (m != l) {
                if (iter++ == 90) throw("Too many iterations in tqli"); //used to be 30
                g=((*d)(l+1,l+1)-(*d)(l,l))/(2.0*(*e)(l,l));
                r=pythag(g,1.0);
                g=(*d)(m,m)-(*d)(l,l)+(*e)(l,l)/(g+SIGN(r,g));
                s=c=1.0;
                p=0.0;
                for (i=m-1;i>=l;i--) {
                    f=s*(*e)(i,i);
                    b=c*(*e)(i,i);
                    (*e)(i+1,i+1)=(r=pythag(f,g));
                    if (r == 0.0) {
                        (*d)(i+1,i+1) -= p;
                        (*e)(m,m)=0.0;
                        break;
                    }
                    s=f/r;
                    c=g/r;
                    g=(*d)(i+1,i+1)-p;
                    r=((*d)(i,i)-g)*s+2.0*c*b;
                    (*d)(i+1,i+1)=g+(p=s*r);
                    g=c*r-b;
                    if (yesvecs) {
                        for (k=0;k<n;k++) {
                            f=(*z)(k,i+1);
                            (*z)(k,i+1)=s*(*z)(k,i)+c*f;
                            (*z)(k,i)=c*(*z)(k,i)-s*f;
                        }
                    }
                }
                if (r == 0.0 && i >= l) continue;
                (*d)(l,l) -= p;
                (*e)(l,l)=g;
                (*e)(m,m)=0.0;
            }
        } while (m != l);
    }
}




template<class dataType> void C_matrix<dataType>::sortEig(C_matrix<dataType> *d, C_matrix<dataType> *z, bool yesvecs)
{ //sort descending order
    if (yesvecs)
        eigsrt(d,z);
    else
        eigsrt(d,NULL);
}
template<class dataType> void C_matrix<dataType>::eigsrt(C_matrix<dataType> *d, C_matrix<dataType> *z)
{
    int k;
    int n=getNbColumn();

    for (int i=0;i<n-1;i++) {
        k=i;
        double p=(*d)(k,k);
        for (int j=i;j<n;j++)
            if ((*d)(j,j) >= p)
            {
                k=j;
                p=(*d)(k,k);
            }
        if (k != i) {
            (*d)(k,k)=(*d)(i,i);
            (*d)(i,i)=p;
            if (z != NULL)
                for (int j=0;j<n;j++) {
                    p=(*z)(j,i);
                    (*z)(j,i)=(*z)(j,k);
                    (*z)(j,k)=p;
                }
        }
    }
}
#endif // C_MATRIX_H
