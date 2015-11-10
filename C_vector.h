//good to go!

#ifndef C_VECTOR_H
#define C_VECTOR_H

#include <stdlib.h>
#include <iostream>

//this is a template, meaning that you need to specify the type when you instanciate a object.
template<class dataType> class C_vector
{
public:
    C_vector(unsigned short _N); //constructor
    virtual ~C_vector(); //destructor
    C_vector(C_vector const& X); // constructeur par copie

    //overload usefull operators
    C_vector& operator= (C_vector const& c);
    C_vector& operator= (dataType const& x);
    C_vector& operator+ (C_vector const& c);
    C_vector& operator+ (dataType const& x);
    C_vector& operator- (C_vector const& c);
    C_vector& operator- ( dataType const& x);
    dataType operator[](unsigned short i)const; //acces array
    bool set (unsigned short i,  dataType x);

    //tool that you could use
    unsigned short length() const; // get length
    void show();
private:
    //size and datas
    unsigned short m_N;
    dataType *m_V;

};


template<class dataType> C_vector<dataType>::C_vector(unsigned short _N):m_N(_N)
{
    /*cout << "create Vector" << endl;
    m_V = new dataType[m_N];
    for(int i=0; i<m_N; i++)
            m_V[i]=0;*/

    //std::cout << "create vector" << std::endl;
    m_V = NULL;
    m_V = new dataType[_N];
    if(m_V==NULL)
    {
        m_N=0;
    }
    else
    {
        m_N=_N;
    }

}

template<class dataType> C_vector<dataType>::C_vector(C_vector const& X)
{
    m_V = NULL;
    m_V = new dataType[X.length()];
    if(m_V==NULL)
    {
        m_N=0;
    }
    else
    {
        m_N=X.length();
        for(unsigned short i=0 ; i<m_N ; i++)
            this->m_V[i] = X[i];
    }

}

template<class dataType> C_vector<dataType>::~C_vector()
{
    //std::cout << "delete vector" << std::endl;
    if(m_V!=NULL)
    {
        delete m_V;
    }
}

template<class dataType> void C_vector<dataType>::show()
{
    for(unsigned short i=0 ; i<m_N ; i++)
        std::cout << "V[" << i << "] = " << (*this)[i] << std::endl;
}

template<class dataType> unsigned short C_vector<dataType>::length() const
{
    return m_N;
}

template<class dataType> dataType C_vector<dataType>::operator[] (unsigned short i) const
{
    if(i<m_N)
    {
        return m_V[i];
    }
    else
    {
        throw "index must be positive and smaller than max length of vector";
    }

}
/*template<class dataType> C_vector<dataType>& C_vector<dataType>::operator= (C_vector const& c)
{
        if(this==&c)
                return *this;

        delete m_V; // nouvelle taille
        m_N = c.length();
        m_V = new dataType[m_N];
        for(int i=0; i<m_N; i++)
                m_V[i]=c[i];
        return *this;
}*/
template<class dataType> C_vector<dataType>& C_vector<dataType>::operator= (C_vector const& c)
{
    if(this==&c) return *this;
    if(c.length()!=this->length()) throw "dimension vector must agree";

    for(unsigned short i=0 ; i<this->m_N ; i++)
        this->m_V[i] = c[i];

    return *this;
}

template<class dataType> C_vector<dataType>& C_vector<dataType>::operator= (dataType const& x)
{
    for(unsigned short i=0 ; i<this->m_N ; i++)
        this->m_V[i] = x;

    return *this;
}

template<class dataType> C_vector<dataType>& C_vector<dataType>::operator+ (const C_vector& c)
{
    if(c.length()!=this->length()) throw "dimension vector must agree";

    C_vector<dataType>* Y = new C_vector<dataType>(this->m_N);
    for(unsigned short i=0 ; i<this->m_N ; i++)
        Y->set(i,this->m_V[i] + c[i]);

    return *Y;
}

template<class dataType> C_vector<dataType>& C_vector<dataType>::operator+ (const dataType& x)
{
    C_vector<dataType>* Y = new C_vector<dataType>(this->m_N);
    for(unsigned short i=0 ; i<this->m_N ; i++)
        Y->set(i,this->m_V[i] + x);

    return *Y;
}

template<class dataType> C_vector<dataType>& C_vector<dataType>::operator- (C_vector const& c)
{
    if(c.length()!=this->length()) throw "dimension vector must agree";

    C_vector<dataType>* Y = new C_vector<dataType>(this->m_N);
    for(unsigned short i=0 ; i<this->m_N ; i++)
        Y->set(i,this->m_V[i] - c[i]);

    return *Y;
}

template<class dataType> C_vector<dataType>& C_vector<dataType>::operator- (dataType const& x)
{
    C_vector<dataType>* Y = new C_vector<dataType>(this->m_N);
    for(unsigned short i=0 ; i<this->m_N ; i++)
        Y->set(i,this->m_V[i] - x);

    return *Y;
}

template<class dataType> bool C_vector<dataType>::set (unsigned short i,  dataType x)
{
    if(i>=m_N) return false;
    this->m_V[i] = x;
    return true;
}

#endif // C_VECTOR_H
