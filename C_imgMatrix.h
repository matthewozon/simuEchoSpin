#ifndef C_IMGMATRIX_H
#define C_IMGMATRIX_H

#include <C_matrix.h> //the class inherit all method from this one

#include <CImg.h> //load and display png

#define NO_NORMALIZATION 0 //it is assumed that data to be displayed are in the range [0,255]
#define NORMALIZE 1 //the displayed values are normalized in the range [0,255]

template<class dataType> class C_imgMatrix : public C_matrix<dataType>
{
public:
    C_imgMatrix();
    C_imgMatrix(std::string fileNamePNG);
    C_imgMatrix(int L, int C);
    C_imgMatrix(const C_imgMatrix &M);
    C_imgMatrix(const C_matrix<dataType> &M);
    virtual ~C_imgMatrix();


    //overload operators (I guess I was wrong, those operator are not really inherited)
    C_imgMatrix<dataType> operator= (C_imgMatrix const& c)
    {
        C_matrix<dataType>::operator =(c);
        return (*this);
    }
    C_imgMatrix<dataType> operator= (dataType const& x)
    {
        C_matrix<dataType>::operator =(x);
        return (*this);
    }
    C_imgMatrix<dataType> operator+ (C_imgMatrix const& c)
    {
        C_imgMatrix<dataType> B = C_matrix<dataType>::operator +(c);
        return B;
    }
    C_imgMatrix<dataType> operator+ (dataType const& x)
    {
        C_imgMatrix<dataType> B = C_matrix<dataType>::operator +(x);
        return B;
    }
    C_imgMatrix<dataType> operator- (C_imgMatrix const& c)
    {
        C_imgMatrix<dataType> B = C_matrix<dataType>::operator -(c);
        return B;
    }
    C_imgMatrix<dataType> operator- (dataType const& x)
    {
        C_imgMatrix<dataType> B = C_matrix<dataType>::operator -(x);
        return B;
    }
    C_imgMatrix<dataType> operator* (C_imgMatrix const& c)
    {
        C_imgMatrix<dataType> B = C_matrix<dataType>::operator *(c);
        return B;
    }
    //C_vector<dataType>& operator* (C_vector<dataType> const& c); //obselete
    C_imgMatrix operator* (const dataType& x)
    {
        C_imgMatrix<dataType> B = C_matrix<dataType>::operator *(x);
        return B;
    }

    C_imgMatrix<dataType> operator== (C_imgMatrix const& c)
    {
        C_imgMatrix<dataType> B = C_matrix<dataType>::operator ==(c);
        return B;
    }
    C_imgMatrix<dataType> operator== (dataType const& x)
    {
        C_imgMatrix<dataType> B = C_matrix<dataType>::operator ==(x);
        return B;
    }

    C_imgMatrix<dataType> operator> (C_imgMatrix const& c)
    {
        C_imgMatrix<dataType> B = C_matrix<dataType>::operator >(c);
        return B;
    }
    C_imgMatrix<dataType> operator> (dataType const& x)
    {
        C_imgMatrix<dataType> B = C_matrix<dataType>::operator >(x);
        return B;
    }
    C_imgMatrix<dataType> operator>= (C_imgMatrix const& c)
    {
        C_imgMatrix<dataType> B = C_matrix<dataType>::operator >=(c);
        return B;
    }
    C_imgMatrix<dataType> operator>= (dataType const& x)
    {
        C_imgMatrix<dataType> B = C_matrix<dataType>::operator >=(x);
        return B;
    }

    C_imgMatrix<dataType> operator< (C_imgMatrix const& c)
    {
        C_imgMatrix<dataType> B = C_matrix<dataType>::operator <(c);
        return B;
    }
    C_imgMatrix<dataType> operator< (dataType const& x)
    {
        C_imgMatrix<dataType> B = C_matrix<dataType>::operator <(x);
        return B;
    }
    C_imgMatrix<dataType> operator<= (C_imgMatrix const& c)
    {
        C_imgMatrix<dataType> B = C_matrix<dataType>::operator <=(c);
        return B;
    }
    C_imgMatrix<dataType> operator<= (dataType const& x)
    {
        C_imgMatrix<dataType> B = C_matrix<dataType>::operator <=(x);
        return B;
    }

    void load(std::string fileNamePNG);
    void savepng(std::string fileNamePNG);
    void savetiff(std::string fileNameTIFF);

    void display(unsigned short NORMALIZATION);
    void displayHisto(int nbBins);
};

template<class dataType> C_imgMatrix<dataType>::C_imgMatrix():C_matrix<dataType>()
{
    //
}

template<class dataType> C_imgMatrix<dataType>::C_imgMatrix(std::string fileNamePNG):C_matrix<dataType>(0,0)
{
    //load image
    cimg_library::CImg<float> image(fileNamePNG.data());

    //resize container
    this->resize(image.width(), image.height());

    //store data in pointer m_A
    for(unsigned short l=0 ; l<this->m_L ; l++)
    {
        //this->m_A[l] = &((imgData->data())[l*C]);
        for(unsigned short c=0 ; c<this->m_C ; c++)
        {
            this->m_A[l][c] = (dataType) image(l,c);//,0,0
        }
    }
}


template<class dataType> C_imgMatrix<dataType>::C_imgMatrix(int L, int C):C_matrix<dataType>(L,C)
{
    //
}

template<class dataType> C_imgMatrix<dataType>::C_imgMatrix(const C_imgMatrix &M):C_matrix<dataType>(M)
{
    //offsetL = M.offsetL;
    //offsetC = M.offsetC;
    //m_fileNamePNG = M.m_fileNamePNG;
}
template<class dataType> C_imgMatrix<dataType>::C_imgMatrix(const C_matrix<dataType> &M):C_matrix<dataType>(M)
{
    //
    //offsetL = 0;
    //offsetC = 0;
}

template<class dataType> C_imgMatrix<dataType>::~C_imgMatrix()
{
}


template<class dataType> void C_imgMatrix<dataType>::display(unsigned short NORMALIZATION)
{
    cimg_library::CImg<dataType> image(this->m_L,this->m_C);
    if(image.spectrum()==3) image.RGBtoBayer();
    ///fill in object
    for(unsigned short l=0 ; l<this->m_L ; l++)
    {
        for(unsigned short c=0 ; c<this->m_C ; c++)
        {
            image(l,c) = (*this)(l,c);
        }
    }

    //add lines
    //const dataType red[] = { 255,0,0 }, green[] = { 0,255,0 }, blue[] = { 0,0,255 };
    //double x0=0.0;
    //double y0=0.0;
    //double x1=10.0;
    //double y1=10.0;
    //image.draw_line(x0,y0,x1,y1,red,1);
    //image.draw_arrow(10,10,25,20,red);

    cimg_library::CImgDisplay main_disp(image,"My image",1);
    while (!main_disp.is_closed())
    {
        main_disp.wait();
        //main_disp.resize(true);
    }
}

template<class dataType> void C_imgMatrix<dataType>::displayHisto(int nbBins)
{
    cimg_library::CImg<dataType> image(this->m_L,this->m_C);
    cimg_library::CImgDisplay main_disp;
    for(unsigned short l=0 ; l<this->m_L ; l++)
    {
        for(unsigned short c=0 ; c<this->m_C ; c++)
        {
            image(l,c) = (*this)(l,c);
        }
    }
    const cimg_library::CImg< dataType > histogram = image.histogram( nbBins );

    cimg_library::CImg<unsigned char> visu(800,600,1,3,0);
    const unsigned char blue[] = { 0,0,255 };
    cimg_library::CImgDisplay draw_disp(visu,"histogram");
    double maxV = (double) histogram(0,0);
    for(int i=1 ; i<nbBins ; i++)
    {
	if( ((double) histogram(i,0))>maxV) maxV = (double) histogram(i,0);
    }
    visu.draw_graph(histogram,blue,1,1,0.0,1.1*maxV,0).display(draw_disp);
    std::cin.ignore();
    return;
}


template<class dataType> void C_imgMatrix<dataType>::load(std::string fileNamePNG)
{
    //load image
    cimg_library::CImg<float> image(fileNamePNG.data());

    //resize container
    this->resize(image.width(), image.height());

    //store data in pointer m_A
    for(unsigned short l=0 ; l<this->m_L ; l++)
    {
        //this->m_A[l] = &((imgData->data())[l*C]);
        for(unsigned short c=0 ; c<this->m_C ; c++)
        {
            this->m_A[l][c] = image(l,c);//,0,0
        }
    }
}

template<class dataType> void C_imgMatrix<dataType>::savepng(std::string fileNamePNG)
{
    //load image
    cimg_library::CImg<float> image(this->m_L,this->m_C);

    //store data in pointer m_A
    for(unsigned short l=0 ; l<this->m_L ; l++)
    {
        //this->m_A[l] = &((imgData->data())[l*C]);
        for(unsigned short c=0 ; c<this->m_C ; c++)
        {
            image(l,c) = (float) this->m_A[l][c];//,0,0
        }
    }
    image.save_png(fileNamePNG.data());
}

template<class dataType> void C_imgMatrix<dataType>::savetiff(std::string fileNameTIFF)
{
    //load image
    cimg_library::CImg<float> image(this->m_L,this->m_C);

    //store data in pointer m_A
    for(unsigned short l=0 ; l<this->m_L ; l++)
    {
        //this->m_A[l] = &((imgData->data())[l*C]);
        for(unsigned short c=0 ; c<this->m_C ; c++)
        {
            image(l,c) = (float) this->m_A[l][c];//,0,0
        }
    }
    image.save_tiff(fileNameTIFF.data());
}

#endif // C_IMGMATRIX_H
