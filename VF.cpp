#include "VF.h"

using namespace std;

Probleme::Probleme(int NbLignes, int NbCol, int choix)
{
  _choix=choix;
  _NbCol=NbCol;
  _NbLignes=NbLignes;
  _U.resize(8);
  _Uapres.resize(8);
  _f.resize(8);
  _g.resize(8);
  _ff.resize(8);
  _gg.resize(8);
  _UE.resize(8);
  _UN.resize(8);
  _US.resize(8);
  _UW.resize(8);
  _Uy.resize(8);
  _Ux.resize(8);
  _f1.resize(8);
  _f2.resize(8);
  _g1.resize(8);
  _g2.resize(8);
  _F2.resize(8);
  _G2.resize(8);
  for (int j=0; j<8; j++)
  {
    _U[j].resize(_NbLignes*_NbCol);
    _Uapres[j].resize(_NbLignes*_NbCol);
    _f[j].resize(_NbLignes*_NbCol);
    _g[j].resize(_NbLignes*_NbCol);
    _ff[j].resize(_NbLignes*_NbCol);
    _gg[j].resize(_NbLignes*_NbCol);
    _UN[j].resize(_NbLignes*_NbCol);
    _US[j].resize(_NbLignes*_NbCol);
    _UE[j].resize(_NbLignes*_NbCol);
    _UW[j].resize(_NbLignes*_NbCol);
    _Ux[j].resize(_NbLignes*_NbCol);
    _Uy[j].resize(_NbLignes*_NbCol);
    _f1[j].resize(_NbLignes*_NbCol);
    _f2[j].resize(_NbLignes*_NbCol);
    _g1[j].resize(_NbLignes*_NbCol);
    _g2[j].resize(_NbLignes*_NbCol);
    _F2[j].resize(_NbLignes*_NbCol);
    _G2[j].resize(_NbLignes*_NbCol);

  }

  _b.resize(3);
  for (int j=0; j<3; j++)
  {
    _b[j].resize(_NbLignes*_NbCol);
  }

  _a2.resize(_NbLignes*_NbCol);
  _c1car.resize(_NbLignes*_NbCol);
  _c2car.resize(_NbLignes*_NbCol);
  _alpha.resize(_NbLignes*_NbCol);
  _beta.resize(_NbLignes*_NbCol);
  _p.resize(_NbLignes*_NbCol);
  _ptilde.resize(_NbLignes*_NbCol);
  _pscal.resize(_NbLignes*_NbCol);

  _Lx=2*acos(-1);
  _Ly=2*acos(-1);
  //_Lx=1.0;
  //_Ly=1.0;
  _Dx=_Lx/(_NbCol);
  _Dy=_Ly/(_NbLignes);
  _gamma=1.4;
  //_tmax=0.125;
  _tmax=acos(-1);


  for (int j=0; j<_NbLignes; j++)
  {
    for (int i=0; i<_NbCol; i++)
    {
      for(int indice=0; indice<8; indice++)
      {
        _Uapres[indice][j*_NbCol+i]=0.;
        _U[indice][j*_NbCol+i]=0.;
        _f[indice][j*_NbCol+i]=0.;
        _g[indice][j*_NbCol+i]=0.;
        _ff[indice][j*_NbCol+i]=0.;
        _gg[indice][j*_NbCol+i]=0.;
        _UN[indice][j*_NbCol+i]=0.;
        _US[indice][j*_NbCol+i]=0.;
        _UE[indice][j*_NbCol+i]=0.;
        _UW[indice][j*_NbCol+i]=0.;
        _Ux[indice][j*_NbCol+i]=0.;
        _Uy[indice][j*_NbCol+i]=0.;
        _f1[indice][j*_NbCol+i]=0.;
        _f2[indice][j*_NbCol+i]=0.;
        _g1[indice][j*_NbCol+i]=0.;
        _g2[indice][j*_NbCol+i]=0.;
        _F2[indice][j*_NbCol+i]=0.;
        _G2[indice][j*_NbCol+i]=0.;
      }

      for (int k=0; k<3; k++)
      {
        _b[k][j*_NbCol+i]=0.;
      }

      _a2[j*_NbCol+i]=0.;
      _c1car[j*_NbCol+i]=0.;
      _c2car[j*_NbCol+i]=0.;
      _alpha[j*_NbCol+i]=0.;
      _beta[j*_NbCol+i]=0.;
      _p[j*_NbCol+i]=0.;
      _ptilde[j*_NbCol+i]=0.;
      _pscal[j*_NbCol+i]=0.;

    }
  }
}

Probleme::~Probleme()
{}

double Probleme::rho(double x, double y)
{
  if(_choix==1)
  {
    return _gamma*_gamma;
  }
  else if (_choix==2)
  {
    double distance;
    distance= pow(x-0.25,2)+pow(y-0.5,2);
    if(x<0.05)
    {
      return 3.86859;
    }
    if (distance<pow(0.15,2))
    {
      return 10.;
    }
    else
    {
      return 1.0;
    }
  }
  else if (_choix==3)
  {
    return 1.0-0.5*tanh(y-0.5);
  }

}

double Probleme::u1(double x, double y)
{
  double pi;
  pi=acos(-1);
  if(_choix==1)
  {
    return -sin(y);
  }
  else if (_choix==2)
  {
    if(x<0.05)
    {
      return 11.2536;
    }
    else
    {
      return 0.;
    }
  }
  else if (_choix==3)
  {
    return 2*pow(sin(pi*x),2)*sin(pi*y)*cos(pi*y);
  }
}

double Probleme::u2(double x, double y)
{
  double pi;
  pi=acos(-1);
  if(_choix==1)
  {
    return sin(x);
  }
  else if (_choix==2)
  {
    if(x<0.05)
    {
      return 0.;
    }
    else
    {
      return 0.;
    }
  }
  else if (_choix==3)
  {
      return -2*sin(pi*x)*cos(pi*x)*pow(sin(pi*y),2);
  }
}

double Probleme::u3(double x, double y)
{
  if(_choix==1)
  {
    return 0.;
  }
  else if (_choix==2)
  {
    if(x<0.05)
    {
      return 0.;
    }
    else
    {
      return 0.;
    }
  }
  else if (_choix==3)
  {
    return 0.;
  }
}

double Probleme::B1(double x, double y)
{
  if(_choix==1)
  {
    return -sin(y);
  }
  else if (_choix==2)
  {
    if(x<0.05)
    {
      return 0.;
    }
    else
    {
      return 0.;
    }
  }
  else if (_choix==3)
  {
    return 0.;
  }
}

double Probleme::B2(double x, double y)
{
  if(_choix==1)
  {
    return sin(2.*x);
  }
  else if (_choix==2)
  {
    if(x<0.05)
    {
      return 2.1826182;
    }
    else
    {
      return 0.56418958;
    }
  }
  else if (_choix==3)
  {
    return 0.;
  }
}

double Probleme::B3(double x, double y)
{
  if(_choix==1)
  {
    return 0.;
  }
  else if (_choix==2)
  {
    if(x<0.05)
    {
      return -2.1826182;
    }
    else
    {
      return 0.56418958;
    }
  }
  else if (_choix==3)
  {
    return 0.;
  }
}

double Probleme::p(double x, double y)
{
  if(_choix==1)
  {
    return _gamma;
  }
  else if (_choix==2)
  {
    if(x<0.05)
    {
      return 167.345;
    }
    else
    {
      return 1.0;
    }
  }
  else if (_choix==3)
  {
    return 0.001;
  }
}

double Probleme::E(double x, double y)
{
  double pscal, calculE;
  pscal= Probleme::rho(x,y)*(Probleme::u1(x,y)*Probleme::u1(x,y)+Probleme::u2(x,y)*Probleme::u2(x,y)+Probleme::u3(x,y)*Probleme::u3(x,y))+Probleme::B1(x,y)*Probleme::B1(x,y)+Probleme::B2(x,y)*Probleme::B2(x,y)+Probleme::B3(x,y)*Probleme::B3(x,y);
  calculE= Probleme::p(x,y)*(1./(_gamma-1.))+0.5*pscal;
  return calculE;
}

void Probleme::initialize_u()
{
  for (int i=0; i<_NbLignes; i++)
  {
    for (int j=0; j<_NbCol; j++)
    {
      _U[0][j*_NbCol+i]=Probleme::rho(i*_Dx,j*_Dy);
      _U[1][j*_NbCol+i]=Probleme::rho(i*_Dx,j*_Dy)*Probleme::u1(i*_Dx,j*_Dy);
      _U[2][j*_NbCol+i]=Probleme::rho(i*_Dx,j*_Dy)*Probleme::u2(i*_Dx,j*_Dy);
      _U[3][j*_NbCol+i]=Probleme::rho(i*_Dx,j*_Dy)*Probleme::u3(i*_Dx,j*_Dy);
      _U[4][j*_NbCol+i]=Probleme::B1(i*_Dx,j*_Dy);
      _U[5][j*_NbCol+i]=Probleme::B2(i*_Dx,j*_Dy);
      _U[6][j*_NbCol+i]=Probleme::B3(i*_Dx,j*_Dy);
      _U[7][j*_NbCol+i]=Probleme::E(i*_Dx,j*_Dy);
    }
  }
}

void Probleme::calcul_p(double x, double y, int i, int j)
{
  _p[j*_NbCol+i]=(_gamma-1.)*(_U[7][j*_NbCol+i]-(0.5/_U[0][j*_NbCol+i])*(pow(_U[1][j*_NbCol+i],2)+pow(_U[2][j*_NbCol+i],2)+pow(_U[3][j*_NbCol+i],2))-0.5*(pow(_U[4][j*_NbCol+i],2)+pow(_U[5][j*_NbCol+i],2)+pow(_U[6][j*_NbCol+i],2)));
  // if (_p[j*_NbCol+i]<0)
  // {
  //   _p[j*_NbCol+i]=0.;
  // }
}

void Probleme::calcul_p_tilde(double x, double y, int i, int j)
{
  _ptilde[j*_NbCol+i]= _p[j*_NbCol+i]+0.5*(pow(_U[4][j*_NbCol+i],2)+pow(_U[5][j*_NbCol+i],2)+pow(_U[6][j*_NbCol+i],2));
}

void Probleme::calcul_pscal(double x, double y, int i, int j)
{
  _pscal[j*_NbCol+i]=(_U[1][j*_NbCol+i]*_U[4][j*_NbCol+i]+_U[2][j*_NbCol+i]*_U[5][j*_NbCol+i]+_U[3][j*_NbCol+i]*_U[6][j*_NbCol+i])/_U[0][j*_NbCol+i];
}

void Probleme::calcul_f(double x, double y, int i, int j)
{
  _f[0][j*_NbCol+i]=_U[1][j*_NbCol+i];
  _f[1][j*_NbCol+i]=pow(_U[1][j*_NbCol+i],2)/_U[0][j*_NbCol+i]+_ptilde[j*_NbCol+i]-0.5*pow(_U[4][j*_NbCol+i],2);
  _f[2][j*_NbCol+i]=(_U[1][j*_NbCol+i]*_U[2][j*_NbCol+i])/_U[0][j*_NbCol+i]-_U[4][j*_NbCol+i]*_U[5][j*_NbCol+i];
  _f[3][j*_NbCol+i]=(_U[1][j*_NbCol+i]*_U[3][j*_NbCol+i])/_U[0][j*_NbCol+i]-_U[4][j*_NbCol+i]*_U[6][j*_NbCol+i];
  _f[4][j*_NbCol+i]=0.;
  _f[5][j*_NbCol+i]=-(_U[2][j*_NbCol+i]*_U[4][j*_NbCol+i]-_U[1][j*_NbCol+i]*_U[5][j*_NbCol+i])/_U[0][j*_NbCol+i];
  _f[6][j*_NbCol+i]=(_U[1][j*_NbCol+i]*_U[6][j*_NbCol+i]-_U[3][j*_NbCol+i]*_U[4][j*_NbCol+i])/_U[0][j*_NbCol+i];
  _f[7][j*_NbCol+i]=(_U[7][j*_NbCol+i]+_ptilde[j*_NbCol+i])*(_U[1][j*_NbCol+i]/_U[0][j*_NbCol+i]) - _pscal[j*_NbCol+i]*_U[4][j*_NbCol+i];
}

void Probleme::calcul_g(double x, double y, int i, int j)
{
  _g[0][j*_NbCol+i]=_U[2][j*_NbCol+i];
  _g[1][j*_NbCol+i]=(_U[1][j*_NbCol+i]*_U[2][j*_NbCol+i])/_U[0][j*_NbCol+i]-_U[4][j*_NbCol+i]*_U[5][j*_NbCol+i];
  _g[2][j*_NbCol+i]=pow(_U[2][j*_NbCol+i],2)/_U[0][j*_NbCol+i] + _ptilde[j*_NbCol+i] - 0.5*(pow(_U[5][j*_NbCol+i],2));
  _g[3][j*_NbCol+i]=(_U[1][j*_NbCol+i]*_U[3][j*_NbCol+i])/_U[0][j*_NbCol+i]-_U[4][j*_NbCol+i]*_U[6][j*_NbCol+i];
  _g[4][j*_NbCol+i]=(_U[2][j*_NbCol+i]*_U[4][j*_NbCol+i]-_U[1][j*_NbCol+i]*_U[5][j*_NbCol+i])/_U[0][j*_NbCol+i];
  _g[5][j*_NbCol+i]=0.;
  _g[6][j*_NbCol+i]=(_U[2][j*_NbCol+i]*_U[6][j*_NbCol+i]-_U[3][j*_NbCol+i]*_U[5][j*_NbCol+i])/_U[0][j*_NbCol+i];
  _g[7][j*_NbCol+i]=(_U[7][j*_NbCol+i] + _ptilde[j*_NbCol+i])*(_U[2][j*_NbCol+i]/_U[0][j*_NbCol+i]) - _pscal[j*_NbCol+i]*_U[5][j*_NbCol+i];
}

void Probleme::calcul_a2(double x, double y, int i, int j)
{
  _a2[j*_NbCol+i]=(_gamma*_p[j*_NbCol+i])/_U[0][j*_NbCol+i];
}

void Probleme::calcul_b(double x, double y, int i, int j)
{
  _b[0][j*_NbCol+i]=_U[4][j*_NbCol+i]/(sqrt(_U[0][j*_NbCol+i]));
  _b[1][j*_NbCol+i]=_U[5][j*_NbCol+i]/(sqrt(_U[0][j*_NbCol+i]));
  _b[2][j*_NbCol+i]=_U[6][j*_NbCol+i]/(sqrt(_U[0][j*_NbCol+i]));
}

void Probleme::calcul_ck(double x, double y, int i, int j, int k)
{
  double bkcar,c,bcarre;
  bcarre=pow(_b[0][j*_NbCol+i],2)+pow(_b[1][j*_NbCol+i],2)+pow(_b[2][j*_NbCol+i],2);
  c=0.5*(_a2[j*_NbCol+i]+bcarre);
  if (k==1)
  {
    bkcar=pow(_b[0][j*_NbCol+i],2);
    _c1car[j*_NbCol+i]=c+0.5*sqrt((pow(_a2[j*_NbCol+i]+bcarre,2))-4*_a2[j*_NbCol+i]*bkcar);
  }
  else if (k==2)
  {
    bkcar=pow(_b[1][j*_NbCol+i],2);
    _c2car[j*_NbCol+i]=c+0.5*sqrt((pow(_a2[j*_NbCol+i]+bcarre,2))-4*_a2[j*_NbCol+i]*bkcar);
  }
}

void Probleme::calcul_alpha_beta(double x, double y, int i, int j)
{
  _alpha[j*_NbCol+i]=abs(_U[1][j*_NbCol+i]/_U[0][j*_NbCol+i])+abs(sqrt(_c1car[j*_NbCol+i]));
  _beta[j*_NbCol+i]=abs(_U[2][j*_NbCol+i]/_U[0][j*_NbCol+i])+abs(sqrt(_c2car[j*_NbCol+i]));
}

void Probleme::flux_f(double x, double y, int i, int j)
{
  for (int k=0; k<5; k++)
  {
    _ff[k][j*_NbCol+i]=0.5*(_f[k][j*_NbCol+i]+_f[k][j*_NbCol+i+1]-max(abs(_alpha[j*_NbCol+i]),abs(_alpha[j*_NbCol+i+1]))*(_U[k][j*_NbCol+i+1]-_U[k][j*_NbCol+i]));
  }

  _ff[5][j*_NbCol+i]=0.5*(_f[5][j*_NbCol+i]+_f[5][j*_NbCol+i+1]+max(abs(_alpha[j*_NbCol+i]),abs(_alpha[j*_NbCol+i+1]))*(_U[5][j*_NbCol+i+1]-_U[5][j*_NbCol+i]));
  _ff[6][j*_NbCol+i]=0.5*(_f[6][j*_NbCol+i]+_f[6][j*_NbCol+i+1]-max(abs(_alpha[j*_NbCol+i]),abs(_alpha[j*_NbCol+i+1]))*(_U[6][j*_NbCol+i+1]-_U[6][j*_NbCol+i]));
  _ff[7][j*_NbCol+i]=0.5*(_f[7][j*_NbCol+i]+_f[7][j*_NbCol+i+1]-max(abs(_alpha[j*_NbCol+i]),abs(_alpha[j*_NbCol+i+1]))*(_U[7][j*_NbCol+i+1]-_U[7][j*_NbCol+i]));
}

void Probleme::flux_g(double x, double y, int i, int j)
{
  for (int k=0; k<8; k++)
  {
    _gg[k][j*_NbCol+i]=0.5*(_g[k][j*_NbCol+i]+_g[k][(j+1)*_NbCol+i]-max(abs(_beta[j*_NbCol+i]),abs(_beta[(j+1)*_NbCol+i]))*(_U[k][(j+1)*_NbCol+i]-_U[k][j*_NbCol+i]));
  }
}

void Probleme::TimeIteration_ordre1()
{
  int compteur;
  compteur =0;
  _t=0.;
  while (_t<=_tmax)
  {
    Probleme::SaveIeration("result" + std::to_string(compteur) );

    for (int j=0; j<_NbLignes; j++)
    {
      for (int i=0; i<_NbCol; i++)
      {
        Probleme::calcul_p(i*_Dx, j*_Dy, i, j);
        Probleme::calcul_p_tilde(i*_Dx, j*_Dy, i, j);
        Probleme::calcul_pscal(i*_Dx, j*_Dy, i, j);
        Probleme::calcul_f(i*_Dx, j*_Dy, i, j);
        Probleme::calcul_g(i*_Dx, j*_Dy, i, j);
        Probleme::calcul_a2(i*_Dx, j*_Dy, i, j);
        Probleme::calcul_b(i*_Dx, j*_Dy, i, j);
        Probleme::calcul_ck(i*_Dx, j*_Dy, i, j, 1);
        Probleme::calcul_ck(i*_Dx, j*_Dy, i, j, 2);
        Probleme::calcul_alpha_beta(i*_Dx, j*_Dy, i, j);
      }
    }
    for (int j=0; j<_NbLignes-1; j++)
    {
      for (int i=0; i<_NbCol-1; i++)
      {
        Probleme::flux_f(i*_Dx, j*_Dy, i, j);
        Probleme::flux_g(i*_Dx, j*_Dy, i, j);
      }
    }
    for (int j=0; j<_NbLignes-1; j++)
    {
      Probleme::flux_g((_NbCol-1)*_Dx, j*_Dy, _NbCol-1, j);
    }
    for (int i=0; i<_NbCol-1; i++)
    {
      Probleme::flux_f(i*_Dx, (_NbLignes-1)*_Dy, i, _NbLignes-1);
    }

    // conditions aux limites
    for (int i=0; i<_NbCol; i++)
    {
      for (int k=0; k<8; k++)
      {
        if(_choix==1)
        {
          _gg[k][(_NbLignes-1)*_NbCol+i]=0.5*(_g[k][(_NbLignes-1)*_NbCol+i]+_g[k][i]-max(abs(_beta[(_NbLignes-1)*_NbCol+i]),abs(_beta[i]))*(_U[k][i]-_U[k][(_NbLignes-1)*_NbCol+i]));
        }
        else if (_choix==2 ||_choix==3)
        {
          _gg[k][(_NbLignes-1)*_NbCol+i]=_gg[k][(_NbLignes-2)*_NbCol+i];
        }
      }

    }

    for (int j=0; j<_NbLignes;j++)
    {
      if (_choix==1)
      {
        for (int k=0; k<5; k++)
        {
          _ff[k][j*_NbCol+_NbCol-1]=0.5*(_f[k][j*_NbCol+_NbCol-1]+_f[k][j*_NbCol]-max(abs(_alpha[j*_NbCol+_NbCol-1]),abs(_alpha[j*_NbCol]))*(_U[k][j*_NbCol]-_U[k][j*_NbCol+_NbCol-1]));
        }

        _ff[5][j*_NbCol+_NbCol-1]=0.5*(_f[5][j*_NbCol+_NbCol-1]+_f[5][j*_NbCol]+max(abs(_alpha[j*_NbCol+_NbCol-1]),abs(_alpha[j*_NbCol]))*(_U[5][j*_NbCol]-_U[5][j*_NbCol+_NbCol-1]));
        _ff[6][j*_NbCol+_NbCol-1]=0.5*(_f[6][j*_NbCol+_NbCol-1]+_f[6][j*_NbCol]-max(abs(_alpha[j*_NbCol+_NbCol-1]),abs(_alpha[j*_NbCol]))*(_U[6][j*_NbCol]-_U[6][j*_NbCol+_NbCol-1]));
        _ff[7][j*_NbCol+_NbCol-1]=0.5*(_f[7][j*_NbCol+_NbCol-1]+_f[7][j*_NbCol]-max(abs(_alpha[j*_NbCol+_NbCol-1]),abs(_alpha[j*_NbCol]))*(_U[7][j*_NbCol]-_U[7][j*_NbCol+_NbCol-1]));
      }
      else if (_choix==2 ||_choix==3)
      {
        for (int k=0; k<8; k++)
        {
          _ff[k][j*_NbCol+_NbCol-1]=_ff[k][j*_NbCol+_NbCol-2];
        }
      }

    }


    double maxeigenvalues,test;
    maxeigenvalues=0;

    for (int j=0; j<_NbLignes; j++)
    {
      for (int i=0; i<_NbCol; i++)
      {
        test=max(abs(_alpha[j*_NbCol+i]),abs(_beta[j*_NbCol+i]));
        if (test>maxeigenvalues)
        {
          maxeigenvalues=test;
        }
      }
    }

  //  std::cout<< "max B1 = " << maxB1 <<std::endl;
    _Dt=0.45*_Dx/(2*maxeigenvalues);

    for (int j=1; j<_NbLignes; j++)
    {
      for (int i=1; i<_NbCol; i++)
      {
        for(int indice=0; indice<8; indice++)
        {
          _Uapres[indice][j*_NbCol+i]=_U[indice][j*_NbCol+i]-_Dt*((_ff[indice][j*_NbCol+i]-_ff[indice][j*_NbCol+i-1])/_Dx+(_gg[indice][j*_NbCol+i]-_gg[indice][(j-1)*_NbCol+i])/_Dy);
        }
      }
    }

    for (int i=1; i<_NbCol; i++)
    {
      for(int indice=0; indice<8; indice++)
      {
        if(_choix==1)
        {
          _Uapres[indice][i]=_U[indice][i]-_Dt*((_ff[indice][i]-_ff[indice][i-1])/_Dx+(_gg[indice][i]-_gg[indice][(_NbLignes-1)*_NbCol+i])/_Dy);
        }
        else if (_choix==2 ||_choix==3)
        {
          _Uapres[indice][i]=_Uapres[indice][_NbCol+i];
        }
      }
    }


    for (int j=1; j<_NbLignes; j++)
    {
      for(int indice=0; indice<8; indice++)
      {
        if (_choix==1)
        {
          _Uapres[indice][j*_NbCol]=_U[indice][j*_NbCol]-_Dt*((_ff[indice][j*_NbCol]-_ff[indice][j*_NbCol+_NbCol-1])/_Dx+(_gg[indice][j*_NbCol]-_gg[indice][(j-1)*_NbCol])/_Dy);
        }
        else if (_choix==2 ||_choix==3)
        {
          _Uapres[indice][j*_NbCol]=_Uapres[indice][j*_NbCol+1];
        }
      }
    }

    for(int indice=0; indice<8; indice++)
    {
      if (_choix==1)
      {
        _Uapres[indice][0]=_U[indice][0]-_Dt*((_ff[indice][0]-_ff[indice][_NbCol-1])/_Dx+(_gg[indice][0]-_gg[indice][(_NbLignes-1)*_NbCol])/_Dy);
      }
      else if (_choix==2 ||_choix==3)
      {
        _Uapres[indice][0]=_Uapres[indice][1];
      }
    }

    _t=_t+_Dt;

    compteur +=1;
   //  std::cout<< compteur <<std::endl;

    _U=_Uapres;

  }
    //  Probleme::SaveIeration("result" + std::to_string(_t));
}

double Probleme::minmod(double a, double b, double c)
{
  double test1, test2;
  if (a*b<0 || b*c <0 || a*c<0)
  {
    return 0.;
  }
  else
  {
    test1=min(abs(a), abs(b));
    test2=min(test1, abs(c));
    if (a<0)
    {
      return -test2;
    }
    else
    {
      return test2;
    }
  }
}

void Probleme::calcul_Ux(double x, double y, int i, int j)
{
  double a,b,c;
  for (int indice=0; indice<8; indice++)
  {
    a=_U[indice][j*_NbCol+i+1]-_U[indice][j*_NbCol+i];
    b=0.5*(_U[indice][j*_NbCol+i+1]-_U[indice][j*_NbCol+i-1]);
    c=_U[indice][j*_NbCol+i]-_U[indice][j*_NbCol+i-1];
    _Ux[indice][j*_NbCol+i]=Probleme::minmod(a,b,c);

  }
}

void Probleme::calcul_Uy(double x, double y, int i, int j)
{
  double q,r,s;
  for (int indice=0; indice<8; indice++)
  {

    q=_U[indice][(j+1)*_NbCol+i]-_U[indice][j*_NbCol+i];
    r=0.5*(_U[indice][(j+1)*_NbCol+i]-_U[indice][(j-1)*_NbCol+i]);
    s=_U[indice][j*_NbCol+i]-_U[indice][(j-1)*_NbCol+i];
    _Uy[indice][j*_NbCol+i]=Probleme::minmod(q,r,s);
  }
}

std::vector<double> Probleme::calcul_pij(double x, double y, int i, int j)
{
  std::vector<double> pij;
  pij.resize(8);

  for(int k=0; k<8; k++)
  {
    pij[k]=_U[k][j*_NbCol+i]+(x-i*_Dx)*_Ux[k][j*_NbCol+i]/_Dx+(y-j*_Dy)*_Uy[k][j*_NbCol+i]/_Dy;
  }
  return pij;
}

void Probleme::calcul_g1_g2(double x, double y, int i, int j)
{
  _g1[0][j*_NbCol+i]=_UN[2][j*_NbCol+i];
  _g1[1][j*_NbCol+i]=(_UN[1][j*_NbCol+i]*_UN[2][j*_NbCol+i])/_UN[0][j*_NbCol+i]-_UN[4][j*_NbCol+i]*_UN[5][j*_NbCol+i];
  _g1[2][j*_NbCol+i]=pow(_UN[2][j*_NbCol+i],2)/_UN[0][j*_NbCol+i] + _ptilde[j*_NbCol+i] - 0.5*(pow(_UN[5][j*_NbCol+i],2));
  _g1[3][j*_NbCol+i]=(_UN[1][j*_NbCol+i]*_UN[3][j*_NbCol+i])/_UN[0][j*_NbCol+i]-_UN[4][j*_NbCol+i]*_UN[6][j*_NbCol+i];
  _g1[4][j*_NbCol+i]=(1./_UN[0][j*_NbCol+i])*(_UN[2][j*_NbCol+i]*_UN[4][j*_NbCol+i]-_UN[1][j*_NbCol+i]*_UN[5][j*_NbCol+i]);
  _g1[5][j*_NbCol+i]=0.;
  _g1[6][j*_NbCol+i]=(1./_UN[0][j*_NbCol+i])*(_UN[2][j*_NbCol+i]*_UN[6][j*_NbCol+i]-_UN[3][j*_NbCol+i]*_UN[5][j*_NbCol+i]);
  _g1[7][j*_NbCol+i]=(_UN[7][j*_NbCol+i] + _ptilde[j*_NbCol+i])*(_UN[2][j*_NbCol+i]/_UN[0][j*_NbCol+i]) - _pscal[j*_NbCol+i]*_UN[5][j*_NbCol+i];

  _g2[0][j*_NbCol+i]=_US[2][j*_NbCol+i];
  _g2[1][j*_NbCol+i]=(_US[1][j*_NbCol+i]*_US[2][j*_NbCol+i])/_US[0][j*_NbCol+i]-_US[4][j*_NbCol+i]*_US[5][j*_NbCol+i];
  _g2[2][j*_NbCol+i]=pow(_US[2][j*_NbCol+i],2)/_US[0][j*_NbCol+i] + _ptilde[j*_NbCol+i] - 0.5*(pow(_US[5][j*_NbCol+i],2));
  _g2[3][j*_NbCol+i]=(_US[1][j*_NbCol+i]*_US[3][j*_NbCol+i])/_US[0][j*_NbCol+i]-_US[4][j*_NbCol+i]*_US[6][j*_NbCol+i];
  _g2[4][j*_NbCol+i]=(1./_US[0][j*_NbCol+i])*(_US[2][j*_NbCol+i]*_US[4][j*_NbCol+i]-_US[1][j*_NbCol+i]*_US[5][j*_NbCol+i]);
  _g2[5][j*_NbCol+i]=0.;
  _g2[6][j*_NbCol+i]=(1./_US[0][j*_NbCol+i])*(_US[2][j*_NbCol+i]*_US[6][j*_NbCol+i]-_US[3][j*_NbCol+i]*_US[5][j*_NbCol+i]);
  _g2[7][j*_NbCol+i]=(_US[7][j*_NbCol+i] + _ptilde[j*_NbCol+i])*(_US[2][j*_NbCol+i]/_US[0][j*_NbCol+i]) - _pscal[j*_NbCol+i]*_US[5][j*_NbCol+i];

}
//f pour l'ordre 2 MUSCL
void Probleme::calcul_f1_f2(double x, double y, int i, int j)
{
  _f1[0][j*_NbCol+i]=_UE[1][j*_NbCol+i];
  _f1[1][j*_NbCol+i]=pow(_UE[1][j*_NbCol+i],2)/_UE[0][j*_NbCol+i]+_ptilde[j*_NbCol+i]-0.5*pow(_UE[4][j*_NbCol+i],2);
  _f1[2][j*_NbCol+i]=(_UE[1][j*_NbCol+i]*_UE[2][j*_NbCol+i])/_UE[0][j*_NbCol+i]-_UE[4][j*_NbCol+i]*_UE[5][j*_NbCol+i];
  _f1[3][j*_NbCol+i]=(_UE[1][j*_NbCol+i]*_UE[3][j*_NbCol+i])/_UE[0][j*_NbCol+i]-_UE[4][j*_NbCol+i]*_UE[6][j*_NbCol+i];
  _f1[4][j*_NbCol+i]=0.;
  _f1[5][j*_NbCol+i]=-(1./_UE[0][j*_NbCol+i])*(_UE[2][j*_NbCol+i]*_UE[4][j*_NbCol+i]-_UE[1][j*_NbCol+i]*_UE[5][j*_NbCol+i]);
  _f1[6][j*_NbCol+i]=(1./_UE[0][j*_NbCol+i])*(_UE[1][j*_NbCol+i]*_UE[6][j*_NbCol+i]-_UE[3][j*_NbCol+i]*_UE[4][j*_NbCol+i]);
  _f1[7][j*_NbCol+i]=(_UE[7][j*_NbCol+i]+_ptilde[j*_NbCol+i])*(_UE[1][j*_NbCol+i]/_UE[0][j*_NbCol+i]) - _pscal[j*_NbCol+i]*_UE[4][j*_NbCol+i];

  _f2[0][j*_NbCol+i]=_UW[1][j*_NbCol+i];
  _f2[1][j*_NbCol+i]=pow(_UW[1][j*_NbCol+i],2)/_UW[0][j*_NbCol+i]+_ptilde[j*_NbCol+i]-0.5*pow(_UW[4][j*_NbCol+i],2);
  _f2[2][j*_NbCol+i]=(_UW[1][j*_NbCol+i]*_UW[2][j*_NbCol+i])/_UW[0][j*_NbCol+i]-_UW[4][j*_NbCol+i]*_UW[5][j*_NbCol+i];
  _f2[3][j*_NbCol+i]=(_UW[1][j*_NbCol+i]*_UW[3][j*_NbCol+i])/_UW[0][j*_NbCol+i]-_UW[4][j*_NbCol+i]*_UW[6][j*_NbCol+i];
  _f2[4][j*_NbCol+i]=0.;
  _f2[5][j*_NbCol+i]=-(1./_UW[0][j*_NbCol+i])*(_UW[2][j*_NbCol+i]*_UW[4][j*_NbCol+i]-_UW[1][j*_NbCol+i]*_UW[5][j*_NbCol+i]);
  _f2[6][j*_NbCol+i]=(1./_UW[0][j*_NbCol+i])*(_UW[1][j*_NbCol+i]*_UW[6][j*_NbCol+i]-_UW[3][j*_NbCol+i]*_UW[4][j*_NbCol+i]);
  _f2[7][j*_NbCol+i]=(_UW[7][j*_NbCol+i]+_ptilde[j*_NbCol+i])*(_UW[1][j*_NbCol+i]/_UW[0][j*_NbCol+i]) - _pscal[j*_NbCol+i]*_UW[4][j*_NbCol+i];

}

void Probleme::calcul_flux_ordre2F(double x, double y, int i, int j)
{
  for (int k=0; k<5; k++)
  {
    _F2[k][j*_NbCol+i]=0.5*(_f1[k][j*_NbCol+i]+_f2[k][j*_NbCol+i+1])-max(abs(_alpha[j*_NbCol+i]),abs(_alpha[j*_NbCol+i+1]))*(_UW[k][j*_NbCol+i+1]-_UE[k][j*_NbCol+i]);
  }

  _F2[5][j*_NbCol+i]=0.5*(_f1[5][j*_NbCol+i]+_f2[5][j*_NbCol+i+1])+max(abs(_alpha[j*_NbCol+i]),abs(_alpha[j*_NbCol+i+1]))*(_UW[5][j*_NbCol+i+1]-_UE[5][j*_NbCol+i]);
  _F2[6][j*_NbCol+i]=0.5*(_f1[6][j*_NbCol+i]+_f2[6][j*_NbCol+i+1])-max(abs(_alpha[j*_NbCol+i]),abs(_alpha[j*_NbCol+i+1]))*(_UW[6][j*_NbCol+i+1]-_UE[6][j*_NbCol+i]);
  _F2[7][j*_NbCol+i]=0.5*(_f1[7][j*_NbCol+i]+_f2[7][j*_NbCol+i+1])-max(abs(_alpha[j*_NbCol+i]),abs(_alpha[j*_NbCol+i+1]))*(_UW[7][j*_NbCol+i+1]-_UE[7][j*_NbCol+i]);
}

void Probleme::calcul_flux_ordre2G(double x, double y, int i, int j)
{
  for (int k=0; k<8; k++)
  {
    _G2[k][j*_NbCol+i]=0.5*(_g1[k][j*_NbCol+i]+_g2[k][(j+1)*_NbCol+i])-max(abs(_beta[j*_NbCol+i]),abs(_beta[(j+1)*_NbCol+i]))*(_US[k][(j+1)*_NbCol+i]-_UN[k][j*_NbCol+i]);
  }
}

void Probleme::calcul_u_directions(double x, double y, int i, int j)
{
  std::vector<double> a,b,c,d;
  a=Probleme::calcul_pij((i+0.5)*_Dx, j*_Dy, i,j);
  b=Probleme::calcul_pij((i-0.5)*_Dx, j*_Dy, i,j);
  c=Probleme::calcul_pij(i*_Dx, (j+0.5)*_Dy, i,j);
  d=Probleme::calcul_pij(i*_Dx, (j-0.5)*_Dy, i,j);
  for (int k=0; k<8; k++)
  {
    _UE[k][j*_NbCol+i]=a[k];
    _UW[k][j*_NbCol+i]=b[k];
    _UN[k][j*_NbCol+i]=c[k];
    _US[k][j*_NbCol+i]=d[k];
  }
}

void Probleme::Time_iteration_MUSCL()
{
  int compteur;
  compteur =0;
  _t=0.;
  while (_t<=_tmax)
  {
    Probleme::SaveIeration("resultMUSCL" + std::to_string(compteur) );

    //CALCUL DE UX ET UY
    double a,b,c,q,r,s, a1, b1, c1, q1,r1,s1;
    for (int indice=0; indice<8; indice++)
    {
      for (int j=0;j<_NbLignes; j++)
      {
        //i=0
        a=_U[indice][j*_NbCol+1]-_U[indice][j*_NbCol];
        b=0.5*(_U[indice][j*_NbCol+1]-_U[indice][j*_NbCol+_NbCol-1]);
        c=_U[indice][j*_NbCol]-_U[indice][j*_NbCol+_NbCol-1];
        _Ux[indice][j*_NbCol]=Probleme::minmod(a,b,c);

        //i=Nb_Col-1
        a1=_U[indice][j*_NbCol]-_U[indice][j*_NbCol+_NbCol-1];
        b1=0.5*(_U[indice][j*_NbCol]-_U[indice][j*_NbCol+_NbCol-2]);
        c1=_U[indice][j*_NbCol+_NbCol-1]-_U[indice][j*_NbCol+_NbCol-2];
        _Ux[indice][j*_NbCol+_NbCol-1]=Probleme::minmod(a1,b1,c1);
      }

      for (int i=0; i<_NbCol; i++)
      {
        //j=0
        q=_U[indice][_NbCol+i]-_U[indice][i];
        r=0.5*(_U[indice][_NbCol+i]-_U[indice][(_NbLignes-1)*_NbCol+i]);
        s=_U[indice][i]-_U[indice][(_NbLignes-1)*_NbCol+i];
        _Uy[indice][i]=Probleme::minmod(q,r,s);

        //j=_NbLignes-1
        q1=_U[indice][i]-_U[indice][(_NbLignes-1)*_NbCol+i];
        r1=0.5*(_U[indice][i]-_U[indice][(_NbLignes-2)*_NbCol+i]);
        s1=_U[indice][(_NbLignes-1)*_NbCol+i]-_U[indice][(_NbLignes-2)*_NbCol+i];
        _Uy[indice][(_NbLignes-1)*_NbCol+i]=Probleme::minmod(q1,r1,s1);
      }
    }

    for (int j=1; j<_NbLignes-1; j++)
    {
      //i=0
      Probleme::calcul_Uy(0, j*_Dy, 0, j);
      //i=_Ncol-1
      Probleme::calcul_Uy((_NbCol-1)*_Dx, j*_Dy, _NbCol-1, j);
    }

    for (int i=1; i<_NbCol-1; i++)
    {
      //j=0
      Probleme::calcul_Ux(i*_Dx, 0, i, 0);
      //j=_NbLignes-1
      Probleme::calcul_Ux(i*_Dx, (_NbLignes-1)*_Dy, i, _NbLignes-1);
    }

    for (int j=1; j<_NbLignes-1; j++)
    {
      for (int i=1; i<_NbCol-1; i++)
      {
        Probleme::calcul_Ux(i*_Dx, j*_Dy, i, j);
        Probleme::calcul_Uy(i*_Dx, j*_Dy, i, j);
      }
    }
    //--------------------------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------------------------

    //u directions
    for (int j=0; j<_NbLignes; j++)
    {
      for (int i=0; i<_NbCol; i++)
      {
        Probleme::calcul_u_directions(i*_Dx, j*_Dy, i, j);
      }
    }

    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------


    for (int j=0; j<_NbLignes; j++)
    {
      for (int i=0; i<_NbCol; i++)
      {
        Probleme::calcul_p(i*_Dx, j*_Dy, i, j);
        Probleme::calcul_p_tilde(i*_Dx, j*_Dy, i, j);
        Probleme::calcul_pscal(i*_Dx, j*_Dy, i, j);
        Probleme::calcul_f1_f2(i*_Dx, j*_Dy, i, j);
        Probleme::calcul_g1_g2(i*_Dx, j*_Dy, i, j);
        Probleme::calcul_a2(i*_Dx, j*_Dy, i, j);
        Probleme::calcul_b(i*_Dx, j*_Dy, i, j);
        Probleme::calcul_ck(i*_Dx, j*_Dy, i, j, 1);
        Probleme::calcul_ck(i*_Dx, j*_Dy, i, j, 2);
        Probleme::calcul_alpha_beta(i*_Dx, j*_Dy, i, j);
      }
    }
    for (int j=0; j<_NbLignes-1; j++)
    {
      for (int i=0; i<_NbCol-1; i++)
      {
        Probleme::calcul_flux_ordre2F(i*_Dx, j*_Dy, i, j);
        Probleme::calcul_flux_ordre2G(i*_Dx, j*_Dy, i, j);
      }
    }
    for (int j=0; j<_NbLignes-1; j++)
    {
      Probleme::calcul_flux_ordre2G((_NbCol-1)*_Dx, j*_Dy, _NbCol-1, j);
    }
    for (int i=0; i<_NbCol-1; i++)
    {
      Probleme::calcul_flux_ordre2F(i*_Dx, (_NbLignes-1)*_Dy, i, _NbLignes-1);
    }


    // conditions aux limites
    for (int i=0; i<_NbCol; i++)
    {
      for (int k=0; k<8; k++)
      {
        _G2[k][(_NbLignes-1)*_NbCol+i]=0.5*(_g1[k][(_NbLignes-1)*_NbCol+i]+_g2[k][i])-max(abs(_beta[(_NbLignes-1)*_NbCol+i]),abs(_beta[i]))*(_US[k][i]-_UN[k][(_NbLignes-1)*_NbCol+i]);
      }
    }

    for (int j=0; j<_NbLignes;j++)
    {
      for (int k=0; k<5; k++)
      {
        _F2[k][j*_NbCol+_NbCol-1]=0.5*(_f1[k][j*_NbCol+_NbCol-1]+_f2[k][j*_NbCol])-max(abs(_alpha[j*_NbCol+_NbCol-1]),abs(_alpha[j*_NbCol]))*(_UW[k][j*_NbCol]-_UE[k][j*_NbCol+_NbCol-1]);
      }

      _F2[5][j*_NbCol+_NbCol-1]=0.5*(_f1[5][j*_NbCol+_NbCol-1]+_f2[5][j*_NbCol])+max(abs(_alpha[j*_NbCol+_NbCol-1]),abs(_alpha[j*_NbCol]))*(_UW[5][j*_NbCol]-_UE[5][j*_NbCol+_NbCol-1]);
      _F2[6][j*_NbCol+_NbCol-1]=0.5*(_f1[6][j*_NbCol+_NbCol-1]+_f2[6][j*_NbCol])-max(abs(_alpha[j*_NbCol+_NbCol-1]),abs(_alpha[j*_NbCol]))*(_UW[6][j*_NbCol]-_UE[6][j*_NbCol+_NbCol-1]);
      _F2[7][j*_NbCol+_NbCol-1]=0.5*(_f1[7][j*_NbCol+_NbCol-1]+_f2[7][j*_NbCol])-max(abs(_alpha[j*_NbCol+_NbCol-1]),abs(_alpha[j*_NbCol]))*(_UW[7][j*_NbCol]-_UE[7][j*_NbCol+_NbCol-1]);    }


    double maxeigenvalues,test;
    maxeigenvalues=0;
    for (int j=0; j<_NbLignes; j++)
    {
      for (int i=0; i<_NbCol; i++)
      {
        test=max(abs(_alpha[j*_NbCol+i]),abs(_beta[j*_NbCol+i]));
        if (test>maxeigenvalues)
        {
          maxeigenvalues=test;
        }
      }
    }

    _Dt=0.45*_Dx/maxeigenvalues;

    for (int j=1; j<_NbLignes; j++)
    {
      for (int i=1; i<_NbCol; i++)
      {
        for(int indice=0; indice<8; indice++)
        {
          _Uapres[indice][j*_NbCol+i]=_U[indice][j*_NbCol+i]-_Dt*((1./_Dx)*(_F2[indice][j*_NbCol+i]-_F2[indice][j*_NbCol+i-1])+(1./_Dy)*(_G2[indice][j*_NbCol+i]-_G2[indice][(j-1)*_NbCol+i]));
        }
      }
    }

    for (int i=1; i<_NbCol; i++)
    {
      for(int indice=0; indice<8; indice++)
      {
        _Uapres[indice][i]=_U[indice][i]-_Dt*((1./_Dx)*(_F2[indice][i]-_F2[indice][i-1])+(1./_Dy)*(_G2[indice][i]-_G2[indice][(_NbLignes-1)*_NbCol+i]));

      }
    }

    for (int j=1; j<_NbLignes; j++)
    {
      for(int indice=0; indice<8; indice++)
      {
        _Uapres[indice][j*_NbCol]=_U[indice][j*_NbCol]-_Dt*((1./_Dx)*(_F2[indice][j*_NbCol]-_F2[indice][j*_NbCol+_NbCol-1])+(1./_Dy)*(_G2[indice][j*_NbCol]-_G2[indice][(j-1)*_NbCol]));
      }
    }

    for(int indice=0; indice<8; indice++)
    {
      _Uapres[indice][0]=_U[indice][0]-_Dt*((1./_Dx)*(_F2[indice][0]-_F2[indice][_NbCol-1])+(1./_Dy)*(_G2[indice][0]-_G2[indice][(_NbLignes-1)*_NbCol]));
    }

    _t=_t+_Dt;

    compteur +=1;
   //  std::cout<< compteur <<std::endl;

    _U=_Uapres;

  }
    //  Probleme::SaveIeration("result" + std::to_string(_t));
}

void Probleme::SaveIeration(std::string fichier)
{
  std::ofstream mon_fluxP;
  mon_fluxP.open(fichier + "_p.txt", std::ios::out);
  mon_fluxP << _t << std::endl;

  for (int i=0; i<_NbLignes; i++)
  {
    for (int j=0; j<_NbCol; j++)
    {
      //mon_fluxP << i*_Dx << "    " << j*_Dy << "     " <<_a2[j*_NbCol+i] <<"    "<< _b[0][j*_NbCol+i]<< "    "<<  _b[1][j*_NbCol+i]<< "    "<<  _b[2][j*_NbCol+i]<<std::endl;//<< "     "<< _c1car[j*_NbCol+i]<< "     " <<_c2car[j*_NbCol+i];
      //mon_fluxP << "     "<< _alpha[j*_NbCol+i]<< "         "<<_beta[j*_NbCol+i]<< std::endl;
      mon_fluxP << i*_Dx << "  " << j*_Dy << "   " << _U[0][j*_NbCol+i]<< "   " << _U[1][j*_NbCol+i]<< "   " << _U[2][j*_NbCol+i];
      mon_fluxP<< "   " << _U[3][j*_NbCol+i]<< "   " << _U[4][j*_NbCol+i]<< "   " << _U[5][j*_NbCol+i]<< "   " << _U[6][j*_NbCol+i];
      mon_fluxP<<"   " << _U[7][j*_NbCol+i]<<"    "<<_p[j*_NbCol+i]<<std::endl; //<< "        " << _p[j*_NbCol+i] << std::endl;
      //mon_fluxP << i*_Dx << " " << j*_Dy << " " << _U[0][j*_NbCol+i] << std::endl;

    }
  }

  mon_fluxP.close();

  // std::ofstream mon_fluxU;
  // mon_fluxU.open(fichier + "_u.txt", std::ios::out);
  // mon_fluxU << "# champ de vitesses sur un maillage carré" << std::endl;
  //
  // std::vector<double> norme;
  // norme.resize(Nx1*Nx1);
  //
  // double normeMax = 0;
  //
  // for (int i = 0; i < Nx1; i++)
  // {
  //   for (int j = 0; j < Nx1; j++)
  //   {
  //     double u1 = U[j + i*Nx1];
  //     double u2 = U[Nx1*Nx1 + j + i*Nx1];
  //     norme[j + i*Nx1] = sqrt(u1*u1 + u2*u2);
  //     if (norme[j+i*Nx1] > normeMax)
  //     {
  //       normeMax = norme[j + i*Nx1];
  //     }
  //   }
  // }
  //
  // for (int i = 0; i < Nx1; i++)
  // {
  //   for (int j = 0; j < Nx1; j++)
  //   {
  //     double u1 = U[j + i*Nx1];
  //     double u2 = U[Nx1*Nx1 + j + i*Nx1];
  //     double coeff = dx1/normeMax;
  //     mon_fluxU << j*dx1 << " " << i*dx1 << " " << u1*coeff << " " << u2*coeff << " " << norme[j + i*Nx1] << std::endl;
  //   }
  // }
  //
  // mon_fluxU.close();
}
