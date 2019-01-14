#include "VF.h"

using namespace std;

Probleme::Probleme(int NbLignes, int NbCol)
{
  _NbCol=NbCol;
  _NbLignes=NbLignes;
  _U.resize(8);
  _Uapres.resize(8);
  _f.resize(8);
  _g.resize(8);
  _ff.resize(8);
  _gg.resize(8);

  for (int j=0; j<8; j++)
  {
    _U[j].resize(_NbLignes*_NbCol);
    _Uapres[j].resize(_NbLignes*_NbCol);
    _f[j].resize(_NbLignes*_NbCol);
    _g[j].resize(_NbLignes*_NbCol);
    _ff[j].resize(_NbLignes*_NbCol);
    _gg[j].resize(_NbLignes*_NbCol);
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

  _Lx=2*3.1415;
  _Ly=2*3.1415;
  _Dx=_Lx/(_NbCol+1);
  _Dy=_Ly/(_NbLignes+1);
  _gamma=1.4;
  _tmax=3.1415;


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

    }
  }
}

Probleme::~Probleme()
{}

  double Probleme::rho(double x, double y)
  {
    return _gamma*_gamma;
  }

  double Probleme::u1(double x, double y)
  {
    return -sin(y);
  }

  double Probleme::u2(double x, double y)
  {
    return sin(x);
  }

  double Probleme::u3(double x, double y)
  {
    return 0.;
  }

  double Probleme::B1(double x, double y)
  {
    return -sin(y);
  }

  double Probleme::B2(double x, double y)
  {
    return sin(2*x);
  }

  double Probleme::B3(double x, double y)
  {
    return 0.;
  }

  double Probleme::p(double x, double y)
  {
    return _gamma;
  }

  double Probleme::E(double x, double y)
  {
    double pscal, calculE;
    pscal= Probleme::rho(x,y)*(Probleme::u1(x,y)*Probleme::u1(x,y)+Probleme::u2(x,y)*Probleme::u2(x,y)+Probleme::u3(x,y)*Probleme::u3(x,y))+Probleme::B1(x,y)*Probleme::B1(x,y)+Probleme::B2(x,y)*Probleme::B2(x,y)+Probleme::B3(x,y)*Probleme::B3(x,y);
    calculE= Probleme::p(x,y)*(1/(_gamma-1))+0.5*pscal;
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

  void Probleme::calcul_f(double x, double y, int i, int j)
  {
    _f[0][j*_NbCol+i]=_U[1][j*_NbCol+i];
    _f[1][j*_NbCol+i]=pow(_U[1][j*_NbCol+i],2)/_U[0][j*_NbCol+i]+(_gamma -1)*(_U[7][j*_NbCol+i]-0.5*((1/(_U[0][j*_NbCol+i]))*(pow(_U[1][j*_NbCol+i],2)+pow(_U[2][j*_NbCol+i],2)+pow(_U[3][j*_NbCol+i],2))+pow(_U[4][j*_NbCol+i],2)+pow(_U[5][j*_NbCol+i],2)+pow(_U[6][j*_NbCol+i],2)))+0.5*(pow(_U[5][j*_NbCol+i],2)+pow(_U[6][j*_NbCol+i],2));
    _f[2][j*_NbCol+i]=(_U[1][j*_NbCol+i]*_U[2][j*_NbCol+i])/_U[0][j*_NbCol+i]-_U[4][j*_NbCol+i]*_U[5][j*_NbCol+i];
    _f[3][j*_NbCol+i]=(_U[1][j*_NbCol+i]*_U[3][j*_NbCol+i])/_U[0][j*_NbCol+i]-_U[4][j*_NbCol+i]*_U[6][j*_NbCol+i];
    _f[4][j*_NbCol+i]=0.;
    _f[5][j*_NbCol+i]=-(1/_U[0][j*_NbCol+i])*(_U[2][j*_NbCol+i]*_U[4][j*_NbCol+i]-_U[1][j*_NbCol+i]*_U[5][j*_NbCol+i]);
    _f[6][j*_NbCol+i]=(1/_U[0][j*_NbCol+i])*(_U[1][j*_NbCol+i]*_U[6][j*_NbCol+i]-_U[3][j*_NbCol+i]*_U[4][j*_NbCol+i]);
    double pscal,ptilde;
    ptilde=(_gamma -1)*(_U[7][j*_NbCol+i]-0.5*((1/(_U[0][j*_NbCol+i]))*(pow(_U[1][j*_NbCol+i],2)+pow(_U[2][j*_NbCol+i],2)+pow(_U[3][j*_NbCol+i],2))+pow(_U[4][j*_NbCol+i],2)+pow(_U[5][j*_NbCol+i],2)+pow(_U[6][j*_NbCol+i],2)))+0.5*(pow(_U[4][j*_NbCol+i],2)+pow(_U[5][j*_NbCol+i],2)+pow(_U[6][j*_NbCol+i],2));
    pscal=(_U[4][j*_NbCol+i]/_U[0][j*_NbCol+i])*(_U[1][j*_NbCol+i]*_U[4][j*_NbCol+i]+_U[2][j*_NbCol+i]*_U[5][j*_NbCol+i]+_U[3][j*_NbCol+i]*_U[6][j*_NbCol+i]);
    _f[7][j*_NbCol+i]=_U[1][j*_NbCol+i]*_U[7][j*_NbCol+i]/_U[0][j*_NbCol+i]+ _U[1][j*_NbCol+i]/_U[0][j*_NbCol+i]*ptilde -pscal;
  }

  void Probleme::calcul_g(double x, double y, int i, int j)
  {
    double ptilde;
    ptilde=(_gamma -1)*(_U[7][j*_NbCol+i]-0.5*((1/(_U[0][j*_NbCol+i]))*(pow(_U[1][j*_NbCol+i],2)+pow(_U[2][j*_NbCol+i],2)+pow(_U[3][j*_NbCol+i],2))+pow(_U[4][j*_NbCol+i],2)+pow(_U[5][j*_NbCol+i],2)+pow(_U[6][j*_NbCol+i],2)))+0.5*(pow(_U[4][j*_NbCol+i],2)+pow(_U[5][j*_NbCol+i],2)+pow(_U[6][j*_NbCol+i],2));

    _g[0][j*_NbCol+i]=_U[2][j*_NbCol+i];
    _g[1][j*_NbCol+i]=(_U[1][j*_NbCol+i]*_U[2][j*_NbCol+i])/_U[0][j*_NbCol+i]-_U[4][j*_NbCol+i]*_U[5][j*_NbCol+i];
    _g[2][j*_NbCol+i]=pow(_U[2][j*_NbCol+i],2)/_U[0][j*_NbCol+i] + ptilde - 0.5*(pow(_U[4][j*_NbCol+i],2));
    _g[3][j*_NbCol+i]=(_U[1][j*_NbCol+i]*_U[3][j*_NbCol+i])/_U[0][j*_NbCol+i]-_U[4][j*_NbCol+i]*_U[6][j*_NbCol+i];
    _g[4][j*_NbCol+i]=(1/_U[0][j*_NbCol+i])*(_U[2][j*_NbCol+i]*_U[4][j*_NbCol+i]-_U[1][j*_NbCol+i]*_U[5][j*_NbCol+i]);
    _g[5][j*_NbCol+i]=0.;
    _g[6][j*_NbCol+i]=(1/_U[0][j*_NbCol+i])*(_U[2][j*_NbCol+i]*_U[6][j*_NbCol+i]-_U[3][j*_NbCol+i]*_U[5][j*_NbCol+i]);
    double pscal;
    pscal=(_U[5][j*_NbCol+i]/_U[0][j*_NbCol+i])*(_U[1][j*_NbCol+i]*_U[4][j*_NbCol+i]+_U[2][j*_NbCol+i]*_U[5][j*_NbCol+i]+_U[3][j*_NbCol+i]*_U[6][j*_NbCol+i]);
    _g[7][j*_NbCol+i]=_U[2][j*_NbCol+i]*_U[7][j*_NbCol+i]/_U[0][j*_NbCol+i]+ _U[2][j*_NbCol+i]*ptilde/_U[0][j*_NbCol+i] -pscal;
  }

  void Probleme::calcul_a2(double x, double y, int i, int j)
  {
    double p;
    p=(_gamma-1)*(_U[7][j*_NbCol+i]-(0.5/_U[0][j*_NbCol+i])*(pow(_U[1][j*_NbCol+i],2)+pow(_U[2][j*_NbCol+i],2)+pow(_U[3][j*_NbCol+i],2))+0.5*(pow(_U[4][j*_NbCol+i],2)+pow(_U[5][j*_NbCol+i],2)+pow(_U[6][j*_NbCol+i],2)));
    _a2[j*_NbCol+i]=(_gamma*p)/_U[0][j*_NbCol+i];
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
      _ff[k][j*_NbCol+i]=0.5*(_f[k][j*_NbCol+i]+_f[k][j*_NbCol+i+1])-max(abs(_alpha[j*_NbCol+i]),abs(_alpha[j*_NbCol+i+1]))*(_U[k][j*_NbCol+i+1]-_U[k][j*_NbCol+i]);
    }

    _ff[5][j*_NbCol+i]=0.5*(_f[5][j*_NbCol+i]+_f[5][j*_NbCol+i+1])+max(abs(_alpha[j*_NbCol+i]),abs(_alpha[j*_NbCol+i+1]))*(_U[5][j*_NbCol+i+1]-_U[5][j*_NbCol+i]);
    _ff[6][j*_NbCol+i]=0.5*(_f[6][j*_NbCol+i]+_f[6][j*_NbCol+i+1])-max(abs(_alpha[j*_NbCol+i]),abs(_alpha[j*_NbCol+i+1]))*(_U[6][j*_NbCol+i+1]-_U[6][j*_NbCol+i]);
    _ff[7][j*_NbCol+i]=0.5*(_f[7][j*_NbCol+i]+_f[7][j*_NbCol+i+1])-max(abs(_alpha[j*_NbCol+i]),abs(_alpha[j*_NbCol+i+1]))*(_U[7][j*_NbCol+i+1]-_U[7][j*_NbCol+i]);
  }

  void Probleme::flux_g(double x, double y, int i, int j)
  {
    for (int k=0; k<8; k++)
    {
      _gg[k][j*_NbCol+i]=0.5*(_g[k][j*_NbCol+i]+_g[k][(j+1)*_NbCol+i])-max(abs(_beta[j*_NbCol+i]),abs(_beta[(j+1)*_NbCol+i]))*(_U[k][(j+1)*_NbCol+i]-_U[k][j*_NbCol+i]);
    }
  }

  void Probleme::TimeIteration()
  {
    _t=0.;
    while (_t<=_tmax)
    {
      for (int j=0; j<_NbLignes; j++)
      {
        for (int i=0; i<_NbCol; i++)
        {
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
          _gg[k][(_NbLignes-1)*_NbCol+i]=0.5*(_g[k][(_NbLignes-1)*_NbCol+i]+_g[k][i])-max(abs(_beta[(_NbLignes-1)*_NbCol+i]),abs(_beta[i]))*(_U[k][i]-_U[k][(_NbLignes-1)*_NbCol+i]);
        }
      }

      for (int j=0; j<_NbLignes;j++)
      {
        for (int k=0; k<5; k++)
        {
          _ff[k][j*_NbCol+_NbCol-1]=0.5*(_f[k][j*_NbCol+_NbCol-1]+_f[k][j*_NbCol])-max(abs(_alpha[j*_NbCol+_NbCol-1]),abs(_alpha[j*_NbCol]))*(_U[k][j*_NbCol]-_U[k][j*_NbCol+_NbCol-1]);
        }

        _ff[5][j*_NbCol+_NbCol-1]=0.5*(_f[5][j*_NbCol+_NbCol-1]+_f[5][j*_NbCol])+max(abs(_alpha[j*_NbCol+_NbCol-1]),abs(_alpha[j*_NbCol]))*(_U[5][j*_NbCol]-_U[5][j*_NbCol+_NbCol-1]);
        _ff[6][j*_NbCol+_NbCol-1]=0.5*(_f[6][j*_NbCol+_NbCol-1]+_f[6][j*_NbCol])-max(abs(_alpha[j*_NbCol+_NbCol-1]),abs(_alpha[j*_NbCol]))*(_U[6][j*_NbCol]-_U[6][j*_NbCol+_NbCol-1]);
        _ff[7][j*_NbCol+_NbCol-1]=0.5*(_f[7][j*_NbCol+_NbCol-1]+_f[7][j*_NbCol])-max(abs(_alpha[j*_NbCol+_NbCol-1]),abs(_alpha[j*_NbCol]))*(_U[7][j*_NbCol]-_U[7][j*_NbCol+_NbCol-1]);
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

       _Dt=0.45*_Dx/maxeigenvalues;

      for (int j=1; j<_NbLignes; j++)
      {
        for (int i=1; i<_NbCol; i++)
        {
          for(int indice=0; indice<8; indice++)
          {
            _Uapres[indice][j*_NbCol+i]=_U[indice][j*_NbCol+i]-_Dt*((1/_Dx)*(_ff[indice][j*_NbCol+i]-_ff[indice][j*_NbCol+i-1])+(1/_Dy)*(_gg[indice][j*_NbCol+i]-_gg[indice][(j-1)*_NbCol+i]));
          }
        }
      }

      for (int i=1; i<_NbCol; i++)
      {
        for(int indice=0; indice<8; indice++)
        {
          _Uapres[indice][i]=_U[indice][i]-_Dt*((1/_Dx)*(_ff[indice][i]-_ff[indice][i-1])+(1/_Dy)*(_gg[indice][i]-_gg[indice][(_NbLignes-1)*_NbCol+i]));

        }
      }

      for (int j=1; j<_NbLignes; j++)
      {
        for(int indice=0; indice<8; indice++)
        {
          _Uapres[indice][j*_NbCol]=_U[indice][j*_NbCol]-_Dt*((1/_Dx)*(_ff[indice][j*_NbCol]-_ff[indice][j*_NbCol+_NbCol-1])+(1/_Dy)*(_gg[indice][j*_NbCol]-_gg[indice][(j-1)*_NbCol]));
        }
      }

      for(int indice=0; indice<8; indice++)
      {
        _Uapres[indice][0]=_U[indice][0]-_Dt*((1/_Dx)*(_ff[indice][0]-_ff[indice][_NbCol-1])+(1/_Dy)*(_gg[indice][0]-_gg[indice][_NbLignes-1]));
      }

      _t=_t+_Dt;

    //  std::cout<<"t  "<< _t <<endl;

      _U=_Uapres;
      Probleme::SaveIeration("result" + std::to_string(_t));

    }
  //  Probleme::SaveIeration("result" + std::to_string(_t));
  }

  void Probleme::SaveIeration(std::string fichier)
  {
    std::ofstream mon_fluxP;
    mon_fluxP.open(fichier + "_p.txt", std::ios::out);
    mon_fluxP << "# champ de pressions sur un maillage carré" << std::endl;

    for (int i=0; i<_NbLignes; i++)
    {
      for (int j=0; j<_NbCol; j++)
      {
          mon_fluxP << i*_Dx << "    " << j*_Dy << "     " <<_a2[j*_NbCol+i] <<"      "<< _b[0][j*_NbCol+i]<< "      "<<  _b[1][j*_NbCol+i]<< "      "<<  _b[2][j*_NbCol+i]<< "     "<< _c1car[j*_NbCol+i]<< "     " <<_c2car[j*_NbCol+i];
          mon_fluxP << "     "<< _alpha[j*_NbCol+i]<< "         "<<_beta[j*_NbCol+i]<<endl;

        //mon_fluxP << i*_Dx << " " << j*_Dy << " " <<(_gamma-1)*(_U[7][j*_NbCol+i]-(0.5/_U[0][j*_NbCol+i])*(pow(_U[1][j*_NbCol+i],2)+pow(_U[2][j*_NbCol+i],2)+pow(_U[3][j*_NbCol+i],2))+0.5*(pow(_U[4][j*_NbCol+i],2)+pow(_U[5][j*_NbCol+i],2)+pow(_U[6][j*_NbCol+i],2))) << std::endl;
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