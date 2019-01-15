#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>

class Probleme
{
private:
  int _NbLignes, _NbCol;
  double _Dx, _Dy, _Dt, _tmax, _t,_Lx,_Ly, _gamma;
  std::vector<std::vector<double>> _f, _g, _U, _b, _ff, _gg, _Uapres;
  std::vector<double> _a2, _c1car, _c2car, _alpha, _beta, _p, _ptilde, _pscal;

public:

  Probleme(int NbLignes, int NbCol);

  ~Probleme();

  double rho(double x, double y);

  double u1(double x, double y);

  double u2(double x, double y);

  double u3(double x, double y);

  double B1(double x, double y);

  double B2(double x, double y);

  double B3(double x, double y);

  double p(double x, double y);

  double E(double x, double y);

  void initialize_u();

  void calcul_p(double x, double y, int i, int j);

  void calcul_p_tilde(double x, double y, int i, int j);

  void calcul_pscal(double x, double y, int i, int j);

  void calcul_f(double x, double y, int i, int j);

  void calcul_g(double x, double y, int i, int j);

  void calcul_a2(double x, double y, int i, int j);

  void calcul_b(double x, double y, int i, int j);

  void calcul_ck(double x, double y, int i, int j, int k);

  void calcul_alpha_beta(double x, double y, int i, int j);

  void flux_f(double x, double y, int i, int j);

  void flux_g(double x, double y, int i, int j);

  void TimeIteration();

  void SaveIeration(std::string fichier);

};
