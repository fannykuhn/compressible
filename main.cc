#include "VF.h"

int main()
{
  Probleme Pb1(200,200,1);
  Pb1.initialize_u();
  Pb1.Time_iteration_MUSCL();

  return 0;
}
