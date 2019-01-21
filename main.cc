#include "VF.h"

int main()
{
  Probleme Pb1(50,50,1);
  Pb1.initialize_u();
  Pb1.TimeIteration_ordre1();

  return 0;
}
