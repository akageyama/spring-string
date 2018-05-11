
class Springs
{
  final int N_SPRINGS = N_PAIRS + 2*(N_PAIRS-1);

            //
            //   o---------o  upper layer
            //   |         |
            //   |         |
            //   |         |
            //   o---------o  self layer
            //   |         |
            //   |         |
            //   |         |
            //   o---------o  lower layer

  SpringElement[] element = new SpringElement[N_SPRINGS];

  Springs(float characteristicPeriod)
  {
    float omega = PI*2 / characteristicPeriod;
    float spc = PARTICLE_MASS * omega * omega;
              // spc = spring constant:  omega^2 = spc / mass

    int sCtr = 0; // spring counter
    for (int l=0; l<N_PAIRS; l++) { // ladder
      int p0 = particles.id(l,0); // 1st vertex in the triangle
      int p1 = particles.id(l,1); // 2nd
      register(spc, sCtr++, p0, p1);
    }
    for (int l=1; l<N_PAIRS; l++) { // skip the lowest layer.
      for (int j=0; j<2; j++) {
        int p0 = particles.id(l,j);
        int p2 = particles.id(l-1,j);
        register(spc, sCtr++, p0, p2);
      }
    }
  }


  void register(float springConst, int springId,
                        int alpha, int beta)
  {
    //
    // ids of particles on the both ends
    //           alpha         beta
    //             \           /
    //              O=========O
    //
    element[springId] = new SpringElement(springConst,alpha,beta);

    particles.connectedSpringsAppend(alpha, springId);
    particles.connectedSpringsAppend(beta,  springId);
  }


  void display(float[] posx,
               float[] posy,
               float[] posz)
  {
    stroke(150, 100, 70);

    for (int s=0; s<N_SPRINGS; s++) {
      element[s].display(posx, posy, posz);
    }
  }

  float energy(float[] posx,
               float[] posy,
               float[] posz)
  {
    float sum = 0.0;
    for (int s=0; s<N_SPRINGS; s++) {
      sum += element[s].energy(posx,posy,posz);
    }
    return sum;
  }


}
