
class Springs
{
  final int N_SPRINGS = 3*N_TRIANGLES + 6*(N_TRIANGLES-1);
    // 3*N_TRIANGLES: In each "hirozntal" layer, a triangle has three edges.
    // 6*(N_TRIANGLES-1): Each vertex in the triangle has two springs
    //                    connected to two vertices in the lower layer.

  SpringElement[] element = new SpringElement[N_SPRINGS];
    //
    //
    //                H
    //              x   x  "self triangle layer"
    //            x       x
    //          H  x  x  x  H
    //           \         /
    //            \       /
    //          L .\. . ./. L
    //            . \   / .
    //              .\ /.   "lower triangle layer"
    //                L
    //
    //
    //             U     U     U
    //        .   / \   / \   /
    //         \ /   \ /   \ /
    //        . H x x H x x H .
    //         / \   / \   / \
    //        .   \ /   \ /   \
    //             L     L     L
    //
    // each edge in the layer 'H' has two
    // springs connecting to a vertex
    // in the lower layer 'L'.

  Springs(float characteristicPeriod)
  {
    float omega = PI*2 / characteristicPeriod;
    float spc = PARTICLE_MASS * omega * omega;
              // spc = spring constant:  omega^2 = spc / mass

    int sCtr = 0; // spring counter
    for (int tl=0; tl<N_TRIANGLES; tl++) { // for "hirizontal" triangles.
      int pId0 = particles.id(tl,0); // 1st vertex in the triangle
      int pId1 = particles.id(tl,1); // 2nd
      int pId2 = particles.id(tl,2); // 3rd

      register(spc, sCtr++, pId0, pId1);
      register(spc, sCtr++, pId1, pId2);
      register(spc, sCtr++, pId2, pId0);
    }
    for (int tl=1; tl<N_TRIANGLES; tl++) { // skip the lowest layer.
      for (int me=0; me<3; me++) {
        //
        // when tl=even
        //
        // each vertex in the layer 'S' has two
        // springs connecting to two vertices
        // in the upper and lower layers 'U' and 'L'.
        //
        //    vertexId (0, 1, 2) in each layer.
        //
        //             0     1     2    upper layer
        //        .   / \   / \   /
        //         \ /   \ /   \ /
        //        . 0 x x 1 x x 2 .    tl (even layer)
        //         / \   / \   / \
        //        .   \ /   \ /   \
        //             0     1     2    lower layer
        //
        // Connection table. me and its counterparts.
        //
        //       same layer    lower layer
        //            /  \     /   \
        //   me  |  m1   m2   l1   l2
        //   ----+--------------------
        //    0  |   1    2    0    2
        //    1  |   2    0    1    0
        //    2  |   0    1    2    1
        //       +---------------------
        //       |  k1   k2   me   k2
        //
        //
        // when tl=odd
        //
        //       same layer    lower layer
        //            /  \     /   \
        //   me  |  m1   m2   l1   l2
        //   ----+--------------------
        //    0  |   1    2    0    1
        //    1  |   2    0    1    2
        //    2  |   0    1    2    0
        //       +--------------------
        //       |  k1   k2   me   k1

        int k1 = (me+1) % 3;
        int k2 = (me+2) % 3;

        int myPid = particles.id(tl,me);

        if ( tl%2==0 ) {
          register(spc, sCtr++, myPid, particles.id(tl-1,me)); // lower layer
          register(spc, sCtr++, myPid, particles.id(tl-1,k2)); // lower layer
        }
        else {
          register(spc, sCtr++, myPid, particles.id(tl-1,me)); // lower layer
          register(spc, sCtr++, myPid, particles.id(tl-1,k1)); // lower layer
        }
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

    particles.sixSpringsAppend(alpha, springId);
    particles.sixSpringsAppend(beta,  springId);
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

}
