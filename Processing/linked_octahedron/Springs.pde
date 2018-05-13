
class Springs
{
  final int N_SPRINGS = 3*N_TRIANGLES + 6*(N_TRIANGLES-1);
    // 3*N_TRIANGLES: In each "hirozntal" layer, a triangle has three edges.
    // 6*(N_TRIANGLES-1): Each vertex in the triangle has two springs
    //                    connected to two vertices in the lower layer.
    // 2: Tension force applying to the two ends.
    //
    //           triangle
    //              tl=0 tl=1 tl=2            tl=N_TRIANGLES-1
    //               /   /   /               /
    //      x - - - o===o===o= ... =o===o===o - - - x
    //

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

  Springs()
  {
    float spc = SPRING_CONST;

    int sCtr = 0; // spring counter
    for (int l=0; l<N_TRIANGLES; l++) { // for "hirizontal" triangles.
      int pId0 = particles.id(l,0); // 1st vertex in the triangle
      int pId1 = particles.id(l,1); // 2nd
      int pId2 = particles.id(l,2); // 3rd

      register(spc, sCtr++, pId0, pId1);
      register(spc, sCtr++, pId1, pId2);
      register(spc, sCtr++, pId2, pId0);
    }
    for (int l=1; l<N_TRIANGLES; l++) { // skip the lowest layer.
      for (int me=0; me<3; me++) {
        //
        // when l=even
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
        //        . 0 x x 1 x x 2 .    l (even layer)
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
        // when l=odd
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

        int myPid = particles.id(l,me);

        if ( l%2==0 ) {
          register(spc, sCtr++, myPid, particles.id(l-1,me)); // lower layer
          register(spc, sCtr++, myPid, particles.id(l-1,k2)); // lower layer
        }
        else {
          register(spc, sCtr++, myPid, particles.id(l-1,me)); // lower layer
          register(spc, sCtr++, myPid, particles.id(l-1,k1)); // lower layer
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
