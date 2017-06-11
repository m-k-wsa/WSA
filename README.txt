++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++                                                      ++
++       WSA: Weakly-hard Schedulability Analyzer       ++
++                                                      ++
++            Author: xxxxxxxxxxxx                      ++
++          xxxxxxxxxxxxxxxxxxxxxxxx                    ++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This is the weakly hard analysis solution for the
industrial challenge of WATERS 2015 is proposed by
Thales Research & Technology France.

  https://ecrts.eit.uni-kl.de/forum/viewtopic.php?f=32&t=86

For the 2nd challenge, we provide formal analysis for the
frequency of possible temporal violations in the system.

As a result, we list in the following table the maximum
number (m*) of temporal violations of task T2, within an
arbitrary sequence of K its activations.

      m*  |  K
    ++++++++++++
      1   |  2
    ++++++++++++
      2   |  5
    ++++++++++++
      4   |  10
    ++++++++++++
      6   |  15
    ++++++++++++

To check the m-K for this system:
./src/jitter-demo m K
