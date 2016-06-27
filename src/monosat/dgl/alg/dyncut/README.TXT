################################################################################
#                                                                              #
#           Software for Efficiently Solving Multi-Label MRFs/CRFs             #
#                               Version 1.0                                    #
#                                                                              #
#      Karteek Alahari          Pushmeet Kohli           Philip H. S. Torr     #
#    (karteek.alahari@*)    (pkohli@microsoft.com)        (philiptorr@*)       #
#                                                                              #
# * brookes.ac.uk                                                              #
################################################################################


1. INTRODUCTION

   Our software uses the following publicly available implementations:
   a. Dynamic Graph Cuts, version 2
   b. Tree-reweighted max-product message passing algorithm, version 1.3.

   It has been tested on the Fedora Core 4 Linux distribution, using g++ version
   4.0. Please report any issues to Karteek Alahari.

   If you use this software for research purposes, you should cite the following
   papers in any resulting publication.

   "Reduce, Reuse & Recycle: Efficiently Solving Multi-Label MRFs"
   Karteek Alahari, Pushmeet Kohli, Philip H. S. Torr
   In IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2008.

   "Efficiently Solving Dynamic Markov Random Fields Using Graph Cuts"
   Pushmeet Kohli, Philip H. S. Torr
   In International Conference on Computer Vision (ICCV), 2005.

   "Convergent Tree-reweighted Message Passing for Energy Minimization"
   Vladimir Kolmogorov
   In IEEE Trans. Pattern Analysis and Machine Intelligence (PAMI), Oct 2006.

   "An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy
   Minimization in Vision"
   Yuri Boykov and Vladimir Kolmogorov
   In IEEE Trans. Pattern Analysis and Machine Intelligence (PAMI), Sep 2004.

################################################################################

2. LICENSE & DISCLAIMER

   Various copyrights apply to this release.

   ------------------------------- PART I --------------------------------------
   Copyright 2009 Karteek Alahari, Pushmeet Kohli, Philip H. S. Torr
   Oxford Brookes University

   This software can be used for research purposes only.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
   -----------------------------------------------------------------------------

   ------------------------------- PART II -------------------------------------
   Dynamic Graph Cuts, version 2
   Copyright 2005 Pushmeet Kohli (pushmeet.kohli@brookes.ac.uk), Philip HS Torr
   (philiptorr@brookes.ac.uk).
   [For files block.h graph.h]

   This software library implements the dynamic maxflow algorithm described in:

   Efficiently Solving Dynamic Markov Random Fields using Graph Cuts
   Pushmeet Kohli and Philip H. S. Torr
   In the Tenth IEEE International Conference on Computer Vision (ICCV 2005).

   The algorithm uses the maxflow algorithm code described in: 

   An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy
   Minimization in Vision, Yuri Boykov and Vladimir Kolmogorov.
   In IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI),
   September 2004.

   The source code also comes under the following license:

   Copyright 2001 Vladimir Kolmogorov (vnk@adastral.ucl.ac.uk), Yuri Boykov
   (yuri@csd.uwo.ca).

   This program is free software; you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2 of the License, or (at your option) any later
   version.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
   details.

   You should have received a copy of the GNU General Public License along with
   this program; if not, write to the Free Software Foundation, Inc., 59 Temple
   Place, Suite 330, Boston, MA  02111-1307  USA
   -----------------------------------------------------------------------------

   ------------------------------ PART III -------------------------------------
   Tree-reweighted max-product message passing algorithm (TRW-S), version 1.3
   Written by Vladimir Kolmogorov (vnk@microsoft.com), 2005.
   (c) Microsoft Corporation. All rights reserved.
   [For files MRFEnergy.h typePotts.h TRWBP.h]

   This software implements two algorithms for minimizing energy functions of
   discrete variables with unary and pairwise terms. They are max-product belief
   propagation (BP, Pearl'88) and sequential tree-reweighted max-product message
   passing (TRW-S, Kolmogorov'05). The sofware is covered by the Microsoft
   Research Shared Source license agreement (MSR-SSLA), and is available for
   non-commercial purposes, subject to the restrictions in MSR-SSLA. Please see
   the terms in MSR-SSLA.TXT
   -----------------------------------------------------------------------------

################################################################################

3. EXAMPLE

   We explain the usage with the stereo computation problem as an example.
   Tsukuba stereo pair (courtesy of the University of Tsukuba) and the
   corresponding unary costs are available in images/

   a. Write code to initialize the energy parameters (unary and pairwise terms)
      in image::image(char*, int). See the comments in image.cpp for the
      variables to be initialised.

      // Example code to add in image::image(char*, int)
      int truncation=1, lambda=20;
      int i, j, k, l, ti, tj;
      int pair_id;
      FILE* fp;

      eng = new Energy(num_labels, truncation);
      eng->nvar = num_nodes;
      eng->npair= 2*width*height - (width+height);
      eng->allocateMemory();

      fp = fopen("images/datacost.txt", "r");
      for(k=0; k<num_labels; k++)
      {
           for(i=0; i<height; i++)
           for(j=0; j<width; j++)
		fscanf(fp, "%d", &(eng->unaryCost[i*width+j][k]));
      }
      fclose(fp);

      pair_id = 0;
      for(i=0; i<height; i++)
      for(j=1; j<width; j++)
      {
           k = i*width + j;
           eng->pairCost[pair_id] = lambda;
           eng->pairIndex[pair_id][0] = k;
           eng->pairIndex[pair_id][1] = k-1;
           pair_id++;
      }
      for(i=1; i<height; i++)
      for(j=0; j<width; j++)
      {
           k = i*width + j;
           eng->pairCost[pair_id] = lambda;
           eng->pairIndex[pair_id][0] = k;
           eng->pairIndex[pair_id][1] = k-width;
           pair_id++;
      }

   b. Create an image object using a ppm file:
      image* Image = image("file.ppm", number_of_labels);
      pgm files are also supported. Change image::image(char*, int) accordingly.

      // Example
      Image = new image("images/tsukuba_l.ppm", num_labels);

   c. Solve the energy function
      Image -> kovtun();     // Compute the partially optimal solution only
      Image -> kovtun(true); // Compute the partially optimal solution and 
                             // solve the remaining nodes using alpha expansion.
      Image -> trw();        // Compute standard TRW/BP solutions.
      Image -> trw(true);    // Compute the partially optimal solution and 
                             // solve the remaining nodes using TRW/BP.

   d. Write your own method to save the result, eg. image::save_disparity();

      // Example
      void image::save_disparity(const char *file)
      {
           int i, j;
           unsigned char lbl;
           unsigned char disp;
           for(i=0; i<height; i++)
           for(j=0; j<width; j++)
           {
                lbl = label_map[i][j];
                if(lbl != NOLABEL)
                {
                     disp = (unsigned char)(lbl * 254.0/float((num_labels-1)));
                     buf[i][3*j]  = disp;
                     buf[i][3*j+1]= disp;
                     buf[i][3*j+2]= disp;
                }
                else
                {
                     buf[i][3*j]  = 255;
                     buf[i][3*j+1]= 0;
                     buf[i][3*j+2]= 0;
                }
           }
           writePPM(file);
      }

################################################################################

