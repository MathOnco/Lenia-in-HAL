package Lenia;

import HAL.Rand;
import HAL.GridsAndAgents.Grid2Ddouble;
import static HAL.Util.*;
import FFT.FFTGrid;

public class Stochastic1Player extends Lenia1Player {

    public static double Rstar = 2;
    public static double u0 = 1; // initial density of cells
    private int N = 1;

    public Stochastic1Player(String filename, int sideLength, double deltaT, int scalefactor, int ClipMax) {
        super(filename,sideLength,deltaT,scalefactor, new boolean[]{false,false,false,false},ClipMax);
        this.KernelFunction = this::K;
        this.GrowthFunction = this::G;
        RecalcKernels();
        // Ks = kSums[0][0];
        this.C = ClipMax;
        SetClip(0, C);
    }

    public Stochastic1Player(String filename, int sideLength, double deltaT, int scalefactor, boolean[] vis_options, int ClipMax) {
        super(filename,sideLength,deltaT,scalefactor, vis_options,ClipMax);
        this.KernelFunction = this::K;
        this.GrowthFunction = this::G;
        rn = new Rand();
        RecalcKernels();
        // Ks = kSums[0][0];
        C = ClipMax;
        SetClip(0, C);
    }

    // default K
    @Override
    public double K(double r) {
        if (r<Rstar) {
            return 1.0;
        } else {
            return 0.0;
        }
    }

    // default G:
    @Override
    public double G(double u) {
        return 0.1*u*(C - u);
    }

    public double U(int x, int y) {
        return Uij(0,0,x,y);
    }

    public double U(int k) {
        return Uij(0,0,k);
    }

    @Override
    public void Update() {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                Grid2Ddouble convField = this.A[j];
                FFTGrid currFftField = this.fftFields[i][j];
                currFftField.SetGrid(convField);
                currFftField.fft2();
                currFftField.ElementwiseMultiplication(fftKernels[i][j]);
                currFftField.ifft2();
            }
        }
        // update each pixel in A[0]:
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < A[i].length; k++) {
                
                // use binomial to add number of cells:
                double rate = this.GrowthFunction.Eval(U(k))*deltaT;
                if (rate < 0 ){
                    double rateNeg = -rate;

                    int Ncells = rn.Binomial((int)(A[i].Get(k)), rateNeg);
                    A[i].Add(k, -Ncells);
                }else{
                    int Ncells = rn.Binomial((int)(C-A[i].Get(k)), rate);
                    A[i].Add(k, Ncells);
                }

//                int Ncells = rn.Binomial((int)A[i].Get(k), rate); // how many cells to add here?



                // Ncells \in [0, 1, 2, ... ClipMax]


                
                // previous implementation (using Poisson probability) is now depreciated:
                // add cells at poisson rate:
                //A[i].Add(k, NPoissonEvents(this.GrowthFunction.Eval(U(k)), deltaT));

                // clip:
                A[i].Set(k, Bound(A[i].Get(k),CLIP[0],CLIP[1]));
            }
        }
        this.tick++;
        this.t += deltaT;
    }

    public void SetGrowthBars() {
        double min = 0;
        double max = 0;
        for (int u = 0; u < CLIP[1]; u+=1) {
            double g = this.GrowthFunction.Eval((double)u);
            min = (min < g) ? min : g;
            max = (max > g) ? max : g;
        }
        SetGrowthScale(min, max);
    }
}

