package Lenia;
import FFT.FFTGrid;

import HAL.Rand;
import HAL.GridsAndAgents.Grid2Ddouble;
import static HAL.Util.*;

public class StochasticNPlayer extends LeniaNPlayer {

    public static double u0 = 1; // initial density of cells
    public static final double[][] P = new double[][]{{1.0,0.5},{0.9,0.4}}; // payoff matrix

    public StochasticNPlayer(String filename, int sideLength, double deltaT, int nPlayers, int scalefactor, int ClipMax) {
        super(filename,sideLength,deltaT,nPlayers,scalefactor, new boolean[]{false,false,false,false},ClipMax);
        this.KernelFunction = this::K;
        this.GrowthFunction = this::G;
        this.N = nPlayers;
        RecalcKernels();
        this.C = ClipMax;
        SetClip(0, C);
    }

    public StochasticNPlayer(String filename, int sideLength, double deltaT, int nPlayers, int scalefactor, boolean[] vis_options, int ClipMax) {
        super(filename,sideLength,deltaT,nPlayers,scalefactor, vis_options,ClipMax);
        this.KernelFunction = this::K;
        this.GrowthFunction = this::G;
        this.N = nPlayers;
        rn = new Rand();
        RecalcKernels();
        // Ks = kSums[0][0];
        C = ClipMax;
        SetClip(0, C);
    }

    // design your kernal function
    @Override
    public double K(int i,int j,double r) {
        if (r <= Rstar[i][j]) {
            return 1;
        } else {
            return 0;
        }
    }

    @Override
    public double G(int iPlayer, double[]U) {
        double[] pi = new double[N];
        double phi = 0.0;
        double[] growth = new double[N];

        for (int j = 0; j < N; j++) {
            pi[j] = 0.0;
            for (int k = 0; k < N; k++) {
                pi[j] += P[j][k]*U[k];
            }
            phi += U[j]/C * pi[j];
        }

        for (int j = 0; j < N; j++) {
            growth[j] = U[j]/C * (pi[j] - phi);
        }
        return growth[iPlayer]*C;
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
        for (int iPlayer = 0; iPlayer < N; iPlayer++) {
            for (int k = 0; k < A[iPlayer].length; k++) {

                SetUiScratch(iPlayer,k);

                // Use binomial distribution to determine the number of cells to add or remove
                double rate = this.GrowthFunction.Eval(iPlayer,Uscratch)*deltaT;
                if (rate < 0) {
                    // If the growth rate is negative, remove cells
                    double rateNeg = -rate;
                    int Ncells = rn.Binomial((int) A[iPlayer].Get(k), rateNeg);
                    A[iPlayer].Add(k, -Ncells); // Remove cells
                } else {
                    // If the growth rate is positive, add cells
                    int Ncells = rn.Binomial((int) (C - A[iPlayer].Get(k)), rate);
                    A[iPlayer].Add(k, Ncells); // Add cells
                }

                // clip:
                A[iPlayer].Set(k, Bound(A[iPlayer].Get(k),CLIP[0],CLIP[1]));
            }
        }
        this.tick++;
        this.t += deltaT;
    }

    // loop over every location in the grid, set to a random cell type (A, B, ...)
    public void SetInitialCondition(int CarryingCapacity) {
        Rand rn = new Rand();
        for (int x = 0; x < this.xDim; x++) {
            for (int y = 0; y < this.yDim; y++) {
                // random cell type (between 0 and the number of cell types)
                int i = rn.Int(N);

                double x1 = (i == 0) ? 1 : 0;
                double x2 = 1 - x1;

                if (CarryingCapacity > 1) {
                    x1 = rn.Int(CarryingCapacity);
                    x2 = CarryingCapacity - x1;
                }

                for (int j = 0; j < N; j++) {
                    this.Set(0,x,y,x1);
                    this.Set(1,x,y,x2);
                }
            }
        }
    }

    public void SetGrowthBars(){
        for (int i = 0; i < N; i++) {
            double min = 0;
            double max = 0;
            for (double u1 = 0; u1 <= 1; u1+=0.001) {
                double u2 = 1.0-u1;
                double[]u = {u1,u2};

                min = Math.min(min, this.G(i,u));
                max = Math.max(max, this.G(i,u));
            }
            SetGrowthScale(i, min, max);
        }
    }

    public static void main(String[] args) {

        int side_length = 64; // length of the side of simulation domain
        int scalefactor = 5; // to scale the drawing size

        double deltaT = 1; // time step

        // double[] Rvec = new double[]{1000};

        int CarryingCapacity = 1;
        int Nplayers = 2;


        String filename = "data/Stochastic/stochasticNPlayer";

        StochasticNPlayer model=new StochasticNPlayer(filename,side_length,deltaT, Nplayers,scalefactor,new boolean[]{true,true,true,true}, CarryingCapacity);
        model.SetupGifSaving();
        model.SetInitialCondition(CarryingCapacity);

        model.SetGrowthBars();
        model.Draw();

        // model.SetPause(0.05); // seconds

        // set up an initial condition (1 cell in the center of the domain)



        while (model.GetTime() < 300.0) {
            model.Output(); // output must be called before Update()
            model.Update(); // increment model.time
            if (model.GetTick() % (int)(1.0/deltaT) == 0) {
                model.Draw(); // draw the current state of the model
            }
        }
        model.Close();


        return;

    }


}

