package Lenia;

import HAL.Rand;
import HAL.Util;

public class ExampleStochastic2Player extends StochasticNPlayer {

    
    public static final double[][] P = new double[][]{{1.0,0.5},{0.9,0.4}};
    //kernel size
    public static final double[][] Rstar = new double[][]{{45.0,5.0},{5.0,45.0}};

    // Colors
    public static final int ImmuneColor = Util.RGB256(153,0,153);
    public static final int TumorColor = Util.RGB256(32,32,32);


    // Parameters
    public double gamma = 5.0;
    public double L = 0.08; //Allee threshold
    public double b = 12.0; //predation rate
    public double C = 1.0; //carrying capacity
    public double d = 1.0, g = 1.5;

    // Kernel function
    @Override
    public double K(int i, int j, double r) {
        if (r <= Rstar[i][j]) {
            return 1;
        } else {
            return 0;
        }
    }

    // Growth function
    @Override
    public double G(int iPlayer, double[] U) {
        double[] growth = new double[N];
        growth[0] = gamma * U[0] * (U[0] - L) * (C - U[0]) - (b * U[0] * U[1]);
        growth[1] = g * (b * U[0] * U[1]) - d * U[1];
        return growth[iPlayer];
    }

    // Constructor
    public ExampleStochastic2Player(int sideLength, int nPlayers, double dt, String filename, int scalefactor, boolean[] vis_options, int ClipMax) {
        super(filename, sideLength, dt, nPlayers, scalefactor, vis_options, ClipMax);
        this.colors = new int[] {TumorColor, ImmuneColor};
        SetGrowthBars(ClipMax);
    }

    // Set initial conditions
    public void SetInitialCondition(int CarryingCapacity) {
        Rand rn = new Rand();
        int[] hood = Util.CircleHood(true, 12);
        int H = MapHood(hood, xDim/2, yDim/2);

        for (int i = 0; i < H; i++) {
            this.Set(0, hood[i], 1);
        }

        double dummy = 0;
        for (int i = 0; i < length; i++) {
            if (this.Get(0, i) == 0) {
                double a = rn.Double();
                if (a < 0.021) {
                    this.Set(1, i, 0.5);
                    dummy = dummy + 1;
                }
            }
        }
        System.out.println(dummy);
    }

    // Set growth scale bars
    public void SetGrowthBars(int ClipMax) {
        for (int i = 0; i < N; i++) {
            double min = 0;
            double max = 0;
            for (double u1 = 0; u1 <= ClipMax; u1 += 0.001) {
                double u2 = ClipMax - u1;
                double[] u = {u1, u2};
                min = Math.min(min, this.G(i, u));
                max = Math.max(max, this.G(i, u));
            }
            SetGrowthScale(i, min, max);
        }
    }

    public static void main(String[] args) {
        // Define the side length of the grid
        int side_length = 64;
        // Define the scale factor for visualization
        int scalefactor = 4;
        // Define the number of players
        int N = 2;
        // Define the maximum value for clipping
        int ClipMax = 1;
        // Define the time step for simulation
        double deltaT = 0.02;
        // Define the time step for drawing
        double drawDeltaT = 0.5;
        // Define the name of the file for saving data
        String name = "data/ExampleStochastic2Player";

        // Create an instance of Example2Player
        ExampleStochastic2Player model = new ExampleStochastic2Player(side_length, N, deltaT, name, scalefactor, new boolean[]{true,true,true,true}, ClipMax);
        // Disable displaying average on plot
        model.DISPLAY_AVG_ON_PLOT = false;
        // Disable displaying time on plot
        model.DISPLAY_TIME_ON_PLOT = false;
        // Disable displaying scale bars on plot
        model.DISPLAY_SCALE_BARS_ON_PLOT = false;

        // Setup GIF saving
        model.SetupGifSaving();
        // Set the initial condition
        model.SetInitialCondition(ClipMax);
        // Draw the initial state
        model.Draw();

        // Run the simulation until the specified time
        while (model.GetTime() < 20.0) {
            // Output data
            model.Output();
            // Update the simulation
            model.Update();
            // Draw the state at regular intervals
            if (model.GetTick() % (int) (drawDeltaT / deltaT) == 0) {
                model.Draw();
            }
        }
        // Close the model
        model.Close();
    }
}
