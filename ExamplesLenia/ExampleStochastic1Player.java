package ExamplesLenia;

import HAL.Util;
import Lenia.Stochastic1Player;

public class ExampleStochastic1Player extends Stochastic1Player {

    // Constants
    public static double Rstar = 5; //Neighborhood radius
    public static double gamma = 5.0; //growth rate
    public static double L = 0.08; //Allee threshold
    public static double C = 1.0; //Carrying capacity

    // Kernel function
    @Override
    public double K(double r) {
        if (r <= Rstar) {
            return 1.0;
        } else {
            return 0.0;
        }
    }

    // Growth function
    @Override
    public double G(double u) {
        return gamma * u * (u - L) * (C - u);
    }

    // Constructor
    public ExampleStochastic1Player(String filename, int side_length, double dt, int scalefactor, int ClipMax) {
        super(filename, side_length, dt, scalefactor, new boolean[]{true, true, true, true}, ClipMax);
    }

    // Main method
    public static void main(String[] args) {
        int side_length = 64; // Length of the simulation domain side
        int scalefactor = 5; // Scale factor for drawing size
        int ClipMax = 1; // Carrying capacity

        double deltaT = 0.02; // Timestep
        double drawDeltaT = 0.5; // Timesteps to draw grid at
        String name = "data/test"; // Folder for saving data

        // Create an instance of ExampleStochastic1Player with the specified parameters
        ExampleStochastic1Player model = new ExampleStochastic1Player(name, side_length, deltaT, scalefactor, ClipMax);
        model.SetupGifSaving(); // Set up saving GIF files
        model.DISPLAY_SCALE_BARS_ON_PLOT = false; // Disable displaying scale bars on the plot
        model.DISPLAY_AVG_ON_PLOT = false; // Disable displaying average on the plot
        model.DISPLAY_TIME_ON_PLOT = false; // Disable displaying time on the plot

        model.SetInitialCondition(ClipMax); // Set the initial condition with the specified carrying capacity
        model.Draw(); // Draw the initial state

        while (model.GetTime() < 20.0) { // Loop until the time reaches 20.0
            model.Output(); // Output the current state
            model.Update(); // Update the state
            if (model.GetTick() % (int) (drawDeltaT / deltaT) == 0) { // Check if it's time to draw the grid
                model.Draw(); // Draw the grid
            }
        }
        model.Close(); // Close the model
    }

    // Set initial condition
    public void SetInitialCondition(double CarryingCapacity) {
        int[] hood = Util.CircleHood(true, 12);
        int H = MapHood(hood, xDim / 2, yDim / 2);

        for (int i = 0; i < H; i++) {
            this.Set(hood[i], CarryingCapacity);
        }
    }
}
