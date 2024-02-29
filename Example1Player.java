package Lenia;

import HAL.Rand;
import HAL.Util;

public class Example1Player extends Lenia1Player {

    //Lenia1Player model;
    public static double Rstar = 45;
    public static double gamma = 1.0;
    public static double ClipMax = 1;
    public double a = 5;
    public double K = ClipMax, L = 0.23*K;
//    public double b = 10.0, alpha = 1.0, h = 0.0;
    // design your kernal function
    @Override
    public double K(double r) {
        if (r<Rstar) {
            return 1.0;
        } else {
            return 0.0;
        }
    }

    // G is the growth function of the tumor, which is a function of density
    @Override
    public double G(double u) {
//        return a*u*(K-u);
            return a * u * (u - L) * (K - u);
    }

    // create a constructor for the Gompertz class, but just inherit the constructor from Lenia1Player
    public Example1Player(String filename, int side_length, double dt, int scalefactor, double ClipMax) {
        super(filename,side_length, dt, scalefactor, new boolean[]{false, false, false, false}, ClipMax);
//        super(filename,side_length, dt, scalefactor, new boolean[]{true, true, true, true}, ClipMax);
    }


    public static void main(String[] args) {

        int side_length = 64; // length of the side of simulation domain
        int scalefactor = 5; // to scale the drawing size

        double deltaT = 0.02; // timestep
        double drawDeltaT = 0.1; // how many timesteps to draw grid at
//        double ClipMax = 100;
//        int ClipMax = 100;

//        for (int R = 1; R < 20; R++) {

//            String filename = "data/test";
//            System.out.println(R);

        String name = "data/test"; // folder must exist
//            String name = filenameTi.concat(Integer.toString(R));

        Example1Player model = new Example1Player(name, side_length, deltaT, scalefactor, ClipMax);
//        model.CLIP = new double[]{0,1};
//        model.SetupGifSaving();
////        model.DISPLAY_SCALE_BARS_ON_PLOT = false;
//        model.DISPLAY_AVG_ON_PLOT = false;
//        model.DISPLAY_TIME_ON_PLOT = false;

//        model.SetPause(0.1); // seconds

            // set up an initial condition (1 cell in the center of the domain)
//        model.Set(model.xDim/2,model.yDim/2,1);
        model.SetInitialCondition(ClipMax);

//        model.Draw();
//        model.SetupGifSaving();

        while (model.GetTime() < 20.0) {
            model.Output(); // output must be called before Update()
            model.Update(); // increment model.time
//            if (model.GetTick() % (int) (drawDeltaT / deltaT) == 0) {
//                model.Draw(); // draw the current state of the model
//                System.out.println(model.A[0].GetAvg());
//            }
        }
//        model.Close();

    }
//    }

    public void SetInitialCondition(double CarryingCapacity) {

        int[] hood = Util.CircleHood(true, 12);
        int H = MapHood(hood, xDim / 2, yDim / 2);

//        for (int i = 0; i < H; i++) {
//            this.Set(hood[i], rn.Double(CarryingCapacity));
//        }
        for (int i = 0; i < H; i++) {
                this.Set(hood[i], CarryingCapacity);
            }
    }
}
