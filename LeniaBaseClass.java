package Lenia;
import HAL.Rand;
import HAL.Util;
import HAL.GridsAndAgents.Grid2Ddouble;
import HAL.Gui.GifMaker;
import HAL.Gui.TickTimer;
import HAL.Gui.UIGrid;
import HAL.Gui.UIWindow;
import HAL.Interfaces.Grid2D;
import HAL.Tools.FileIO;

import java.io.Serializable;
import java.util.ArrayList;

import static HAL.Util.*;

public class LeniaBaseClass implements Grid2D, Serializable {
    public int N; // number of cell types

    // colors:
    final public static int RED = RGB(0.84,0.18,0.14); // average marker color
    final public static int Kmin_RED = RGB(1.0,0.961,0.941), Kmax_RED = RGB(0.404,0.0,0.051);
    final public static int Umin_GREEN = RGB(0.969,0.988,0.961), Umax_GREEN = RGB(0.0,0.267,0.106);
    final public static int Gmin_BLUE = RGB(0.969,0.984,1.0), Gmax_BLUE = RGB(0.031,0.188,0.42);
    public boolean KERNEL_DRAWN = false;

    // enums
    public static final int Ax = 0, Kr = 1, Ux = 2, Gu = 3;

    public final int xDim, yDim, length;
    public int sf = 1; // scale-factor (for drawing)

    // public boolean CONTOURS = false;
    public boolean DISPLAY_TIME_ON_PLOT = true; // print the time on each grid
    public boolean DISPLAY_AVG_ON_PLOT = true; // print the time on each grid
    public boolean DISPLAY_SCALE_BARS_ON_PLOT = true;
    public boolean[] vis_options;
    public int[] colors;
    public int SCALE_BAND_SIZE; // for the colormap scale on the right of each map
    public static final double NUMBER_OF_CONTOURS = 40.0;

    TickTimer tt = new TickTimer();
    public double pause = 0.0; // pause between ticks (seconds)
    String filepath; // used for saving a csv file, and/or gif.
    public final double deltaT;
    public int tick;
    public double t;

    // merging N and 1-player:
    UIGrid[] aVis; // A (state of each Ai world, plus a combined map)
    UIWindow win;
    UIGrid[][] kVis; // K
    UIGrid[][] uVis; // U(x)
    UIGrid[] gVis; // G(U(x)))

    // gif save:
    public boolean toSaveGif = false; // defaults off
    ArrayList<GifMaker> AgifList; // array of gifs for world, Ai
    ArrayList<GifMaker> KgifList; // array of gifs for kernel, Kij
    ArrayList<GifMaker> UgifList; // array of gifs for world, Uij
    ArrayList<GifMaker> GgifList; // array of gifs for world, Gi
    int gif_time_step; // defaults to saving a gif every 1.0 seconds (converted to nearest tick)

    public double C = 1; // max number cells on grid, 1 by default
    double[] CLIP = new double[]{0.0, 1.0}; // lower, upper bound for both U and Aij
//    double[] CLIP = new double[]{0.0, 10.0}; // lower, upper bound for both U and Aij

    // stochastic methods:
    Rand rn;

    public LeniaBaseClass(String filename, int side_length, double deltaT, int nPlayers, int scalefactor) {
        int sideLen = SetSideLength(side_length);
        xDim = sideLen;
        yDim = sideLen;
        length = xDim * yDim;
        this.deltaT = deltaT;
        this.filepath = filename;
        this.sf = scalefactor;
        this.N = nPlayers; // default to 1 player
        rn = new Rand();
    }

    private int SetSideLength(int sideLength) {
        // 2^x = i, solve for x
        double x = (Math.log(sideLength) / Math.log(2));
        int sideLen = (int) Math.pow(2, x);
        if (Math.pow(2,x) > Math.pow(2,Math.floor(x))) {
            sideLen = (int) Math.pow(2, Math.ceil(x));
            System.out.println("For fast computation, must use side length that is a power of 2. Using side length of " + sideLen + " instead of " + sideLength);
        }
        return sideLen;
    }

    // round for output (precision probably isn't very important here)
    public static double Round ( double val){
        return ((double) Math.round(val * 100000d) / 100000d);
    }

    public static double Round3( double val){
        return ((double) Math.round(val * 1000d) / 1000d);
    }

    @Override
    public int Xdim() {
        return xDim;
    }

    @Override
    public int Ydim() {
        return yDim;
    }

    @Override
    public int Length() {
        return length;
    }

    @Override
    public boolean IsWrapX() {
        return true;
    }

    @Override
    public boolean IsWrapY() {
        return true;
    }

    // set the clipped lower, upper bounds
    public void SetClip(double lower, double upper) {
        // lower, upper bound for both U and Aij
        this.CLIP[0] = lower;
        this.CLIP[1] = upper; 
        this.C = (int)CLIP[1];
    }

    public int GetTick() {
        return tick;
    }

    public double GetTime() {
        return t;
    }

    public void ResetTick() {
        this.tick = 0;
        this.t = 0.0;
    }

    public void Pause(double seconds) {
        int milliseconds = (int) (seconds * 1000);
        tt.TickPause(milliseconds);
    }

    // set the default pause length, which pauses the simulation after the Draw() step.
    public void SetPause(double seconds) {
        this.pause = seconds;
    }


    // output fnc for multiple players
    public void Output(Grid2Ddouble[] fields){
        boolean appendBool = (t < deltaT) ? false : true;
        StringBuilder sb = new StringBuilder();
        FileIO fileIO;
        if (appendBool) {
            fileIO = new FileIO((this.filepath+".csv"), "a");
        } else {
            fileIO = new FileIO((this.filepath+".csv"), "w");

            sb.append("Time");
            for (int i = 0; i < fields.length; i++) {
                sb.append(",x" + i + "(t)");
            }
            sb.append("\n");
        }


        sb.append(Round(t));
        for (int i = 0; i < fields.length; i++) {
            sb.append("," + Round(fields[i].GetAvg()));
        }
        sb.append("\n");
        fileIO.Write(sb.toString());
        fileIO.Close();
    }

    // output fnc for single player
    public void Output(Grid2Ddouble field){
        boolean appendBool = (t < deltaT) ? false : true;
        StringBuilder sb = new StringBuilder();
        FileIO fileIO;
        if (appendBool) {
            fileIO = new FileIO((this.filepath+".csv"), "a");
        } else {
            fileIO = new FileIO((this.filepath+".csv"), "w");
            sb.append("Time,x(t)");
            sb.append("\n");
        }

        sb.append(Round(t) + "," + Round(field.GetAvg()));
        sb.append("\n");
        fileIO.Write(sb.toString());
        fileIO.Close();
    }

    

    // val \in [0, 1]:
    public static int ZeroToOneMap(double val) {
        // val = (double)((int)(val*10))/10.0;
        return Util.GreyScale(1-val);
    }

    public static int ZeroToOneColorMap(double val, int color1, int color2){
        return Util.ColorMap(val, color1, color2);
    }

    public static int ZeroToMaxColorMap(double val, double max, int color1, int color2){
        val = Scale0to1(val, 0, max);
        return Util.ColorMap(val, color1, color2);
    }

    // used for G(x) scale bar
    public static int MinToZeroToMaxColorMap(double val, int color1, int color2, int color3){
        if (val >= 0.5) {
            val = Scale0to1(val, 0.5, 1.0);
            return Util.ColorMap(val, color2, color3);
        } else {
            val = Scale0to1(val, 0, 0.5);
            return Util.ColorMap(val, color1, color2);
        }
        
    }

    // A(x) colormap:
    public int AColorMap(double val) {
        val = (double)((int)(val*NUMBER_OF_CONTOURS))/NUMBER_OF_CONTOURS;
        return ZeroToMaxColorMap(val,CLIP[1],WHITE,BLACK);
//        this.colors = new int[] {TumorColor,ImmuneColor};
//        return ZeroToMaxColorMap(val,CLIP[1],WHITE,Util.RGB256(100, 175, 120));
    }

    // A(x) colormap:
    public int AColorMap(double val, int iPlayer) {
        double smallest_interval = (double)CLIP[1] / (double)NUMBER_OF_CONTOURS;

        val = Scale0to1(val, smallest_interval, CLIP[1]); // normalize by max CLIP val.
        val = (double)((int)(val*NUMBER_OF_CONTOURS))/NUMBER_OF_CONTOURS;

        int background = Util.RGB(0.95,0.95,0.95);
        int dark_red = Util.RGB256(210,130,131);
        int light_red = Util.RGB256(210/3,130/3,131/3);

        if (val > smallest_interval) {
            // int color1 = Util.RGB(0.3,0.3,0.3);
            int color1 = Util.RGB(0.8,0.8,0.8);
//            return ColorMap(val, dark_red,light_red); // colors[iPlayer]
//            return ColorMap(val, dark_red,light_red); // colors[iPlayer]

             return RGB(1-(1-GetRed(colors[iPlayer]))*val, 1-(1-GetGreen(colors[iPlayer]))*val, 1-(1-GetBlue(colors[iPlayer]))*val);
        } else {
            return Util.WHITE;
        }

        
    }
    // val = [0, CLIPMAX]
    public int UColorMap(double val, int i) {
        if (N > 1) {
            return AColorMap(val, i);
        } else {
            val = Scale0to1(val, 0, CLIP[1]); // normalize by max CLIP val.
            val = (double)((int)(val*NUMBER_OF_CONTOURS))/NUMBER_OF_CONTOURS;
            return ZeroToMaxColorMap(val,1,Umin_GREEN,Umax_GREEN);
        }
    }

    // K(x) colormap:
    public int KColorMap(double val, double max) {
        val = (double)((int)(val*NUMBER_OF_CONTOURS))/NUMBER_OF_CONTOURS;
        return ZeroToMaxColorMap(val,max,Kmin_RED,Kmax_RED);
    }

    // U(x) colormap:
    public int UColorMap(double val) {
        val = (double)((int)(val*NUMBER_OF_CONTOURS))/NUMBER_OF_CONTOURS;
        return ZeroToMaxColorMap(val,CLIP[1],Umin_GREEN,Umax_GREEN);
    }

    //G(x) colormap:
    // public static int GColorMap(double val, double[] min_max) {
    //     if  (min_max[0] >= -0.000000001) {
    //         // white to blue colormap:
    //         val = Scale0to1(val, min_max[0], min_max[1]);
    //         val = (double)((int)(val*NUMBER_OF_CONTOURS))/NUMBER_OF_CONTOURS;
    //         return ColorMap(val, Gmin_BLUE, Gmax_BLUE);
    //     } else {
    //         // blue to white to red colormap:
    //         double absolute_value_of_max = Math.max(Math.abs(min_max[0]), Math.abs(min_max[1]));
    //         val = Scale0to1(val, -absolute_value_of_max, absolute_value_of_max);
    //         if (val > 0.5) {
    //             val = Scale0to1(val, 0.5, 1);
    //             //val = (double)((int)(val*(NUMBER_OF_CONTOURS/2)))/(NUMBER_OF_CONTOURS/2);
    //             return ColorMap(val, Gmin_BLUE, Gmax_BLUE);
    //         } else {
    //             val = Scale0to1(val, 0, 0.5);
    //             //val = (double)((int)(val*(NUMBER_OF_CONTOURS/2)))/(NUMBER_OF_CONTOURS/2);
    //             return ColorMap(val, BLACK, Gmin_BLUE);
    //         }
    //     }
    // }

    public static int GColorMap(double val, double min, double max) {
        if  (min >= -0.000000001) {
            // white to blue colormap:
            val = Scale0to1(val, min, max);
            val = (double)((int)(val*NUMBER_OF_CONTOURS))/NUMBER_OF_CONTOURS;
            return ColorMap(val, Gmin_BLUE, Gmax_BLUE);
        } else {

            // val is between [min, max]

            // color bar is between [-abs_max, abs_max]

            // blue to white to red colormap:
            double abs_max = Math.max(Math.abs(min), Math.abs(max));
            val = Scale0to1(val, -abs_max, abs_max);

            if (val > 0.5) {
                val = Scale0to1(val, 0.5, 1);
                val = (double)((int)(val*(NUMBER_OF_CONTOURS/2)))/(NUMBER_OF_CONTOURS/2);
                return ColorMap(val, Gmin_BLUE, Gmax_BLUE);
            } else {
                val = Scale0to1(val, 0, 0.5);
                val = (double)((int)(val*(NUMBER_OF_CONTOURS/2)))/(NUMBER_OF_CONTOURS/2);
                return ColorMap(val, BLACK, Gmin_BLUE);
            }
        }
    }


    public void AddPause(double seconds) {
        this.Pause(seconds);
    }






    ////////////////////////////////////////////////////////////////////////////////////////
    /// drawing functions:


    // saves a gifs, every gif_time_step seconds
    public void SetupGifSaving(){
        this.toSaveGif = true;

        this.AgifList = (vis_options[Ax]) ? new ArrayList<>() : null;
        this.KgifList = (vis_options[Kr]) ? new ArrayList<>() : null;
        this.UgifList = (vis_options[Ux]) ? new ArrayList<>() : null;
        this.GgifList = (vis_options[Gu]) ? new ArrayList<>() : null;

        for (int i = 0; i < N; i++) {
            if (vis_options[Ax]) { this.AgifList.add(new GifMaker(this.filepath+"_A"+((N>1)?(i+1):"")+".gif", 150, true)); }
            if (vis_options[Gu]) { this.GgifList.add(new GifMaker(this.filepath+"_G"+((N>1)?(i+1):"")+".gif", 150, true)); }

            for (int j = 0; j < N; j++) {
                if (vis_options[Kr]) { this.KgifList.add(new GifMaker(this.filepath+"_K"+((N>1)?Integer.toString(i+1)+Integer.toString(j+1):"")+".gif", 150, true)); }
                if (vis_options[Ux]) { this.UgifList.add(new GifMaker(this.filepath+"_U"+((N>1)?Integer.toString(i+1)+Integer.toString(j+1):"")+".gif", 150, true)); }
            }
        }
        this.gif_time_step = (int)Math.round(1.0/this.deltaT); // assumes dt < 1
    }

    public void SetupWindows(boolean[] vis_options) {

        // all UIGrids are in the same window:
        boolean window_drawn = false;
        

        // A:
        if (vis_options[Ax]) {
            this.aVis = new UIGrid[N+1]; // include also a combined map
            for (int i = 0; i < N; i++) {
                aVis[i] = new UIGrid(xDim, yDim, sf);
                if (!window_drawn) { 
                    this.win = new UIWindow("Lenia",true);
                    window_drawn = true;
                 }
                win.AddCol(i, "A"+ ((N>1)?(i+1):"") + "(x)");
                win.AddCol(i, aVis[i]);
            }
        }

        // K:
        if (vis_options[Kr]) {
            this.kVis = new UIGrid[N][N];
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    kVis[i][j] = new UIGrid(xDim, yDim, sf);
                    if (!window_drawn) { 
                        this.win = new UIWindow("Lenia",true);
                        window_drawn = true;
                    }
                    win.AddCol(j+N, "K"+ ((N>1)?(i+1):"") + ((N>1)?(j+1):"") + "(x)");
                    win.AddCol(j+N,kVis[i][j]);
                }
            }
        }

        // U:
        if (vis_options[Ux]) {
            this.uVis = new UIGrid[N][N];
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    if (!window_drawn) { 
                        this.win = new UIWindow("Lenia",true);
                        window_drawn = true;
                    }
                    uVis[i][j] = new UIGrid(xDim, yDim, sf);
                    win.AddCol(j+2*N, "U"+ ((N>1)?(i+1):"") + ((N>1)?(j+1):"") + "(x)");
                    win.AddCol(j+2*N,uVis[i][j]);
                }
            }
        }
        

        // G:
        int Gcol = (N>1) ? 0 : 3;
        if (vis_options[Gu]) {
            this.gVis = new UIGrid[N];
            for (int i = 0; i < N; i++) {
                if (!window_drawn) { 
                    this.win = new UIWindow("Lenia",true);
                    window_drawn = true;
                }
                gVis[i] = new UIGrid(xDim,yDim,sf);
                win.AddCol(Gcol+i, "G"+ ((N>1)?(i+1):"") + "(x)");
                win.AddCol(Gcol+i,gVis[i]);
            }
        }

        if (window_drawn) { win.RunGui(); }

        if (N==1) {
            this.colors = new int[] {BLACK};
        } else if (N<=3) {
            this.colors = new int[] {RED,BLUE,GREEN};
        } else {
            this.colors = Util.CategoricalColors(0, N);
        }
    }

    // L is lambda (rate of events) and deltaT is the time interval
    // public int NPoissonEvents(double L, double deltaT) {
    //     // if rate is negative, then returns negative # events (integer values only)
    //     // if rate is positive, then returns positive # events (integer values only)
    //     int sign = (L < 0) ? -1 : 1;
    //     double draw = rn.Double();


    //     int n = 0;
    //     int nFactorial = 1;
    //     double poisson = PoissonRate(L, deltaT, n, nFactorial);

    //     if (rn.Double() < poisson) {
    //         return n; // 0 events
    //     } else {
    //         while (draw > poisson) {
    //             // try n events:
    //             n = n + 1;
    //             nFactorial = nFactorial*n;
    //             poisson += PoissonRate(L, deltaT, n, nFactorial);
    //         }
    //         return n*sign;
    //     }
    // }

    public static double PoissonRate(double L, double deltaT, int n, int nFactorial) {
        return Math.pow(Math.abs(L),n) / nFactorial * Math.exp(-Math.abs(L)*deltaT);
    }

}

