package Lenia;
import HAL.Interfaces.DoubleToDouble;
import static HAL.Util.*;

public class Lenia1Player extends LeniaNPlayer {

    // also used in stochastic lenia:
    public DoubleToDouble GrowthFunction;
    public DoubleToDouble KernelFunction; // K(r)
    
    
    private double kSum;
    public static double Rstar = 2;

    // use the default kernal, growth function.
    Lenia1Player(String filename, int sideLength, double deltaT, int scalefactor, double ClipMax) {
        this(filename,sideLength,deltaT,scalefactor,new boolean[]{true,true,true,true},ClipMax);
        this.C = ClipMax;
        SetClip(0, C);
    }

    // use the default Kc and G function:
    public Lenia1Player(String filename, int sideLength, double deltaT, int scalefactor, boolean[] vis_options, double ClipMax) {
        super(filename,sideLength,deltaT,1,scalefactor,vis_options,ClipMax);
        this.KernelFunction = this::K;
        this.GrowthFunction = this::G;
        RecalcKernels();
        if (toSaveGif) { SetupGifSaving(); }
        SetupWindows(vis_options); // set up only the desired windows
        this.SCALE_BAND_SIZE = Math.max(2, (int)Math.round(0.03 * sideLength));
        this.C = ClipMax;
        SetClip(0, C);
    }

    // single-player G (over ride this one in Lenia1Player extensions)
    public double G(double u) {
        return u*(1.0 - Math.pow(u,1.0));
    }

    // single-player G (over ride this one in Lenia1Player extensions)
    public double K(double r) {
        return (r<Rstar) ? 1.0 : 0.0;
    }

    // do not change:
    @Override 
    public double G(int i, double[] U) {
        return G(U[0]);
    }

    // do not change:
    @Override
    public double K(int i, int j, double r) {
        return K(r);
    }
    public void Update() {
        // FFT2d, elementwise multiply, then IFFT2d
        fftFields[0][0].SetGrid(A[0]);
        this.fftFields[0][0].fft2();
        this.fftFields[0][0].ElementwiseMultiplication(this.fftKernels[0][0]);
        this.fftFields[0][0].ifft2();

        // update every pixel in the field:
        for (int i = 0; i < A[0].length; i++) {
            if (vis_options[Gu]) {
                GrowthScratch.get(0).Set(i, this.GrowthFunction.Eval(U(i)));
            }
            A[0].Set(i, Bound(A[0].Get(i) + this.GrowthFunction.Eval(U(i)) * deltaT, CLIP[0],CLIP[1]));
        }
        this.tick++;
        this.t += deltaT;
    }
    // U(x): ith location in the field
    public double U(int i){
        return Bound(fftFields[0][0].REAL.Get(i),CLIP[0],CLIP[1]);
    }
    // U(x): (x,y)th location in the field
    public double U(int x, int y){
        return Bound(fftFields[0][0].REAL.Get(x,y),CLIP[0],CLIP[1]);
    }

    // get average of potential field:
    public double GetUAvg(){
        return fftFields[0][0].REAL.GetAvg();
    }
    // get average of potential field:
    public double GetUMax(){
        return fftFields[0][0].REAL.GetMax();
    }
    //used for visualization
    public double GetKernelVal(int i) {
        return kernelFields[0][0].Get(i)*kSum;
    }
    public double GetKernelVal(int x, int y) {
        return kernelFields[0][0].Get(x,y)*kSum;
    }
    public double GetKernelMax() {
        return kernelFields[0][0].GetMax()*kSum;
    }
    public void Set(int i, double val) {
        this.A[0].Set(i, val);
    }

    public void Set(int x, int y, double val) {
        this.A[0].Set(x, y, val);
    }

    public double Get(int i) {
        return this.A[0].Get(i);
    }

    public double Get(int x, int y) {
        return this.A[0].Get(x, y);
    }
    public void Output() {
        this.Output(this.A[0]);
    }
    public void DrawKernel() {
        for (int i = 0; i < this.length; i++) {
            kVis[0][0].SetPix(i,KColorMap(this.GetKernelVal(i), GetKernelMax()));
        }
        KERNEL_DRAWN=true;
    }
    public void SetGrowthScale(double min, double max) {
        this.SetGrowthScale(0, min, max);
    }
}

