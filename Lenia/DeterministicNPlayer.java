package Lenia;

import HAL.Util;
import HAL.GridsAndAgents.Grid2Ddouble;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import static HAL.Util.*;
import FFT.FFTGrid;
import FFT.Coords1DDoubleArrToDouble;
import FFT.Coords2DDoubleToDouble;


public class DeterministicNPlayer extends LeniaBaseClass implements Iterable<Grid2Ddouble> {
    
    public Coords1DDoubleArrToDouble GrowthFunction; // G(u)
    public final Grid2Ddouble[] A; // A(x)
    public final FFTGrid[][]fftFields; // FFT of Ai (changes every time step)
    public FFTGrid[][]fftKernels; // FFT of Kij (never changes)

    // private HashMap<String, KernelInfo>kernelStorage=new HashMap<>();
    public Coords2DDoubleToDouble KernelFunction; // K(r)
    public Grid2Ddouble[][]kernelFields;
    public double[][]kSums;
    // private int[][]fftFieldRefs;
    

    public boolean KERNEL_DRAWN = false;
    public boolean GROWTH_MAX_CALCULATED = false;
    public double[][] Gmin_max;
    public double[] Uscratch,UNormscratch;
    public double USumScratch;
    
    public ArrayList<Grid2Ddouble> GrowthScratch; // only used when drawing

    public double[][] Rstar = new double[][]{{50.0,50.0},{50.0,50.0}}; // payoff matrix


    // use the default kernal, growth function.
    DeterministicNPlayer(String filename, int sideLength, double deltaT, int nPlayers, int scalefactor, int ClipMax) {
        this(filename,sideLength,deltaT,nPlayers,scalefactor,new boolean[]{true,true,true,true},ClipMax);
    }

    // String filename, int sideLength, double deltaT, int scalefactor, boolean[] vis_options
    DeterministicNPlayer(String filename, int sideLength, double deltaT, int nPlayers, int scalefactor, boolean[] vis_options, double ClipMax) {
        super(filename,sideLength,deltaT,nPlayers,scalefactor);
        this.KernelFunction = this::K;
        this.GrowthFunction = this::G;
        this.A = new Grid2Ddouble[N];
        this.fftFields = new FFTGrid[N][N];
        // this.fftFieldRefs=new int[N][N];
        this.fftKernels = new FFTGrid[N][N];
        this.kernelFields = new Grid2Ddouble[N][N];
        this.kSums = new double[N][N];
        this.Uscratch = new double[N];
        this.vis_options = vis_options;
        this.GrowthScratch = (vis_options[Gu]) ? new ArrayList<>() : null;
        this.SCALE_BAND_SIZE = Math.max(2, (int)Math.round(0.03 * sideLength));

        for (int i = 0; i < N; i++) {
            this.A[i] = new Grid2Ddouble(xDim, yDim,true,true);
            if (vis_options[Gu]) { this.GrowthScratch.add(new Grid2Ddouble(sideLength, sideLength,true,true) ); }
        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                fftFields[i][j] = new FFTGrid(A[i]);
                fftKernels[i][j] = new FFTGrid(xDim, yDim);
                kernelFields[i][j] = new Grid2Ddouble(xDim, yDim);
            }
        }
        RecalcKernels();
        SetupWindows(vis_options); // set up only the desired windows
        CalculateGrowthMax();

        this.C = ClipMax;
        SetClip(0, C);
    }
    
    // default kernal
    // this is for the (i,j)th kernal, at r distance from the center of the kernal
    public double K(int i,int j,double r){
        return (r<1.8) ? 1.0 : 0.0;
    }

    // default growth function
    // this is for the (i)th world, where i \in [0,nPlayers-1]
    public double G(int i,double[] U){
        return 0.0;
    }
    
    public void RecalcKernels(){
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                Grid2Ddouble newKernel = kernelFields[i][j];
                double kSum=0;
                for (int x = 0; x < xDim; x++) {
                    for (int y = 0; y < yDim; y++) {
                        double val=KernelFunction.Eval(i,j,Dist(x,y,xDim/2.0,yDim/2.0));
                        kSum+=val;
                        newKernel.Set(x,y,val);
                    }
                }
                kSums[i][j]=kSum;
                for (int k = 0; k < newKernel.length; k++) {
                    newKernel.Set(k,newKernel.Get(k)/kSum);
                }
                fftKernels[i][j].SetGrid(newKernel);
                fftKernels[i][j].fftshift();
                fftKernels[i][j].fft2();
            }
        }
        // FindKernelDuplicates();
    }
    public void WriteGrowthField(int i){
        if (vis_options[Gu]) {
            Grid2Ddouble currField = this.A[i];
            for (int k = 0; k < currField.length; k++) {
                SetUiScratch(i,k);
                GrowthScratch.get(i).Set(k,this.GrowthFunction.Eval(i, Uscratch));
            }
        }
        
    }

    // naive way to calculate growth max:
    public void CalculateGrowthMax() {
        if (!GROWTH_MAX_CALCULATED) {
            this.Gmin_max = GetGrowthFunctionMinMax();
            for (int i = 0; i < this.Gmin_max.length; i++) {
//                System.out.println("Growth function " + i + " min: " + this.Gmin_max[i][0] + " max: " + this.Gmin_max[i][1]);
            }
            GROWTH_MAX_CALCULATED=true;
        }
    }

    public void Update() {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                //if (fftFieldRefs[i][j] == i) {
                Grid2Ddouble convField = this.A[j];
                FFTGrid currFftField = this.fftFields[i][j];
                currFftField.SetGrid(convField);
                currFftField.fft2();
                currFftField.ElementwiseMultiplication(fftKernels[i][j]);
                currFftField.ifft2();
                //}
            }
        }
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < A[i].length; k++) {
                SetUiScratch(i,k);
                A[i].Set(k, Bound(A[i].Get(k) + this.GrowthFunction.Eval(i, Uscratch) * deltaT, CLIP[0],CLIP[1]));
            }
        }
        //update Ais with Gis
        this.tick++;
        this.t += deltaT;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // this is U = K * A
    public double Uij(int i, int j, int index) {
        return Bound(fftFields[i][j].REAL.Get(index),CLIP[0],CLIP[1]);
    }
    public double Uij(int i, int j, int x, int y) {
        return Bound(fftFields[i][j].REAL.Get(x,y),CLIP[0],CLIP[1]);
    }
    // sets Uscrath vector of Ui
    public void SetUiScratch(int iPlayer, int k) {
        for (int j = 0; j < N; j++) {
            Uscratch[j]=Uij(iPlayer,j, k);
        }
    }
    public void SetUiScratch(int iPlayer, int x, int y) {
        for (int j = 0; j < N; j++) {
            Uscratch[j]=Uij(iPlayer,j, x, y);
        }
    }

    // // returns normalized vector of Unorm = U / sum(U)
    // public void SetUiNormScratch(int iPlayer, int k) {
    //     USumScratch=0;
    //     for (int j = 0; j < N; j++) { USumScratch+=Uij(iPlayer,j, k); }
    //     for (int j = 0; j < N; j++) { UNormscratch[j]=Uij(iPlayer,j, k)/USumScratch; }
    // }

    // // returns normalized vector of Unorm = U / sum(U)
    // public void SetUiNormScratch(int iPlayer, int x, int y) {
    //     USumScratch=0;
    //     for (int j = 0; j < N; j++) { USumScratch+=Uij(iPlayer,j, x, y); }
    //     for (int j = 0; j < N; j++) { UNormscratch[j]=Uij(iPlayer,j, x, y)/USumScratch; }
    // }

    // get average of potential field:
    public double UijGetAvg(int i, int j){
        return Bound(fftFields[i][j].REAL.GetAvg(),CLIP[0],CLIP[1]);
    }
    public double GetKernelVal(int iPlayer, int jPlayer, int i){
        return kernelFields[iPlayer][jPlayer].Get(i)*kSums[iPlayer][jPlayer];
    }
    public double GetKernelVal(int iPlayer, int jPlayer, int x, int y){
        return kernelFields[iPlayer][jPlayer].Get(x,y)*kSums[iPlayer][jPlayer];
    }
    public double GetKernelMax(int i, int j) {
        return kernelFields[i][j].GetMax()*kSums[i][j];
    }
    public double GetKernelAvg(int i, int j){
        return kernelFields[i][i].GetAvg()*kSums[i][i];
    }
    public void ResetTick(){
        this.tick=0;
    }
    public int GetTick() {
        return tick;
    }
    public double GetTime() {
        return t;
    }
    public double Get(int iPlayer,int i){
        return this.A[iPlayer].Get(i);
    }
    public double Get(int iPlayer,int x,int y){
        return this.A[iPlayer].Get(x,y);
    }
    public void Set(int iPlayer,int i,double v){
        this.A[iPlayer].Set(i,Bound(v,CLIP[0],CLIP[1]));
    }
    public void Set(int iPlayer,int x,int y,double v){
        this.A[iPlayer].Set(x,y,Bound(v,CLIP[0],CLIP[1]));
    }
    // // add N cells to location:
    // public void AddNCells(int iPlayer, int x, int y, int events) {
    //     this.A[iPlayer].Add(x,y,events);
    // }
    

    @Override
    public Iterator<Grid2Ddouble> iterator() {
        return new FieldIterator();
    }
    private class FieldIterator implements Iterator<Grid2Ddouble>,Iterable<Grid2Ddouble>,Serializable{
        int i;
        @Override
        public Iterator<Grid2Ddouble> iterator() {
            return this;
        }

        @Override
        public boolean hasNext() {
            if(i<N){
                return true;
            }
            return false;
        }

        @Override
        public Grid2Ddouble next() {
            return A[i++];
        }
    }

    

    // output average domain density to csv:
    public void Output() {
        this.Output(this.A);
    }
    

    // public void SetGrowthScaleBars(double[][] min_max) {
    //     this.Gmin_max = min_max;
    // }

    public void SetGrowthScale(int iPlayer, double min, double max) {
        this.Gmin_max[iPlayer][0] = min;
        this.Gmin_max[iPlayer][1] = max;
    }

    public double GetGrowthScaleMin(int iPlayer) {
        return this.Gmin_max[iPlayer][0];
    }
    public double GetGrowthScaleMax(int iPlayer) {
        return this.Gmin_max[iPlayer][1];
    }

    // determine max & min of the GrowthFunction
    // should change this to nested for loop (N-loops)
    public double[][] GetGrowthFunctionMinMax() {
        double[][] all_min_max = new double[N][2];
        if (N == 1) {
            int i = 0;
            for (double u1 = 0; u1 <= CLIP[1]; u1+=(CLIP[1]/1000)) {
                double[]u = {u1};
                double g = GrowthFunction.Eval(i,u);

                all_min_max[i][0] = Math.min(all_min_max[i][0], g);
                all_min_max[i][1] = Math.max(all_min_max[i][1], g);
            }
        } else if (N == 2) {
            for (int i = 0; i < N; i++) {
                for (double u1 = 0; u1 <= CLIP[1]; u1+=(CLIP[1]/1000)) {
                    for (double u2 = 0; u2 <= CLIP[1]; u2+=(CLIP[1]/1000)) {
                        double[]u = {u1,u2};
                        double g = GrowthFunction.Eval(i,u);

                        all_min_max[i][0] = Math.min(all_min_max[i][0], g);
                        all_min_max[i][1] = Math.max(all_min_max[i][1], g);
                    }
                }
            }
        } else if (N == 3) {
            for (int i = 0; i < N; i++) {
                for (double u1 = 0; u1 <= CLIP[1]; u1+=(CLIP[1]/1000)) {
                    for (double u2 = 0; u2 <= CLIP[1]; u2+=(CLIP[1]/1000)) {
                        for (double u3 = 0; u3 <= CLIP[1]; u3+=0.001) {
                            double[]u = {u1,u2,u3};
                            double g = GrowthFunction.Eval(i,u);

                            all_min_max[i][0] = Math.min(all_min_max[i][0], g);
                            all_min_max[i][1] = Math.max(all_min_max[i][1], g);
                        }
                    }
                }
            }
        }
        return all_min_max;
    }

    // these things could be moved into a base class

    public void DrawKernels() {
        if (!KERNEL_DRAWN) {
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    for (int k = 0; k < this.length; k++) {
                        if (vis_options[Kr]) { kVis[i][j].SetPix(k, KColorMap(GetKernelVal(i,j,k),GetKernelMax(i,j))); }
                    }
                }
            }
            if (DISPLAY_AVG_ON_PLOT){
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        if (vis_options[Kr]) { kVis[i][j].SetString("avg="+Round3(GetKernelAvg(i, j)),1,8,Kmax_RED,Kmin_RED); }
                    }
                }
            }
            // these gif's are only one frame:
            if (toSaveGif){
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        if (vis_options[Kr]) { KgifList.get(i*N+j).AddFrame(kVis[i][j]); }
                    }
                }
            }
        }
        KERNEL_DRAWN=true;
    }
    public void Draw(){
        DrawKernels();
        CalculateGrowthMax();
    
        for (int i = 0; i < N; i++) {
            if (vis_options[Gu]) {
                WriteGrowthField(i);
            }
            

            for (int index = 0; index < this.length; index++) {
                if (vis_options[Ax]) { 
                    aVis[i].SetPix(index,AColorMap(A[i].Get(index),i)); 
                }

                for (int j = 0; j < N; j++) {
                    if (vis_options[Ux]) { uVis[i][j].SetPix(index,UColorMap(Uij(i, j, index),j)); }
                }

                if (vis_options[Gu]) { 
                    SetUiScratch(i, index);
                    gVis[i].SetPix(index,GColorMap(this.GrowthFunction.Eval(i, Uscratch), GetGrowthScaleMin(i), GetGrowthScaleMax(i)));
                }
            }

            if (this.DISPLAY_TIME_ON_PLOT){
                if (vis_options[Ax]) { aVis[i].SetString("t="+Round(t),1,yDim-1,Util.BLACK,Util.WHITE); }
                if (vis_options[Gu]) { gVis[i].SetString("t="+Round(t),1,yDim-1,Gmax_BLUE,Gmin_BLUE); }
                for (int j = 0; j < N; j++) {
                    if (vis_options[Ux]) { uVis[i][j].SetString("t="+Round(t),1,yDim-1,Util.BLACK,Util.WHITE); }
                }
            }
        }        

        if (this.DISPLAY_SCALE_BARS_ON_PLOT) {
            AddScaleBars();
            AddAverageMarkers();
        }
        
        SaveAllGifs();

        System.out.println("Time: " + Round(t));
        Pause(pause); // pause for "pause" # of seconds
    }
    public void SaveAllGifs() {
        for (int i = 0; i < N; i++) {
            // add gifs:
            if (this.toSaveGif) {
                if (GetTick() % gif_time_step == 0 ) {
                    if (vis_options[Ax]) { AgifList.get(i).AddFrame(aVis[i]); }
                    if (vis_options[Gu]) { GgifList.get(i).AddFrame(gVis[i]); }
                    for (int j = 0; j < N; j++) {
                        if (vis_options[Ux]) { UgifList.get(i*N+j).AddFrame(uVis[i][j]); }
                    }
                }
            }
        }
    }
    public void AddScaleBars() {
        double p;
        for (int i = 0; i < N; i++) { // type i = [0,1, ..., N]
            for (int y = 0; y < this.yDim; y++) {
                for (int x = this.xDim-SCALE_BAND_SIZE; x < this.xDim; x++) {
                    // A:
                    p = (double)(y)/(double)(this.yDim);
                    if (vis_options[Ax]) { aVis[i].SetPix(x,y,AColorMap(p*this.CLIP[1],i)); }
                    if (vis_options[Gu]) { gVis[i].SetPix(x,y,GColorMap(p,GetGrowthScaleMin(i), GetGrowthScaleMax(i))); }

                    for (int j = 0; j < N; j++) { // type j = [0,1, ..., N]
                        if (vis_options[Ux]) { uVis[i][j].SetPix(x,y,UColorMap(p*this.CLIP[1],j)); }
                        if (vis_options[Kr]) { kVis[i][j].SetPix(x,y,KColorMap(p,GetKernelMax(i, j))); }
                    }

                    // scale p from gmin to gmax
                    if (vis_options[Gu]) { 
                        // double abs_max = Math.max(Math.abs(GetGrowthScaleMin(i)), Math.abs(GetGrowthScaleMax(i)));
                        // // p = p*(GetGrowthScaleMax(i) - GetGrowthScaleMin(i)) + GetGrowthScaleMin(i);
                        // p = Scale0to1(p, -abs_max/2, abs_max/2);

                        // if (GetGrowthScaleMin(i) > 0) {
                        //     gVis[i].SetPix(x,y,GColorMap(p, -1, 1));
                        // }

                        // gVis[i].SetPix(x,y,GColorMap(p, -1, 1));
                        
                        gVis[i].SetPix(x,y,MinToZeroToMaxColorMap(p, BLACK, Gmin_BLUE, Gmax_BLUE));
                    }
                }
            }
        }
    }
    public void AddAverageMarkers() {
        int y;
        if (this.DISPLAY_AVG_ON_PLOT){
            for (int x = this.xDim-SCALE_BAND_SIZE; x < this.xDim; x++) {
                for (int i = 0; i < N; i++) {
                    // add average marker bar to A(x) field:
                    if (vis_options[Ax]) { 
                        aVis[i].SetPix(x,(int)((this.A[i].GetAvg()/CLIP[1])*(yDim-1)),BLACK); 
                    }
                    for (int j = 0; j < N; j++) {
                        // add average marker bar to U(x) field:
                        if (vis_options[Ux]) { 
                            uVis[i][j].SetPix(x,(int)((UijGetAvg(i,j)/CLIP[1])*(yDim-1)),BLACK); 
                        }
                        // add average marker bar to K(x) field:
                        if (vis_options[Kr]) { 
                            kVis[i][j].SetPix(x,(int)(GetKernelAvg(i,j)*(yDim-1)),BLACK); 
                        }
                    }

                    // scale p from gmin to gmax
                    if (vis_options[Gu]) { 
                        // double val = 


                        y = (int)(Bound((yDim-1)*(this.GrowthScratch.get(i).GetAvg() - this.Gmin_max[i][0])/(this.Gmin_max[i][1] - this.Gmin_max[i][0]),0,yDim-1));
                        gVis[i].SetPix(x,y,RED); 
                    }
                }
            }
        
            for (int i = 0; i < N; i++) {
                if (vis_options[Ax]) { aVis[i].SetString("avg="+Round3(A[i].GetAvg()),1,8,Util.BLACK,Util.WHITE); }
                if (vis_options[Gu]) { 
                    gVis[i].SetString("avg="+Round3(GrowthScratch.get(i).GetAvg()),1,8,Gmax_BLUE,Gmin_BLUE); 
                }
                for (int j = 0; j < N; j++) {
                    if (vis_options[Ux]) { uVis[i][j].SetString("avg="+Round3(UijGetAvg(i, j)),1,8,Umax_GREEN,Umin_GREEN); }
                }
            }
        }

        
    }
    public void Close() {
        CloseGifs();
        win.Close();
    }

    public void CloseGifs() {
        for (int i = 0; i < N; i++) {
            if (vis_options[Ax]) { AgifList.get(i).Close(); }
            if (vis_options[Gu]) { GgifList.get(i).Close(); }
                        
            for (int j = 0; j < N; j++) {
                if (vis_options[Ux]) { UgifList.get(i*N+j).Close(); }
                if (vis_options[Kr]) { KgifList.get(i*N+j).Close(); }
            }
        }
    }
    


    /////////////////////////////////////////////////////
    // may add this functionality back in later:


    // public void StoreKernel(String label){
    //     this.kernelStorage.put(label, SaveKernelInfo());
    // }
    // public void LoadKernel(String label){
    //     LoadKernelInfo(this.kernelStorage.get(label));
    // }
    
    

    // private void FindKernelDuplicates(){
    //     for (int j = 0; j < N; j++) {
    //         for (int i = 0; i < N; i++) {
    //             Grid2Ddouble g1Real=fftKernels[i][j].REAL;
    //             Grid2Ddouble g1Imag=fftKernels[i][j].IMAG;
    //             for (int k = 0; k < i; k++) {
    //                 Grid2Ddouble g2Real=fftKernels[k][j].REAL;
    //                 Grid2Ddouble g2Imag=fftKernels[k][j].IMAG;
    //                 boolean identical=true;
    //                 for (int l = 0; l < g2Real.length; l++) {
    //                     if(g1Real.Get(l)!=g2Real.Get(l)){
    //                         identical=false;
    //                         break;
    //                     }
    //                     if(g1Imag.Get(l)!=g2Imag.Get(l)){
    //                         identical=false;
    //                         break;
    //                     }
    //                 }
    //                 if(identical){
    //                     fftFieldRefs[i][j]=k;
    //                 }
    //                 else{
    //                     fftFieldRefs[i][j]=i;
    //                 }
    //             }
    //         }
    //     }
    // }

    // private KernelInfo SaveKernelInfo(){
    //     return new KernelInfo(KernelFunction,kernelFields,kSums,fftFieldRefs,fftKernels);
    // }
    // private void LoadKernelInfo(KernelInfo from){
    //     this.KernelFunction=from.Kernel;
    //     this.kernelFields=from.kernelFields;
    //     this.kSums=from.kSums;
    //     this.fftFieldRefs=from.fftFieldRefs;
    //     this.fftKernels=from.fftKernels;
    // }

    // private class KernelInfo {
    //     public final Coords2DDoubleToDouble Kernel;
    //     public final Grid2Ddouble[][]kernelFields;
    //     public final double[][]kSums;
    //     public final int[][]fftFieldRefs;
    //     public final FFTGrid[][]fftKernels;
    //     KernelInfo(Coords2DDoubleToDouble Kernel,Grid2Ddouble[][] kernelFields, double[][] kSums, int[][] fftFieldRefs, FFTGrid[][] fftKernels){
    //         this.Kernel=Kernel;
    //         this.kernelFields = kernelFields;
    //         this.kSums = kSums;
    //         this.fftFieldRefs = fftFieldRefs;
    //         this.fftKernels = fftKernels;
    //     }
    // }
}
