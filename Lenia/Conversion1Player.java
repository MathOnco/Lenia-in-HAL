package Lenia;
import FFT.FFTGrid.*;

import HAL.Rand;
import HAL.Util;
import HAL.GridsAndAgents.Grid2Ddouble;
import HAL.Tools.FileIO;

public class Conversion1Player extends Lenia1Player {

    //Lenia1Player model;
    public static double Rstar = 3;
    public static double gamma = 10.0;
    public static double u0 = 0.1;
    public int[] full_domain_hood;
    public Rand rn;
    
    
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
        return gamma*u*(1.0-u);
    }

    // create a constructor for the Gompertz class, but just inherit the constructor from Lenia1Player
    public Conversion1Player(String filename, int side_length, double dt, int scalefactor, double ClipMax) {
        super(filename,side_length, dt, scalefactor, new boolean[]{false,false,false,false}, ClipMax);
        this.rn = new Rand();
        this.full_domain_hood = new int[length];

        for (int i = 0; i < length; i++) {
            full_domain_hood[i] = i;
        }

        // full_domain_hood = Util.RectangleHood(true,this.xDim,this.yDim);
    }

    public void SetNfCells(int nf) {
        // shuffle my array
        rn.Shuffle(full_domain_hood);

        for (int i = 0; i < nf; i++) {
            this.Set(full_domain_hood[i], 1);
        }
    }

    public void SetZero() {
        for (int i = 0; i < length; i++) {
            this.Set(i, 0);
        }
    }

    public double Avg() {
        return A[0].GetAvg();
    }

    public static double Round ( double val){
        return ((double) Math.round(val * 100000d) / 100000d);
    }

    public void RunAndWrite() {

        this.RecalcKernels();

        String filepath = "Conversion/L="+Integer.toString(this.xDim) + "/R-" + Double.toString(Rstar);
        StringBuilder sb = new StringBuilder();
        FileIO fileIO = new FileIO((filepath+".csv"), "w");
        double avg0,avgF,Gestimate;


        int TOTAL_SIMS = 10;

        sb.append(Round(0));
        for (int nf = 1; nf <= (this.length); nf++) {
            Gestimate = 0.0;

            for (int nSims = 0; nSims < TOTAL_SIMS; nSims++) {
                this.SetNfCells(nf);
                avg0=this.Avg();
                this.Update();
                avgF=this.Avg();
                Gestimate+=(avgF-avg0)/this.deltaT;
                this.SetZero();
            }

            Gestimate /= (double)TOTAL_SIMS;

            sb.append("," + Round(Gestimate));

            
        }

        fileIO.Write(sb.toString());
        fileIO.Close();

    }

    public static void main(String[] args) {

        int side_length = 64*2; // length of the side of simulation domain
        int scalefactor = 5; // to scale the drawing size

        double deltaT = 0.01; // timestep
        double ClipMax = 1;

        Conversion1Player model=new Conversion1Player("data/test",side_length,deltaT, scalefactor, ClipMax);
        

        double[] R = new double[]{1,2,4,6,10,50};
        // K = [5,13,49,113,317,L2];

        

        for (int i = 0; i < R.length; i++) {
            Rstar = R[i];

            System.out.println(model.Rstar);
            model.RunAndWrite();
        }

        model.Close();
        


    }
}
