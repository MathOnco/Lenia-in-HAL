package FFT;
import FFT.FFTGrid.*;

@FunctionalInterface
public interface Coords1DDoubleArrToDouble {
    double Eval(int i,double[] vals);
}
