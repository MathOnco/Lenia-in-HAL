package FFT;
import FFT.FFTGrid.*;

@FunctionalInterface
public interface Coords2DDoubleToDouble {
    double Eval(int i,int j,double v);
}
