package FFT;

import java.text.DecimalFormat;

public class Main {
    public static final double  pi = Math.atan2(1, 1) * 4;
    public static final int samples = 64;//Number of samples in the buffer
    public static final double time = 0.5;//total time recorded in the buffer
    public static final double sine_freq = 4;// just a generated sine frequency as an example
    public static final double sine_freq_period = pi / (((double)samples / time) / sine_freq);//calculating sine wave sample period

    private static double[]  x1, x2, y2, res;

    public static void main(String[] args)
    {
        x1 = new double[samples];
        x2 = new double[samples];
        y2 = new double[samples];
        res = new double[samples];

        for(int i=0;i<samples;i++)
        {
            x1[i]= Math.sin((double)(2)*(double)i*sine_freq_period) + 2 * Math.sin((double)(2)*(double)i*sine_freq_period * 2) + 0.5* Math.sin((double)(2)*(double)i*sine_freq_period * 4);

            x2[i]=x1[i];
            y2[i]=0;
        }

        //FFT
        if (fft(samples, x2, y2)==1) System.exit(-9);


        System.out.println("\n Frequency\t\t\tFFT\n");
        for (int i = 0; i < samples/2; i++)
        {
            res[i]=Math.sqrt(x2[i]*x2[i]+y2[i]*y2[i]);
            System.out.println(" " + formatDoubleToString(((samples / time) / samples)*i) + " Hz:\t\t\t\t" + formatDoubleToString(res[i]));
        }
    }

    static void makeSineTable(int n, double[] sintbl)
    {
        int i, n2, n4, n8;
        double c, s, dc, ds, t;
        n2 = n / 2;  n4 = n / 4;  n8 = n / 8;
        t = Math.sin(pi / n);
        dc = 2 * t * t;  ds = Math.sqrt(dc * (2 - dc));
        t = 2 * dc;  c = sintbl[n4] = 1;  s = sintbl[0] = 0;
        for (i = 1; i < n8; i++) {
            c -= dc;  dc += t * c;
            s += ds;  ds -= t * s;
            sintbl[i] = s;
            sintbl[n4 - i] = c;
        }
        if (n8 != 0) sintbl[n8] = Math.sqrt(0.5);
        for (i = 0; i < n4; i++)
            sintbl[n2 - i] = sintbl[i];
        for (i = 0; i < n2 + n4; i++)
            sintbl[i + n2] = - sintbl[i];
    }

    static void makeBitReverse(int n, int[] bitrev)
    {
        int i, j, k, n2;
        n2 = n / 2;  i = j = 0;
        for ( ; ; ) {
            bitrev[i] = j;
            if (++i >= n) break;
            k = n2;
            while (k <= j) {  j -= k;  k /= 2;  }
            j += k;
        }
    }

    static int fft(int n, double[] x, double[] y)
    {
        int    last_n = 0;
        int[]   bitrev = null;
        double[] sintbl = null;
        int i, j, k, ik, h, d, k2, n4, inverse;
        double t, s, c, dx, dy;

        if (n < 0) {
            n = -n;  inverse = 1;
        } else inverse = 0;
        n4 = n / 4;
        if (n != last_n || n == 0) {
            last_n = n;
            if (sintbl != null) sintbl = null;
            if (bitrev != null) bitrev = null;
            if (n == 0) return 0;
            sintbl = new double[n + n4];
            bitrev = new int[n];
            if (sintbl == null || bitrev == null) {
                return 1;
            }
            makeSineTable(n, sintbl);
            makeBitReverse(n, bitrev);
        }
        for (i = 0; i < n; i++) {
            j = bitrev[i];
            if (i < j) {
                t = x[i];  x[i] = x[j];  x[j] = t;
                t = y[i];  y[i] = y[j];  y[j] = t;
            }
        }
        for (k = 1; k < n; k = k2) {
            h = 0;  k2 = k + k;  d = n / k2;
            for (j = 0; j < k; j++) {
                c = sintbl[h + n4];
                if (inverse!=0) s = - sintbl[h];
                else         s =   sintbl[h];
                for (i = j; i < n; i += k2) {
                    ik = i + k;
                    dx = s * y[ik] + c * x[ik];
                    dy = c * y[ik] - s * x[ik];
                    x[ik] = x[i] - dx;  x[i] += dx;
                    y[ik] = y[i] - dy;  y[i] += dy;
                }
                h += d;
            }
        }

        if (inverse==0)
            for (i = 0; i < n; i++)
            {
                x[i] /= n;
                y[i] /= n;
            }

        return 0;
    }

    private static String formatDoubleToString(double number)
    {
        return new DecimalFormat("##########.########").format(number).toString();
    }
}
