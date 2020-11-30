package FFT;

import java.text.DecimalFormat;

public class Main {
    public static final double  pi = Math.atan2(1, 1) * 4.0;//calculated pi
    public static final int samples = 64;//Number of samples in the buffer
    public static final double time = 1.0;//total time recorded in the buffer
    public static final double sine_freq_period = pi / ((double)samples / time);//calculating sine wave sample period

    private static double[] in, real, img, res;

    public static void main(String[] args)
    {
        in = new double[samples];
        real = new double[samples];
        img = new double[samples];
        res = new double[samples];

        for(int i=0;i<samples;i++)
        {
            in[i]= Math.sin(2.0*(double)i*sine_freq_period)
                    + 3 * Math.sin(2.0*(double)i*sine_freq_period * 2)
                    + 0.5 * Math.sin(2.0*(double)i*sine_freq_period * 4);

            real[i]= in[i];
            img[i]=0;
            System.out.println(real[i]);
        }

        //FFT
        fft(samples, real, img);


        System.out.println("\n Frequency\t\t\tFFT magnitude");
        for (int i = 0; i < samples/2; i++)
        {
            res[i]=Math.sqrt(real[i]* real[i]+ img[i]* img[i])*2.0;
            System.out.println(" " + formatDoubleToString(((samples / time) / samples)*i) + " Hz:\t\t\t\t" + formatDoubleToString(res[i]));
        }

        System.out.println("\n Frequency\t\t\tFFT real+img");

        for (int i = 0; i < samples/2; i++)
        {
            System.out.println(" " + formatDoubleToString(((samples / time) / samples)*i) + " Hz:\t\t\t\tr: " + formatDoubleToString(real[i]) + "\t\ti: " + formatDoubleToString(img[i]));
        }
    }

    static void makeSineTable(int n, double[] sinTable)
    {
        int i, n2, n4, n8;
        double c, s, dc, ds, t;
        n2 = n / 2;  n4 = n / 4;  n8 = n / 8;
        t = Math.sin(pi / n);
        dc = 2 * t * t;  ds = Math.sqrt(dc * (2 - dc));
        t = 2 * dc;  c = sinTable[n4] = 1;  s = sinTable[0] = 0;
        for (i = 1; i < n8; i++) {
            c -= dc;  dc += t * c;
            s += ds;  ds -= t * s;
            sinTable[i] = s;
            sinTable[n4 - i] = c;
        }
        if (n8 != 0) sinTable[n8] = Math.sqrt(0.5);
        for (i = 0; i < n4; i++)
            sinTable[n2 - i] = sinTable[i];
        for (i = 0; i < n2 + n4; i++)
            sinTable[i + n2] = - sinTable[i];
    }

    static void makeBitReverse(int n, int[] bitrev)
    {
        int i, j, k, n2;
        n2 = n / 2;  i = j = 0;
        while(true)
        {
            bitrev[i] = j;
            if (++i >= n) break;
            k = n2;
            while (k <= j) {  j -= k;  k /= 2;  }
            j += k;
        }
    }

    static void fft(int n, double[] x, double[] y)
    {
        int[]   bitrev;
        double[] sintbl;
        int i, j, k, ik, h, d, k2, n4, inverse;
        double t, s, c, dx, dy;

        if(n < 0)
        {
            n = -n;  inverse = 1;
        } else inverse = 0;
        n4 = n / 4;
        if (n == 0) return;
        sintbl = new double[n + n4];
        bitrev = new int[n];

        makeSineTable(n, sintbl);
        makeBitReverse(n, bitrev);
        for(i = 0; i < n; i++)
        {
            j = bitrev[i];
            if (i < j)
            {
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
                for(i = j; i < n; i += k2)
                {
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
    }

    private static String formatDoubleToString(double number)
    {
        return new DecimalFormat("##########.########").format(number);
    }
}
