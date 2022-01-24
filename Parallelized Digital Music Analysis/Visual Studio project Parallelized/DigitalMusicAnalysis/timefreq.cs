using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading;

namespace DigitalMusicAnalysis
{
    public class timefreq
    {
        public float[][] timeFreqData;
        public int wSamp;
        public Complex[] twiddles;
        private delegate float parDeg(int start, int size);
        public timefreq(float[] x, int windowSamp)
        {
            int ii;
            double pi = 3.14159265;
            Complex i = Complex.ImaginaryOne;
            this.wSamp = windowSamp;
            twiddles = new Complex[wSamp];
            for (ii = 0; ii < wSamp; ii++)
            {
                double a = 2 * pi * ii / (double)wSamp;
                twiddles[ii] = Complex.Pow(Complex.Exp(-i), (float)a);
            }

            timeFreqData = new float[wSamp / 2][];

            int nearest = (int)Math.Ceiling((double)x.Length / (double)wSamp);
            nearest = nearest * wSamp;

            Complex[] compX = new Complex[nearest];
            for (int kk = 0; kk < nearest; kk++)
            {
                if (kk < x.Length)
                {
                    compX[kk] = x[kk];
                }
                else
                {
                    compX[kk] = Complex.Zero;
                }
            }


            int cols = 2 * nearest / wSamp;

            for (int jj = 0; jj < wSamp / 2; jj++)
            {
                timeFreqData[jj] = new float[cols];
            }

            timeFreqData = stft(compX, wSamp);

        }

        float[][] stft(Complex[] x, int wSamp)
        {

            int core_count =4;
            int N = x.Length;
            float fftMax = 0;
            float[][] Y = new float[wSamp / 2][];
            float[] fftMax_array = new float[core_count];
            int max_size = (int)(2 * Math.Floor((double)N / (double)wSamp) - 1);
            int split = (int)Math.Ceiling((double)max_size / (double)core_count);
            //parDeg pardeg = par;
            for (int ll = 0; ll < wSamp / 2; ll++)
            {
                Y[ll] = new float[2 * (int)Math.Floor((double)N / (double)wSamp)];
            }
            using (var countdownEvent = new CountdownEvent(core_count))
            {
                for (int i = 0; i < core_count; i++)
                {
                    int core_n = i;
                    int start_loc = split * core_n;
                    ThreadPool.QueueUserWorkItem(
                         a =>
                         {
                             fftMax_array[core_n] = par(start_loc, split);
                             countdownEvent.Signal();
                             Console.WriteLine(core_n);
                         });
                }
                countdownEvent.Wait();
            }

            for (int i = 0; i < core_count; i++)
            {
                if (fftMax_array[i] > fftMax)
                {
                    fftMax = fftMax_array[i];
                }
            }
            for (int ii = 0; ii < 2 * Math.Floor((double)N / (double)wSamp) - 1; ii++)
            {
                for (int kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] /= fftMax;
                }
            }

            return Y;
            float par(int start, int size)
            {
                float local_max = 0;
                int stop = Math.Min(start + size, max_size - 1);
                Complex[] temp = new Complex[wSamp];
                Complex[] tempFFT = new Complex[wSamp];
                //Console.WriteLine("start: {0}, end: {1}", start, stop);
                for (int ii = start; ii < stop; ii++)
                {

                    for (int jj = 0; jj < wSamp; jj++)
                    {
                        temp[jj] = x[ii * (wSamp / 2) + jj];
                    }
                    //Console.WriteLine("{0}", ii);
                    tempFFT = fft(temp);

                    for (int kk = 0; kk < ((wSamp / 2) - 1); kk++)
                    {
                        Y[kk][ii] = (float)Complex.Abs(tempFFT[kk]);

                        if (Y[kk][ii] > local_max)
                        {
                            //Console.WriteLine(Y[kk][ii]);
                            local_max = Y[kk][ii];
                        }
                    }


                }
                return local_max;
            }
        }
        //new fft function created removing recursion 
        Complex[] fft(Complex[] x)
        {
            int N = x.Length;
            int aSize = N;
            int gap = N / 2;
            int subArrays = N / aSize;
            int subSize = 2;
            Complex E;
            Complex O;
            Complex[] Y = new Complex[N];
            int mul1;
            int mul2;
            int[] mulArray = new int[N];
            for (int i = 0; i < N; i++)
            {
                mulArray[i] = 0;
            }
            for (int split = N / 2; split > 0; split /= 2)
            {
                for (int s = 0; s < subArrays; s++)
                {
                    for (int p = 0; p < gap; p++)
                    {
                        E = x[p + s * aSize];
                        O = x[p + gap + s * aSize];
                        mulArray[p + gap + s * aSize] += subArrays;
                        mul1 = mulArray[p + s * aSize];
                        mul2 = mulArray[p + gap + s * aSize];
                        if (split == 1)
                        {
                            Y[mul1] = E + O * twiddles[mul1 * wSamp / subSize];
                            Y[mul2] = E + O * twiddles[mul2 * wSamp / subSize];
                        }
                        else
                        {
                            x[p + s * aSize] = E + O * twiddles[mul1 * wSamp / subSize];
                            x[p + gap + s * aSize] = E + O * twiddles[mul2 * wSamp / subSize];
                        }
                    }
                }
                subSize *= 2;
                aSize /= 2;
                gap /= 2;
                subArrays = N / aSize;

            }
            return Y;
        }

    }
}

