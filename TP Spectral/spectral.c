#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

#include <sndfile.h>

#include <math.h>
#include <complex.h>

#include <time.h>
/* --------- Fourrier rapide --------- */
#include <complex.h>
#include <fftw3.h>

#include "gnuplot_i.h"


/* taille de la fenetre */
#define	FRAME_SIZE 1024
#define	ZERO_PADDING 0
#define	FFT_SIZE (FRAME_SIZE+ZERO_PADDING)
/* avancement */
#define HOP_SIZE 1024
/* réglage autocorrélation */
#define JUMP 2
/*octave*/
#define OCTAVE 0

static gnuplot_ctrl *h;

static void print_usage (char *progname)
{	printf ("\nUsage : %s <input file> \n", progname) ;
    puts ("\n"
          ) ;

} 

static void fill_buffer(double *buffer, double *new_buffer)
{
    int i;
    double tmp[FRAME_SIZE-HOP_SIZE];

    /* save */
    for (i=0;i<FRAME_SIZE-HOP_SIZE;i++)
        tmp[i] = buffer[i+HOP_SIZE];

    /* save offset */
    for (i=0;i<(FRAME_SIZE-HOP_SIZE);i++)
    {
        buffer[i] = tmp[i];
    }

    for (i=0;i<HOP_SIZE;i++)
    {
        buffer[FRAME_SIZE-HOP_SIZE+i] = new_buffer[i];
    }
}

static int read_n_samples (SNDFILE * infile, double * buffer, int channels, int n)
{

    if (channels == 1)
    {
        /* MONO */
        int readcount ;

        readcount = sf_readf_double (infile, buffer, n);

        return readcount==n;
    }
    else if (channels == 2)
    {
        /* STEREO */
        double buf [2 * n] ;
        int readcount, k ;
        readcount = sf_readf_double (infile, buf, n);
        for (k = 0 ; k < readcount ; k++)
            buffer[k] = (buf [k * 2]+buf [k * 2+1])/2.0 ;

        return readcount==n;
    }
    else
    {
        /* FORMAT ERROR */
        printf ("Channel format error.\n");
    }

    return 0;
} 

static int read_samples (SNDFILE * infile, double * buffer, int channels)
{
    return read_n_samples (infile, buffer, channels, HOP_SIZE);
}

void dft(double* sample, complex * buffer, size_t size , bool inverse)
{
    int op =-1;
    if (inverse)
    {
        op = 1;
    }

    for (int j=0; j< size; j++)
    {
        buffer[j] = 0;
        for (int i=0; i< size; i++)
        {
            buffer[j] += sample[i] * cexp(op* I * (2 *M_PI / size) * j * i);
        }
    }
}

void cartesian_to_polar (complex * Sf, size_t size, double* amp, double* phs)
{
    for (int i=0; i< size; i++)
    {
        amp[i] = cabsf(Sf[i]);
        phs[i] = cargf(Sf[i]);
    }
}

double estimateFreqMax(double* amp, size_t size, int freqIndex, double* d)
{
    (void)size;
    double al = 20. * log10(amp[freqIndex-1]);
    double ac = 20. * log10(amp[freqIndex]);
    double ar = 20. * log10(amp[freqIndex+1]);
    *d = 0.5 * (al - ar) / (al - 2.0 * ac +  ar);
    return (ac-(al-ar) * 0.25*(*d));
}

double maxFreq(double* amp, size_t size, double* freqMax, int sampleRate)
{
    double max =-1;
    int freqIndex = -1;
    for (int i=1; i< (size/2)-1; i++)
    {
        if(amp[i] > max && amp[i] > amp[i-1] && amp[i] > amp[i+1])
        {
            max = amp[i];
            freqIndex = i;
        }
    }
    double d;
    max = estimateFreqMax(amp, size, freqIndex, &d);
    *freqMax = (freqIndex+d)  * sampleRate  / (size*1.0);
    return max;
}

double maxFreqBis(double* amp, size_t size, double* freqMax, int sampleRate)
{
    double max =-1;
    int freqIndex = -1;
    double maxBis = -1;
    int freqIndexBis = -1;
    for (int i=1; i< (size/2)-1; i++)
    {
        
        if(amp[i] > max && amp[i] > amp[i-1] && amp[i] > amp[i+1])
        {
            maxBis = max;
            freqIndexBis = freqIndex;
            max = amp[i];
            freqIndex = i;
        }else{
            if(amp[i] > maxBis && amp[i] > amp[i-1] && amp[i] > amp[i+1]){
                maxBis = amp[i];
                freqIndexBis = i;
            }
        }
        
        
    }
    double d;
    maxBis = estimateFreqMax(amp, size, freqIndexBis, &d);
    *freqMax = (freqIndexBis+d)  * sampleRate  / (size*1.0);
    return maxBis;
}


double fenetreHann(int index, size_t size)
{
    return 0.5 - 0.5* cos(2 *M_PI *index /(size*1.0));
}
int max(double* r1,int trunc_index){
    double max = 0;
    int max_index = 0;
    for(int i=trunc_index;i<FRAME_SIZE;i++)
        if(r1[i]>max){ max = r1[i]; max_index=i;}
    return max_index;
}

int max_with_jump(double* r1){
    int pics_visited = 0;
    for(int i=0;i<FRAME_SIZE;i++){
        if(r1[i]>0 && r1[i+1]<0) pics_visited+=1;
        if(pics_visited == JUMP) return max(r1,i);
    }
    return 0;
}
int notes(double* fourier_samples,char** notes){
    int notes_found = -1;
    for(int freq=15;freq<FFT_SIZE/2-1;freq++){
        printf("Sample %d\n",freq);
        if(fourier_samples[freq]>0.02 && fourier_samples[freq]>fourier_samples[freq-1] && fourier_samples[freq]>fourier_samples[freq+1]){
            printf("Testing notes\n");
            if(freq>15*pow(2,OCTAVE) && freq<17*pow(2,OCTAVE)){
                notes[notes_found]="Do";
            }
            if(freq>17*pow(2,OCTAVE) && freq<18*pow(2,OCTAVE) ){
                notes[notes_found]="Do#";
            }
            if(freq>18*pow(2,OCTAVE) && freq<19*pow(2,OCTAVE)){
                notes[notes_found]="Ré";
            }
            if(freq>19*pow(2,OCTAVE) && freq<20*pow(2,OCTAVE)){
                notes[notes_found]="Ré#";
            }
            if(freq>20*pow(2,OCTAVE) && freq<21*pow(2,OCTAVE) ){
                notes[notes_found]="Mi";
            }
            if(freq>21*pow(2,OCTAVE) && freq<22.5*pow(2,OCTAVE) ){
                notes[notes_found]="Fa";
            }
            if(freq>22.5*pow(2,OCTAVE) && freq<24*pow(2,OCTAVE) ){
                notes[notes_found]="Fa#";
            }
            if(freq>24*pow(2,OCTAVE) && freq<25.5*pow(2,OCTAVE) ){
                notes[notes_found]="Sol";
            }
            if(freq>25.5*pow(2,OCTAVE) && freq<26.8*pow(2,OCTAVE) ){
                notes[notes_found]="Sol#";
            }
            if(freq>26.8*pow(2,OCTAVE) && freq<28*pow(2,OCTAVE) ){
                notes[notes_found]="La";
            }
            if(freq>28*pow(2,OCTAVE) && freq<29.5*pow(2,OCTAVE) ){
                notes[notes_found]="La#";
            }
            if(freq>29.5*pow(2,OCTAVE) && freq<32*pow(2,OCTAVE) ){
                notes[notes_found]="Si";
            }
            notes_found++;
        }
    }
    return notes_found;
}


int main (int argc, char * argv [])
{
    char 		*progname, *infilename;
    SNDFILE	 	*infile = NULL ;
    SF_INFO	 	sfinfo ;

    progname = strrchr (argv [0], '/') ;
    progname = progname ? progname + 1 : argv [0] ;

    if (argc != 2)
    {	print_usage (progname) ;
        return 1 ;
    } ;

    infilename = argv [1] ;

    if ((infile = sf_open (infilename, SFM_READ, &sfinfo)) == NULL)
    {	printf ("Not able to open input file %s.\n", infilename) ;
        puts (sf_strerror (NULL)) ;
        return 1 ;
    } ;

    /* Read WAV */
    int nb_frames = 0;
    double new_buffer[HOP_SIZE];
    double buffer[FRAME_SIZE];

    /* Plot Init */
    h=gnuplot_init();
    gnuplot_setstyle(h, "lines");

    int i;
    for (i=0;i<(FRAME_SIZE/HOP_SIZE-1);i++)
    {
        if (read_samples (infile, new_buffer, sfinfo.channels)==1)
            fill_buffer(buffer, new_buffer);
        else
        {
            printf("not enough samples !!\n");
            return 1;
        }
    }

    /* Info file */
    printf("sample rate %d\n", sfinfo.samplerate);
    printf("channels %d\n", sfinfo.channels);
    printf("size %d\n", (int)sfinfo.frames);

    /*FFTW Init*/
    fftw_complex data_in[FFT_SIZE];
    fftw_complex fourrier[FFT_SIZE];

    double ampSpectre[FFT_SIZE];
    double phaseSpectre[FFT_SIZE];

    fftw_complex data_out[FFT_SIZE];
    fftw_plan tofreq= fftw_plan_dft_1d(FFT_SIZE, data_in, fourrier, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan toreal =fftw_plan_dft_1d(FFT_SIZE, fourrier, data_out, FFTW_BACKWARD, FFTW_ESTIMATE);

    /*Clock*/
    clock_t t1,t2;
    double delta_t=0;

    /*double numbers[12];
    double numbersbis[12];*/
    while (read_samples (infile, new_buffer, sfinfo.channels)==1)
    {
        /* Process Samples */
        printf("Processing frame %d\n",nb_frames);

        /* hop size */
        fill_buffer(buffer, new_buffer);

        /* Fourrier */
        
        t1=clock();
        dft(buffer, data_in, FRAME_SIZE, false);
        t2=clock();
        delta_t += t2-t1;
        
       //Auto-corrélation
        /*double r1[FRAME_SIZE];
        for(int to=0;to<FRAME_SIZE;to++){
            double tmp = 0;
            for(int n=0;n<FRAME_SIZE-to;n++){
                tmp+=buffer[n] * buffer[n+to];
            }
            r1[to] = (1.0/FRAME_SIZE) * tmp;
        }
        int max = max_with_jump(r1);
        printf("Max index: %d\n",max);
        double freq = 44100.0/max;*/
        
        /* Fourrier rapide */

        int i;
        for(i=0;i<FRAME_SIZE;i++){
            data_in[i]=buffer[i] * fenetreHann(i, FRAME_SIZE);
        }
        //Padding zero
        for(; i<FFT_SIZE; i++){
            data_in[i]= 0;
        }
        
        t1=clock();
        fftw_execute(tofreq);
        t2=clock();
        delta_t += t2-t1;
        

        /* Visu amplitude */

        cartesian_to_polar(fourrier, FFT_SIZE, ampSpectre, phaseSpectre);

        /* Frequence max */

        double freqMax;
        double ampMax = maxFreq(ampSpectre, FFT_SIZE, &freqMax, sfinfo.samplerate);
        double freqMaxBis;
        double ampMax_bis = maxFreqBis(ampSpectre, FFT_SIZE, &freqMaxBis, sfinfo.samplerate);
        printf("Frequence max : %lf Hertz, value %lf \n", freqMax, ampMax);
        printf("Frequence max bis : %lf Hertz, value %lf\n",freqMaxBis, ampMax_bis);
        
        //Normalisation amplitude
        
        for (int i=0; i< FFT_SIZE; i++)
        {
            ampSpectre[i] = ampSpectre[i]/(FFT_SIZE/2);
        }
        char** notes_char;
        int notes_found = notes(ampSpectre,notes_char);
        printf("Accord : ");
        for(int i=0;i<notes_found;i++) printf("%s ",notes_char[i]);
        printf("\n");
        
        /* FFTW to real */
        /*
        fftw_execute(toreal);
        for(int i=0;i<FRAME_SIZE;i++){
            buffer[i]=creal(data_in[i]/FRAME_SIZE);

        }
        */

        /* PLOT */
        double x_axis[FRAME_SIZE];
        // Temporal
        /*
        for (int i=0 ; i<N ; i++) {
            x_axis[i] = i * 1.0 / SAMPLING_RATE ;
        }
        */

        // Frequential
        /*for (int i=0 ; i<FFT_SIZE ; i++) {
            x_axis[i] = i * sfinfo.samplerate / FFT_SIZE;
        }*/
        for (int i=0 ; i<FFT_SIZE/2 ; i++) {
            x_axis[i] = i;
        }
        gnuplot_resetplot(h);
        //gnuplot_plot_xy(h, x_axis, buffer,FRAME_SIZE,"temporal frame");
        gnuplot_plot_xy(h, x_axis, ampSpectre, FRAME_SIZE, "Autocorrélation frame");
        usleep(200*(10e3));

        nb_frames++;
    }
    
    fftw_destroy_plan(tofreq);
    fftw_destroy_plan(toreal);

    sf_close (infile) ;
    return 0 ;
} /* main */

