// LMS-spec~ - spectral envelope estimation via LMS filtering
//           - enables cross-synthesis by using obtained spectral estimate
//             as filter frequency response
//
// (c) 2014 greg surges
// surgesg@gmail.com
// http://www.gregsurges.com/ 

#include <flext.h>
#include <stdlib.h>
#include <math.h>
 
#if !defined(FLEXT_VERSION) || (FLEXT_VERSION < 401)
#error You need at least flext version 0.4.1
#endif


class allpasscascade: 
	public flext_dsp
{
	FLEXT_HEADER(allpasscascade, flext_dsp)

	public:
		allpasscascade(int argc, t_atom *argv)
		{
            // arguments [allpasscascade~ u L]
            // or just [allpasscascade~ u] 

			AddInSignal("filter audio in");       // left audio in
			AddInSignal("center freq");      // right audio in
			AddInSignal("bandwidth");      // right audio in
			AddOutSignal("output");          // 1 audio out 
            
            c = 0.0;
            d = 0.0;

            // handle variable argument numbers
            switch(argc){
                case 3:
                    n_filters = GetInt(argv[0]);
                    f_c = GetFloat(argv[1]);
                    f_b = GetFloat(argv[2]);
                    break;
                case 2:
                    n_filters = GetInt(argv[0]);
                    f_c = GetFloat(argv[1]);
                    f_b = 500;
                    break;
                case 1:
                    n_filters = GetInt(argv[0]);
                    f_c = 1000;
                    f_b = 500;
                    break;
                case 0:
                    n_filters = 10;
                    f_c = 1000;
                    f_b = 500;
                    break;
                default:
                    break;
            }

            post("n filters %d", n_filters);
       
            z_n1 = new double[n_filters];
            z_n2 = new double[n_filters];
            for(int i = 0; i < n_filters; i++){
                z_n1[i] = z_n2[i] = 0;
            }

            fs = Samplerate();
			
			FLEXT_ADDMETHOD_(2, "init", m_init);
            		
			// We're done constructing:
			post("-- allpass-cascade~ -- \n 2014 greg surges \n surgesg@gmail.com");
			
		}
	protected:
		virtual void m_signal(int n, float *const *in, float *const *out);
        void m_init();
        ~allpasscascade(); // destructor
	private:	
        float c, d; // filter coefficients
        int n_filters; // number of filters in cascade
        float f_c, f_b; // filter center freq and bandwidth params
        double xmat[3];
        double ymat[3];
        double rmat1[3][3];
        double rmat2[3][3];
        double tmp[3];
        double y;
        double *z_n1;
        double *z_n2;
        int fs; 
		FLEXT_CALLBACK(m_init)
}; 

FLEXT_NEW_DSP_V("allpass-cascade~ allpasscascade~", allpasscascade)

allpasscascade::~allpasscascade(){
    delete[] z_n1;
    delete[] z_n2;
}

void allpasscascade::m_init(){

}


void allpasscascade::m_signal(int n, float *const *in, float *const *out)
{
	float *insig    =  InSig(0);
	float *fc_mod =  InSig(1);
	float *fb_mod =  InSig(2);
	
	float *outs1          = OutSig(0);
    float x_n, x_n_1, x_n_2;
    float y_n, y_n_1, y_n_2;
    double r1, r2;
    r1 = 0;
    r2 = 0;
    int s = 0; // sample counter
	
	while (n--)
	{
        // update coefficients
        float fc = *fc_mod++;
        float fb = *fb_mod++;
        if(fc <= 0){
            fc = 1;
        }
        if(fb <= 0){
            fb = 1;
        }
        //r1 = fb / (fs/2) * M_PI;        
        //r2 = fc / (fs/2) * M_PI;        
        d = -cos(2*M_PI*fc/fs);
        c = (tan(M_PI*fb/fs) - 1.0) / (tan(M_PI*fb/fs)+1.0);
        r1 = acos(-c);
        r2 = acos(-d);
        //post("fb: %f, fc: %f, c: %f, d: %f", fb, fc, c, d);
        //post("%f", r1);
        //post("%f", r2);

        // compute first filter
        // and update state memory
        rmat1[0][0] = cos(r1);
        rmat1[0][1] = -sin(r1);
        rmat1[0][2] = 0;
        rmat1[1][0] = sin(r1);
        rmat1[1][1] = cos(r1);
        rmat1[1][2] = 0;
        rmat1[2][0] = 0;
        rmat1[2][1] = 0;
        rmat1[2][2] = 1;

        rmat2[0][0] = 1;
        rmat2[0][1] = 0;
        rmat2[0][2] = 0;
        rmat2[1][0] = 0;
        rmat2[1][1] = cos(r2);
        rmat2[1][2] = -sin(r2);
        rmat2[2][0] = 0;
        rmat2[2][1] = sin(r2);
        rmat2[2][2] = cos(r2);

        /*
        r1 = fb;
        r2 = cos(2*M_PI*fc/fs);

        rmat1[0][0] = r1;
        rmat1[0][1] = sqrt(1.0 - r1*r1);
        rmat1[0][2] = 0;
        rmat1[1][0] = sqrt(1.0 - r1*r1);
        rmat1[1][1] = r1;
        rmat1[1][2] = 0;
        rmat1[2][0] = 0;
        rmat1[2][1] = 0;
        rmat1[2][2] = 1;

        rmat2[0][0] = 1;
        rmat2[0][1] = 0;
        rmat2[0][2] = 0;
        rmat2[1][0] = 0;
        rmat2[1][1] = r2;
        rmat2[1][2] = sqrt(1.0 - r2*r2);
        rmat2[2][0] = 0;
        rmat2[2][1] = sqrt(1.0 - r2*r2);
        rmat2[2][2] = r2;
        */






        /*
        post("-----------------");
        post("%f %f %f", rmat1[0][0], rmat1[0][1], rmat1[0][2]);
        post("%f %f %f", rmat1[1][0], rmat1[1][1], rmat1[1][2]);
        post("%f %f %f", rmat1[2][0], rmat1[2][1], rmat1[2][2]);
        */

        // compute consecutive filters
        // and update state memory

        xmat[0] = insig[s];
        for(int i = 0; i < n_filters; i++){
            xmat[1] = z_n1[i];
            xmat[2] = z_n2[i];
            // tmp = xmat * rmat1
            tmp[0] = rmat1[0][0]*xmat[0] + rmat1[0][1]*xmat[1] + rmat1[0][2]*xmat[2];
            tmp[1] = rmat1[1][0]*xmat[0] + rmat1[1][1]*xmat[1] + rmat1[1][2]*xmat[2];
            tmp[2] = rmat1[2][0]*xmat[0] + rmat1[2][1]*xmat[1] + rmat1[2][2]*xmat[2];

            // ymat = tmp * rmat2
            ymat[0] = rmat2[0][0]*tmp[0] + rmat2[0][1]*tmp[1] + rmat2[0][2]*tmp[2];
            ymat[1] = rmat2[1][0]*tmp[0] + rmat2[1][1]*tmp[1] + rmat2[1][2]*tmp[2];
            ymat[2] = rmat2[2][0]*tmp[0] + rmat2[2][1]*tmp[1] + rmat2[2][2]*tmp[2];

            z_n1[i] = ymat[1];
            z_n2[i] = ymat[2];

            xmat[0] = ymat[0];
        }
        //
        //post("%f, %f", xmat[0], ymat[0]);

        outs1[s] = ymat[0];
        s++;
	}
}
