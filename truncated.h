#define TAIL_LIMIT  2.983851594898812
#define EXPTLOVTL  0.00390733345803262
#define LOG_TAIL_LIMIT  1.09321494749176
#define SQ_TAIL_LIMIT  4.45168517019009
#define MAX_PTAIL  0.00390733345803262
#define basep   0.0078125
#define infinity (std::numeric_limits<double>::infinity() )
#define NAN  (std::numeric_limits<double>::quiet_NaN())
#define M_PI 3.141592653589793 

#define RTAIL_BIN  162
#define LTAIL_BIN  (-163)

class truncated
{

private:
	double *xhist;
	double *yhist;
	double *dxhist;
	int *yhistratio;
	int *whichbin;
	gsl_rng * rr;

	inline double rnd()
	{
		return (rand()/(RAND_MAX+1.0));
	}
	inline double itaillimit( double b)
	{
		return EXPTLOVTL - exp(SQ_TAIL_LIMIT - TAIL_LIMIT * b - LOG_TAIL_LIMIT );
	}

	inline double sample_tail( double a, double b)	
	{
		double sa =  1 ;
		double sb = exp(a*(a-b));
		if ( sa == sb )
			return (a*a<b*b?a:b);

		double u;
		do {
			u = rnd() * (sb-sa) + sa;
		} while( u == 0 );		
		u = a - log(u)/a;

		double v = rnd() * exp(0.5*(a-u)*(a-u));
		if ( v <= 1 )
			return u;
		else
			return NAN;
	}
	double truncate2(double a, double b)
	{
		double z;
		int y,ia,ib;
		double ah,bh;
		ah = a; bh = b;
		
		double ap = basep;
		double bp = basep;
		
		double p_tail_left = 0;
		if ( a <= -TAIL_LIMIT ) {
			ia = LTAIL_BIN;
			p_tail_left = MAX_PTAIL;
			ah = -TAIL_LIMIT;
			ap = 0;
		}
		else if ( a >= TAIL_LIMIT ) {
			do { 
				z = sample_tail( a, b );
			} while ( z != z );
			return z;
		}
		else {
			double ax = fabs(a);
			ia = whichbin[ (int)(ax*128.0) ];
			if (xhist[ia+1] <= ax )
				++ia;
			ia = ((a>=0)?ia:-ia-1);
		}
		
		double p_tail_right = 0;
		if (b <= -TAIL_LIMIT ) {
			do {
				z = -sample_tail(-b,-a);
			} while (z != z);
			return z;
		}
		else if ( b >= TAIL_LIMIT ) {
			p_tail_right = MAX_PTAIL;
			ib = RTAIL_BIN;
			bh = TAIL_LIMIT;
			bp = 0;
		}
		else {
			double bx = fabs(b);
			ib = whichbin[ (int)(bx*128.0) ];
			if (xhist[ib+1] <= bx )
				++ib;
			ib = ((b>=0)?ib:-ib-1);
		}
		
		if (ia == ib ) {
			do {
				z = rnd() * (bh-ah) + ah;
				y = rand();
			} while ( y >= yhistratio[ia] && y * yhist[ia] >= exp(-0.5*z*z)*(RAND_MAX+1.0) );
			return(z);
		}
		
		double middlep = (ib-ia-1)*basep;		
		do {			
			double u = rnd() * (middlep + ap + bp + p_tail_left + p_tail_right );
			if ( u < middlep ) {
				u /= basep;
				int ui = (int)u;
				u = modf(u,&z);
				int col = ui + ia  + 1;
				z = u * dxhist[col] + xhist[col];
				y = rand();
				
				if ( (y>=yhistratio[col] ) && y * yhist[col] >= exp(-0.5*z*z) * (RAND_MAX+1.0) )
					z = NAN;
				continue;
			}
			if (ap == basep )
				ap = (xhist[ia+1]-ah)*yhist[ia];
			if (u < middlep + ap ) {
				z = (u-middlep)/ap * (xhist[ia+1]-ah) + ah;
				y = rand();
				if ( y >= yhistratio[ia] && y * yhist[ia] >= exp(-0.5*z*z) * (RAND_MAX+1.0) )
					z = NAN;
				continue;
			}			
			if (bp == basep)
				bp = (bh-xhist[ib])*yhist[ib];
			if ( u < middlep + ap + bp) {
				z = (u-middlep-ap)/bp * (bh-xhist[ib])+xhist[ib];
				y = rand();
				if ( y >= yhistratio[ib] && y * yhist[ib] >= exp(-0.5*z*z) * (RAND_MAX+1.0) )
					z = NAN;
				continue;
			}			
			if ( p_tail_left == MAX_PTAIL )
				p_tail_left = itaillimit(-a);
			if ( u < middlep + ap + bp + p_tail_left ) {
				z = -sample_tail( TAIL_LIMIT, -a );
				continue;
			}
			if ( p_tail_right == MAX_PTAIL )
				p_tail_right = itaillimit(b);			
			if ( u < middlep + ap + bp + p_tail_left + p_tail_right ) {
				z = sample_tail( TAIL_LIMIT, b );
				continue;
			}			
			z = NAN;
		} while ( z != z );	
		return z;					
	}
	
public:	

	double draw( double mu, double sigma, double aa, double bb )
	{
		
		double a = (aa-mu)/sigma;
		double b = (bb-mu)/sigma;

		if ( a == -infinity && b == infinity )
			return NAN;		
		else {
			double z = sigma * truncate2(a,b) + mu;
			if ( z < aa - 0.0000001|| z > bb +  0.0000001 )
				std::cout << "truncated failed big time! " << z << " not in [" << aa << ", " << bb << "] mu=" << mu << " sigma=" << sigma <<  std::endl;
			return z;
		}		
	}
			
	// creates the object by initializing the arrays we will need
	truncated()
	{
	
		rr = gsl_rng_alloc( gsl_rng_rand );
	
		double x = 0;
		int i = 0, c = 0;
		while( x < 3.5 ) {
			x += basep * exp(0.5*x*x);
			++c;
		}
		xhist = (double*) malloc((2*c+1)*sizeof(double)) + c;
		dxhist =(double*) malloc(2*c*sizeof(double)) + c;
		yhist = (double*)malloc(2*(c+1)*sizeof(double)) + (c+1);

		yhistratio = (int*)malloc(2*c*sizeof(int)) + c;		
		whichbin = (int*)malloc( ((int)(TAIL_LIMIT * 128.0 )+1) * sizeof(int));
		
		x = 0; xhist[0] = 0; yhist[0] = 1; yhist[-1] = 1;
		
		for( i = 1; i <= c; ++i )  {
			x += basep * exp(0.5*x*x);	
			xhist[i] = x;
			xhist[-i] = -x;
			dxhist[i-1] = xhist[i] - xhist[i-1];
			dxhist[-i] = dxhist[i-1];			
			yhist[i] = exp(-0.5*x*x);
			yhist[-i-1] = yhist[i];
			yhistratio[i-1] = (int)ceil( yhist[i]/yhist[i-1] * (RAND_MAX+1.0) );
			yhistratio[-i] = yhistratio[i-1];
		}
		c = 0;
		for( i = 0; i <= (int)(TAIL_LIMIT / basep); ++i ) {
			while( (int)(xhist[c] / basep ) < i )
				++c;
			whichbin[i] = c-1;
		}
	}	
};
