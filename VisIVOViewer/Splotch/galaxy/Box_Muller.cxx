# include <math.h>
# include <stdlib.h>
# include <stdio.h>


float box_muller(float m, float s)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
	float x1, x2, w, y1;
	static float y2;
	static int use_last = 0;
        float x_rand_max = (float) RAND_MAX;
        float ran_result;
        float x1ran, x2ran;


	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1ran = ((float)rand())/x_rand_max;
			x2ran = ((float)rand())/x_rand_max;
			x1 = 2.0 * x1ran - 1.0;
			x2 = 2.0 * x2ran - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	ran_result = ( m + y1 * s );

	return( ran_result );
}

float box_uniform(float m, float s)	/* normal random variate generator */
{				        /* mean m, uniform within s */
  float y1 = ((float)rand())/(float)RAND_MAX-0.5;
  float ran_result = ( m + y1 * s );

  return( ran_result );
}
