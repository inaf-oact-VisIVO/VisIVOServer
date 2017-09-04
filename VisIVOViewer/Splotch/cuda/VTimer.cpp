/*
Copyright things go here.
*/
#include "VTimer.h"
//#include "windows.h"

VTimer::VTimer(void)
{
	QueryPerformanceFrequency( (LARGE_INTEGER *)&performance_frequency );
	one_over_frequency = 1.0/((double)performance_frequency.QuadPart);
	time0 = elapsed = 0.0;
	running = false;
}

VTimer::~VTimer(void)
{
}
		
double VTimer::fetchSysTime()
{
	QueryPerformanceCounter( (LARGE_INTEGER *)&performance_counter );

/*uncomment this to debug 
	TRACE("\nget time:low:%d, high:%d, whole:%f, time:%f", 
		  t->performance_counter.LowPart,
		  t->performance_counter.HighPart,
		  (double) t->performance_counter.QuadPart, 
		  (double) t->performance_counter.QuadPart*t->one_over_frequency);
*/
	return ((double)performance_counter.QuadPart) * one_over_frequency;
}

void	VTimer::start()
{
    running = true;
    time0 = fetchSysTime();
}

void	VTimer::stop()
{
    running = false;
    elapsed += fetchSysTime() - time0;
//TRACE("\nint stoptimer() time0:%f",t->time0);
}

void	VTimer::reset()
{
    running = false;
    elapsed = 0.0;
}

double	VTimer::getTime()
{
	if ( running ) {
		stop();
		//TRACE("\nafter stop timer:%f", t->time0) ;
		start();
		//TRACE("\nafter start timer%f", t->time0) ;
	}
	return elapsed;
}