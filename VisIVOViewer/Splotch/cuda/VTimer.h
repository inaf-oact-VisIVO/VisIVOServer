/*
Copyright things go here.
*/
#pragma once
#include <windows.h>
#if !defined(_M_IX86)
 typedef __int64 LONGLONG; 
#else
// typedef double LONGLONG;
#endif

typedef union _LARGE_INTEGER1 {  
	struct {    
	unsigned long LowPart;    
	long HighPart;  
	};  
	struct {    
	unsigned long LowPart;    
	long HighPart;  
	} u;  
	LONGLONG QuadPart;
} LARGE_INTEGER1; 

class VTimer
{
private:
	    double	time0, elapsed;
        bool	running;
        LARGE_INTEGER1	performance_counter;
		LARGE_INTEGER1	performance_frequency;
        double	one_over_frequency;
		
		double fetchSysTime();
public:
		VTimer(void);
		~VTimer(void);
	
		void	start();
		void	stop();
		void	reset();
		double	getTime();
};