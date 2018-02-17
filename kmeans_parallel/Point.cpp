#include "Point.h"
#include <math.h>

// Calculates distance between two points
double getPointsDistance(Position* p1, Position* p2){
	double x = p2->x - p1->x;
	double y = p2->y - p1->y;
	return sqrt(x*x + y*y);
}