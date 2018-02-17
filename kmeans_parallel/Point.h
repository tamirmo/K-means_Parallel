#pragma once

#define INITIAL_DISTANCE -1

struct Cluster;

struct Position {
	double x;
	double y;
};

struct Velocity {
	double vx;
	double vy;
};

struct Point {
	Position position;
	Cluster* cluster;
	Velocity velocity;
};

double getPointsDistance(Position* p1, Position* p2);