#pragma once
#include <iostream>
#include <algorithm>
using namespace std;

#define MINP -1
#define MAXP 4

unsigned int expandBits(unsigned int v);

unsigned int morton3D(float x, float y, float z);

unsigned int getMortonCode(float x, float y, float z);
